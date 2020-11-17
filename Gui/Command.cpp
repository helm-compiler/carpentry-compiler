/***************************************************************************
 *   Copyright (c) 2008 Jürgen Riegel (juergen.riegel@web.de)              *
 *                                                                         *
 *   This file is part of the FreeCAD CAx development system.              *
 *                                                                         *
 *   This library is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Library General Public           *
 *   License as published by the Free Software Foundation; either          *
 *   version 2 of the License, or (at your option) any later version.      *
 *                                                                         *
 *   This library  is distributed in the hope that it will be useful,      *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with this library; see the file COPYING.LIB. If not,    *
 *   write to the Free Software Foundation, Inc., 59 Temple Place,         *
 *   Suite 330, Boston, MA  02111-1307, USA                                *
 *                                                                         *
 ***************************************************************************/

#include "PreCompiled.h"

#ifndef _PreComp_
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <GeomLib_IsPlanarSurface.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <QMessageBox>
#include <QTextStream>
#include <Inventor/nodes/SoCamera.h>
#endif
#include <QProgressDialog>
#include <sstream>
#include <algorithm>

#include <App/DocumentObjectGroup.h>
#include <App/Origin.h>
#include <App/OriginFeature.h>
#include <App/Part.h>
#include <Base/Interpreter.h>
#include <Gui/Application.h>
#include <Gui/Command.h>
#include <Gui/Control.h>
#include <Gui/Selection.h>
#include <Gui/MainWindow.h>
#include <Gui/Document.h>
#include <Gui/View3DInventor.h>
#include <Gui/View3DInventorViewer.h>
#include <Gui/WorkbenchManager.h>
#include <Mod/Sketcher/App/SketchObject.h>

#include <Mod/PartDesign/App/Body.h>

#include <Mod/PartDesign/App/FeaturePrimitive.h>
#include <Mod/PartDesign/App/FeatureHole.h>
#include <Mod/PartDesign/App/FeaturePocket.h>
#include <Mod/PartDesign/Gui/PHELM.h>
#include <Mod/PartDesign/Gui/ToolsConfiguration.h>
#include <Mod/PartDesign/Gui/CodeEditor.h>
#include <Mod/PartDesign/Gui/GeomCommonFunctions.h>
#include <Mod/PartDesign/Gui/HLParser.h>
#include <Mod/PartDesign/Gui/CompilerConfig.h>

#include "TaskFeaturePick.h"
#include "ReferenceSelection.h"
#include "Utils.h"
#include "WorkflowManager.h"
#include "Workbench.h"
#include "Compiler.h"
#include "opt/optimizer.hpp"

 // TODO Remove this header after fixing code so it won;t be needed here (2015-10-20, Fat-Zer)
#include "ui_DlgReference.h"
#include <QInputDialog>
#include "GC_MakeArcOfCircle.hxx"
#include "GC_MakeSegment.hxx"
#include "BRepPrimAPI_MakePrism.hxx"
#include <boost/lexical_cast.hpp>

using namespace OptCompiler;
//#include <boost/graph/isomorphism.hpp>
//
//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
//typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;


FC_LOG_LEVEL_INIT("PartDesign", true, true);

using namespace std;
using namespace Attacher;

TopoDS_Shape MakePart(const Standard_Real& part_length, const Standard_Integer& end_angle_0, const Standard_Integer& end_angle_1);

//===========================================================================
// PartDesign_Datum
//===========================================================================

/**
 * @brief UnifiedDatumCommand is a common routine called by datum plane, line and point commands
 * @param cmd (i/o) command, to have shortcuts to doCommand, etc.
 * @param type (input)
 * @param name (input). Is used to generate new name for an object, and to fill undo messages.
 *
 */
void UnifiedDatumCommand(Gui::Command& cmd, Base::Type type, std::string name)
{
	try
	{
		std::string fullTypeName(type.getName());

		App::PropertyLinkSubList support;
		cmd.getSelection().getAsPropertyLinkSubList(support);

		bool bEditSelected = false;
		if (support.getSize() == 1 && support.getValue())
		{
			if (support.getValue()->isDerivedFrom(type))
				bEditSelected = true;
		}

		PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ false);

		if (bEditSelected)
		{
			std::string tmp = std::string("Edit ") + name;
			cmd.openCommand(tmp.c_str());
			PartDesignGui::setEdit(support.getValue(), pcActiveBody);
		}
		else if (pcActiveBody)
		{
			// TODO Check how this will work outside of a body (2015-10-20, Fat-Zer)
			std::string FeatName = cmd.getUniqueObjectName(name.c_str(), pcActiveBody);

			std::string tmp = std::string("Create ") + name;

			cmd.openCommand(tmp.c_str());
			FCMD_OBJ_CMD(pcActiveBody, "newObject('" << fullTypeName << "','" << FeatName << "')");

			// remove the body from links in case it's selected as
			// otherwise a cyclic dependency will be created
			support.removeValue(pcActiveBody);

			auto Feat = pcActiveBody->getDocument()->getObject(FeatName.c_str());
			if (!Feat)
				return;

			//test if current selection fits a mode.
			if (support.getSize() > 0)
			{
				Part::AttachExtension* pcDatum = Feat->getExtensionByType<Part::AttachExtension>();
				pcDatum->attacher().setReferences(support);
				SuggestResult sugr;
				pcDatum->attacher().suggestMapModes(sugr);
				if (sugr.message == Attacher::SuggestResult::srOK)
				{
					//fits some mode. Populate support property.
					FCMD_OBJ_CMD(Feat, "Support = " << support.getPyReprString());
					FCMD_OBJ_CMD(Feat, "MapMode = '" << AttachEngine::getModeName(sugr.bestFitMode) << "'");
				}
				else
				{
					QMessageBox::information(Gui::getMainWindow(), QObject::tr("Invalid selection"), QObject::tr("There are no attachment modes that fit selected objects. Select something else."));
				}
			}
			cmd.doCommand(Gui::Command::Doc, "App.activeDocument().recompute()"); // recompute the feature based on its references
			PartDesignGui::setEdit(Feat, pcActiveBody);
		}
		else
		{
			QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Error"), QObject::tr("There is no active body. Please make a body active before inserting a datum entity."));
		}
	}
	catch (Base::Exception & e)
	{
		QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Error"), QString::fromLatin1(e.what()));
	}
	catch (Standard_Failure & e)
	{
		QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Error"), QString::fromLatin1(e.GetMessageString()));
	}
}

/* Variable command*/

DEF_STD_CMD_A(CmdPartDesignVariable);

CmdPartDesignVariable::CmdPartDesignVariable()
	: Command("PartDesign_Variable")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Create a new variable");
	sToolTipText = QT_TR_NOOP("Create a new variable");
	sWhatsThis = "PartDesign_Variable";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_Variable";
}

void CmdPartDesignVariable::activated(int iMsg)
{
	Q_UNUSED(iMsg);
	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ true);
	if (pcActiveBody == 0)
		return;

	Gui::Command* cmd = this;
	std::string FeatName = cmd->getUniqueObjectName("Variable", pcActiveBody);
	std::cerr << "Create a feature as " << FeatName << std::endl;
	std::string tmp = std::string("Create ") + "Variable";
	cmd->openCommand(tmp.c_str());
	FCMD_OBJ_CMD(pcActiveBody, "newObject('PartDesign::Variable" << "','" << FeatName << "')");
}

bool CmdPartDesignVariable::isActive(void)
{
	if (getActiveGuiDocument())
		return true;
	else
		return false;
}

//===========================================================================
// PartDesign_Sketch
//===========================================================================

/* Sketch commands =======================================================*/
DEF_STD_CMD_A(CmdPartDesignNewSketch);

CmdPartDesignNewSketch::CmdPartDesignNewSketch()
	: Command("PartDesign_NewSketch")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Create sketch");
	sToolTipText = QT_TR_NOOP("Create a new sketch");
	sWhatsThis = "PartDesign_NewSketch";
	sStatusTip = sToolTipText;
	sPixmap = "Sketcher_NewSketch";
}

void CmdPartDesignNewSketch::activated(int iMsg)
{
}

bool CmdPartDesignNewSketch::isActive(void)
{
	if (getActiveGuiDocument())
		return true;
	else
		return false;
}

//===========================================================================
// Common utility functions for all features creating solids
//===========================================================================

void finishFeature(const Gui::Command* cmd, App::DocumentObject* Feat,
	App::DocumentObject* prevSolidFeature = nullptr,
	const bool hidePrevSolid = true,
	const bool updateDocument = true)
{
	PartDesign::Body* pcActiveBody;

	if (prevSolidFeature)
	{
		pcActiveBody = PartDesignGui::getBodyFor(prevSolidFeature, /*messageIfNot = */ false);
	}
	else
	{ // insert into the same body as the given previous one
		pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ false);
	}

	if (hidePrevSolid && prevSolidFeature && (prevSolidFeature != NULL))
		FCMD_OBJ_HIDE(prevSolidFeature);

	if (updateDocument)
		cmd->updateActive();

	auto base = dynamic_cast<PartDesign::Feature*>(Feat);
	if (base)
		base = dynamic_cast<PartDesign::Feature*>(base->getBaseObject(true));
	App::DocumentObject* obj = base;
	if (!obj)
		obj = pcActiveBody;

	// Do this before calling setEdit to avoid to override the 'Shape preview' mode (#0003621)
	if (obj)
	{
		cmd->copyVisual(Feat, "ShapeColor", obj);
		cmd->copyVisual(Feat, "LineColor", obj);
		cmd->copyVisual(Feat, "PointColor", obj);
		cmd->copyVisual(Feat, "Transparency", obj);
		cmd->copyVisual(Feat, "DisplayMode", obj);
	}

	// #0001721: use '0' as edit value to avoid switching off selection in
	// ViewProviderGeometryObject::setEditViewer
	PartDesignGui::setEdit(Feat, pcActiveBody);
	cmd->doCommand(cmd->Gui, "Gui.Selection.clearSelection()");
	//cmd->doCommand(cmd->Gui,"Gui.Selection.addSelection(App.ActiveDocument.ActiveObject)");
}

//===========================================================================
// Common utility functions for ProfileBased features
//===========================================================================

// Take a list of Part2DObjects and classify them for creating a
// ProfileBased feature. FirstFreeSketch is the first free sketch in the same body
// or sketches.end() if non available. The returned number is the amount of free sketches
unsigned validateSketches(std::vector<App::DocumentObject*>& sketches,
	std::vector<PartDesignGui::TaskFeaturePick::featureStatus>& status,
	std::vector<App::DocumentObject*>::iterator& firstFreeSketch)
{
	// TODO Review the function for non-part bodies (2015-09-04, Fat-Zer)
	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(false);
	App::Part* pcActivePart = PartDesignGui::getPartFor(pcActiveBody, false);

	// TODO: If the user previously opted to allow multiple use of sketches or use of sketches from other bodies,
	// then count these as valid sketches!
	unsigned freeSketches = 0;
	firstFreeSketch = sketches.end();

	for (std::vector<App::DocumentObject*>::iterator s = sketches.begin(); s != sketches.end(); s++)
	{
		if (!pcActiveBody)
		{
			// We work in the old style outside any body
			if (PartDesign::Body::findBodyOf(*s))
			{
				status.push_back(PartDesignGui::TaskFeaturePick::otherPart);
				continue;
			}
		}
		else if (!pcActiveBody->hasObject(*s))
		{
			// Check whether this plane belongs to a body of the same part
			PartDesign::Body* b = PartDesign::Body::findBodyOf(*s);
			if (!b)
				status.push_back(PartDesignGui::TaskFeaturePick::notInBody);
			else if (pcActivePart && pcActivePart->hasObject(b, true))
				status.push_back(PartDesignGui::TaskFeaturePick::otherBody);
			else
				status.push_back(PartDesignGui::TaskFeaturePick::otherPart);

			continue;
		}

		//Base::Console().Error("Checking sketch %s\n", (*s)->getNameInDocument());
		// Check whether this sketch is already being used by another feature
		// Body features don't count...
		std::vector<App::DocumentObject*> inList = (*s)->getInList();
		std::vector<App::DocumentObject*>::iterator o = inList.begin();
		while (o != inList.end())
		{
			//Base::Console().Error("Inlist: %s\n", (*o)->getNameInDocument());
			if ((*o)->getTypeId().isDerivedFrom(PartDesign::Body::getClassTypeId()))
				o = inList.erase(o); //ignore bodies
			else if (!((*o)->getTypeId().isDerivedFrom(PartDesign::Feature::getClassTypeId())))
				o = inList.erase(o); //ignore non-partDesign
			else
				++o;
		}
		if (inList.size() > 0)
		{
			status.push_back(PartDesignGui::TaskFeaturePick::isUsed);
			continue;
		}

		if (pcActiveBody && pcActiveBody->isAfterInsertPoint(*s))
		{
			status.push_back(PartDesignGui::TaskFeaturePick::afterTip);
			continue;
		}

		// Check whether the sketch shape is valid
		Part::Part2DObject* sketch = static_cast<Part::Part2DObject*>(*s);
		const TopoDS_Shape& shape = sketch->Shape.getValue();
		if (shape.IsNull())
		{
			status.push_back(PartDesignGui::TaskFeaturePick::invalidShape);
			continue;
		}

		// count free wires
		int ctWires = 0;
		TopExp_Explorer ex;
		for (ex.Init(shape, TopAbs_WIRE); ex.More(); ex.Next())
		{
			ctWires++;
		}
		if (ctWires == 0)
		{
			status.push_back(PartDesignGui::TaskFeaturePick::noWire);
			continue;
		}

		// All checks passed - found a valid sketch
		if (firstFreeSketch == sketches.end())
			firstFreeSketch = s;
		freeSketches++;
		status.push_back(PartDesignGui::TaskFeaturePick::validFeature);
	}

	return freeSketches;
}

bool base_worker(PartDesign::Body* pcActiveBody, Gui::Command* cmd, const std::string& which,
	const std::string& specifiedName, App::DocumentObject* feature, std::string sub,
	boost::function<void(Part::Feature*, App::DocumentObject*)> func)
{
	if (!feature || !feature->isDerivedFrom(Part::Feature::getClassTypeId()))
		return false;

	// Related to #0002760: when an operation can't be performed due to a broken
	// profile then make sure that it is recomputed when cancelling the operation
	// otherwise it might be impossible to see that it's broken.
	if (feature->isTouched())
		feature->recomputeFeature();

	std::string FeatName = specifiedName;
	if (specifiedName.empty())
	{
		FeatName = cmd->getUniqueObjectName(which.c_str(), pcActiveBody);
		if (which == "Pocket")
			FeatName = cmd->getUniqueObjectName("Cut");
		else
			FeatName = cmd->getUniqueObjectName(which.c_str());
	}
	else
	{
		if (cmd->getObject(FeatName.c_str()) != nullptr) return false;
	}

	Gui::Command::openCommand((std::string("Make ") + which).c_str());

	FCMD_OBJ_CMD(pcActiveBody, "newObject('PartDesign::" << which << "','" << FeatName << "')");
	auto Feat = pcActiveBody->getDocument()->getObject(FeatName.c_str());

	// We would like to add the new features into the workbench and check them whenever it is re-activated
	Gui::WorkbenchManager* pcWbMan = Gui::WorkbenchManager::instance();
	auto cdWb = dynamic_cast<PartDesignGui::Workbench*>(pcWbMan->active());

	if (Feat->isDerivedFrom(PartDesign::Feature::getClassTypeId()))
	{
		auto dFeat = static_cast<PartDesign::Feature*>(Feat);
		cdWb->AddFeatureToCheck(dFeat);
	}

	auto objCmd = Gui::Command::getObjectCmd(feature);
	if (feature->isDerivedFrom(Part::Part2DObject::getClassTypeId()))
	{
		FCMD_OBJ_CMD(Feat, "Profile = " << objCmd);
	}
	else
	{
		FCMD_OBJ_CMD(Feat, "Profile = (" << objCmd << ", ['" << sub << "'])");
	}

	func(static_cast<Part::Feature*>(feature), Feat);
	return true;
}

bool base_worker(PartDesign::Body* pcActiveBody, Gui::Command* cmd, const std::string& which,
	App::DocumentObject* feature, std::string sub,
	boost::function<void(Part::Feature*, App::DocumentObject*)> func)
{
	return base_worker(pcActiveBody, cmd, which, "", feature, sub, func);
}

void prepareProfileBased(PartDesign::Body* pcActiveBody, Gui::Command* cmd, const std::string& which,
	boost::function<void(Part::Feature*, App::DocumentObject*)> func)
{
	//if a profile is selected we can make our life easy and fast
	std::vector<Gui::SelectionObject> selection = cmd->getSelection().getSelectionEx();
	if (!selection.empty() && selection.front().hasSubNames())
	{
		base_worker(pcActiveBody, cmd, which, selection.front().getObject(), selection.front().getSubNames().front(), func);
		return;
	}

	//no face profile was selected, do the extended sketch logic

	bool bNoSketchWasSelected = false;
	// Get a valid sketch from the user
	// First check selections
	std::vector<App::DocumentObject*> sketches = cmd->getSelection().getObjectsOfType(Part::Part2DObject::getClassTypeId());
	if (sketches.empty())
	{ //no sketches were selected. Let user pick an object from valid ones available in document
		sketches = cmd->getDocument()->getObjectsOfType(Part::Part2DObject::getClassTypeId());
		bNoSketchWasSelected = true;
	}

	if (sketches.empty())
	{
		QMessageBox::warning(Gui::getMainWindow(), QObject::tr("No sketch to work on"),
			QObject::tr("No sketch is available in the document"));
		return;
	}

	std::vector<PartDesignGui::TaskFeaturePick::featureStatus> status;
	std::vector<App::DocumentObject*>::iterator firstFreeSketch;
	int freeSketches = validateSketches(sketches, status, firstFreeSketch);

	auto accepter = [=](const std::vector<App::DocumentObject*>& features) -> bool {
		if (features.empty())
			return false;

		return true;
	};

	auto sketch_worker = [&](std::vector<App::DocumentObject*> features) {
		base_worker(pcActiveBody, cmd, which, features.front(), "", func);
	};

	//if there is a sketch selected which is from another body or part we need to bring up the
	//pick task dialog to decide how those are handled
	bool ext = std::find_if(status.begin(), status.end(),
		[](const PartDesignGui::TaskFeaturePick::featureStatus& s) {
		return s == PartDesignGui::TaskFeaturePick::otherBody ||
			s == PartDesignGui::TaskFeaturePick::otherPart ||
			s == PartDesignGui::TaskFeaturePick::notInBody;
	}) != status.end();

	// TODO Clean this up (2015-10-20, Fat-Zer)
	if (pcActiveBody && !bNoSketchWasSelected && ext)
	{
		auto* pcActivePart = PartDesignGui::getPartFor(pcActiveBody, true);
		// getPartFor() already has reported an error
		if (!pcActivePart)
			return;

		QDialog* dia = new QDialog;
		Ui_Dialog dlg;
		dlg.setupUi(dia);
		dia->setModal(true);
		int result = dia->exec();
		if (result == QDialog::DialogCode::Rejected)
			return;
		else if (!dlg.radioXRef->isChecked())
		{
			auto copy = PartDesignGui::TaskFeaturePick::makeCopy(sketches[0], "", dlg.radioIndependent->isChecked());
			auto oBody = PartDesignGui::getBodyFor(sketches[0], false);
			if (oBody)
				pcActiveBody->addObject(copy);
			else
				pcActivePart->addObject(copy);

			sketches[0] = copy;
			firstFreeSketch = sketches.begin();
		}
	}

	// Show sketch choose dialog and let user pick sketch if no sketch was selected and no free one available or
	// multiple free ones are available
	if (bNoSketchWasSelected && (freeSketches != 1))
	{
		Gui::TaskView::TaskDialog* dlg = Gui::Control().activeDialog();
		PartDesignGui::TaskDlgFeaturePick* pickDlg = qobject_cast<PartDesignGui::TaskDlgFeaturePick*>(dlg);
		if (dlg && !pickDlg)
		{
			QMessageBox msgBox;
			msgBox.setText(QObject::tr("A dialog is already open in the task panel"));
			msgBox.setInformativeText(QObject::tr("Do you want to close this dialog?"));
			msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
			msgBox.setDefaultButton(QMessageBox::Yes);
			int ret = msgBox.exec();
			if (ret == QMessageBox::Yes)
				Gui::Control().closeDialog();
			else
				return;
		}

		if (dlg)
			Gui::Control().closeDialog();

		Gui::Selection().clearSelection();
		pickDlg = new PartDesignGui::TaskDlgFeaturePick(sketches, status, accepter, sketch_worker);
		if (!bNoSketchWasSelected && ext)
			pickDlg->showExternal(true);

		Gui::Control().showDialog(pickDlg);
	}
	else
	{
		std::vector<App::DocumentObject*> theSketch;
		if (!bNoSketchWasSelected)
			theSketch.push_back(sketches[0]);
		else
			theSketch.push_back(*firstFreeSketch);

		sketch_worker(theSketch);
	}
}

void finishProfileBased(const Gui::Command* cmd, const Part::Feature* sketch, App::DocumentObject* Feat)
{
	if (sketch && sketch->isDerivedFrom(Part::Part2DObject::getClassTypeId()))
		FCMD_OBJ_HIDE(sketch);
	finishFeature(cmd, Feat);
}

//===========================================================================
// PartDesign_Pocket
//===========================================================================
DEF_STD_CMD_A(CmdPartDesignPocket);

CmdPartDesignPocket::CmdPartDesignPocket()
	: Command("PartDesign_Pocket")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Cut");
	sToolTipText = QT_TR_NOOP("Create a cutting sketch on a face");
	sWhatsThis = "PartDesign_Pocket";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_Pocket";
}

void CmdPartDesignPocket::activated(int iMsg)
{
	Q_UNUSED(iMsg);
	App::Document* doc = getDocument();
	std::vector<int> vecLinks;

	if (!PartDesignGui::assureModernWorkflow(doc))
		return;

	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(true);

	if (!pcActiveBody)
		return;

	Gui::SelectionFilter FaceFilter("SELECT Part::Feature SUBELEMENT Face COUNT 1");

	if (FaceFilter.match())
	{
		// get the selected object
		std::string supportString;
		App::DocumentObject* obj = FaceFilter.Result[0][0].getObject();

		if (!obj->isDerivedFrom(Part::Feature::getClassTypeId()))
			return;

		Part::Feature* feat = static_cast<Part::Feature*>(obj);

		const std::vector<std::string>& sub = FaceFilter.Result[0][0].getSubNames();
		if (sub.size() > 1)
		{
			// No assert for wrong user input!
			QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Several sub-elements selected"),
				QObject::tr("You have to select a single face as support for a sketch!"));
			return;
		}

		// get the selected sub shape (a Face)
		const Part::TopoShape& shape = feat->Shape.getValue();
		TopoDS_Shape sh = shape.getSubShape(sub[0].c_str());
		const TopoDS_Face& face = TopoDS::Face(sh);
		if (face.IsNull())
		{
			// No assert for wrong user input!
			QMessageBox::warning(Gui::getMainWindow(), QObject::tr("No support face selected"),
				QObject::tr("You have to select a face as support for a sketch!"));
			return;
		}

		BRepAdaptor_Surface adapt(face);
		if (adapt.GetType() != GeomAbs_Plane)
		{
			TopLoc_Location loc;
			Handle(Geom_Surface) surf = BRep_Tool::Surface(face, loc);
			if (surf.IsNull() || !GeomLib_IsPlanarSurface(surf).IsPlanar())
			{
				QMessageBox::warning(Gui::getMainWindow(), QObject::tr("No planar support"),
					QObject::tr("You need a planar face as support for a sketch!"));
				return;
			}
		}

		supportString = FaceFilter.Result[0][0].getAsPropertyLinkSubString();

		// Add all edges in the selected faces as external links
		for (TopExp_Explorer edgeExplorer(face, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next())
		{
			TopoDS_Edge edgeForward = TopoDS::Edge(edgeExplorer.Current());
			TopoDS_Edge edgeReversed = TopoDS::Edge(edgeExplorer.Current().Reversed());
			int nEdges = shape.countSubShapes("Edge");
			for (int i = 1; i <= nEdges; ++i)
			{
				TopoDS_Edge sh = TopoDS::Edge(shape.getSubShape(("Edge" + std::to_string(i)).c_str()));
				if (sh.IsEqual(edgeForward))
				{
					vecLinks.emplace_back(i);
					break;
				}
				else if (sh.IsEqual(edgeReversed))
				{
					vecLinks.emplace_back(i);
					break;
				}
			}
		}

		if (!pcActiveBody->hasObject(obj))
		{
			if (!obj->isDerivedFrom(App::Plane::getClassTypeId()))
			{
				// TODO check here if the plane associated with right part/body (2015-09-01, Fat-Zer)

				auto pcActivePart = PartDesignGui::getPartFor(pcActiveBody, false);

				//check the prerequisites for the selected objects
				//the user has to decide which option we should take if external references are used
				// TODO share this with UnifiedDatumCommand() (2015-10-20, Fat-Zer)
				QDialog* dia = new QDialog;
				Ui_Dialog dlg;
				dlg.setupUi(dia);
				dia->setModal(true);
				int result = dia->exec();
				if (result == QDialog::DialogCode::Rejected)
					return;
				else if (!dlg.radioXRef->isChecked())
				{
					std::string sub;
					if (FaceFilter.match())
						sub = FaceFilter.Result[0][0].getSubNames()[0];
					auto copy = PartDesignGui::TaskFeaturePick::makeCopy(obj, sub, dlg.radioIndependent->isChecked());

					if (pcActiveBody)
						pcActiveBody->addObject(copy);
					else if (pcActivePart)
						pcActivePart->addObject(copy);

					//it is ensured that only a single face is selected, hence it must always be Face1 of the shapebinder
					supportString = getObjectCmd(copy, "(", ",'Face1')");
				}
			}
		}

		// create Sketch on Face or Plane
		std::string FeatName = getUniqueObjectName("Sketch", pcActiveBody);

		openCommand("Create a Sketch on Face");
		FCMD_OBJ_CMD(pcActiveBody, "newObject('Sketcher::SketchObject','" << FeatName << "')");
		auto Feat = pcActiveBody->getDocument()->getObject(FeatName.c_str());

		Gui::Command* cmd = this;
		auto worker = [cmd](Part::Feature* sketch, App::DocumentObject* Feat) {
			if (!Feat)
				return;

			FCMD_OBJ_CMD(Feat, "Length = 5.0");
			finishProfileBased(cmd, sketch, Feat);
			cmd->adjustCameraPosition();
		};

		// prepareProfileBased(pcActiveBody, this, "Pocket", worker);
		base_worker(pcActiveBody, this, "Pocket", Feat, FeatName, worker);

		FCMD_OBJ_CMD(Feat, "Support = " << supportString);
		FCMD_OBJ_CMD(Feat, "MapMode = '" << Attacher::AttachEngine::getModeName(Attacher::mmFlatFace) << "'");

		for (auto& edge : vecLinks)
		{
			doCommand(Gui::Command::Doc, "App.ActiveDocument.%s.addExternal(\"%s\",\"%s\")",
				FeatName.c_str(),
				obj->getNameInDocument(), ("Edge" + std::to_string(edge)).c_str());
		}

		updateActive();
		PartDesignGui::setEdit(Feat, pcActiveBody);
	}
	else
	{
		QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Error"),
			QObject::tr("You should specify a face on which the cutting sketch lies!"));
	}
}

bool CmdPartDesignPocket::isActive(void)
{
	return hasActiveDocument();
}

//===========================================================================
// PartDesign_Hole
//===========================================================================
DEF_STD_CMD_A(CmdPartDesignHole);

CmdPartDesignHole::CmdPartDesignHole()
	: Command("PartDesign_Hole")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Hole");
	sToolTipText = QT_TR_NOOP("Create a hole with the selected sketch");
	sWhatsThis = "PartDesign_Hole";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_Hole";
}

void CmdPartDesignHole::activated(int iMsg)
{
	Q_UNUSED(iMsg);
	App::Document* doc = getDocument();
	std::vector<int> vecLinks;

	if (!PartDesignGui::assureModernWorkflow(doc))
		return;

	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(true);

	if (!pcActiveBody)
		return;

	Gui::SelectionFilter FaceFilter("SELECT Part::Feature SUBELEMENT Face COUNT 1");

	if (FaceFilter.match())
	{
		// get the selected object
		std::string supportString;
		App::DocumentObject* obj = FaceFilter.Result[0][0].getObject();

		if (!obj->isDerivedFrom(Part::Feature::getClassTypeId()))
			return;

		Part::Feature* feat = static_cast<Part::Feature*>(obj);

		const std::vector<std::string>& sub = FaceFilter.Result[0][0].getSubNames();
		if (sub.size() > 1)
		{
			// No assert for wrong user input!
			QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Several sub-elements selected"),
				QObject::tr("You have to select a single face as support for a sketch!"));
			return;
		}

		// get the selected sub shape (a Face)
		const Part::TopoShape& shape = feat->Shape.getValue();
		TopoDS_Shape sh = shape.getSubShape(sub[0].c_str());
		const TopoDS_Face& face = TopoDS::Face(sh);
		if (face.IsNull())
		{
			// No assert for wrong user input!
			QMessageBox::warning(Gui::getMainWindow(), QObject::tr("No support face selected"),
				QObject::tr("You have to select a face as support for a sketch!"));
			return;
		}

		BRepAdaptor_Surface adapt(face);
		if (adapt.GetType() != GeomAbs_Plane)
		{
			TopLoc_Location loc;
			Handle(Geom_Surface) surf = BRep_Tool::Surface(face, loc);
			if (surf.IsNull() || !GeomLib_IsPlanarSurface(surf).IsPlanar())
			{
				QMessageBox::warning(Gui::getMainWindow(), QObject::tr("No planar support"),
					QObject::tr("You need a planar face as support for a sketch!"));
				return;
			}
		}

		supportString = FaceFilter.Result[0][0].getAsPropertyLinkSubString();

		// Add all edges in the selected faces as external links
		for (TopExp_Explorer edgeExplorer(face, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next())
		{
			TopoDS_Edge edgeForward = TopoDS::Edge(edgeExplorer.Current());
			TopoDS_Edge edgeReversed = TopoDS::Edge(edgeExplorer.Current().Reversed());
			int nEdges = shape.countSubShapes("Edge");
			for (int i = 1; i <= nEdges; ++i)
			{
				TopoDS_Edge sh = TopoDS::Edge(shape.getSubShape(("Edge" + std::to_string(i)).c_str()));
				if (sh.IsEqual(edgeForward))
				{
					vecLinks.emplace_back(i);
					break;
				}
				else if (sh.IsEqual(edgeReversed))
				{
					vecLinks.emplace_back(i);
					break;
				}
			}
		}

		if (!pcActiveBody->hasObject(obj))
		{
			if (!obj->isDerivedFrom(App::Plane::getClassTypeId()))
			{
				// TODO check here if the plane associated with right part/body (2015-09-01, Fat-Zer)

				auto pcActivePart = PartDesignGui::getPartFor(pcActiveBody, false);

				//check the prerequisites for the selected objects
				//the user has to decide which option we should take if external references are used
				// TODO share this with UnifiedDatumCommand() (2015-10-20, Fat-Zer)
				QDialog* dia = new QDialog;
				Ui_Dialog dlg;
				dlg.setupUi(dia);
				dia->setModal(true);
				int result = dia->exec();
				if (result == QDialog::DialogCode::Rejected)
					return;
				else if (!dlg.radioXRef->isChecked())
				{
					std::string sub;
					if (FaceFilter.match())
						sub = FaceFilter.Result[0][0].getSubNames()[0];
					auto copy = PartDesignGui::TaskFeaturePick::makeCopy(obj, sub, dlg.radioIndependent->isChecked());

					if (pcActiveBody)
						pcActiveBody->addObject(copy);
					else if (pcActivePart)
						pcActivePart->addObject(copy);

					//it is ensured that only a single face is selected, hence it must always be Face1 of the shapebinder
					supportString = getObjectCmd(copy, "(", ",'Face1')");
				}
			}
		}

		// create Sketch on Face or Plane
		std::string FeatName = getUniqueObjectName("Sketch", pcActiveBody);

		openCommand("Create a Sketch on Face");
		FCMD_OBJ_CMD(pcActiveBody, "newObject('Sketcher::SketchObject','" << FeatName << "')");
		auto Feat = pcActiveBody->getDocument()->getObject(FeatName.c_str());

		Gui::Command* cmd = this;
		auto worker = [cmd](Part::Feature* sketch, App::DocumentObject* Feat) {
			if (!Feat)
				return;

			finishProfileBased(cmd, sketch, Feat);
			cmd->adjustCameraPosition();
		};

		// prepareProfileBased(pcActiveBody, this, "Pocket", worker);
		base_worker(pcActiveBody, this, "Hole", Feat, FeatName, worker);

		FCMD_OBJ_CMD(Feat, "Support = " << supportString);
		FCMD_OBJ_CMD(Feat, "MapMode = '" << Attacher::AttachEngine::getModeName(Attacher::mmFlatFace) << "'");

		for (auto& edge : vecLinks)
		{
			doCommand(Gui::Command::Doc, "App.ActiveDocument.%s.addExternal(\"%s\",\"%s\")",
				FeatName.c_str(),
				obj->getNameInDocument(), ("Edge" + std::to_string(edge)).c_str());
		}

		updateActive();
		PartDesignGui::setEdit(Feat, pcActiveBody);
	}
	else
	{
		QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Error"),
			QObject::tr("A cutting should be started on an existing face!"));
	}
}

bool CmdPartDesignHole::isActive(void)
{
	return hasActiveDocument();
}

/* Boolean commands =======================================================*/
DEF_STD_CMD_A(CmdPartDesignBoolean);

CmdPartDesignBoolean::CmdPartDesignBoolean()
	: Command("PartDesign_Boolean")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Boolean operation");
	sToolTipText = QT_TR_NOOP("Compile to HL-HELM and LL-HELM");
	sWhatsThis = "PartDesign_Boolean";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_Boolean";
}

void CmdPartDesignBoolean::activated(int iMsg)
{
	Q_UNUSED(iMsg);
	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ true);
	if (!pcActiveBody)
		return;



	App::Document* doc = getDocument();
	std::vector<App::DocumentObject*> sortedObjects = doc->topologicalSort();

	Gui::WorkbenchManager* pcWbMan = Gui::WorkbenchManager::instance();
	auto cdWb = dynamic_cast<PartDesignGui::Workbench*>(pcWbMan->active());

	//cdWb->langCompiler->CreateVisualizer();

	Compiler cpl_;
	cpl_.CompileDown(sortedObjects);

	if (cdWb != nullptr)
	{
		cdWb->hlEditor->setPlainText(cpl_.GetHLHELM());
		cdWb->llEditor->setPlainText(cpl_.GetLLHELM());
	}
	return;
}

bool CmdPartDesignBoolean::isActive(void)
{
	if (getActiveGuiDocument())
		return true;
	else
		return false;
}

//===========================================================================
// Carpentry_Packing_Visualization for 2d
//===========================================================================

TopoDS_Shape MakePart(const Standard_Real& part_length, const Standard_Integer& end_angle_0, const Standard_Integer& end_angle_1)
{
	auto MakeWires = [](const std::vector<Vector3d>& vecs) {
		TopoDS_Wire wire;
		for (int iter = 0; iter < vecs.size(); iter++) {
			Vector3d start = vecs[iter];
			Vector3d end = vecs[(iter + 1) % vecs.size()];
			TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(gp_Pnt(start[0], start[1], start[2]), gp_Pnt(end[0], end[1], end[2]));
			wire = BRepBuilderAPI_MakeWire(wire, edge);
		}
		return wire;
	};

	const Standard_Real part_height = 88.9;

	std::vector<Vector3d> contours;
	contours.emplace_back(Vector3d(0.0, 0.0, 0.0));
	contours.emplace_back(Vector3d(part_length, 0.0, 0.0));
	contours.emplace_back(Vector3d(part_length + end_angle_1 * 38.1, 38.1, 0.0));
	contours.emplace_back(Vector3d(0.0 + end_angle_0 * 38.1, 38.1, 0.0));

	if (part_length + end_angle_1 * 38.1 < end_angle_0 * 38.1)
	{
		std::cerr << "if (part_length + end_angle_1 * 38.1 < end_angle_0 * 38.1)" << std::endl;
		system("pause");
	}

	TopoDS_Wire aWire = MakeWires(contours);
	TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(aWire);
	gp_Vec aPrismVec(0, 0, part_height);
	TopoDS_Shape myBody = BRepPrimAPI_MakePrism(myFaceProfile, aPrismVec);

	return myBody;
}

/* Boolean commands =======================================================*/
DEF_STD_CMD_A(CarpentryPacking);

CarpentryPacking::CarpentryPacking()
	: Command("Carpentry_Packing")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("Carpentry_Packing operation");
	sToolTipText = QT_TR_NOOP("Carpentry Packing");
	sWhatsThis = "Carpentry_Packing";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_Mirrored";
}

void RunParserScript(const char* folder)
{
#ifndef MAINOBJ
	Base::PyGILStateLocker lock;

	string chdir_cmd = string("sys.path.append('../../Mod/Partdesign/Scripts')");
	const char* cstr_cmd = chdir_cmd.c_str();

	int nRet = PyRun_SimpleString("import sys");
	nRet = PyRun_SimpleString(cstr_cmd);

	PyObject* pModule = PyImport_ImportModule("Parser");
	if (!pModule)
	{
		QMessageBox::warning(NULL, QObject::tr("Warning"), QObject::tr("Cannot load python scripts: Parser.py"), QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	auto pFunc = PyObject_GetAttrString(pModule, "parse");

	PyObject* temp = Py_BuildValue("(s)", folder);
	PyObject_CallObject(pFunc, temp);

#else

	std::string cmd = "python ../../../../../src/Mod/PartDesign/Scripts/Parser.py " + std::string(folder);
	system(cmd.c_str());

#endif

}

void CarpentryHELM(PartDesignGui::PCE *pce, PHELM *langCompiler, PartDesignGui::GLMATS& glmats)
{
	const auto start = clock();
	
	#ifndef MAINOBJ
	auto progressDlg = new QProgressDialog(Gui::getMainWindow());
	//Set progress dialog in undeterminate state: 
	progressDlg->setMaximum(0);
	progressDlg->setMinimum(0);
	progressDlg->setMinimumWidth(400);
	progressDlg->setMinimumHeight(200);
	progressDlg->setWindowModality(Qt::WindowModal);
	progressDlg->setCancelButton(nullptr);
	progressDlg->setWindowFlag(Qt::FramelessWindowHint);
	progressDlg->setLabelText(QString::fromUtf8("Optimizing... Please wait patiently!"));
	progressDlg->show();
	#endif

	std::atomic<bool> worker_canceled(false), worker_stopped(false);
	std::exception_ptr worker_exception;
	std::thread worker_thread([&]() {
		try {
			const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

			std::ofstream log_file(outputFolder + "\\log.txt");

			if (CompilerConfig::Instance().GetRunPack())
			{
				Math::Functs::ClearFolder(outputFolder);

				// Enumerating
				PartDesignGui::PCE::PartsRegistration(glmats, CompilerConfig::Instance().GetVecShape());
				std::vector<int> parts_encodes;
				PartDesignGui::PCE::BuildShapes(CompilerConfig::Instance().GetVecWLC(),0, glmats,
					pce->GetPackingFolder(), CompilerConfig::Instance().GetVecShape(), parts_encodes);
				pce->EnumerateOrders(CompilerConfig::Instance().GetVecWLC(),parts_encodes);

				//cutting order
				if (CompilerConfig::Instance().GetRunOrder())
				{
					// Compilation
					std::unordered_set<int> compiledNodes;
					std::unordered_set<int> nonCompiledNodes;
					std::unordered_set<int> nonCompiledArranges;
					auto resCompiledEgraph = langCompiler->CompileEgraphEnhanced(pce,-1, compiledNodes, nonCompiledNodes, nonCompiledArranges, outputFolder);

					//pce->SetPHELM(*langCompiler);

					if (!nonCompiledNodes.empty())
					{
						std::cerr << "These e-classes have been failed to compile: ";
						for (auto& node : nonCompiledNodes)std::cerr << node << " "; std::cerr << "\n";
						std::cerr << "These arranges have been failed to compile: ";
						for (auto& arrange : nonCompiledArranges)std::cerr << arrange << " "; std::cerr << "\n";
						//system("pause");

						std::ofstream file(CompilerConfig::Instance().GetPackingOutputFolder()+"\\nonCompiledNodes.txt",ios::app);
						file << "E_Class: ";
						for (auto& node : nonCompiledNodes)file << node << " "; std::cerr << "\n";
						file << "Arrange: ";
						for (auto& arrange : nonCompiledArranges)file << arrange << " "; std::cerr << "\n";
						file.close();
					}

					langCompiler->PostprocessEGraph(pce,-1, compiledNodes, outputFolder);

					//std::vector<int> nodes;
					//cdWb->langCompiler->TestPost(cdWb->pce, cdWb->pce->GetEClass(-1),nodes);
					//TestPost(PartDesignGui::PCE* pce, PCEEClass& ec,std::vector<int>& nodes)
					//std::cerr << nodes.size() << std::endl;
					langCompiler->ExportResultEGraphEnhanced(pce, -1, compiledNodes, outputFolder);

					// Python parser
					RunParserScript(outputFolder.c_str());

					std::cerr << "MapShapeToUuid Size: " << langCompiler->GetMapShapeToUuid().Size() << std::endl;
				}
			}

			if (CompilerConfig::Instance().GetRunOpt())
			{
				// Optimization
				auto optimizer = std::make_shared<OptCompiler::OptimizingCompiler>();
				optimizer->execute(outputFolder.c_str());

				//visualization_read
				OptCompiler::FabInfo fab_info= OptCompiler::OptimizingCompiler::ReadParetoFronts(outputFolder + "\\collapse\\result.xml");

				for (int i = 0; i < fab_info.fabs.size(); i++)
				{
					std::cerr << "Visualization: " << i << std::endl;
					auto& fab = fab_info.fabs[i];
					auto score = fab.score;

					std::vector<std::pair<Vector3d, Vector3d>> refer_edges;
					Vector3d2 refer_faces;
					Vector3d1 refer_values;
					for (int j = 0; j < fab.referFE.size(); j++)
					{
						//edge
						auto edge = langCompiler->GetShapeDS(boost::lexical_cast<boost::uuids::uuid>(std::get<1>(fab.referFE[j])));
						const TopoDS_Edge& curEdge = TopoDS::Edge(edge);
						gp_Pnt st, nd;
						GeomFunc::GetEdgeEndPnt(curEdge, st, nd);
						refer_edges.emplace_back(Vector3d(st.X(), st.Y(), st.Z()), Vector3d(nd.X(), nd.Y(), nd.Z()));

						//face
						auto face = langCompiler->GetShapeDS(boost::lexical_cast<boost::uuids::uuid>(std::get<0>(fab.referFE[j])));
						auto contour = GeomFunc::LoadTopoDSGeometry(TopoDS::Face(face));
						refer_faces.emplace_back(contour[0]);

						refer_values.emplace_back(std::get<2>(fab.referFE[j]), std::get<3>(fab.referFE[j]), std::get<4>(fab.referFE[j]));
					}

					std::string str = std::to_string(score[0]) + "_" + std::to_string(score[1]) + "_" + std::to_string(score[2]);

					pce->Visualization(
						score,
						CompilerConfig::Instance().GetVecWLC(),
						fab.order,
						fab.order_stack,
						fab.order_merge,
						fab.ArrIDs,
						fab.ENodeIDs,
						refer_edges,
						refer_faces,
						refer_values,
						outputFolder+"\\vis_"+std::to_string(i)+"_"+ str +".obj");
				}
			}

			const auto end = clock();
			const auto endtime = static_cast<double>(end - start) / CLOCKS_PER_SEC;
			std::cerr << "Total time:" << Math::Functs::DoubleString(endtime) << "s" << std::endl;

			log_file << "Total time: " << Math::Functs::DoubleString(endtime) << "s" << std::endl;
		}
		catch (...) {
			worker_exception = std::current_exception();
		}
		worker_stopped = true;
		});

	try {
		while (!worker_stopped) {
			QApplication::processEvents(QEventLoop::WaitForMoreEvents);

#ifndef MAINOBJ
			if (progressDlg->wasCanceled()) {
				worker_canceled = true;
				break;
			}
#endif
		}
	}
	catch (...) {
		worker_canceled = true;
		worker_thread.join();
		throw;
	}
	worker_thread.join();
	if (worker_exception) {
		std::rethrow_exception(worker_exception);
	}
#ifndef MAINOBJ
	progressDlg->cancel();
#endif
}

void CarpentryPacking::activated(int iMsg)
{
	Q_UNUSED(iMsg);

	//PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ true);
	//if (!pcActiveBody) return;

	//get doc and workbench
	auto doc = getDocument();
	auto pcWbMan = Gui::WorkbenchManager::instance();
	auto cdWb = dynamic_cast<PartDesignGui::Workbench*>(pcWbMan->active());

	if (!doc)
	{
		std::cerr << "if (!cdWb)" << std::endl;
		system("pause");
		return;
	}

	PartDesignGui::GLMATS glmats;

	//get object
	for (auto& feat : doc->getObjects())
	{
		if (feat->isDerivedFrom(PartDesign::Body::getClassTypeId()))
		{
			const auto& features = dynamic_cast<PartDesign::Body*>(feat)->Group.getValue();
			auto g_m = dynamic_cast<PartDesign::Body*>(feat)->globalGroupPlacement().toMatrix();
			//glm::dmat4 M = dynamic_cast<glm::dmat4>&g_m;

			for (auto rot = features.crbegin(); rot != features.crend(); ++rot)
			{
				if ((*rot)->isDerivedFrom(Part::Feature::getClassTypeId()) && PartDesign::Body::isSolidFeature(*rot))
				{
					const TopoDS_Shape& local_ashape = Part::Feature::getShape(*rot, 0, false, 0, 0, true, false);
					const TopoDS_Shape& global_ashape = Part::Feature::getShape(*rot, 0, false, &g_m, 0, true, true);

					glm::dmat4 M;
					for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) M[i][j] = g_m[j][i];

					//random
					//srand((unsigned)time(NULL));
					auto random_M =  Math::Functs::TranslationMatrix(Vector3d(50.0 * rand() / RAND_MAX, 50.0 * rand() / RAND_MAX, 50.0 * rand() / RAND_MAX));
					const auto trsf = GeomFunc::ConvertTogpTrsf(random_M);
					TopoDS_Shape curShape = local_ashape;
					TopoDS_Shape local_ashape_= BRepBuilderAPI_GTransform(curShape, trsf).Shape();
					auto mat = M * glm::inverse(random_M);//random

					glmats.emplace_back(PartDesignGui::GLMAT(global_ashape, local_ashape_, mat));
					break;
				}
			}
		}
	}

	if (glmats.size() == 0)
	{
		std::cerr << "if (part_nb == 0)" << std::endl;
		system("pause");
		return;
	}

	CarpentryHELM(cdWb->pce,cdWb->langCompiler, glmats);
	
}

bool CarpentryPacking::isActive(void)
{
	return true;
	//return getActiveGuiDocument() != nullptr;
}

//===========================================================================
// PartsCombination_Emulating for 2d
//===========================================================================

/* Boolean commands =======================================================*/
DEF_STD_CMD_A(CmdParseHLHelmProgram);

CmdParseHLHelmProgram::CmdParseHLHelmProgram()
	: Command("CmdParseHLHelmProgram")
{
	sAppModule = "PartDesign";
	sGroup = QT_TR_NOOP("PartDesign");
	sMenuText = QT_TR_NOOP("CmdParseHLHelmProgram");
	sToolTipText = QT_TR_NOOP("CmdParseHLHelmProgram");
	sWhatsThis = "CmdParseHLHelmProgram";
	sStatusTip = sToolTipText;
	sPixmap = "PartDesign_MultiTransform";
}

std::string QueryFaceByClosestPoint(PartDesign::Feature* prm, const gp_Pnt& testPoint)
{
	using namespace HighLevelParser;
	std::string selectedFace = "";
	auto shape = Part::Feature::getTopoShape(prm, 0, false, 0, 0, true, false);
	int nFaces = shape.countSubShapes("Face");

	auto minDistPntToFace = DBL_MAX;
	for (int i = 1; i <= nFaces; ++i)
	{
		std::string tmpFaceName = std::string("Face") + std::to_string(i);
		TopoDS_Face sh = TopoDS::Face(shape.getSubShape(tmpFaceName.c_str()));

		// only select planar faces
		Handle(Geom_Surface) aSurface = BRep_Tool::Surface(sh);
		BRepAdaptor_Surface adapt(sh);
		if (adapt.GetType() != GeomAbs_Plane)
		{
			TopLoc_Location loc;
			if (aSurface.IsNull() || !GeomLib_IsPlanarSurface(aSurface).IsPlanar())
				continue;
		}

		auto pps = GeomAPI_ProjectPointOnSurf(testPoint, aSurface);

		const auto tmpDist = pps.LowerDistance();
		//std::cerr << tmpDist << " " << tmpFaceName << std::endl;
		if (tmpDist < minDistPntToFace)
		{
			minDistPntToFace = tmpDist;
			selectedFace = tmpFaceName;
		}
	}
	std::cerr << "Selecte: " << selectedFace << std::endl;
	return selectedFace;
}

std::string QueryEdgeByClosestPoint(PartDesign::Feature* prm, const gp_Pnt& testPoint)
{
	using namespace HighLevelParser;
	std::string selectedEdge = "";
	auto shape = Part::Feature::getTopoShape(prm, 0, false, 0, 0, true, false);
	int nEdges = shape.countSubShapes("Edge");

	auto minDistPntToFace = DBL_MAX;
	for (int i = 1; i <= nEdges; ++i)
	{
		std::string tmpFaceName = std::string("Edge") + std::to_string(i);
		TopoDS_Edge sh = TopoDS::Edge(shape.getSubShape(tmpFaceName.c_str()));

		Standard_Real first, last;
		Handle(Geom_Curve) aCurve = BRep_Tool::Curve(sh, first, last);

		auto pps = GeomAPI_ProjectPointOnCurve(testPoint, aCurve);

		const auto tmpDist = pps.LowerDistance();
		//std::cerr << tmpDist << " " << tmpFaceName << std::endl;
		if (tmpDist < minDistPntToFace)
		{
			minDistPntToFace = tmpDist;
			selectedEdge = tmpFaceName;
		}
	}
	std::cerr << "Selecte: " << selectedEdge << std::endl;
	return selectedEdge;
}

std::string QueryArcByClosestCenterAndRadius(PartDesign::Feature* prm, const gp_Pnt& testPoint, const double& rad)
{
	using namespace HighLevelParser;
	std::string selectedEdge = "";
	auto shape = Part::Feature::getTopoShape(prm, 0, false, 0, 0, true, false);
	int nEdges = shape.countSubShapes("Edge");

	auto minDistPntToFace = DBL_MAX;
	for (int i = 1; i <= nEdges; ++i)
	{
		std::string tmpFaceName = std::string("Edge") + std::to_string(i);
		TopoDS_Edge sh = TopoDS::Edge(shape.getSubShape(tmpFaceName.c_str()));
		BRepAdaptor_Curve curve(sh);
		if (curve.GetType() == GeomAbs_Circle)
		{
			gp_Circ circle = curve.Circle();
			gp_Pnt cnt = circle.Location();

			const gp_XYZ diff = cnt.XYZ() - testPoint.XYZ();
			const auto tmpDist = diff.Dot(diff) + std::abs(rad - circle.Radius());
			//std::cerr << tmpDist << " " << tmpFaceName << std::endl;
			if (tmpDist < minDistPntToFace)
			{
				minDistPntToFace = tmpDist;
				selectedEdge = tmpFaceName;
			}
		}
	}
	std::cerr << "Selected: " << selectedEdge << std::endl;
	return selectedEdge;
}

std::vector<HighLevelParser::GVertex> GetOrderedChildren(const HighLevelParser::Graph& g, const HighLevelParser::GVertex tv)
{
	using namespace HighLevelParser;
	std::vector<GVertex> vecChildren;
	boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
	tie(vi, vi_end) = boost::out_edges(tv, g);
	for (next = vi; vi != vi_end; vi = next)
	{
		++next;
		auto vtv = target(*vi, g);
		vecChildren.push_back(vtv);
	}
	return vecChildren;
}

int GetGeometryIdOfExternalGeom(const HighLevelParser::Graph& g, const HighLevelParser::GVertex& v,
	const std::unordered_map<std::string, void*>& mapNameToObj,
	Sketcher::SketchObject* sketchObj)
{
	using namespace HighLevelParser;
	std::vector<double> testPointVals;
	std::string queryBody;
	for (auto vv : GetOrderedChildren(g, v))
	{
		if (g[vv]->GetType() == TYPE_Float)
		{
			auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(g[vv]);
			testPointVals.push_back(tmpFloat->value);
		}
		else
		{
			queryBody = g[vv]->str;
		}
	}

	if (g[v]->GetType() == TYPE_Query_Edge_By_Closest_Point)
	{
		std::cerr << "Query_Edge_By_Closest_Point" << std::endl;
		auto* prm = static_cast<PartDesign::Feature*>(mapNameToObj.find(queryBody)->second);
		gp_Pnt testPoint(testPointVals[0], testPointVals[1], testPointVals[2]);
		auto selectedEdge = QueryEdgeByClosestPoint(prm, testPoint);
		auto& externalGeometry = sketchObj->ExternalGeometry.getValues();
		auto& externalSub = sketchObj->ExternalGeometry.getSubValues();
		auto& externalGeo = sketchObj->ExternalGeo.getValues();

		if (externalGeo.size() - 2 != externalGeometry.size())
		{
			std::cerr << "insame length " << std::endl;
			system("pause");
		}

		for (auto i = 0; i < externalGeometry.size(); ++i)
		{
			if (externalGeometry[i] == prm && externalSub[i] == selectedEdge)
			{
				return (-3 - i);
			}
		}

		int eId = sketchObj->addExternal(prm, selectedEdge.c_str());

		return (-3 - eId);
	}
	else if (g[v]->GetType() == TYPE_Query_Arc_By_Closest_Center_And_Radius)
	{
		std::cerr << "Query_Arc_By_Closest_Center_And_Radius" << std::endl;
		auto* prm = static_cast<PartDesign::Feature*>(mapNameToObj.find(queryBody)->second);
		gp_Pnt testPoint(testPointVals[0], testPointVals[1], testPointVals[2]);
		auto selectedEdge = QueryArcByClosestCenterAndRadius(prm, testPoint, testPointVals[3]);

		auto& externalGeometry = sketchObj->ExternalGeometry.getValues();
		auto& externalSub = sketchObj->ExternalGeometry.getSubValues();
		auto& externalGeo = sketchObj->ExternalGeo.getValues();

		if (externalGeo.size() - 2 != externalGeometry.size())
		{
			std::cerr << "insame length " << std::endl;
			system("pause");
		}

		for (auto i = 0; i < externalGeometry.size(); ++i)
		{
			if (externalGeometry[i] == prm && externalSub[i] == selectedEdge)
			{
				return (-3 - i);
			}
		}

		int eId = sketchObj->addExternal(prm, selectedEdge.c_str());

		return (-3 - eId);
	}
	return 1000;
}

void SetConstraint(const HighLevelParser::Graph& g, const HighLevelParser::GVertex& v,
	const std::unordered_map<std::string, void*>& mapNameToObj,
	int order, Sketcher::Constraint* cstr,
	Sketcher::SketchObject* sketchObj)
{
	using namespace HighLevelParser;
	if (g[v]->GetType() == TYPE_Float)
	{
		const auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(g[v]);
		cstr->setValue(tmpFloat->value);
		return;
	}

	auto* setGeom = &(cstr->First);
	auto* setGeomPos = &(cstr->FirstPos);
	if (order == 1)
	{
		setGeom = &(cstr->Second);
		setGeomPos = &(cstr->SecondPos);
	}
	else if (order == 2)
	{
		setGeom = &(cstr->Third);
		setGeomPos = &(cstr->ThirdPos);
	}

	if (g[v]->GetType() == TYPE_Start)
	{
		std::cerr << "Start" << std::endl;
		auto firstChild = GetOrderedChildren(g, v)[0];
		*setGeomPos = Sketcher::start;

		if (g[firstChild]->GetType() == TYPE_Query_Edge_By_Closest_Point ||
			g[firstChild]->GetType() == TYPE_Query_Arc_By_Closest_Center_And_Radius)
		{
			*setGeom = GetGeometryIdOfExternalGeom(g, firstChild, mapNameToObj, sketchObj);
		}
		else
		{
			if (g[firstChild]->GetType() == TYPE_Horizontal)*setGeom = -1;
			else if (g[firstChild]->GetType() == TYPE_Vertical)*setGeom = -2;
			else *setGeom = static_cast<Part::Geometry*>(mapNameToObj.find(g[firstChild]->str)->second)->Id;
		}
	}
	else if (g[v]->GetType() == TYPE_End)
	{
		std::cerr << "End" << std::endl;
		const auto firstChild = GetOrderedChildren(g, v)[0];
		*setGeomPos = Sketcher::end;

		if (g[firstChild]->GetType() == TYPE_Query_Edge_By_Closest_Point ||
			g[firstChild]->GetType() == TYPE_Query_Arc_By_Closest_Center_And_Radius)
		{
			*setGeom = GetGeometryIdOfExternalGeom(g, firstChild, mapNameToObj, sketchObj);
		}
		else
		{
			if (g[firstChild]->GetType() == TYPE_Horizontal)*setGeom = -1;
			else if (g[firstChild]->GetType() == TYPE_Vertical)*setGeom = -2;
			else *setGeom = static_cast<Part::Geometry*>(mapNameToObj.find(g[firstChild]->str)->second)->Id;
		}
	}
	else if (g[v]->GetType() == TYPE_Query_Edge_By_Closest_Point ||
		g[v]->GetType() == TYPE_Query_Arc_By_Closest_Center_And_Radius)	// External relationship, none type
	{
		*setGeom = GetGeometryIdOfExternalGeom(g, v, mapNameToObj, sketchObj);
		*setGeomPos = Sketcher::none;
	}
	else // None type
	{
		if (g[v]->GetType() == TYPE_Horizontal)*setGeom = -1;
		else if (g[v]->GetType() == TYPE_Vertical)*setGeom = -2;
		else *setGeom = static_cast<Part::Geometry*>(mapNameToObj.find(g[v]->str)->second)->Id;
		*setGeomPos = Sketcher::none;
	}
}

void CmdParseHLHelmProgram::activated(int iMsg)
{
	Q_UNUSED(iMsg);
	using namespace HighLevelParser;
	HLParser hlp;

	std::ifstream fin;
	fin.open("D:\\test.hlhelm");
	std::string line;
	while (!fin.eof())
	{
		getline(fin, line);
		if (!line.empty())
			hlp.Parse(line);
	}
	fin.close();

	auto& vecAst = hlp.GetAst();

	if (getActiveGuiDocument() == nullptr)
	{
		App::GetApplication().newDocument();
	}

	auto doc = getDocument();

	// TODO: Two choices, one is only create a body, that's enough for compilation; the other is to create individually bodies
	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/* messageIfNot = */ false);
	auto shouldMakeBody(false);
	if (pcActiveBody == nullptr)
	{
		if (doc->getObjectsOfType(PartDesign::Body::getClassTypeId()).empty())
		{
			shouldMakeBody = true;
		}
		else
		{
			PartDesignGui::needActiveBodyError();
		}
	}

	std::unordered_map<std::string, PartDesign::Body*> mapNameToBody;
	std::unordered_map<std::string, void*> mapNameToObj;

	for (auto& ast : vecAst)
	{
		std::cerr << ast.g[ast.assign]->str << std::endl;
		std::string assignValue = ast.g[ast.assign]->str;

		if (ast.g[ast.operation]->GetType() == TYPE_Make_Sketch)
		{
			std::string selectedFace;
			Sketcher::SketchObject* sketchObj = nullptr;
			boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
			tie(vi, vi_end) = boost::out_edges(ast.operation, ast.g);
			for (next = vi; vi != vi_end; vi = next)
			{
				++next;
				auto tv = target(*vi, ast.g);

				if (ast.g[tv]->GetType() == TYPE_Query_Face_By_Closest_Point)
				{
					std::cerr << "TYPE_Query_Face_By_Closest_Point" << std::endl;
					boost::graph_traits<Graph>::out_edge_iterator vvi, vvi_end, vnext;
					tie(vvi, vvi_end) = boost::out_edges(tv, ast.g);

					std::vector<double> testPointVals;
					std::string queryBody;

					for (vnext = vvi; vvi != vvi_end; vvi = vnext)
					{
						++vnext;
						auto vtv = target(*vvi, ast.g);

						if (ast.g[vtv]->GetType() == TYPE_Float)
						{
							auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(ast.g[vtv]);
							testPointVals.push_back(tmpFloat->value);
						}
						else
						{
							std::cerr << ast.g[vtv]->str << std::endl;
							queryBody = ast.g[vtv]->str;
						}
					}

					pcActiveBody = mapNameToBody[queryBody];
					std::cerr << "query body = " << queryBody << std::endl;
					auto* prm = dynamic_cast<PartDesign::Feature*>(pcActiveBody->getDocument()->getObject(queryBody.c_str()));

					gp_Pnt testPoint(testPointVals[0], testPointVals[1], testPointVals[2]);
					selectedFace = QueryFaceByClosestPoint(prm, testPoint);

					auto featName = assignValue;
					//openCommand("Create a Sketch on Face");
					FCMD_OBJ_CMD(pcActiveBody, "newObject('Sketcher::SketchObject','" << featName << "')");
					auto Feat = pcActiveBody->getDocument()->getObject(featName.c_str());

					auto supportObject = std::string("(App.getDocument('") + std::string(doc->getName()) + "').getObject('" + queryBody + "'), ['" + selectedFace + "', ])";
					FCMD_OBJ_CMD(Feat, "Support = " << supportObject);
					FCMD_OBJ_CMD(Feat, "MapMode = '" << Attacher::AttachEngine::getModeName(Attacher::mmFlatFace) << "'");

					sketchObj = dynamic_cast<Sketcher::SketchObject*>(Feat);
					mapNameToObj[assignValue] = sketchObj;
				}
				else if (ast.g[tv]->GetType() == TYPE_Geometry)
				{
					std::cerr << "TYPE_Geometry" << std::endl;
					boost::graph_traits<Graph>::out_edge_iterator vvi, vvi_end, vnext;
					tie(vvi, vvi_end) = boost::out_edges(tv, ast.g);

					std::vector<std::string> vecGeoms;
					for (vnext = vvi; vvi != vvi_end; vvi = vnext)
					{
						++vnext;
						auto vtv = target(*vvi, ast.g);
						std::cerr << ast.g[vtv]->str << std::endl;
						vecGeoms.push_back(ast.g[vtv]->str);
					}

					for (auto& geomName : vecGeoms)
					{
						auto geom = static_cast<Part::Geometry*>(mapNameToObj[geomName]);
						auto gId = sketchObj->addGeometry(geom);
						geom->Id = gId;
					}

					Gui::Command::updateActive();
				}
				else if (ast.g[tv]->GetType() == TYPE_Constraint)
				{
					std::cerr << "TYPE_Constraint" << std::endl;

					// Get external
					boost::graph_traits<Graph>::out_edge_iterator vvi, vvi_end, vnext;
					tie(vvi, vvi_end) = boost::out_edges(tv, ast.g);
					for (vnext = vvi; vvi != vvi_end; vvi = vnext)
					{
						++vnext;
						auto vtv = target(*vvi, ast.g);

						auto cstr = new Sketcher::Constraint();

						int order = 0;
						for (auto v : GetOrderedChildren(ast.g, vtv))
						{
							SetConstraint(ast.g, v, mapNameToObj, order++, cstr, sketchObj);
						}

						if (ast.g[vtv]->GetType() == TYPE_PointOnObject)
						{
							std::cerr << "TYPE_PointOnObject" << std::endl;

							cstr->Type = Sketcher::PointOnObject;
						}
						else if (ast.g[vtv]->GetType() == TYPE_Angle)
						{
							std::cerr << "TYPE_Angle" << std::endl;
							cstr->setValue(cstr->getValue() / 180.0 * M_PI);
							cstr->Type = Sketcher::Angle;
						}
						else if (ast.g[vtv]->GetType() == TYPE_Coincident)
						{
							std::cerr << "TYPE_Coincident" << std::endl;
							cstr->Type = Sketcher::Coincident;
						}
						else if (ast.g[vtv]->GetType() == TYPE_Distance)
						{
							std::cerr << "TYPE_Distance" << std::endl;
							cstr->Type = Sketcher::Distance;
						}
						else if (ast.g[vtv]->GetType() == TYPE_DistanceX)
						{
							std::cerr << "TYPE_DistanceX" << std::endl;
							cstr->Type = Sketcher::DistanceX;
						}
						else if (ast.g[vtv]->GetType() == TYPE_DistanceY)
						{
							std::cerr << "TYPE_DistanceY" << std::endl;
							cstr->Type = Sketcher::DistanceY;
						}
						else if (ast.g[vtv]->GetType() == TYPE_Equal)
						{
							std::cerr << "TYPE_Equal" << std::endl;
							cstr->Type = Sketcher::Equal;
						}
						else if (ast.g[vtv]->GetType() == TYPE_Parallel)
						{
							std::cerr << "TYPE_Parallel" << std::endl;
							cstr->Type = Sketcher::Parallel;
						}

						sketchObj->addConstraint(cstr);
					}
				}
			}
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Make_Cut)
		{
			std::cerr << "TYPE_Make_Cut" << std::endl;
			// Add cut operation
			Gui::Command* cmd = this;
			auto worker = [](Part::Feature* sketch, App::DocumentObject* feat) {
				if (!feat)
					return;
			};

			auto children = GetOrderedChildren(ast.g, ast.operation);

			auto queryBody = ast.g[children[0]]->str;
			auto featName = ast.g[children[1]]->str;
			auto feat = static_cast<App::DocumentObject*>(mapNameToObj[featName]);

			pcActiveBody = mapNameToBody[queryBody];
			auto success = base_worker(pcActiveBody, this, "Pocket", assignValue, feat, featName, worker);

			auto* prm = pcActiveBody->getDocument()->getObject(assignValue.c_str());

			auto* featPocket = dynamic_cast<PartDesign::Pocket*>(prm);
			auto partId = std::lround(dynamic_cast<HighLevelParser::Float*>(ast.g[children.back()])->value);

			mapNameToObj[assignValue] = prm;
			updateActive();

			featPocket->Part.setValue(partId);
			updateActive();
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Make_Hole)
		{
			std::cerr << "TYPE_Make_Cut" << std::endl;
			// Add cut operation
			Gui::Command* cmd = this;
			auto worker = [](Part::Feature* sketch, App::DocumentObject* feat) {
				if (!feat)
					return;
			};

			auto children = GetOrderedChildren(ast.g, ast.operation);

			auto queryBody = ast.g[children[0]]->str;
			auto featName = ast.g[children[1]]->str;
			auto feat = static_cast<App::DocumentObject*>(mapNameToObj[featName]);

			pcActiveBody = mapNameToBody[queryBody];
			auto success = base_worker(pcActiveBody, this, "Hole", assignValue, feat, featName, worker);

			auto* prm = pcActiveBody->getDocument()->getObject(assignValue.c_str());

			auto* featHole = dynamic_cast<PartDesign::Hole*>(prm);
			auto diameter = dynamic_cast<HighLevelParser::Float*>(ast.g[children.back()])->value;
			featHole->Diameter.setValue(diameter);
			featHole->DepthType.setValue(1);
			mapNameToObj[assignValue] = prm;
			updateActive();
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Make_Stock)
		{
			std::vector<double> vecLwh;

			for (auto v : GetOrderedChildren(ast.g, ast.operation))
			{
				if (ast.g[v]->GetType() == TYPE_Float)
				{
					auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(ast.g[v]);
					vecLwh.push_back(tmpFloat->value);
					std::cerr << "pushing " << tmpFloat->value << std::endl;
				}
			}

			pcActiveBody = PartDesignGui::makeBody(doc);

			const auto& featName = assignValue;
			Gui::Command::openCommand((std::string("Make additive Box")).c_str());
			FCMD_OBJ_DOC_CMD(pcActiveBody, "addObject('PartDesign::Additive" << "Box" << "','" << featName << "')");

			auto* prm = dynamic_cast<PartDesign::Box*>(
				pcActiveBody->getDocument()->getObject(featName.c_str()));

			if (!prm) return;
			FCMD_OBJ_CMD(pcActiveBody, "addObject(" << getObjectCmd(prm) << ")");

			prm->SetLWH(vecLwh[0], vecLwh[1], vecLwh[2]);
			mapNameToObj[featName] = prm;

			updateActive();
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Line)
		{
			std::vector<double> vecValues;
			for (auto v : GetOrderedChildren(ast.g, ast.operation))
			{
				if (ast.g[v]->GetType() == TYPE_Float)
				{
					auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(ast.g[v]);
					std::cerr << tmpFloat->value << std::endl;
					vecValues.push_back(tmpFloat->value);
				}
			}

			auto geom = new Part::GeomLineSegment();
			geom->setPoints(Base::Vector3d(vecValues[0], vecValues[1], 0.0), Base::Vector3d(vecValues[2], vecValues[3], 0.0));
			mapNameToObj[assignValue] = geom;
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Circle)
		{
			std::vector<double> vecValues;
			for (auto v : GetOrderedChildren(ast.g, ast.operation))
			{
				if (ast.g[v]->GetType() == TYPE_Float)
				{
					auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(ast.g[v]);
					std::cerr << tmpFloat->value << std::endl;
					vecValues.push_back(tmpFloat->value);
				}
			}

			auto geom = new Part::GeomCircle();
			geom->setCenter(Base::Vector3d(vecValues[0], vecValues[1], 0.0));
			geom->setRadius(vecValues[2]);
			mapNameToObj[assignValue] = geom;
		}
		else if (ast.g[ast.operation]->GetType() == TYPE_Arc)
		{
			std::vector<double> vecValues;
			for (auto v : GetOrderedChildren(ast.g, ast.operation))
			{
				if (ast.g[v]->GetType() == TYPE_Float)
				{
					auto tmpFloat = dynamic_cast<HighLevelParser::Float*>(ast.g[v]);
					std::cerr << tmpFloat->value << std::endl;
					vecValues.push_back(tmpFloat->value);
				}
			}

			auto geom = new Part::GeomArcOfCircle();
			geom->setCenter(Base::Vector3d(vecValues[0], vecValues[1], 0.0));
			geom->setRadius(vecValues[2]);
			geom->setRange(vecValues[3], vecValues[4], false);
			mapNameToObj[assignValue] = geom;
		}

		// Record active body pointer everytime you create a new feature
		mapNameToBody[assignValue] = pcActiveBody;
	}
}

bool CmdParseHLHelmProgram::isActive(void)
{
	// This is always true
	return true;
	//return getActiveGuiDocument() != nullptr;
}

//===========================================================================
// Initialization
//===========================================================================

void CreatePartDesignCommands(void)
{
	Gui::CommandManager& rcCmdMgr = Gui::Application::Instance->commandManager();
	rcCmdMgr.addCommand(new CmdPartDesignVariable());
	rcCmdMgr.addCommand(new CmdPartDesignNewSketch());
	rcCmdMgr.addCommand(new CmdPartDesignPocket());
	rcCmdMgr.addCommand(new CmdPartDesignHole());
	rcCmdMgr.addCommand(new CarpentryPacking());
	rcCmdMgr.addCommand(new CmdParseHLHelmProgram());
	rcCmdMgr.addCommand(new CmdPartDesignBoolean());

}


