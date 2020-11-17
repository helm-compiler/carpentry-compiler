/***************************************************************************
 *   Copyright (c) 2008 Werner Mayer <wmayer[at]users.sourceforge.net>     *
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
#include <boost/bind.hpp>
#include <QMessageBox>
#endif

#include <Gui/Application.h>
#include <Gui/Command.h>
#include <Gui/Control.h>
#include <Gui/MDIView.h>
#include <Gui/MenuManager.h>
#include <Gui/DockWindowManager.h>
#include <Gui/CombiView.h>
#include <Gui/MainWindow.h>
#include <Gui/Command.h>

#include <Mod/Sketcher/Gui/Workbench.h>
#include <Mod/PartDesign/App/Body.h>
#include <Mod/PartDesign/App/Feature.h>
#include <Mod/PartDesign/App/FeatureMultiTransform.h>

#include "Utils.h"

#include "Workbench.h"

#include "WorkflowManager.h"

#include <Mod/PartDesign/Gui/CodeEditor.h>
#include <Mod/PartDesign/Gui/ToolsConfiguration.h>
#include <Mod/PartDesign/Gui/PHELM.h>
#include <Mod/PartDesign/Gui/CompilerConfig.h>
#include <Mod/PartDesign/App/FeaturePocket.h>

#include <QMessageBox>
#include <QFileDialog>

using namespace PartDesignGui;

#if 0 // needed for Qt's lupdate utility
qApp->translate("Workbench", "Part Design");
qApp->translate("Gui::TaskView::TaskWatcherCommands", "Face tools");
qApp->translate("Gui::TaskView::TaskWatcherCommands", "Sketch tools");
qApp->translate("Gui::TaskView::TaskWatcherCommands", "Create Geometry");
#endif

/// @namespace PartDesignGui @class Workbench
TYPESYSTEM_SOURCE(PartDesignGui::Workbench, Gui::StdWorkbench)

Workbench::Workbench(const string& config_file /*= ""*/)
{
	std::string helm_config_path = "HELMConfig.xml";

	if (!config_file.empty()) helm_config_path = config_file;

	std::ifstream configFile(helm_config_path);

	if (!configFile)
	{
		QMessageBox msgBox;
		msgBox.setText(QObject::tr("Can't load configuration file."));
		msgBox.setInformativeText(QObject::tr("Do you agree to use default configuration? Or you can load it from your disk (Choose No)."));
		msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
		msgBox.setDefaultButton(QMessageBox::Yes);
		msgBox.setMinimumSize(280, 100);
		int ret = msgBox.exec();

		QString fileName = QString();
		if (ret == QMessageBox::No)
		{
			fileName = QFileDialog::getOpenFileName(Gui::getMainWindow(), QObject::tr("Load Config"),
				QObject::tr(helm_config_path.c_str()), QObject::tr("HELM Config Files (*.xml)"));
		}

		if (fileName.isEmpty())
		{
			CompilerConfig::Instance().WriteDefaultConfig(helm_config_path);
		}
	}

	if (!CompilerConfig::Instance().HasLoaded())
	{
		CompilerConfig::Instance().LoadConfigFile(helm_config_path);
	}

	// Load configuration file before instancing objects
#ifndef MAINOBJ
	hlEditor = new CodeEditor(0);
	llEditor = new CodeEditor(0);
	tools = new ToolsConfiguration(0);
#endif

	langCompiler = new PHELM();
	pce = new PCE();
}

Workbench::~Workbench()
{
	if (langCompiler != nullptr)
		delete langCompiler;
	WorkflowManager::destruct();
}

void Workbench::_switchToDocument(const App::Document* /*doc*/)
{
	// TODO Commented out for thurther remove or rewrite  (2015-09-04, Fat-Zer)
	//    if (doc == NULL) return;
	//
	//    PartDesign::Body* activeBody = NULL;
	//    std::vector<App::DocumentObject*> bodies = doc->getObjectsOfType(PartDesign::Body::getClassTypeId());
	//
	//    // No tip, so build up structure or migrate
	//    if (!doc->Tip.getValue())
	//    {
	//        ;/*if (doc->countObjects() == 0){
	//            buildDefaultPartAndBody(doc);
	//            activeBody = Gui::Application::Instance->activeView()->getActiveObject<PartDesign::Body*>(PDBODYKEY);
	//            assert(activeBody);
	//        } else {
	//            // empty document with no tip, so do migration
	//            _doMigration(doc);
	//                        activeBody = Gui::Application::Instance->activeView()->getActiveObject<PartDesign::Body*>(PDBODYKEY);
	//                        assert(activeBody);
	//                }
	//                */
	//    }
	//    else
	//    {
	//      App::Part *docPart = dynamic_cast<App::Part *>(doc->Tip.getValue());
	//      if (docPart) {
	//          App::Part *viewPart = Gui::Application::Instance->activeView()->getActiveObject<App::Part *>("Part");
	//          if (viewPart != docPart)
	//            Gui::Application::Instance->activeView()->setActiveObject(docPart, "Part");
	//          //if (docPart->countObjectsOfType(PartDesign::Body::getClassTypeId()) < 1)
	//          //  setUpPart(docPart);
	//          PartDesign::Body *tempBody = dynamic_cast<PartDesign::Body *> (docPart->getObjectsOfType(PartDesign::Body::getClassTypeId()).front());
	//          if (tempBody) {
	//              PartDesign::Body *viewBody = Gui::Application::Instance->activeView()->getActiveObject<PartDesign::Body*>(PDBODYKEY);
	//              activeBody = viewBody;
	//              if (!viewBody)
	//                activeBody = tempBody;
	//              else if (!docPart->hasObject(viewBody))
	//                activeBody = tempBody;
	//
	//              if (activeBody != viewBody)
	//                Gui::Application::Instance->activeView()->setActiveObject(activeBody, PDBODYKEY);
	//          }
	//      }
	//    }
	//
	//    /*if (activeBody == NULL) {
	//        QMessageBox::critical(Gui::getMainWindow(), QObject::tr("Could not create body"),
	//            QObject::tr("No body was found in this document, and none could be created. Please report this bug."
	//                        "We recommend you do not use this document with the PartDesign workbench until the bug has been fixed."
	//                        ));
	//    }*/
}

void Workbench::slotActiveDocument(const Gui::Document& /*Doc*/)
{
	//     _switchToDocument(Doc.getDocument());
}

void Workbench::slotNewDocument(const App::Document& /*Doc*/)
{
	//     _switchToDocument(&Doc);
}

void Workbench::slotFinishRestoreDocument(const App::Document& /*Doc*/)
{
	//     _switchToDocument(&Doc);
}

void Workbench::slotDeleteDocument(const App::Document&)
{
	//ActivePartObject = 0;
	//ActiveGuiDoc = 0;
	//ActiveAppDoc = 0;
	//ActiveVp = 0;
}
/*
  This does not work for Std_DuplicateSelection:
  Tree.cpp gives: "Cannot reparent unknown object", probably because the signalNewObject is emitted
  before the duplication of the object has been completely finished

void Workbench::slotNewObject(const App::DocumentObject& obj)
{
	if ((obj.getDocument() == ActiveAppDoc) && (ActivePartObject != NULL)) {
		// Add the new object to the active Body
		// Note: Will this break Undo? But how else can we catch Edit->Duplicate selection?
		Gui::Command::doCommand(Gui::Command::Doc,"App.activeDocument().%s.addObject(App.activeDocument().%s)",
								ActivePartObject->getNameInDocument(), obj.getNameInDocument());
	}
}
*/

void Workbench::setupContextMenu(const char* recipient, Gui::MenuItem* item) const
{
	auto selection = Gui::Selection().getSelection();
	// Add move Tip Command
	if (selection.size() >= 1)
	{
		App::DocumentObject* feature = selection.front().pObject;
		PartDesign::Body* body = nullptr;

		// if PD workflow is not new-style then add a command to the context-menu
		bool assertModern = true;
		if (feature && !isModernWorkflow(feature->getDocument()))
		{
			assertModern = false;
			*item << "PartDesign_Migrate";
		}

		body = PartDesignGui::getBodyFor(feature, false, false, assertModern);
		// lote of assertion so feature should be marked as a tip
		if (selection.size() == 1 && feature && (feature->isDerivedFrom(PartDesign::Body::getClassTypeId()) || (feature->isDerivedFrom(PartDesign::Feature::getClassTypeId()) && body) || (feature->isDerivedFrom(Part::Feature::getClassTypeId()) && body && body->BaseFeature.getValue() == feature)))
		{
			*item << "PartDesign_MoveTip";
		}

		if (strcmp(recipient, "Tree") == 0)
		{
			Gui::MDIView* activeView = Gui::Application::Instance->activeView();

			if (selection.size() > 0 && activeView)
			{
				bool docHaveBodies = activeView->getAppDocument()->countObjectsOfType(
					PartDesign::Body::getClassTypeId()) > 0;

				if (docHaveBodies)
				{
					bool addMoveFeature = true;
					bool addMoveFeatureInTree = (body != nullptr);
					for (auto sel : selection)
					{
						// if at least one selected feature cannot be moved to a body
						// disable the entry
						if (addMoveFeature && !PartDesign::Body::isAllowed(sel.pObject))
						{
							addMoveFeature = false;
						}
						// if all at least one selected feature doesn't belong to the same body
						// disable the menu entry
						if (addMoveFeatureInTree && !body->hasObject(sel.pObject))
						{
							addMoveFeatureInTree = false;
						}

						if (!addMoveFeatureInTree && !addMoveFeature)
						{
							break;
						}
					}

					if (addMoveFeature)
					{
						*item << "PartDesign_MoveFeature";
					}

					if (addMoveFeatureInTree)
					{
						*item << "PartDesign_MoveFeatureInTree";
					}
				}
			}

			if (Gui::Selection().countObjectsOfType(PartDesign::Transformed::getClassTypeId()) -
				Gui::Selection().countObjectsOfType(PartDesign::MultiTransform::getClassTypeId()) ==
				1)
				*item << "PartDesign_MultiTransform";

			if (Gui::Selection().countObjectsOfType(App::DocumentObject::getClassTypeId()) > 0)
			{
				*item << "Std_SetAppearance"
					<< "Std_RandomColor";
			}
		}
	}
}

void Workbench::CheckFeatures()
{
	PartDesign::Body* pcActiveBody = PartDesignGui::getBody(/*messageIfNot = */ true);

	if (!toCheck.empty())
	{
		for (auto feat : toCheck)
		{
			if (feat->isDerivedFrom(PartDesign::Pocket::getClassTypeId()))
			{
				auto cut = static_cast<PartDesign::Pocket*>(feat);
				auto part2dObj = cut->getVerifiedSketch();
				if (part2dObj->isDerivedFrom(Sketcher::SketchObject::getClassTypeId()))
				{
					auto sketch = static_cast<Sketcher::SketchObject*>(part2dObj);
					const auto& profileWires = cut->getProfileWires();

					if (profileWires.empty())
					{
						QMessageBox::warning(Gui::getMainWindow(), QObject::tr("Empty error"),
							QObject::tr("You basically drew nothing in sketch, won't be effective."));
						continue;
					}
				}
			}
		}
	}

	toCheck.clear();
}

void Workbench::activated()
{
	// Before activating the workbench, check the features first
	CheckFeatures();

	Gui::Workbench::activated();

	WorkflowManager::init();

	std::vector<Gui::TaskView::TaskWatcher*> Watcher;

	const char* Vertex[] = {
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT Part::Feature SUBELEMENT Vertex COUNT 1..",
		Vertex,
		"Vertex tools",
		"Part_Box"));

	const char* Edge[] = {
		"PartDesign_Fillet",
		"PartDesign_Chamfer",
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT Part::Feature SUBELEMENT Edge COUNT 1..",
		Edge,
		"Edge tools",
		"Part_Box"));

	const char* Face[] = {
		"PartDesign_NewSketch",
		"PartDesign_Fillet",
		"PartDesign_Chamfer",
		"PartDesign_Draft",
		"PartDesign_Thickness",
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT Part::Feature SUBELEMENT Face COUNT 1",
		Face,
		"Face tools",
		"Part_Box"));

	const char* Body[] = {
		"PartDesign_NewSketch",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::Body COUNT 1",
		Body,
		"Start Body",
		"Part_Box"));

	const char* Body2[] = {
		"PartDesign_Boolean",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::Body COUNT 1..",
		Body2,
		"Start Boolean",
		"Part_Box"));

	const char* Plane1[] = {
		"PartDesign_NewSketch",
		"PartDesign_Plane",
		"PartDesign_Line",
		"PartDesign_Point",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT App::Plane COUNT 1",
		Plane1,
		"Start Part",
		"Part_Box"));
	const char* Plane2[] = {
		"PartDesign_NewSketch",
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::Plane COUNT 1",
		Plane2,
		"Start Part",
		"Part_Box"));

	const char* Line[] = {
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::Line COUNT 1",
		Line,
		"Start Part",
		"Part_Box"));

	const char* Point[] = {
		"PartDesign_Point",
		"PartDesign_Line",
		"PartDesign_Plane",
		"PartDesign_CoordinateSystem",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::Point COUNT 1",
		Point,
		"Start Part",
		"Part_Box"));

	const char* NoSel[] = {
		"PartDesign_Body",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommandsEmptySelection(
		NoSel,
		"Start Part",
		"Part_Box"));

	const char* Faces[] = {
		"PartDesign_Fillet",
		"PartDesign_Chamfer",
		"PartDesign_Draft",
		"PartDesign_Thickness",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT Part::Feature SUBELEMENT Face COUNT 2..",
		Faces,
		"Face tools",
		"Part_Box"));

	const char* Sketch[] = {
		"PartDesign_NewSketch",
		"PartDesign_Pad",
		"PartDesign_Pocket",
		"PartDesign_Hole",
		"PartDesign_Revolution",
		"PartDesign_Groove",
		"PartDesign_AdditivePipe",
		"PartDesign_SubtractivePipe",
		"PartDesign_AdditiveLoft",
		"PartDesign_SubtractiveLoft",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT Sketcher::SketchObject COUNT 1",
		Sketch,
		"Sketch tools",
		"Part_Box"));

	const char* Transformed[] = {
		"PartDesign_Mirrored",
		"PartDesign_LinearPattern",
		"PartDesign_PolarPattern",
		//        "PartDesign_Scaled",
		"PartDesign_MultiTransform",
		0 };
	Watcher.push_back(new Gui::TaskView::TaskWatcherCommands(
		"SELECT PartDesign::SketchBased",
		Transformed,
		"Transformation tools",
		"PartDesign_MultiTransform"));

	Gui::DockWindowManager* pDockMgr = Gui::DockWindowManager::instance();
	Gui::DockWnd::CombiView* combiView = qobject_cast<Gui::DockWnd::CombiView*>(pDockMgr->getDockWindow("Combo View"));
	QTabWidget* tabs = combiView->getTabPanel();

	const int tabToolIdx = tabs->addTab(tools, QObject::trUtf8("Tools"));
	tools->setTabIndex(tabToolIdx);

	const int tabHLIdx = tabs->addTab(hlEditor, QObject::trUtf8("HL-HELM"));
	hlEditor->setTabIndex(tabHLIdx);

	const int tabLLIdx = tabs->addTab(llEditor, QObject::trUtf8("LL-HELM"));
	llEditor->setTabIndex(tabLLIdx);

	// make the previously used active Body active again
	//PartDesignGui::ActivePartObject = NULL;
	_switchToDocument(App::GetApplication().getActiveDocument());

	addTaskWatcher(Watcher);
	Gui::Control().showTaskView();

	// Let us be notified when a document is activated, so that we can update the ActivePartObject
	Gui::Application::Instance->signalActiveDocument.connect(boost::bind(&Workbench::slotActiveDocument, this, _1));
	App::GetApplication().signalNewDocument.connect(boost::bind(&Workbench::slotNewDocument, this, _1));
	App::GetApplication().signalFinishRestoreDocument.connect(boost::bind(&Workbench::slotFinishRestoreDocument, this, _1));
	App::GetApplication().signalDeleteDocument.connect(boost::bind(&Workbench::slotDeleteDocument, this, _1));
	// Watch out for objects being added to the active document, so that we can add them to the body
	//App::GetApplication().signalNewObject.connect(boost::bind(&Workbench::slotNewObject, this, _1));
}

void Workbench::deactivated()
{
	// Let us be notified when a document is activated, so that we can update the ActivePartObject
	Gui::Application::Instance->signalActiveDocument.disconnect(boost::bind(&Workbench::slotActiveDocument, this, _1));
	App::GetApplication().signalNewDocument.disconnect(boost::bind(&Workbench::slotNewDocument, this, _1));
	App::GetApplication().signalFinishRestoreDocument.disconnect(boost::bind(&Workbench::slotFinishRestoreDocument, this, _1));
	App::GetApplication().signalDeleteDocument.disconnect(boost::bind(&Workbench::slotDeleteDocument, this, _1));
	//App::GetApplication().signalNewObject.disconnect(boost::bind(&Workbench::slotNewObject, this, _1));

	// remove code editor
	Gui::DockWindowManager* pDockMgr = Gui::DockWindowManager::instance();
	Gui::DockWnd::CombiView* combiView = qobject_cast<Gui::DockWnd::CombiView*>(pDockMgr->getDockWindow("Combo View"));
	QTabWidget* tabs = combiView->getTabPanel();
	tabs->removeTab(tools->getTabIndex());
	tabs->removeTab(hlEditor->getTabIndex());
	tabs->removeTab(llEditor->getTabIndex());

	removeTaskWatcher();
	// reset the active Body
	Gui::Command::doCommand(Gui::Command::Doc, "import PartDesignGui");

	Gui::Workbench::deactivated();
}

void Workbench::AddFeatureToCheck(PartDesign::Feature* f)
{
	toCheck.push_back(f);
}

CodeEditor* PartDesignGui::Workbench::GetCodeEditor()
{
	return nullptr;
}

ToolsConfiguration* PartDesignGui::Workbench::GetToolsConfig()
{
	return nullptr;
}

PHELM* PartDesignGui::Workbench::GetPHELM()
{
	return nullptr;
}

Gui::MenuItem* Workbench::setupMenuBar() const
{
	Gui::MenuItem* root = StdWorkbench::setupMenuBar();
	Gui::MenuItem* item = root->findItem("&Windows");

	Gui::MenuItem* part = new Gui::MenuItem;
	root->insertItem(item, part);
	part->setCommand("&Part Design");
	*part << "PartDesign_Body"
		<< "PartDesign_NewSketch"
		<< "Sketcher_LeaveSketch"
		<< "Sketcher_ViewSketch"
		<< "Sketcher_MapSketch"
		<< "Sketcher_ReorientSketch"
		<< "Sketcher_ValidateSketch"
		<< "Separator"
		<< "PartDesign_Point"
		<< "PartDesign_Line"
		<< "PartDesign_Plane"
		<< "PartDesign_CoordinateSystem"
		<< "PartDesign_ShapeBinder"
		<< "PartDesign_SubShapeBinder"
		<< "PartDesign_Clone"
		<< "Separator"
		<< "PartDesign_Pad"
		<< "PartDesign_Revolution"
		<< "PartDesign_AdditiveLoft"
		<< "PartDesign_AdditivePipe"
		<< "PartDesign_CompPrimitiveAdditive"
		<< "Separator"
		<< "PartDesign_Pocket"
		<< "PartDesign_Hole"
		<< "PartDesign_Groove"
		<< "PartDesign_SubtractiveLoft"
		<< "PartDesign_SubtractivePipe"
		<< "PartDesign_CompPrimitiveSubtractive"
		<< "Separator"
		<< "PartDesign_Mirrored"
		<< "PartDesign_LinearPattern"
		<< "PartDesign_PolarPattern"
		//          << "PartDesign_Scaled"
		<< "PartDesign_MultiTransform"
		<< "Separator"
		<< "PartDesign_Fillet"
		<< "PartDesign_Chamfer"
		<< "PartDesign_Draft"
		<< "PartDesign_Thickness"
		<< "Separator"
		<< "PartDesign_Boolean"
		<< "Separator"
		//<< "PartDesign_Hole"
		<< "PartDesign_InvoluteGear"
		<< "Separator"
		<< "PartDesign_Migrate";

	// For 0.13 a couple of python packages like numpy, matplotlib and others
	// are not deployed with the installer on Windows. Thus, the WizardShaft is
	// not deployed either hence the check for the existence of the command.
	if (Gui::Application::Instance->commandManager().getCommandByName("PartDesign_InvoluteGear"))
	{
		*part << "PartDesign_InvoluteGear";
	}
	if (Gui::Application::Instance->commandManager().getCommandByName("PartDesign_WizardShaft"))
	{
		*part << "Separator"
			<< "PartDesign_WizardShaft";
	}

	// Replace the "Duplicate selection" menu item with a replacement that is compatible with Body
	item = root->findItem("&Edit");
	Gui::MenuItem* dup = item->findItem("Std_DuplicateSelection");
	dup->setCommand("PartDesign_DuplicateSelection");

	return root;
}

Gui::ToolBarItem* Workbench::setupToolBars() const
{
	Gui::ToolBarItem* root = StdWorkbench::setupToolBars();
	Gui::ToolBarItem* part = new Gui::ToolBarItem(root);
	part->setCommand("Part Design Helper");
	*part << "PartDesign_Body"
		//<< "PartDesign_NewSketch"
		//<< "Sketcher_EditSketch"
		//<< "Sketcher_MapSketch"
		<< "Separator"
		// 		<< "PartDesign_Point"
		// 		<< "PartDesign_Line"
		// 		<< "PartDesign_Plane"
		// 		<< "PartDesign_CoordinateSystem"
		<< "PartDesign_Variable";
	//<< "PartDesign_ShapeBinder"
	//<< "PartDesign_Clone";

	part = new Gui::ToolBarItem(root);
	part->setCommand("Part Design Modeling");
	*part //<< "PartDesign_Pad"
		  //<< "PartDesign_Revolution"
		  //<< "PartDesign_AdditiveLoft"
		  //           << "PartDesign_AdditivePipe"
		<< "PartDesign_CompPrimitiveAdditive"
		//           << "Separator"
		<< "PartDesign_Pocket"
		<< "PartDesign_Hole"
		//           << "PartDesign_Groove"
		//           << "PartDesign_SubtractiveLoft"
		//           << "PartDesign_SubtractivePipe"
		//           << "PartDesign_CompPrimitiveSubtractive"
		//           << "Separator"
		//           << "PartDesign_Mirrored"
		//           << "PartDesign_LinearPattern"
		//           << "PartDesign_PolarPattern"
		// //          << "PartDesign_Scaled"
		//           << "PartDesign_MultiTransform"
		//           << "Separator"
		//           << "PartDesign_Fillet"
		//           << "PartDesign_Chamfer"
		//           << "PartDesign_Draft"
		//           << "PartDesign_Thickness"
		<< "Separator"
		<< "PartDesign_Boolean"
		<< "Carpentry_Packing_Visualization"
		<< "Carpentry_Packing"
		<< "CmdParseHLHelmProgram";
	return root;
}

Gui::ToolBarItem* Workbench::setupCommandBars() const
{
	// Part tools
	Gui::ToolBarItem* root = new Gui::ToolBarItem;
	return root;
}