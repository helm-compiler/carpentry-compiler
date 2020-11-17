/***************************************************************************
 *   Copyright (c) 2010 Juergen Riegel <FreeCAD@juergen-riegel.net>        *
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
#include <Bnd_Box.hxx>
#include <gp_Dir.hxx>
#include <BRep_Builder.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepFeat_MakePrism.hxx>
#include <Geom_Plane.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#endif

#include <Base/Console.h>
#include <Base/Exception.h>
#include <Base/Placement.h>
#include <App/Document.h>
#include <App/Application.h>
#include <Mod/PartDesign/Gui/cgal.h>

#include "FeaturePocket.h"
#include "GEOMAlgo_Splitter.hxx"
#include "../Gui/GeomCommonFunctions.h"

using namespace PartDesign;

const char* Pocket::TypeEnums[] = { "Length", "ThroughAll", "UpToFirst", "UpToFace", "TwoLengths", NULL };

PROPERTY_SOURCE(PartDesign::Pocket, PartDesign::ProfileBased)

Pocket::Pocket()
{
	addSubType = FeatureAddSub::Subtractive;

	ADD_PROPERTY_TYPE(Type, ((long)1), "Pocket", App::Prop_None, "Pocket type");
	Type.setEnums(TypeEnums);
	ADD_PROPERTY_TYPE(Length, (100.0), "Pocket", App::Prop_None, "Pocket length");
	ADD_PROPERTY_TYPE(Length2, (100.0), "Pocket", App::Prop_None, "P");
	ADD_PROPERTY_TYPE(UpToFace, (0), "Pocket", App::Prop_None, "Face where pocket will end");
	ADD_PROPERTY_TYPE(Offset, (0.0), "Pocket", App::Prop_None, "Offset from face in which pocket will end");
	ADD_PROPERTY_TYPE(Part, (long(0)), "Pocket", App::Prop_None, "Numer of parts");
	static const App::PropertyQuantityConstraint::Constraints signedLengthConstraint = { -DBL_MAX, DBL_MAX, 1.0 };
	Offset.setConstraints(&signedLengthConstraint);
}

short Pocket::mustExecute() const
{
	if (Placement.isTouched() ||
		Type.isTouched() ||
		Length.isTouched() ||
		Length2.isTouched() ||
		Offset.isTouched() ||
		UpToFace.isTouched() ||
		Part.isTouched())
		return 1;
	return ProfileBased::mustExecute();
}
static std::vector<gp_Pnt> GetSortedPointsOnWire(const TopoDS_Wire& wire)
{
	std::vector<gp_Pnt> points;
	const auto sortedWire = wire;// GeomFunc::SortWireEdges(wire);
	TopoDS_Edge curEdge;
	for (BRepTools_WireExplorer verExp(sortedWire); verExp.More(); verExp.Next())
	{
		curEdge = TopoDS::Edge(verExp.Current());
		gp_Pnt mp = BRep_Tool::Pnt(TopoDS::Vertex(verExp.CurrentVertex()));
		points.push_back(mp);
	}

	if (!wire.Closed())
	{
		TopTools_IndexedMapOfShape vertices;
		TopExp::MapShapes(curEdge, TopAbs_VERTEX, vertices);

		if (vertices.Extent() != 2)
			throw std::bad_exception();

		const TopoDS_Vertex& vert = TopoDS::Vertex(vertices(2));
		points.emplace_back(BRep_Tool::Pnt(vert));
	}
	return points;
}
App::DocumentObjectExecReturn* Pocket::execute(void)
{
	// Handle legacy features, these typically have Type set to 3 (previously NULL, now UpToFace),
	// empty FaceName (because it didn't exist) and a value for Length
	if (std::string(Type.getValueAsString()) == "UpToFace" &&
		(UpToFace.getValue() == NULL && Length.getValue() > Precision::Confusion()))
		Type.setValue("Length");

	// Validate parameters
	double L = Length.getValue();
	if ((std::string(Type.getValueAsString()) == "Length") && (L < Precision::Confusion()))
		return new App::DocumentObjectExecReturn("Pocket: Length of pocket too small");

	double L2 = Length2.getValue();
	if ((std::string(Type.getValueAsString()) == "TwoLengths") && (L < Precision::Confusion()))
		return new App::DocumentObjectExecReturn("Pocket: Second length of pocket too small");

	Part::Feature* obj = 0;
	Part::TopoShape profileshape;
	try
	{
		obj = getVerifiedObject();
		//profileshape = getVerifiedFace();
		auto wires = getProfileWires();
		if (wires.size() != 1)
			return new App::DocumentObjectExecReturn("Pocket: The size of wires is not equal to 1");
		profileshape = wires[0];
	}
	catch (const Base::Exception & e)
	{
		return new App::DocumentObjectExecReturn(e.what());
	}

	// if the Base property has a valid shape, fuse the prism into it
	TopoShape base;
	try
	{
		base = getBaseShape();
	}
	catch (const Base::Exception&)
	{
		return new App::DocumentObjectExecReturn("No sketch support and no base shape: Please tell me where to remove the material of the pocket!");
	}

	// get the Sketch plane
	Base::Placement SketchPos = obj->Placement.getValue();
	Base::Vector3d SketchVector = getProfileNormal();

	// turn around for pockets
	SketchVector *= -1;

	try
	{
		this->positionByPrevious();
		TopLoc_Location invObjLoc = this->getLocation().Inverted();

		base.move(invObjLoc);

		gp_Dir dir(SketchVector.x, SketchVector.y, SketchVector.z);
		dir.Transform(invObjLoc.Transformation());
		if (profileshape.isNull())
			return new App::DocumentObjectExecReturn("Pocket: Creating a face from sketch failed");
		profileshape.move(invObjLoc);
		
		std::string method(Type.getValueAsString());

		if (base.isNull())
			return new App::DocumentObjectExecReturn("Pocket: Extruding up to a face is only possible if the sketch is located on a face");

		TopoShape surface(0, getDocument()->getStringHasher());
		generatePrism(surface, profileshape, method, dir, L, L2,
			Midplane.getValue(), Reversed.getValue());

		if (surface.isNull())
			return new App::DocumentObjectExecReturn("Pocket: Resulting shape is empty");

		GEOMAlgo_Splitter splitter;
		splitter.AddArgument(base.getShape());
		splitter.AddTool(surface.getShape());
		splitter.Perform();

		TopTools_ListIteratorOfListOfShape iter(splitter.Modified(base.getShape()));

		std::vector<std::string> PartNames;
		int cnt = 0;
		for (; iter.More(); iter.Next(), ++cnt)
		{
			PartNames.push_back(std::to_string(cnt));
		}

		int selectedPart = (Part.getValue() >= cnt || Part.getValue() < 0) ? 0 : Part.getValue();
		PartEnums.setEnums(PartNames);
		Part.setValue(PartEnums);
		Part.setValue(long(selectedPart));

		//std::cerr << "number of components = " << cnt << std::endl;
		//std::cerr << "number of part = " << Part.getValue() << std::endl;

		iter = TopTools_ListIteratorOfListOfShape(splitter.Modified(base.getShape()));

		// extract target part
		if (cnt != 0)
		{
			for (int i = 0; iter.More(); iter.Next(), ++i)
			{
				if (i == selectedPart)
				{
					TopoDS_Shape result = iter.Value();
					auto solRes = this->getSolid(result);
					if (solRes.isNull())
						return new App::DocumentObjectExecReturn("Pocket: Resulting shape is not a solid");

					solRes = refineShapeIfActive(solRes);
					remapSupportShape(solRes.getShape());
					this->Shape.setValue(getSolid(solRes));
				}
			}
		}

		// identify the belongings of the faces
		TopoDS_Face reservedFace;
		if (cnt == 2)
		{
			auto supportFace = getSupportFace();
			supportFace.Move(invObjLoc);
			const TopTools_ListOfShape& splits = splitter.Modified(supportFace);
			std::cerr << "size of splits = " << splits.Size() << std::endl;
			TopTools_ListIteratorOfListOfShape it;
			for (it.Initialize(splits); it.More(); it.Next())
			{
				TopoDS_Face splitface = TopoDS::Face(it.Value());
				gp_Pnt p = GeomFunc::GetCentralFacePoint(splitface);

				// check if point is in shapeToExclude
				BRepClass3d_SolidClassifier classifier;
				classifier.Load(this->Shape.getValue());
				classifier.Perform(p, Precision::Confusion());

				if (classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON)
				{
					reservedFace = splitface;
					std::cerr << "you select the right one" << std::endl;
					break;
				}
				else
				{
					std::cerr << "you select the wrong one" << std::endl;
				}
			}
		}

		return App::DocumentObject::StdReturn;
	}
	catch (Standard_Failure & e)
	{
		if (std::string(e.GetMessageString()) == "TopoDS::Face" &&
			(std::string(Type.getValueAsString()) == "UpToFirst" || std::string(Type.getValueAsString()) == "UpToFace"))
			return new App::DocumentObjectExecReturn("Could not create face from sketch.\n"
				"Intersecting sketch entities or multiple faces in a sketch are not allowed "
				"for making a pocket up to a face.");
		else
			return new App::DocumentObjectExecReturn(e.GetMessageString());
	}
	catch (Base::Exception & e)
	{
		return new App::DocumentObjectExecReturn(e.what());
	}
}