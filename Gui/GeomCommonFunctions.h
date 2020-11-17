#ifndef _GEOM_COMMON_FUNCTIONS
#define _GEOM_COMMON_FUNCTIONS
#include <Standard_TypeDef.hxx>
#include <TopoDS_Edge.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TopTools_ListOfShape.hxx>
#include "TColgp_HArray1OfPnt.hxx"
#include "Geom_Curve.hxx"
#include "Geom_Surface.hxx"
#include "GeomAdaptor_Curve.hxx"
#include "BRep_Tool.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Edge.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopTools_HSequenceOfShape.hxx"
#include "GeomAdaptor_Curve.hxx"
#include "BRepAdaptor_CompCurve.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "BRepTools.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"
#include "BRepExtrema_ExtCF.hxx"
#include "BRepIntCurveSurface_Inter.hxx"
#include "BRepFill.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "BRepClass3d_SolidClassifier.hxx"
#include "BRepClass_FaceClassifier.hxx"
#include "IGESControl_Reader.hxx"

#include <Approx_Curve3d.hxx>
#include <BRepAdaptor_HCompCurve.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include "BRepExtrema_ExtCC.hxx"
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepFill_CompatibleWires.hxx>
#include <BRepFill_Filling.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2dAPI_InterCurveCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <GeomConvert.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeFix_EdgeConnect.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <ShapeFix_Wire.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <Standard_Version.hxx>
#include <StlAPI_Writer.hxx>
#include <StlAPI_Reader.hxx>
#include <IGESControl_Controller.hxx> 
#include <IGESControl_Writer.hxx>
#include <Base/Vector3D.h>
#include <ShapeAnalysis_FreeBounds.hxx>
#include "GEOMAlgo_Splitter.hxx"
#include "BRepPrimAPI_MakePrism.hxx"
#include "BRepOffsetAPI_Sewing.hxx"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopoDS.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include "BRepLib_FindSurface.hxx"
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepExtrema_ShapeProximity.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <Geom_Circle.hxx>


#include <list>
#include <algorithm>
#include <cassert>
#include <limits>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <Base/Exception.h>
#include <unordered_set>
#include <Mod/PartDesign/Gui/math.hpp>
#include <Mod/PartDesign/Gui/cgal.h>

constexpr auto INCH_OVER_16 = 1.5875;

using namespace Math;

namespace
{
	struct IsSame
	{
		IsSame(double tolerance) : _tol(tolerance) {}
		bool operator() (double first, double second)
		{
			return (fabs(first - second) < _tol);
		}

		double _tol;
	};

} // anonymous namespace


class GeomFunc {
public:
	typedef std::unordered_set<std::string> SetArrangement;

	static std::string vector2string(const std::vector<int>& vec)
	{
		std::stringstream ss;
		for (auto& v : vec) ss << v;
		return ss.str();
	}

	static int GetAnglePrecision(double num)
	{
		num = num * 10;
		int inter = Round(num);
		if (IsAlmostZero_Double(num - static_cast<double>(inter), 1e-4))
		{
			if (inter % 450 == 0) return 0;
			else if (inter % 150 == 0) return 1;
			else if (inter % 10 == 0) return 2;
			else return 3;
		}
		else
			return 4;
	}

	static Vector3d GetAnyPerpendicularVector(Vector3d& in)
	{
		auto a = in[0];
		auto b = in[1];
		auto c = in[2];
		int selectIndex = !Math::Functs::IsAlmostZero(c) && !Math::Functs::IsAlmostZero(b + a);
		if (selectIndex)
		{
			return Vector3d(-b - c, a, a);
		}
		else
		{
			return Vector3d(c, c, -a - b);
		}
	}

	static int GetDimensionPrecision(double num)
	{
		double facByInch = num / INCH_OVER_16;
		int closedInch = Round(facByInch);
		if (IsAlmostZero_Double(facByInch - static_cast<double>(closedInch), 1e-4))
		{
			const int inchMod = closedInch % 16;
			switch (inchMod)
			{
			case 0:
				return 0;
				break;
			case 8:
				return 1;
				break;
			default:
				return 2;
				break;
			}
		}
		else
		{
			// has rounding error
			return 4;
		}

		/*
		num = num * 10;
		int inter = Round(num);
		if (IsAlmostZero(num - static_cast<double>(inter)))
		{
			if (inter % 50 == 0) return 1;
			else if (inter % 10 == 0) return 2;
			else return 3;
		}
		else
			return 4;
		*/
	}


	static Standard_Real GetLength(const TopoDS_Wire& wire)
	{
		BRepAdaptor_CompCurve aCompoundCurve(wire, Standard_True);
		return GCPnts_AbscissaPoint::Length(aCompoundCurve);
	}

	static gp_Dir GetFaceNormal(const TopoDS_Face& face)
	{
		Standard_Real umin, umax, vmin, vmax;
		BRepTools::UVBounds(face, umin, umax, vmin, vmax);
		Handle(Geom_Surface) aSurface = BRep_Tool::Surface(face);
		GeomLProp_SLProps props(aSurface, umin, vmin, 1, 0.01);
		gp_Dir normal = props.Normal();

		if (face.Orientation() == TopAbs_REVERSED)
		{
			normal = -normal;
		}
		return normal;
	}

	static Standard_Real GetEdgeLength(const TopoDS_Edge& edge)
	{
		gp_Pnt pts[2];
		int cnt = 0;
		for (TopExp_Explorer vertexExplorer(edge, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next(), ++cnt)
		{
			const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
			pts[cnt] = BRep_Tool::Pnt(vertex);
		}
		return (pts[0].Distance(pts[1]));
	}

	static bool AreTwoShapesSame_Mesh(const TopoDS_Shape& s1, const TopoDS_Shape& s2, const double& part_match_error_)
	{
		Vector3d2 mesh_0 = GeomFunc::LoadTopoDSGeometry(s1);
		Vector3d2 mesh_1 = GeomFunc::LoadTopoDSGeometry(s2);

		double d = Math::Functs::GetDistanceNONE(mesh_0, mesh_1);
		return Math::Functs::IsAlmostZero_Double(d, part_match_error_);
	}

	static bool AreTwoShapesSame(const TopoDS_Shape& s1, const TopoDS_Shape& s2, const double& part_match_error_)
	{
		GProp_GProps system1, system2;
		// 		BRepGProp::LinearProperties(s1, system1);
		// 		BRepGProp::LinearProperties(s2, system2);
		// 		std::cerr << system1.Mass() << " " << system2.Mass() << std::endl;
		// 		if (!Math::Functs::AreAlmostEqual(system1.Mass(), system2.Mass())) return false;

		BRepGProp::SurfaceProperties(s1, system1);
		BRepGProp::SurfaceProperties(s2, system2);
		//std::cerr << system1.Mass() << " " << system2.Mass() << std::endl;
		if (!Math::Functs::AreAlmostEqual(system1.Mass(), system2.Mass())) return false;

		BRepGProp::VolumeProperties(s1, system1);
		BRepGProp::VolumeProperties(s2, system2);
		if (!Math::Functs::AreAlmostEqual(system1.Mass(), system2.Mass())) return false;

		Base::Vector3d c1, c2;
		if (!GetCenterOfGravity(s1, c1)) return false;
		if (!GetCenterOfGravity(s2, c2)) return false;
		if (!Math::Functs::AreAlmostEqual((c1 - c2).Length(), 0.0)) return false;

		return true;
	}

	static void GetEdgeEndPnt(const TopoDS_Edge& edge, gp_Pnt& st, gp_Pnt& nd)
	{
		int cnt = 0;
		for (TopExp_Explorer vertexExplorer(edge, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next(), ++cnt)
		{
			const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
			if (cnt == 0) st = BRep_Tool::Pnt(vertex);
			else nd = BRep_Tool::Pnt(vertex);
		}
	}

	static bool IsOnAndIn(TopAbs_State s)
	{
		return (s == TopAbs_IN || s == TopAbs_ON);
	}

	static bool GetFaceVector(const TopoDS_Face& aFace, std::vector<Vector3d>& poly)
	{
		if (!poly.empty()) poly.clear();
		Handle(Geom_Surface) aSurface = BRep_Tool::Surface(TopoDS::Face(aFace));
		std::vector<std::vector<Vector3d>> segments;
		//std::cerr << "one face://////////////" << std::endl;
		for (TopExp_Explorer edgeExp(aFace, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
		{
			Standard_Real FirstDummy, LastDummy;
			Handle(Geom_Curve) aCurve = BRep_Tool::Curve(TopoDS::Edge(edgeExp.Current()), FirstDummy, LastDummy);
			gp_Pnt _s = aCurve->Value(FirstDummy);
			gp_Pnt _e = aCurve->Value(LastDummy);
			std::vector<Vector3d> segment;
			segment.emplace_back(_s.X(), _s.Y(), _s.Z());
			segment.emplace_back(_e.X(), _e.Y(), _e.Z());
			segments.push_back(segment);
		}
		std::vector<std::vector<Vector3d>> lines;
		Math::Functs::Connecting_Segments(segments, lines);
		poly = lines[0];
		return true;
	}

	static bool GetFaceVector2d(const TopoDS_Face& aFace, std::vector<Vector2d>& poly)
	{
		if (!poly.empty()) poly.clear();
		std::vector<Vector3d> poly3D;
		GetFaceVector(aFace, poly3D);
		poly.reserve(poly3D.size());
		for (auto& p : poly3D)
		{
			poly.emplace_back(p[0], p[1]);
		}
		return true;
	}

	static TopAbs_State StatePointOnFace(const std::vector<Vector3d>& aFace, const gp_Pnt& aPoint, CGAL_3D_Distance_Point_Polygon* funDist)
	{
		// New Function
		/*
		TopoDS_Vertex aVertex = BRepBuilderAPI_MakeVertex(aPoint);
		BRepExtrema_DistShapeShape anExtrema(aFace, aVertex);
		const double dist = anExtrema.Value();
		if ((anExtrema.IsDone() == Standard_True) && (dist < 100 * Precision::Confusion()))
		{
			return TopAbs_IN;
			BRepClass_FaceClassifier classifier = BRepClass_FaceClassifier(aFace, aPoint, 10 * Precision::Confusion());
			return classifier.State();
		}
		else return TopAbs_OUT;
		*/
	};

	static gp_Pnt GetEdgeEndPnt(const TopoDS_Edge& edge, const gp_Dir& dir, int index)
	{
		if (index > 1) index = 1;
		gp_Pnt pts[2];
		double dist[2] = { 0, 0 };
		int cnt = 0;
		for (TopExp_Explorer vertexExplorer(edge, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next(), ++cnt)
		{
			const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
			pts[cnt] = BRep_Tool::Pnt(vertex);
			dist[cnt] = dir.XYZ().Dot(pts[cnt].XYZ());
		}
		if (dist[0] > dist[1]) std::swap(pts[0], pts[1]);
		return pts[index];
	}

	static TopoDS_Wire SortWireEdges(const TopoDS_Wire& wire)
	{
		ShapeFix_Wire fixWire;
		fixWire.Load(wire);
		fixWire.FixReorder();
		TopoDS_Wire result = fixWire.Wire();
		return result;
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
				throw Base::RuntimeError("Error in GeomFunc: vertices.Extent() != 2");

			const TopoDS_Vertex& vert = TopoDS::Vertex(vertices(2));
			points.emplace_back(BRep_Tool::Pnt(vert));
		}
		return points;
	}


	static gp_Dir GetEdgeDirection(const TopoDS_Edge& edge)
	{
		gp_Pnt pts[2];
		int cnt = 0;
		for (TopExp_Explorer vertexExplorer(edge, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next(), ++cnt)
		{
			const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
			pts[cnt] = BRep_Tool::Pnt(vertex);
		}
		return gp_Dir(pts[0].XYZ() - pts[1].XYZ());
	}

	static gp_Pnt GetCentralEdgePoint(const TopoDS_Edge& edge)
	{
		gp_Pnt pts[2];
		int cnt = 0;
		for (TopExp_Explorer vertexExplorer(edge, TopAbs_VERTEX); vertexExplorer.More(); vertexExplorer.Next(), ++cnt)
		{
			const TopoDS_Vertex& vertex = TopoDS::Vertex(vertexExplorer.Current());
			pts[cnt] = BRep_Tool::Pnt(vertex);
		}
		return gp_Pnt((pts[0].XYZ() + pts[1].XYZ()) / 2.0);
	}

	static gp_Pnt GetCentralFacePoint(const TopoDS_Face& face)
	{
		// compute point on face
		Standard_Real umin, umax, vmin, vmax;

		gp_Pnt p;

		Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
		BRepTools::UVBounds(face, umin, umax, vmin, vmax);
		Standard_Real umean = 0.5 * (umin + umax);
		Standard_Real vmean = 0.5 * (vmin + vmax);


		// compute intersection of u-iso line with face boundaries
		Handle(Geom2d_Curve) uiso = new Geom2d_Line(
			gp_Pnt2d(umean, 0.),
			gp_Dir2d(0., 1.)
			);

		TopExp_Explorer exp(face, TopAbs_EDGE);
		std::list<double> intersections;
		for (; exp.More(); exp.Next()) {
			TopoDS_Edge edge = TopoDS::Edge(exp.Current());
			Standard_Real first, last;

			// Get geometric curve from edge
			Handle(Geom2d_Curve) hcurve = BRep_Tool::CurveOnSurface(edge, face, first, last);
			hcurve = new Geom2d_TrimmedCurve(hcurve, first, last);

			Geom2dAPI_InterCurveCurve intersector(uiso, hcurve);
			for (int ipoint = 0; ipoint < intersector.NbPoints(); ++ipoint) {
				gp_Pnt2d p = intersector.Point(ipoint + 1);
				intersections.push_back(p.Y());
			}
		}

		// remove duplicate solutions defined by tolerance
		double tolerance = 1e-5;
		intersections.sort();
		intersections.unique(IsSame((vmax - vmin) * tolerance));

		// normally we should have at least two intersections
		// also the number of sections should be even - else something is really strange
		//assert(intersections.size() % 2 == 0);
		if (intersections.size() >= 2) {
			std::list<double>::iterator it = intersections.begin();
			double int1 = *it++;
			double int2 = *it;
			vmean = (int1 + int2) / 2.;
		}

		surface->D0(umean, vmean, p);

		return p;
	}

	static Standard_Real GetLength(const TopoDS_Edge& edge);
	static unsigned int GetNumberOfEdges(const TopoDS_Shape& shape);
	static unsigned int GetNumberOfFaces(const TopoDS_Shape& shape);
	static unsigned int GetNumberOfSubshapes(const TopoDS_Shape& shape);

	static bool GetCenterOfGravity(const TopoDS_Shape& _Shape, Base::Vector3d& center)
	{
		if (_Shape.IsNull())
			return false;

		// Computing of CentreOfMass
		gp_Pnt pnt;

		if (_Shape.ShapeType() == TopAbs_VERTEX) {
			pnt = BRep_Tool::Pnt(TopoDS::Vertex(_Shape));
		}
		else {
			GProp_GProps prop;
			if (_Shape.ShapeType() == TopAbs_EDGE || _Shape.ShapeType() == TopAbs_WIRE) {
				BRepGProp::LinearProperties(_Shape, prop);
			}
			else if (_Shape.ShapeType() == TopAbs_FACE || _Shape.ShapeType() == TopAbs_SHELL) {
				BRepGProp::SurfaceProperties(_Shape, prop);
			}
			else {
				BRepGProp::VolumeProperties(_Shape, prop);
			}

			pnt = prop.CentreOfMass();
		}

		center.Set(pnt.X(), pnt.Y(), pnt.Z());
		return true;
	}

	static gp_Trsf ConvertTogpTrsf(const glm::dmat4& mtrx)
	{
		gp_Trsf trsf;
		//		const auto mtrx = glm::transpose(m);
		// 		trsf.SetValues(mtrx[0][0], mtrx[0][1], mtrx[0][2], mtrx[0][3],
		// 			mtrx[1][0], mtrx[1][1], mtrx[1][2], mtrx[1][3],
		// 			mtrx[2][0], mtrx[2][1], mtrx[2][2], mtrx[2][3]
		// 		#if OCC_VERSION_HEX < 0x060800
		// 			, 0.00001, 0.00001
		// 		#endif
		//		); //precision was removed in OCCT CR0025194
		trsf.SetValues(mtrx[0][0], mtrx[1][0], mtrx[2][0], mtrx[3][0],
			mtrx[0][1], mtrx[1][1], mtrx[2][1], mtrx[3][1],
			mtrx[0][2], mtrx[1][2], mtrx[2][2], mtrx[3][2]
#if OCC_VERSION_HEX < 0x060800
			, 0.00001, 0.00001
#endif
			); //precision was removed in OCCT CR0025194
		return trsf;
	}


	template<typename _READER_>
	static TopoDS_Shape loadFile(
		const char* fileName)
	{
		TopoDS_Shape result;

		_READER_ reader;
		const int status = reader.ReadFile(const_cast<Standard_CString>(fileName));

		if (status == IFSelect_RetDone) {

			reader.NbRootsForTransfer();
			reader.TransferRoots();
			result = reader.OneShape();

		}
		return result;
	}

	static bool ExportToSTL(const TopoDS_Shape& shape, const char* string)
	{
		StlAPI_Writer STLwriter;
		return STLwriter.Write(shape, string);
	}

	static bool ExportToIGES(const TopoDS_Shape& shape, const char* string)
	{
		IGESControl_Controller::Init();
		IGESControl_Writer ICW("MM", 0);
		ICW.AddShape(shape);
		ICW.ComputeModel();
		Standard_Boolean OK = ICW.Write(string);
		return OK;
	}

	static bool ImportFromIGES(TopoDS_Shape& shape, const char* string)
	{
		shape = loadFile<IGESControl_Reader>(string);
		return true;
	}


	static TopoDS_Shape ImportFromIGES(const char* string)
	{
		return loadFile<IGESControl_Reader>(string);
	}


	static bool ImportFromBrep(TopoDS_Shape& shape, const char* file)
	{
		BRep_Builder brepb;
		BRepTools::Read(shape, file, brepb);
		return true;
	}

	static bool ImportFromStl(TopoDS_Shape& shape, const char* file)
	{
		StlAPI_Reader reader;
		return reader.Read(shape, file);

		//BRepBuilderAPI_MakeSolid make_solid(TopoDS::Shell(shape));
		//TopoDS_Solid topods_solid = make_solid.Solid();
	}

	static bool ExportToBrep(const TopoDS_Shape& shape, const char* file)
	{
		BRepTools::Write(shape, file);
		return true;
	}

	static void GetOutsideContour(const TopoDS_Face& surf, Vector3d& n_vec, Vector3d1& outside_contour)
	{
		Vector3d1 hole_vector;
		double hole_radius;
		bool circle_or_line = false;
		GetOutsideContour(surf, n_vec, outside_contour, hole_vector, hole_radius, circle_or_line);
	}
	//************************************
		// Method:    GetOutsideContour
		// Description: Get outside 3d contour of this surface
		// FullName:  PartDesignGui::PCEPart::GetOutsideContour
		// Access:    private 
		// Returns:   void
		// Parameter: const TopoDS_Face & surf
		// Parameter: Vector3d & n_vec
		// Parameter: Vector3d1 & outside_contour
		// Parameter: Vector3d1 & hole_vector
		// Parameter: double & hole_radius
		// Parameter: bool & circle_or_line
		//************************************
	static void GetOutsideContour(const TopoDS_Face& surf, Vector3d& n_vec, Vector3d1& outside_contour,
		Vector3d1& hole_vector, double& hole_radius, bool& circle_or_line)

	{
		auto CompareFunction = [](const std::pair<int, double>& cn0, const std::pair<int, double>& cn1)
		{
			return cn0.second > cn1.second;
		};

		//surface normal
		gp_Dir n = GeomFunc::GetFaceNormal(surf);
		n_vec = Vector3d(n.X(), n.Y(), n.Z());
		Math::Functs::SetVectorLength(n_vec, 1.0);

		//classification
		std::vector<std::pair<double, gp_Pnt>> circles;
		for (TopExp_Explorer edgeExp(surf, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
		{
			Standard_Real FirstDummy, LastDummy;
			Handle(Geom_Curve) curve = BRep_Tool::Curve(TopoDS::Edge(edgeExp.Current()), FirstDummy, LastDummy);

			if (curve->IsKind("Geom_Circle"))
			{
				// do something with circle
				Handle(Geom_Circle) aCircle = Handle(Geom_Circle)::DownCast(curve);
				double radius = aCircle->Radius();
				gp_Pnt center;
				aCircle->D0(0.0, center);
				circles.emplace_back(radius, center);
			}
		}
		/////////////////////////////
		circle_or_line = circles.size() == 2;

		if (circle_or_line)
		{
			if (!Math::Functs::IsAlmostZero(circles[0].first - circles[1].first))
			{
				std::cerr << "if (!Math::IsAlmostZero(circles[0].first - circles[1].first))" << std::endl;
				system("pause");
			}

			hole_vector.emplace_back(circles[0].second.X(), circles[0].second.Y(), circles[0].second.Z());
			hole_vector.emplace_back(circles[1].second.X(), circles[1].second.Y(), circles[1].second.Z());
			hole_radius = circles[0].first;
		}
		else
		{
			Vector3d2 contours_3d;
			std::vector<std::pair<int, double>> lens;
			for (TopExp_Explorer wireExp(surf, TopAbs_WIRE); wireExp.More(); wireExp.Next())
			{
				TopoDS_Wire  wire = TopoDS::Wire(wireExp.Current());

				//get points
				std::vector<gp_Pnt> vecPolygon = GeomFunc::GetSortedPointsOnWire(TopoDS::Wire(wireExp.Current()));
				if (vecPolygon.size() < 3)continue;
				Vector3d1 contour_3d;
				for (auto& p : vecPolygon) contour_3d.emplace_back(p.X(), p.Y(), p.Z());

				//check normal
				Vector3d c_n = Math::Functs::ComputeNormalFromPolyline(contour_3d);
				if (Math::Functs::GetAngleBetween(c_n, n_vec) < Math::Math_PI / 2.0)
					std::reverse(contour_3d.begin(), contour_3d.end());

				//3d => 2d
				double len = Math::Functs::GetLength(contour_3d);
				if (len > 0)
				{
					lens.emplace_back(contours_3d.size(), len);
					contours_3d.emplace_back(contour_3d);
				}
			}
			if (contours_3d.size() > 0 && lens.size() > 0)
			{
				std::sort(lens.begin(), lens.end(), CompareFunction);
				outside_contour = contours_3d[lens[0].first];
			}
		}
	}


	static void LoadTopoDSGeometry(const TopoDS_Shape& ashape_,
		Vector3d2& surfs, Vector3d1& surf_normals, double& volume,
		std::vector<std::pair<Vector3d1, double>>& holes,
		Vector3d& center)
	{
		//volume
		GProp_GProps gprops;
		BRepGProp::VolumeProperties(ashape_, gprops);
		volume = gprops.Mass();

		int nb = 0;
		for (TopExp_Explorer exp(ashape_, TopAbs_FACE); exp.More(); exp.Next())
		{
			Vector3d n_vec;
			Vector3d1 outside_contour, hole_vector;
			double hole_radius;
			bool circle_or_line = false;
			GetOutsideContour(TopoDS::Face(exp.Current()), n_vec, outside_contour, hole_vector, hole_radius, circle_or_line);

			if (circle_or_line)
				holes.emplace_back(hole_vector, hole_radius);
			else
			{
				//line
				if (outside_contour.size() >= 3)
				{
					surfs.emplace_back(outside_contour);
					surf_normals.emplace_back(n_vec);
				}
			}
			nb++;
		}
		if (surfs.size() < 4)
		{
			std::cerr << "Input surface's number is less than 4..." << std::endl;
			system("pause");
			return;
		}
		center = Math::Functs::GetCenter(surfs);
	}

	static Vector3d2 LoadTopoDSGeometry(const TopoDS_Shape& ashape_, const bool check_4 = true)
	{
		Vector3d2 surfs;

		int nb = 0;
		for (TopExp_Explorer exp(ashape_, TopAbs_FACE); exp.More(); exp.Next())
		{
			Vector3d n_vec;
			Vector3d1 outside_contour, hole_vector;
			double hole_radius;
			bool circle_or_line = false;
			GetOutsideContour(TopoDS::Face(exp.Current()), n_vec, outside_contour, hole_vector, hole_radius, circle_or_line);

			if (!circle_or_line)
			{
				//line
				if (outside_contour.size() >= 3)
				{
					surfs.emplace_back(outside_contour);
				}
			}
			nb++;
		}
		if (surfs.size() < 4 && check_4)
		{
			std::cerr << "Input surface's number is less than 4..." << std::endl;
			system("pause");
		}

		return surfs;
	}

	static Vector3d2 LoadTopoDSGeometry(const TopoDS_Face& aface_)
	{
		Vector3d2 surfs;

		Vector3d n_vec;
		Vector3d1 outside_contour, hole_vector;
		double hole_radius;
		bool circle_or_line = false;
		GetOutsideContour(aface_, n_vec, outside_contour, hole_vector, hole_radius, circle_or_line);

		if (!circle_or_line)
		{
			//line
			if (outside_contour.size() >= 3)
			{
				surfs.emplace_back(outside_contour);
			}
		}
		return surfs;
	}


	static Standard_Real ProjectPointOnWire(const TopoDS_Wire& wire, gp_Pnt p)
	{
		double smallestDist = DBL_MAX;
		double alpha = 0.;
		int edgeIndex = 0;

		// find edge with closest dist to point p
		BRepTools_WireExplorer wireExplorer;
		int iwire = 0;
		for (wireExplorer.Init(wire); wireExplorer.More(); wireExplorer.Next(), iwire++) {
			Standard_Real firstParam, lastParam;
			TopoDS_Edge edge = wireExplorer.Current();
			Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, firstParam, lastParam);

			GeomAPI_ProjectPointOnCurve proj(p, curve, firstParam, lastParam);
			if (proj.NbPoints() > 0 && proj.LowerDistance() < smallestDist) {
				smallestDist = proj.LowerDistance();
				edgeIndex = iwire;
				alpha = proj.LowerDistanceParameter();
			}
		}

		// compute partial length of wire until projection point is reached
		wireExplorer.Init(wire);
		double partLength = 0.;
		for (int iwire = 0; iwire <= edgeIndex; ++iwire) {
			Standard_Real firstParam;
			Standard_Real lastParam;
			TopoDS_Edge edge = wireExplorer.Current();
			Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, firstParam, lastParam);
			GeomAdaptor_Curve adaptorCurve(curve, firstParam, lastParam);
			if (iwire == edgeIndex) {
				auto projPnt = curve->Value(alpha);
				return (projPnt.XYZ() - p.XYZ()).Modulus();
			}
			wireExplorer.Next();
		}

		return DBL_MAX;
	}

	static bool IsPointInsideShape(const TopoDS_Shape& solid, gp_Pnt point)
	{
		// check if solid
		TopoDS_Solid s;
		try {
			s = TopoDS::Solid(solid);
		}
		catch (Standard_Failure) {
			throw std::exception("The shape is not a solid");
		}

		BRepClass3d_SolidClassifier algo(s);

		// test whether a point at infinity lies inside. If yes, then the shape is reversed
		algo.PerformInfinitePoint(1e-3);

		bool shapeIsReversed = (algo.State() == TopAbs_IN);

		algo.Perform(point, 1e-3);

		return ((algo.State() == TopAbs_IN) != shapeIsReversed) || (algo.State() == TopAbs_ON);
	}

	static void DecomposeContinuity(const TopoDS_Wire& inputWire)
	{


	}

	enum CutType
	{
		BandSawForward,
		BandSawReversed,
		OnlyJigSaw,
		ChopOrTrackSaw
	};

	static void EnumerateCuts(const std::vector<std::pair<bool, bool>>& cuts,
		std::vector<std::vector<std::pair<int, CutType>>>& candidates,
		std::vector<std::pair<int, CutType>> currentCandidate,
		int i)
	{
		if (i == cuts.size())
		{
			candidates.push_back(currentCandidate);
			return;
		}

		if (cuts[i].first)
		{
			std::vector<std::pair<int, CutType>> posCan = currentCandidate;
			posCan.emplace_back(i, BandSawForward);
			EnumerateCuts(cuts, candidates, posCan, i + 1);
		}

		if (cuts[i].second)
		{
			std::vector<std::pair<int, CutType>> negCan = currentCandidate;
			negCan.emplace_back(i, BandSawReversed);
			EnumerateCuts(cuts, candidates, negCan, i + 1);
		}

		if (!cuts[i].first && !cuts[i].second)
		{
			std::vector<std::pair<int, CutType>> jigCan = currentCandidate;
			jigCan.emplace_back(i, OnlyJigSaw);
			EnumerateCuts(cuts, candidates, jigCan, i + 1);
		}

		if (cuts[i].first && cuts[i].second)
		{
			std::vector<std::pair<int, CutType>> chopCan = currentCandidate;
			chopCan.emplace_back(i, ChopOrTrackSaw);
			EnumerateCuts(cuts, candidates, chopCan, i + 1);
		}
	}

	static GeomAbs_Shape GetEdgeContinuity(const TopoDS_Edge& edge1, const TopoDS_Edge& edge2)
	{
		// **********************************************************************************
		// Create a wire from both edges
		// **********************************************************************************
		BRepBuilderAPI_MakeWire wireMaker;
		wireMaker.Add(edge1);
		if (!wireMaker.IsDone()) {
			throw Base::RuntimeError("GetEdgeContinuity: Unable to create wire for edge1");
		}
		wireMaker.Add(edge2);
		if (!wireMaker.IsDone()) {
			throw Base::RuntimeError("GetEdgeContinuity: Unable to create wire for edge2");
		}
		TopoDS_Vertex commonVertex;
		TopoDS_Wire wire = wireMaker.Wire();

		// **********************************************************************************
		// Fix connectivity: Create common vertex
		// **********************************************************************************
		ShapeFix_Wire wireFixer;
		wireFixer.Load(wire);
		wireFixer.FixConnected();
		wireFixer.Perform();
		wire = wireFixer.Wire();
		if (wire.IsNull()) {
			throw Base::RuntimeError("GetEdgeContinuity: Unable to create common vertex in connected wire");
		}

		// **********************************************************************************
		// Get new edges with common vertex
		// **********************************************************************************
		TopTools_IndexedMapOfShape edges;
		TopExp::MapShapes(wire, TopAbs_EDGE, edges);
		if (edges.Extent() != 2) {
			throw Base::RuntimeError("checkEdgeContinuity: Unexpected error in connected wire");
		}
		TopoDS_Edge newedge1 = TopoDS::Edge(edges(1));
		TopoDS_Edge newedge2 = TopoDS::Edge(edges(2));

		TopExp::CommonVertex(newedge1, newedge2, commonVertex);
		if (commonVertex.IsNull()) {
			throw Base::RuntimeError("checkEdgeContinuity: Unable to find common vertex");
		}

		// **********************************************************************************
		// Get derivatives at common vertex
		// **********************************************************************************
		BRepAdaptor_Curve firstCurve;
		BRepAdaptor_Curve secondCurve;
		if (TopExp::FirstVertex(newedge1).IsSame(commonVertex) && TopExp::LastVertex(newedge2).IsSame(commonVertex)) {
			secondCurve.Initialize(newedge1);
			firstCurve.Initialize(newedge2);
		}
		else if (TopExp::LastVertex(newedge1).IsSame(commonVertex) && TopExp::FirstVertex(newedge2).IsSame(commonVertex)) {
			firstCurve.Initialize(newedge1);
			secondCurve.Initialize(newedge2);
		}
		else {
			throw Base::RuntimeError("checkEdgeContinuity: Unexpected error: Unable to find common vertex in connected wire");
		}

		gp_Pnt point;
		gp_Vec firstD1;
		gp_Vec firstD2;
		gp_Vec secondD1;
		gp_Vec secondD2;
		firstCurve.D2(firstCurve.LastParameter(), point, firstD1, firstD2);
		secondCurve.D2(secondCurve.FirstParameter(), point, secondD1, secondD2);

		// **********************************************************************************
		// Get continuity depending on derivatives at the common vertex
		// **********************************************************************************
		if ((firstD1.Magnitude() >= gp::Resolution() && secondD1.Magnitude() >= gp::Resolution())
			&& gp_Dir(firstD1).IsEqual(gp_Dir(secondD1), Precision::Confusion())) {
			if ((firstD2.Magnitude() >= gp::Resolution() && secondD2.Magnitude() >= gp::Resolution())
				&& gp_Dir(firstD2).IsEqual(gp_Dir(secondD2), Precision::Confusion())) {
				return GeomAbs_Shape::GeomAbs_C2;
			}
			else {
				return GeomAbs_Shape::GeomAbs_C1;
			}
		}
		else {
			return GeomAbs_Shape::GeomAbs_C0;
		}
	}

	// The following codes are borrower ed from Haisen

	static std::string IntString(int i)
	{
		std::stringstream ss;
		std::string str;
		ss << i;
		ss >> str;
		return str;
	}
	static std::string IntString(std::vector<int> vecs)
	{
		std::string str;
		for (int i = 0; i < vecs.size(); i++)
		{
			str += IntString(vecs[i]);
		}
		return str;
	}



	static std::string DoubleString(double d)
	{
		std::stringstream ss;
		ss.precision(5);
		std::string str;
		ss << d;
		ss >> str;
		return str;
	}

	template <class Type>
	Type StringToNum(const string& str)
	{
		istringstream iss(str);
		Type num;
		iss >> num;
		return num;
	}

	static double GetLength(Vector3d v) {
		return glm::length(v);
	}

	static double GetLength(Vector2d v) {
		return glm::length(v);
	}

	static double Radian2Angle(double radian)
	{
		return radian / Math_PI * 180.0;
	}

	static double Angle2Radian(double angle)
	{
		return angle / 180.0 * Math_PI;
	}

	static glm::dmat4 RotationMatrix(const Vector3d& n, double angle)
	{
		double u = n[0];
		double v = n[1];
		double w = n[2];

		glm::dmat4  rotationMatrix;

		double L = (u * u + v * v + w * w);

		//angle = angle * M_PI / 180.0; //converting to radian value	
		double u2 = u * u;
		double v2 = v * v;
		double w2 = w * w;

		rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
		rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][3] = 0.0;

		rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
		rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][3] = 0.0;

		rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
		rotationMatrix[2][3] = 0.0;

		rotationMatrix[3][0] = 0.0;
		rotationMatrix[3][1] = 0.0;
		rotationMatrix[3][2] = 0.0;
		rotationMatrix[3][3] = 1.0;

		return rotationMatrix;
	}

	static glm::dmat4 TranslationMatrix(const Vector3d& v)
	{
		glm::dmat4  translationMatrix;
		translationMatrix[0][0] = 1.0;
		translationMatrix[0][1] = 0.0;
		translationMatrix[0][2] = 0.0;
		translationMatrix[0][3] = 0.0;

		translationMatrix[1][0] = 0.0;
		translationMatrix[1][1] = 1.0;
		translationMatrix[1][2] = 0.0;
		translationMatrix[1][3] = 0.0;

		translationMatrix[2][0] = 0.0;
		translationMatrix[2][1] = 0.0;
		translationMatrix[2][2] = 1.0;
		translationMatrix[2][3] = 0.0;

		translationMatrix[3][0] = v[0];
		translationMatrix[3][1] = v[1];
		translationMatrix[3][2] = v[2];
		translationMatrix[3][3] = 1.0;

		return translationMatrix;
	}

	static Vector3d RotationAxis(Vector3d p, double angle, Vector3d n)
	{
		glm::dmat4 inputMatrix(0.0);
		inputMatrix[0][0] = p[0];
		inputMatrix[1][0] = p[1];
		inputMatrix[2][0] = p[2];
		inputMatrix[3][0] = 1.0;
		double u = n[0];
		double v = n[1];
		double w = n[2];

		glm::dmat4  rotationMatrix;

		double L = (u * u + v * v + w * w);

		//angle = angle * M_PI / 180.0; //converting to radian value
		double u2 = u * u;
		double v2 = v * v;
		double w2 = w * w;

		rotationMatrix[0][0] = (u2 + (v2 + w2) * cos(angle)) / L;
		rotationMatrix[0][1] = (u * v * (1 - cos(angle)) - w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][2] = (u * w * (1 - cos(angle)) + v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[0][3] = 0.0;

		rotationMatrix[1][0] = (u * v * (1 - cos(angle)) + w * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][1] = (v2 + (u2 + w2) * cos(angle)) / L;
		rotationMatrix[1][2] = (v * w * (1 - cos(angle)) - u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[1][3] = 0.0;

		rotationMatrix[2][0] = (u * w * (1 - cos(angle)) - v * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][1] = (v * w * (1 - cos(angle)) + u * sqrt(L) * sin(angle)) / L;
		rotationMatrix[2][2] = (w2 + (u2 + v2) * cos(angle)) / L;
		rotationMatrix[2][3] = 0.0;

		rotationMatrix[3][0] = 0.0;
		rotationMatrix[3][1] = 0.0;
		rotationMatrix[3][2] = 0.0;
		rotationMatrix[3][3] = 1.0;

		double outputMatrix[4][1] = { 0.0, 0.0, 0.0, 0.0 };

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 1; j++) {
				outputMatrix[i][j] = 0;
				for (int k = 0; k < 4; k++) {
					outputMatrix[i][j] += rotationMatrix[i][k] * inputMatrix[k][j];
				}
			}
		}
		return Vector3d(outputMatrix[0][0], outputMatrix[0][1], outputMatrix[0][2]);
	}

	static Vector2d RotationAxis2d(Vector2d p, double angle, Vector2d center)
	{
		Vector3d r = RotationAxis(Vector3d(p[0] - center[0], 0.0, p[1] - center[1]),
			angle, Vector3d(0.0, 1.0, 0.0)) + Vector3d(center[0], 0.0, center[1]);
		return Vector2d(r[0], r[2]);
	}

	static Vector3d GetCrossproduct(Vector3d& v1, Vector3d& v2) {
		return glm::cross(v1, v2);
	}

	static bool AreAlmostEqual(double value1, double value2) {
		if (value1 == value2) {
			return true;
		}
		double eps = (glm::abs(value1) + glm::abs(value2) + 10.0) * DOUBLE_EPSILON;
		double delta = value1 - value2;
		return (-eps < delta) && (eps > delta);

	}

	static Vector3d SetVectorLength(Vector3d& v, double length)
	{
		double l = GetLength(v);

		v[0] = v[0] / l * length;
		v[1] = v[1] / l * length;
		v[2] = v[2] / l * length;

		return v;
	}

	static bool IsAlmostZero(double value) {
		return (value < 10.0 * DOUBLE_EPSILON) && (value > -10.0 * DOUBLE_EPSILON);
	}

	static bool IsAlmostZero_Double(double value, double EPSILON) {
		return (value < 10.0 * DOUBLE_EPSILON) && (value > -10.0 * EPSILON);
	}

	static double GetAngleBetween(Vector3d v1, Vector3d v2) {

		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));

		if (IsAlmostZero(d - 1.0))
			return 0.0;

		if (IsAlmostZero(d + 1.0))
			return Math_PI;

		return glm::acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
	}

	static double GetAngleBetween(Vector2d v1, Vector2d v2) {
		double d = glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2));
		if (IsAlmostZero(d - 1.0))
			return 0.0;
		if (IsAlmostZero(d + 1.0))
			return Math_PI;
		return glm::acos(glm::dot(v1, v2) / (glm::length(v1) * glm::length(v2)));
	}

	static bool VectorContain(std::vector<int> vecs, int element)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			if (vecs[i] == element)
				return true;
		}
		return false;
	}

	//existing bugs in this function
	static Vector3d Vector3dBase(Vector3d v)
	{
		Vector3d n(1.0, 1.0, 1.0);
		if (!IsAlmostZero(v[0]))
		{
			n[0] = -(v[1] + v[2]) / v[0];
			return n;
		}
		else
		{
			if (!IsAlmostZero(v[1]))
			{
				n[1] = -(v[0] + v[2]) / v[1];
				return n;
			}
			else
			{
				if (!IsAlmostZero(v[2]))
				{
					n[2] = -(v[0] + v[1]) / v[2];
					return n;
				}
				return n;
			}
		}
	}

	static Vector3d GetCenter(const std::vector<Vector3d>& points)
	{
		Vector3d center(0.0, 0.0, 0.0);
		for (int i = 0; i < points.size(); i++)
			center += points[i];
		center = center / (double)points.size();
		return center;
	}
	static Vector3d GetCenter(const std::vector<std::vector<Vector3d>>& points)
	{
		Vector3d center(0.0, 0.0, 0.0);
		int nb = 0;
		for (int i = 0; i < points.size(); i++)
		{
			for (int j = 0; j < points[i].size(); j++)
			{
				center += points[i][j];
			}
			nb += points[i].size();
		}
		center = center / (double)nb;
		return center;
	}

	static void GetBoundingBox(const std::vector<std::vector<Vector3d>>& points, Vector3d& minimal_corner, Vector3d& maximal_corner)
	{
		minimal_corner = Vector3d(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
		maximal_corner = Vector3d(-MAXDOUBLE, -MAXDOUBLE, -MAXDOUBLE);
		for (int i = 0; i < points.size(); i++)
		{
			for (int j = 0; j < points[i].size(); j++)
			{
				minimal_corner[0] = min(minimal_corner[0], points[i][j][0]);
				minimal_corner[1] = min(minimal_corner[1], points[i][j][1]);
				minimal_corner[2] = min(minimal_corner[2], points[i][j][2]);
				maximal_corner[0] = max(maximal_corner[0], points[i][j][0]);
				maximal_corner[1] = max(maximal_corner[1], points[i][j][1]);
				maximal_corner[2] = max(maximal_corner[2], points[i][j][2]);
			}
		}
	}

	static double CircumCircleRaidius(Vector2d v0, Vector2d v1, Vector2d v2)
	{
		double a = GetLength(v0 - v1);
		double b = GetLength(v0 - v2);
		double c = GetLength(v1 - v2);
		double p = (a + b + c) / 2.0;
		double area = (4.0 * pow(p * (p - a) * (p - b) * (p - c), 0.5));
		double radius;

		if (IsAlmostZero(area))
		{
			double max_l = a;
			max_l = max(max_l, b);
			max_l = max(max_l, c);
			radius = 10 * max_l;
		}
		else
			radius = a * b * c / area;
		return radius;
	}



	static double GetLength(Vector3d v0, Vector3d v1) {
		return GetLength(v0 - v1);
	}

	static double GetDistance(Vector3d v0, Vector3d v1) {
		return GetLength(v0 - v1);
	}

	static Vector3d PlaneProject(Vector3d planar_location, Vector3d planar_direction, Vector3d p)
	{
		if (IsAlmostZero(GetLength(planar_location, p)))
			return planar_location;

		double angle = GetAngleBetween(planar_direction, p - planar_location);
		double length = GetLength(planar_location, p);

		if (angle <= Math_PI / 2.0)
			return p - SetVectorLength(planar_direction, length * sin(Math_PI / 2.0 - angle));
		else
			return p + SetVectorLength(planar_direction, length * sin(angle - Math_PI / 2.0));

	}

	static Vector3d RotationAxis(Vector3d p, double angle, Vector3d ray_point, Vector3d ray_vector)
	{
		return RotationAxis(p - ray_point, angle, ray_vector) + ray_point;
	}

	static void RotationAxis(std::vector<Vector3d>& points, double angle, Vector3d ray_point, Vector3d ray_vector)
	{
		for (int i = 0; i < points.size(); i++)
			points[i] = RotationAxis(points[i], angle, ray_point, ray_vector);
	}

	static void RotationAxis(std::vector<std::vector<Vector3d>>& points, double angle, Vector3d ray_point, Vector3d ray_vector)
	{
		for (int i = 0; i < points.size(); i++)
			RotationAxis(points[i], angle, ray_point, ray_vector);
	}
	static void RotationAxis(std::vector<std::vector<std::vector<Vector3d>>>& pointses, double angle, Vector3d ray_point, Vector3d ray_vector)
	{
		for (int i = 0; i < pointses.size(); i++)
			RotationAxis(pointses[i], angle, ray_point, ray_vector);
	}

	static Vector3d Translate(Vector3d p, Vector3d v)
	{
		return p + v;
	}


	static void Connecting_Segments(std::vector<std::vector<Vector3d>>& segments, std::vector<std::vector<Vector3d>>& lines)
	{
		//save connecting relations
		std::vector<bool> used(segments.size(), false);
		std::vector<int> relations;
#pragma region get_relations
		for (int i = 0; i < segments.size(); i++)
		{
			for (int j = 0; j < segments.size(); j++)
			{
				if (i != j && !used[i] && !used[j])
				{
					bool b_0_0 = IsAlmostZero_Double(Math::Functs::GetLength(segments[i][0], segments[j][0]), DOUBLE_EPSILON);
					bool b_0_1 = IsAlmostZero_Double(Math::Functs::GetLength(segments[i][0], segments[j][1]), DOUBLE_EPSILON);
					bool b_1_0 = IsAlmostZero_Double(Math::Functs::GetLength(segments[i][1], segments[j][0]), DOUBLE_EPSILON);
					bool b_1_1 = IsAlmostZero_Double(Math::Functs::GetLength(segments[i][1], segments[j][1]), DOUBLE_EPSILON);

					if ((b_0_0 && b_1_1) || (b_0_1 && b_1_0))
					{
						used[j] = true;
						continue;
					}

					if (b_0_0)
					{
						relations.push_back(i);
						relations.push_back(0);
						relations.push_back(j);
						relations.push_back(0);
						continue;
					}
					if (b_0_1)
					{
						relations.push_back(i);
						relations.push_back(0);
						relations.push_back(j);
						relations.push_back(1);
						continue;
					}
					if (b_1_0)
					{
						relations.push_back(i);
						relations.push_back(1);
						relations.push_back(j);
						relations.push_back(0);
						continue;
					}
					if (b_1_1)
					{
						relations.push_back(i);
						relations.push_back(1);
						relations.push_back(j);
						relations.push_back(1);
						continue;
					}
				}
			}
		}
#pragma endregion

		std::vector<std::vector<int>> ones;


		while (true)
		{
			int index = -1;
			int end = -1;

			for (int i = 0; i < segments.size(); i++)
			{
				if (!used[i]) {
					index = i;
					end = 0;
					used[i] = true;
					break;
				}
			}

			if (index < 0)break;

			std::vector<Vector3d> line(1, segments[index][end]);

			std::vector<int> one(1, index);

			while (true)
			{
				end = 1 - end;
				bool search = false;
				for (int i = 0; i < relations.size(); i = i + 4)
				{
					if (relations[i] == index && relations[i + 1] == end && !used[relations[i + 2]])
					{
						line.push_back(segments[relations[i + 2]][relations[i + 3]]);
						one.push_back(relations[i + 2]);
						index = relations[i + 2];
						end = relations[i + 3];
						used[index] = true;
						search = true;
						break;
					}
					if (relations[i + 2] == index && relations[i + 3] == end && !used[relations[i]])
					{
						line.push_back(segments[relations[i]][relations[i + 1]]);
						one.push_back(relations[i]);
						index = relations[i];
						end = relations[i + 1];
						used[index] = true;
						search = true;
						break;
					}
				}
				if (!search) { break; }
			}

			ones.push_back(one);
			lines.push_back(line);
		}
	}



	static void Translate(std::vector<Vector3d>& points, Vector3d v)
	{
		for (int i = 0; i < points.size(); i++)
			points[i] = Translate(points[i], v);
	}

	static void Translate(std::vector<std::vector<Vector3d>>& points, Vector3d v)
	{
		for (int i = 0; i < points.size(); i++)
			Translate(points[i], v);
	}

	static Vector3d IntersectPointPlane2Ray(Vector3d planar_location, Vector3d planar_direction, Vector3d ray_location, Vector3d ray_vector)
	{
		Vector3d project_point = PlaneProject(planar_location, planar_direction, ray_location);
		double distance = GetDistance(ray_location, project_point);
		if (IsAlmostZero(GetLength(project_point, ray_location)))
			return ray_location;
		double angle = GetAngleBetween(ray_vector, project_point - ray_location);
		double length = distance / cos(angle);
		return ray_location + SetVectorLength(ray_vector, length);
	}
	static Vector3d ComputeNormalFromPolyline(const std::vector<Vector3d>& points)
	{
		Vector3d planar_direction;
		planar_direction = GetCrossproduct(points[0] - points[1], points[2] - points[1]);
		SetVectorLength(planar_direction, 1.0);
		return planar_direction;
	}

	static void  ComputePlanarFromPolyline(Vector3d& planar_location, Vector3d& planar_direction, const std::vector<Vector3d>& points)
	{
		planar_location = points[0];
		planar_direction = GetCrossproduct(points[0] - points[1], points[2] - points[1]);
		SetVectorLength(planar_direction, 1.0);
	}

	static TopoDS_Face Poly2TopoDSFace(const HMODULE& hModule, const Vector3d1& poly)
	{
		auto center = Math::Functs::GetCenter(poly);
		auto normal = Math::Functs::ComputeNormalFromPolyline(poly);
		//CGAL_3D_Plane_Points_Projection
		const auto projection = (CGAL_3D_Plane_Points_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Points_Projection");
		Vector3d1 ps0 = poly;
		Vector3d1 ps;
		projection(center, normal, ps0, ps);

		BRepBuilderAPI_MakePolygon poly_Wire;
		for (int j = 0; j < ps.size(); j++)
		{
			gp_Pnt pnt(ps[j][0], ps[j][1], ps[j][2]);
			poly_Wire.Add(pnt);
		}
		poly_Wire.Close();

		if (!BRep_Tool::IsClosed(poly_Wire))
		{
			int dasd = 0.0;
		}

		TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(poly_Wire, Standard_True);

		return myFaceProfile;
	}

	static TopoDS_Solid Polys2TopoSolid(const HMODULE& hModule, const Vector3d2& polys)
	{
		BRepOffsetAPI_Sewing sew;
		sew.Init(1.0e-01);
		for (int i = 0; i < polys.size(); i++)
		{
			/////////////////
			auto center = Math::Functs::GetCenter(polys[i]);
			auto normal = Math::Functs::ComputeNormalFromPolyline(polys[i]);
			//CGAL_3D_Plane_Points_Projection
			const auto projection = (CGAL_3D_Plane_Points_Projection)GetProcAddress(hModule, "CGAL_3D_Plane_Points_Projection");
			Vector3d1 ps0 = polys[i];
			Vector3d1 ps;
			projection(center, normal, ps0, ps);

			BRepBuilderAPI_MakePolygon poly_Wire;
			for (int j = 0; j < ps.size(); j++)
			{
				gp_Pnt pnt(ps[j][0], ps[j][1], ps[j][2]);
				poly_Wire.Add(pnt);
			}
			poly_Wire.Close();

			if (!BRep_Tool::IsClosed(poly_Wire))
			{
				int dasd = 0.0;
			}

			TopoDS_Face myFaceProfile = BRepBuilderAPI_MakeFace(poly_Wire, Standard_True);
			sew.Add(myFaceProfile);
		}
		sew.Perform();
		TopoDS_Shape aShape = sew.SewedShape();
		TopoDS_Shell aShell = TopoDS::Shell(aShape);
		if (!BRep_Tool::IsClosed(aShell))
		{
			int dasd = 0.0;
		}
		BRepBuilderAPI_MakeSolid make_solid(TopoDS::Shell(aShape));
		return make_solid.Solid();
	};

	static std::vector<TopoDS_Shape> CEOSplit(TopoDS_Shape& prism_0, TopoDS_Shape& prism_1)
	{
		std::vector<TopoDS_Shape> splitting_parts;
		GEOMAlgo_Splitter splitterC;
		splitterC.SetFuzzyValue(0.2);
		splitterC.AddArgument(prism_0);
		splitterC.AddTool(prism_1);
		splitterC.Perform();

		TopTools_ListIteratorOfListOfShape iterC(splitterC.Modified(prism_0));
		int ncParts = 0;
		TopoDS_Shape firstComponent;
		for (; iterC.More(); iterC.Next(), ++ncParts)
		{
			BRepClass3d_SolidClassifier classifier;
			classifier.Load(prism_1);
			Base::Vector3d center;
			GeomFunc::GetCenterOfGravity(iterC.Value(), center);
			classifier.Perform(gp_Pnt(center.x, center.y, center.z), Precision::Confusion());
			if (!(classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON))
			{
				firstComponent = iterC.Value();
				auto dsad = iterC.Value().ShapeType();
				splitting_parts.emplace_back(iterC.Value());
			}
		}

		return splitting_parts;
	};



	static bool CollisionDetection(const HMODULE& hModule, const Vector3d2& polys_1, const Vector3d2& polys_2)
	{
		//bounding box detection
		Vector3d min_corner_1, max_corner_1, min_corner_2, max_corner_2;
		Math::Functs::GetBoundingBox(polys_1, min_corner_1, max_corner_1);
		Math::Functs::GetBoundingBox(polys_2, min_corner_2, max_corner_2);

		if (min_corner_1[0] > max_corner_2[0])return false;
		if (min_corner_2[0] > max_corner_1[0])return false;
		if (min_corner_1[1] > max_corner_2[1])return false;
		if (min_corner_2[1] > max_corner_1[1])return false;
		if (min_corner_1[2] > max_corner_2[2])return false;
		if (min_corner_2[2] > max_corner_1[2])return false;

		//solid 
		auto solid1 = GeomFunc::Polys2TopoSolid(hModule, polys_1);
		auto solid2 = GeomFunc::Polys2TopoSolid(hModule, polys_2);

		const Standard_Real aLinearDeflection = 0.01;
		const Standard_Real anAngularDeflection = 0.01;

		BRepMesh_IncrementalMesh aMesh1(solid1, aLinearDeflection, Standard_False, anAngularDeflection);
		BRepMesh_IncrementalMesh aMesh2(solid2, aLinearDeflection, Standard_False, anAngularDeflection);

		TopoDS_Compound compound;
		BRep_Builder aBuilder;
		aBuilder.MakeCompound(compound);
		aBuilder.Add(compound, solid1);

		TopoDS_Compound compound2;
		BRep_Builder aBuilder2;
		aBuilder2.MakeCompound(compound2);
		aBuilder2.Add(compound2, solid2);

		Standard_Boolean done = 0;

		//Perform proximity test
		BRepExtrema_ShapeProximity proximity(compound, compound2, 0.001);

		proximity.Perform();

		//Initializing compound for store collision faces
		TopoDS_Builder aCompBuilder;
		TopoDS_Compound CollisionFaceCompound;
		aCompBuilder.MakeCompound(CollisionFaceCompound);

		// if proximity test is successful 
		if (proximity.IsDone())
		{
			// gives collision collision face and edges of shape1 
			for (BRepExtrema_MapOfIntegerPackedMapOfInteger::Iterator anIt1(proximity.OverlapSubShapes1()); anIt1.More(); anIt1.Next())
			{
				const TopoDS_Face& aFace1 = proximity.GetSubShape1(anIt1.Key());
				TopoDS_Shape F1 = aFace1;

				// gives collision collision face and as well as edges of shape2, so we use BRepAlgoAPI_Common to //take common faces alone
				for (BRepExtrema_MapOfIntegerPackedMapOfInteger::Iterator anIt2(proximity.OverlapSubShapes2()); anIt2.More(); anIt2.Next())
				{
					const TopoDS_Face& aFace2 = proximity.GetSubShape2(anIt2.Key());

					/// int this step we take 1st shape's collision face is in "aFace1" and compare with all the collision ///faces of second shape
					//and BRepAlgoAPI_Common return exact collision faces 
					TopoDS_Shape S1 = BRepAlgoAPI_Common(aFace1, aFace2);
					aCompBuilder.Add(CollisionFaceCompound, S1);

					return true;
				}
			}
		}
		return false;
	};

};

#endif


