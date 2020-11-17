#include "PHELM.h"
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepExtrema_ExtPC.hxx>
#include <Mod/PartDesign/Gui/Math.hpp>
#include <Mod/Part/App/TopoShape.h>
#include <Geom2dAPI_ProjectPointOnCurve.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

#include "CompilerConfig.h"


void PHELM::ClearAll()
{
	std::vector<std::string>().swap(CADCode);
	std::vector<std::string>().swap(HELMCode);
	std::vector<std::pair<double, double>>().swap(vecDrillTools);
	std::vector<std::tuple<std::list<std::string>, int>>().swap(prog_strs);
}


void PHELM::FaceToVertexList(const TopoDS_Face& face, const gp_Dir& dir, TopoDS_Wire& wire,
	std::vector<gp_Pnt>& vecPolygon, std::vector<Vector2d>& vecPolygon2D)
{
	const auto cgal_3d_2d_point = (CGAL_3D_Plane_3D_to_2D_Point)GetProcAddress(*hModule, "CGAL_3D_Plane_3D_to_2D_Point");
	Vector3d planeNormal(dir.X(), dir.Y(), dir.Z());
	Vector3d planePoint(0.0, 0.0, 0.0);

	for (TopExp_Explorer wireExp(face, TopAbs_WIRE); wireExp.More(); wireExp.Next())
	{
		wire = TopoDS::Wire(wireExp.Current());
		vecPolygon = GeomFunc::GetSortedPointsOnWire(wire);
		for (auto& p : vecPolygon)
		{
			Vector2d projPnt;
			Vector3d oriPnt(p.X(), p.Y(), p.Z());
			cgal_3d_2d_point(planePoint, planeNormal, oriPnt, projPnt);
			vecPolygon2D.push_back(projPnt);
		}
	}
}


std::string PHELM::GetEdgeReference(TopoDS_Face& face, const std::vector<Vector2d>& facePoly, const Vector2d& s)
{
	Geom2dAPI_ProjectPointOnCurve ProjectOnCurve;
	std::vector<RefPnt> vecRefPnts;

	TopTools_IndexedMapOfShape originMap;
	TopExp::MapShapes(face, TopAbs_EDGE, originMap);

	const auto polySize = facePoly.size();
	for (auto i = 0; i < polySize; ++i)
	{
		const auto& p1 = facePoly[i];
		const auto& p2 = facePoly[(i + 1) % polySize];
		Handle(Geom2d_Edge) e1 = new Geom2d_Edge;
		e1->SetPoints(gp_Pnt2d(p1[0], p1[1]), gp_Pnt2d(p2[0], p2[1]));
		ProjectOnCurve.Init(gp_Pnt2d(s[0], s[1]), e1);
		if (ProjectOnCurve.NbPoints() > 0)
		{
			vecRefPnts.emplace_back(ProjectOnCurve.NearestPoint(), ProjectOnCurve.LowerDistance(), e1, i);
		}
	}

	std::pair<int, int> findBestRef;
	bool hasFindBest = false;
	double bestSumDist = DBL_MAX;

	for (int i = 0; i < vecRefPnts.size(); ++i)
	{
		if (vecRefPnts[i].edge->EndParameter() < 10.0) continue;
		const auto& dirEdgeA = vecRefPnts[i].edge->Direction();
		for (int j = i + 1; j < vecRefPnts.size(); ++j)
		{
			if (vecRefPnts[j].edge->EndParameter() < 10.0) continue;
			const auto& dirEdgeB = vecRefPnts[j].edge->Direction();
			if (!dirEdgeA.IsParallel(dirEdgeB, Precision::Confusion()))
			{
				const double sumDist = vecRefPnts[i].dist + vecRefPnts[j].dist;
				if (sumDist < bestSumDist)
				{
					bestSumDist = sumDist;
					hasFindBest = true;
					findBestRef.first = i;
					findBestRef.second = j;
				}
			}
		}
	}

	if (!hasFindBest)
	{
		return "";
	}

	std::stringstream ss;
	auto uid1 = GetShapeUID(TopoDS::Edge(originMap(findBestRef.first + 1)));
	auto uid2 = GetShapeUID(TopoDS::Edge(originMap(findBestRef.second + 1)));

	auto d1 = vecRefPnts[findBestRef.first].dist / 25.4;
	auto d2 = vecRefPnts[findBestRef.second].dist / 25.4;
	if (Math::Functs::IsAlmostZero(d1)) d1 = 0;
	if (Math::Functs::IsAlmostZero(d2)) d2 = 0;
	ss << "(" << uid1 << ", " << d1 << ", " << uid2 << ", " << d2 << ")";
	return ss.str();
}

void PHELM::GetPolygonOfFace(TopoDS_Face& face, std::vector<Vector2d>& points)
{
	gp_Dir dir = GeomFunc::GetFaceNormal(face);
	Vector3d face_normal_dir(dir.X(), dir.Y(), dir.Z());
	Handle(Geom_Surface) aSurface = BRep_Tool::Surface(face);
	Vector3d2 segments;
	//std::cerr << "one face://////////////" << std::endl;
	for (TopExp_Explorer edgeExp(face, TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
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
	Vector3d2 lines;
	Math::Functs::Connecting_Segments(segments, lines);
	Vector3d face_normal = Math::Functs::ComputeNormalFromPolyline(lines[0]);
	auto angle = Math::Functs::GetAngleBetween(face_normal, face_normal_dir);
	if (angle < Math::Math_PI / 2.0)
		std::reverse(lines[0].begin(), lines[0].end());

	//3d => 2d
	const auto base1 = GeomFunc::GetAnyPerpendicularVector(face_normal_dir);
	const auto base2 = glm::cross(face_normal_dir, base1);
	points.reserve(lines[0].size());
	for (int i = 0; i < lines[0].size(); i++)
	{
		const auto px = glm::dot(lines[0][i], base1);
		const auto py = glm::dot(lines[0][i], base2);
		points.emplace_back(px, py);
	}
		
}

void PHELM::GetUpDownFaces(const TopoDS_Shape& aShape, const bool& up_down, std::vector<Vector2d>& points)
{
	TopoDS_Face face;
	GetUpDownFaces(aShape, up_down, points, face);
}

void PHELM::GetUpDownFaces(const TopoDS_Shape& aShape, const bool& up_down, std::vector<Vector2d>& points, TopoDS_Face& face)
{
	for (TopExp_Explorer exp(aShape, TopAbs_FACE); exp.More(); exp.Next())
	{
		gp_Dir dir = GeomFunc::GetFaceNormal(TopoDS::Face(exp.Current()));
		Vector3d face_normal_dir(dir.X(), dir.Y(), dir.Z());
		double angle = Math::Functs::GetAngleBetween(Vector3d(0.0, 0.0, up_down?1:-1), face_normal_dir);

		if (Math::Functs::IsAlmostZero(angle))
		{
			face = TopoDS::Face(exp.Current());

			if (true)
			{
				Vector3d n_vec;
				Vector3d1 outside_contour, hole_vector;
				double hole_radius;
				bool circle_or_line = false;
				GeomFunc::GetOutsideContour(TopoDS::Face(exp.Current()), n_vec, outside_contour, hole_vector, hole_radius, circle_or_line);

				Vector3d face_normal = Math::Functs::ComputeNormalFromPolyline(outside_contour);
				angle = Math::Functs::GetAngleBetween(face_normal, face_normal_dir);
				if (angle < Math::Math_PI / 2.0)
					std::reverse(outside_contour.begin(), outside_contour.end());

				//3d => 2d
				points = Math::Functs::Vector3d2d(outside_contour);
			}
			else
			{
				Handle(Geom_Surface) aSurface = BRep_Tool::Surface(TopoDS::Face(exp.Current()));
				Vector3d2 segments;
				//std::cerr << "one face://////////////" << std::endl;
				for (TopExp_Explorer edgeExp(exp.Current(), TopAbs_EDGE); edgeExp.More(); edgeExp.Next())
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
				Vector3d2 lines;
				Math::Functs::Connecting_Segments(segments, lines);
				Vector3d face_normal = Math::Functs::ComputeNormalFromPolyline(lines[0]);
				angle = Math::Functs::GetAngleBetween(face_normal, face_normal_dir);
				if (angle < Math::Math_PI / 2.0)
					std::reverse(lines[0].begin(), lines[0].end());

				//3d => 2d
				points.reserve(lines[0].size());
				for (int i = 0; i < lines[0].size(); i++)
					points.emplace_back(lines[0][i].x, lines[0][i].y);
			}
			break;
		}
	}
}

void PHELM::GetOffsetUpFaces(const TopoDS_Shape& aShape, std::vector<Vector2d>& points)
{
	auto poly_offset = (CGAL_2D_Polygon_One_Offsets)GetProcAddress(*hModule, "CGAL_2D_Polygon_One_Offsets");
	std::vector<Vector2d> facePnts;
	TopoDS_Face face;
	GetUpDownFaces(aShape,true, facePnts, face);

	std::vector<std::vector<Vector2d>> offset_points;
	poly_offset(facePnts, -kerfDistance2, offset_points);
	points = offset_points[0];

	//	typedef double (*CGAL_2D_Polygon_One_Offsets)(std::vector<Vector2d> & poly,
	//double d, std::vector<std::vector<Vector2d>>& offset_polys);
}

void PHELM::GenerateCuttingPlane(TopoDS_Shape& prism, const Vector2d& s, const Vector2d& e, const gp_Dir& dir)
{
	gp_XYZ st(s.x, s.y, 0.), nd(e.x, e.y, 0.);
	gp_Dir offDir(s.y - e.y, e.x - s.x, 0);
	gp_XYZ offset(offDir.XYZ() * kerfDistance2);

	BRepBuilderAPI_MakePolygon mkWire(gp_Pnt(st + offset), gp_Pnt(nd + offset), gp_Pnt(nd - offset),
		gp_Pnt(st - offset), Standard_True);
	TopoDS_Shape sketchShape = mkWire.Shape();
	Part::TopoShape surface;
	double Ltotal = 1e4;
	double Loffset = Ltotal/2.0;
	gp_Trsf mov;
	mov.SetTranslation(-Loffset * gp_Vec(dir));
	TopLoc_Location loc(mov);
	sketchShape.Move(loc);

	try
	{
		surface.makEPrism(sketchShape, Ltotal * gp_Vec(dir)); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
}

void PHELM::GenerateCuttingPlane(TopoDS_Shape& prism, const Vector3d& s, const Vector3d& e, const gp_Dir& dir)
{
	gp_XYZ st(s.x, s.y, s.z), nd(e.x, e.y, e.z);
	gp_Dir cutDir(st - nd);
	auto offDir = cutDir.Crossed(dir);
	gp_XYZ offset(offDir.XYZ() * kerfDistance2);

	BRepBuilderAPI_MakePolygon mkWire(gp_Pnt(st + offset), gp_Pnt(nd + offset), gp_Pnt(nd - offset),
		gp_Pnt(st - offset), Standard_True);
	TopoDS_Shape sketchShape = mkWire.Shape();
	Part::TopoShape surface;
	double Ltotal = 1e4;
	double Loffset = Ltotal/2.0;
	gp_Trsf mov;
	mov.SetTranslation(-Loffset * gp_Vec(dir));
	TopLoc_Location loc(mov);
	sketchShape.Move(loc);

	try
	{
		surface.makEPrism(sketchShape, Ltotal * gp_Vec(dir)); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
}

void PHELM::GenerateCuttingSlice(TopoDS_Shape& prism, const Vector2d& s, const Vector2d& e, const gp_Dir& dir)
{
	BRepBuilderAPI_MakePolygon mkWire(gp_Pnt(s.x, s.y, 0.), gp_Pnt(e.x, e.y, 0.));
	TopoDS_Shape sketchShape = mkWire.Shape();
	Part::TopoShape surface;
	double Ltotal = 1e4;
	double Loffset = Ltotal/2.0;
	gp_Trsf mov;
	mov.SetTranslation(-Loffset * gp_Vec(dir));
	TopLoc_Location loc(mov);
	sketchShape.Move(loc);

	try
	{
		surface.makEPrism(sketchShape, Ltotal * gp_Vec(dir)); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
}

void PHELM::GenerateCuttingSlice(TopoDS_Shape& prism, const Vector3d& s, const Vector3d& e, const gp_Dir& dir)
{
	BRepBuilderAPI_MakePolygon mkWire(gp_Pnt(s.x, s.y, s.z), gp_Pnt(e.x, e.y, e.z));
	TopoDS_Shape sketchShape = mkWire.Shape();
	Part::TopoShape surface;
	double Ltotal = 1e4;
	double Loffset = Ltotal/2.0;
	gp_Trsf mov;
	mov.SetTranslation(-Loffset * gp_Vec(dir));
	TopLoc_Location loc(mov);
	sketchShape.Move(loc);

	try
	{
		surface.makEPrism(sketchShape, Ltotal * gp_Vec(dir)); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
}

void PHELM::GenerateCuttingPlane(TopoDS_Shape& prism, const std::vector<Vector3d>& poly, const double& offset)
{
	const gp_Dir offDir(1, 0, 0);
	BRepBuilderAPI_MakePolygon mkWire;
	for (const auto& p : poly)
	{
		gp_Pnt pnt(p[0], p[1], p[2]);
		mkWire.Add(pnt);
	}

	mkWire.Close();
	TopoDS_Shape sketchShape = mkWire.Shape();
	gp_Trsf mov;
	mov.SetTranslation(-offset * offDir.XYZ());
	TopLoc_Location loc(mov);
	sketchShape.Move(loc);
	Part::TopoShape surface;
	try
	{
		surface.makEPrism(sketchShape, 2.0 * offset * offDir.XYZ()); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
	//prism = sketchShape;
}

void PHELM::GenerateCuttingPlane(TopoDS_Shape& prism, const std::vector<Vector3d>& poly, TopoDS_Face& sketchShape,
	const double& offset)
{
	const gp_Dir offDir(1, 0, 0);
	BRepBuilderAPI_MakePolygon mkWire;
	for (const auto& p : poly)
	{
		gp_Pnt pnt(p[0], p[1], p[2]);
		mkWire.Add(pnt);
	}
	mkWire.Close();
	TopoDS_Shape wireShape = mkWire.Shape();
	sketchShape = BRepBuilderAPI_MakeFace(TopoDS::Wire(wireShape));
	gp_Trsf mov;
	mov.SetTranslation(-offset * offDir.XYZ());
	TopLoc_Location loc(mov);
	wireShape.Move(loc);
	TopoDS_Face face = BRepBuilderAPI_MakeFace(TopoDS::Wire(wireShape));
	Part::TopoShape surface;
	try
	{
		surface.makEPrism(face, 2.0 * offset * offDir.XYZ()); // finite prism
	}
	catch (Standard_Failure&)
	{
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();
	/*
	const gp_Dir offDir(1, 0, 0);
	Vector3d c1 = ComputeFaceCenter(poly);
	std::vector<Vector3d> poly_;
	GeomFunc::GetFaceVector(sketchShape, poly_);
	Vector3d c2 = ComputeFaceCenter(poly);
	gp_Vec offVec;
	if (c1[0] < c2[0])
	{
		offVec = -2.0 * offset * offDir.XYZ();
	}
	else
	{
		offVec = 2.0 * offset * offDir.XYZ();
	}

	Part::TopoShape surface;
	try {
		//BRepPrimAPI_MakePrism mkPrism(sketchShape, offVec);
		//prism = mkPrism.Shape();
		surface.makEPrism(sketchShape, offVec); // finite prism
	}
	catch (Standard_Failure&) {
		throw Base::RuntimeError("SketchBased: Length: Could not extrude the sketch!");
	}
	prism = surface.getShape();*/
}

bool PHELM::CGALIdentifyPolycutExtend(const std::vector<Vector2d>& polygon, const Vector2d& s, const Vector2d& e,
	Vector2d& ns, Vector2d& ne, Vector2d& refNs, Vector2d& refNe)
{
	auto cgal_2d_distance_point_poly = (CGAL_2D_Distance_Point_Polygon)GetProcAddress(
		*hModule, "CGAL_2D_Distance_Point_Polygon");
	auto cgal_construction_inoutside_poly = (CGAL_Construct_InOutSide_Polygon)GetProcAddress(
		*hModule, "CGAL_Construct_InOutSide_Polygon");
	auto cgal_2d_intersection_ray_seg = (CGAL_2D_Intersection_Ray_Segment)GetProcAddress(
		*hModule, "CGAL_2D_Intersection_Ray_Segment");
	auto get_angle_kerf_offset_tan = (GetAngleKerfOffsetTan)GetProcAddress(
		*hModule, "GetAngleKerfOffsetTan");

	const double eps = 0.1;
	//std::cerr << "input s = " << s[0] << " " << s[1] << std::endl;
	//std::cerr << "input e = " << e[0] << " " << e[1] << std::endl;
	const Vector2d cutDir = normalize(e - s);

	// copy to avoid default returns
	ns = s, ne = e;
	refNs = s, refNe = e;

	//return true;

	const int polySize = polygon.size();
	std::vector<Vector2d> rayD1Int, rayD2Int;
	std::vector<Vector2d> rayD1IntRef, rayD2IntRef;

	const auto pd1 = cgal_2d_distance_point_poly(s, polygon);
	const auto pd2 = cgal_2d_distance_point_poly(e, polygon);

	// if the points are in the outside of polygon, then we intersect them with polygon
	bool isoutside1 = false, isoutside2 = false;
	cgal_construction_inoutside_poly(polygon, s, e, isoutside1, isoutside2);

	if (isoutside1 && isoutside2)
		return true;

	//std::cerr << "inside1 = " << inside1 << " inside2 = " << inside2 << std::endl;

	if (Math::Functs::IsAlmostZero(pd1))
	{
	}
	else if (isoutside1)
	{
		std::vector<Vector2d> raySnap;
		std::vector<Vector2d> raySnapRef;
		for (int i = 0; i < polySize; i++)
		{
			Vector2d inter;
			const auto& pop1 = polygon[i];
			const auto& pop2 = polygon[(i + 1) % polySize];
			const Vector2d segDir = normalize(pop1 - pop2);
			const Vector2d conDir = normalize(pop1 - s);
			if (glm::abs(dot(cutDir, segDir)) > 0.9999 &&
				glm::abs(dot(cutDir, conDir)) > 0.9999)
			{
				raySnap.push_back(pop1);
				raySnap.push_back(pop2);
				raySnapRef.push_back(pop1);
				raySnapRef.push_back(pop2);
			}
			if (cgal_2d_intersection_ray_seg(s - eps * cutDir, cutDir, pop1, pop2, inter))
			{
				auto tanAngle = get_angle_kerf_offset_tan(segDir, cutDir);
				raySnapRef.push_back(inter);

				inter -= (1.5875 / tanAngle * cutDir);
				//std::cerr << "ray snap 1 angle" << tanAngle << std::endl;
				//std::cerr << "inter = " << inter[0] << " " << inter[1] << std::endl;
				raySnap.push_back(inter);
			}
		}
		//std::cerr << "raysnap1 = " << raySnap.size() << std::endl;
		if (raySnap.empty())
			return false;

		if (raySnap.size() != 1)
		{
			double minDist = DBL_MAX;
			for (auto k = 0; k < raySnap.size(); ++k)
			{
				const auto& v = raySnap[k];
				const double tmpDist2 = length2(v - s);
				if (tmpDist2 < minDist)
				{
					ns = v;
					refNs = raySnapRef[k];
					minDist = tmpDist2;
				}
			}
		}
		else
		{
			std::cerr << "raysnap1 % 2 != 0" << " size = " << raySnap.size() << std::endl;
		}
	}
	else
	{
		for (int i = 0; i < polySize; i++)
		{
			Vector2d inter;
			const auto& pop1 = polygon[i];
			const auto& pop2 = polygon[(i + 1) % polySize];
			Vector2d segDir = normalize(pop1 - pop2);

			if (cgal_2d_intersection_ray_seg(s + eps * cutDir, -cutDir, pop1, pop2, inter))
			{
				auto tanAngle = get_angle_kerf_offset_tan(segDir, cutDir);
				rayD1IntRef.push_back(inter);
				//std::cerr << "rayd1int" << tanAngle << std::endl;
				inter -= (1.5875 / tanAngle * cutDir);
				rayD1Int.push_back(inter);
			}
		}
		//std::cerr << "rayd1int = " << rayD1Int.size() << std::endl;
		if (rayD1Int.empty())
		{
			return false;
		}

		double minDist = DBL_MAX;
		for (auto k = 0; k < rayD1Int.size(); ++k)
		{
			const auto& v = rayD1Int[k];
			const double tmpDist2 = length2(v - s);
			if (tmpDist2 < minDist)
			{
				ns = v;
				refNs = rayD1IntRef[k];
				minDist = tmpDist2;
			}
		}
	}

	if (Math::Functs::IsAlmostZero(pd2))
	{
	}
	else if (isoutside2)
	{
		std::vector<Vector2d> raySnap;
		std::vector<Vector2d> raySnapRef;
		Vector2d inter;
		for (int i = 0; i < polySize; i++)
		{
			const auto& pop1 = polygon[i];
			const auto& pop2 = polygon[(i + 1) % polySize];
			const Vector2d segDir = normalize(pop1 - pop2);
			const Vector2d conDir = normalize(pop1 - s);
			if (glm::abs(dot(cutDir, segDir)) > 0.9999 &&
				glm::abs(dot(cutDir, conDir)) > 0.9999)
			{
				raySnap.push_back(pop1);
				raySnap.push_back(pop2);
				raySnapRef.push_back(pop1);
				raySnapRef.push_back(pop2);
			}
			else if (cgal_2d_intersection_ray_seg(e + eps * cutDir, -cutDir, pop1, pop2, inter))
			{
				raySnapRef.push_back(inter);
				auto tanAngle = get_angle_kerf_offset_tan(segDir, cutDir);
				inter += (1.5875 / tanAngle * cutDir);
				raySnap.push_back(inter);
			}
		}
		//std::cerr << "raysnap2 = " << raySnap.size() << std::endl;
		if (raySnap.empty())
			return false;
		if (raySnap.size() % 2 == 0)
		{
			double minDist = DBL_MAX;
			for (auto k = 0; k < raySnap.size(); ++k)
			{
				const auto& v = raySnap[k];
				//std::cerr << "update = " << v[0] << " " << v[1] << std::endl;
				const double tmpDist2 = length2(v - e);
				if (tmpDist2 < minDist)
				{
					ne = v;
					refNe = raySnapRef[k];
					minDist = tmpDist2;
				}
			}
		}
		else
		{
			std::cerr << "raysnap2 % 2 != 0" << " size = " << raySnap.size() << std::endl;
			//std::string filename = "D:\\cgaldebug\\" + std::to_string(rand()).append(".obj");
			//std::cerr << "raysnap2 % 2 != 0 and see " << filename << std::endl;
			// 			std::vector<std::vector<Vector2d> > poly = { polygon };
			// 			OutputRectangle(filename, poly);
			// 			auto fromS = e + eps * cutDir;
			// 			auto toE = -cutDir;
			// 			std::cerr << fromS[0] << " " << fromS[1] << " -> " << toE[0] << " " << toE[1]  << std::endl;
			// 			system("pause");
		}
	}
	else
	{
		for (int i = 0; i < polySize; i++)
		{
			Vector2d inter;
			const auto& pop1 = polygon[i];
			const auto& pop2 = polygon[(i + 1) % polySize];
			if (cgal_2d_intersection_ray_seg(e - eps * cutDir, cutDir, pop1, pop2, inter))
			{
				Vector2d segDir = normalize(pop1 - pop2);
				auto tanAngle = get_angle_kerf_offset_tan(segDir, cutDir);
				//std::cerr << "angle = " << tanAngle << std::endl;
				rayD2IntRef.push_back(inter);
				inter += (1.5875 / tanAngle * cutDir);
				rayD2Int.push_back(inter);
			}
		}
		//std::cerr << "rayd2int = " << rayD2Int.size() << std::endl;
		if (rayD2Int.empty()) return false;

		double minDist = DBL_MAX;
		for (auto k = 0; k < rayD2Int.size(); ++k)
		{
			const auto& v = rayD2Int[k];
			//std::cerr << "update = " << v[0] << " " << v[1] << std::endl;
			const double tmpDist2 = length2(v - e);
			if (tmpDist2 < minDist)
			{
				refNe = rayD2IntRef[k];
				ne = v;
				minDist = tmpDist2;
			}
		}
	}

	//std::cerr << "output s = " << ns[0] << " " << ns[1] << std::endl;
	//std::cerr << "output e = " << ne[0] << " " << ne[1] << std::endl;

	//std::cerr << "pd1 = " << pd1 << " pd2 = " << pd2 << std::endl;
	return true;
}
bool PHELM::ConvertVecUidToBuf(std::vector<boost::uuids::uuid>& vecUids)
{
	if (vecUids.empty()) return false;
	std::stringstream ss;
	ss << "(";

	for (auto i = 0; i < vecUids.size() - 1; ++i)
	{
		ss << vecUids[i] << ",";
	}

	ss << vecUids.back() << ") = ";

	retBuffer = ss.str();
	return true;
}

void PHELM::TestPost(PartDesignGui::PCE* pce, shared_ptr<PCEEClass> ec, std::vector<int>& nodes)
{
	if (ec->valid)
	{
		if (!Math::Functs::CheckContain(nodes, ec->index))
		{
			nodes.emplace_back(ec->index);
			for (auto& edges : ec->egraph_edges)
			{
				for (auto& edge : edges)
				{
					TestPost(pce, pce->GetEGraphCT().ptr_egraph[edge], nodes);
				}
			}
		}
	}
}

bool PHELM::CheckENodeEffectiveness(std::vector<shared_ptr<PCEEClass>> eg, shared_ptr<PCEEClass> ec)
{
	// An ENode is valid iff. 1) it has been compiled, 2) all of its offsprings is valid
	if (ec->valid) return true;

	for (auto edges : ec->egraph_edges)
	{
		bool b = true;
		for (auto edge : edges)
		{
			b &= CheckENodeEffectiveness(eg, eg[edge]);
		}

		if (b)
		{
			ec->valid = true; 
		}
	}
	return ec->valid;
}

bool PHELM::GetAllConnectedEClassId(std::vector<shared_ptr<PCEEClass>> eg, shared_ptr<PCEEClass> ec, std::set<int>& ecIds)
{
	// An ENode is valid iff. 1) it has been compiled, 2) all of its offsprings is valid
	if (ec->valid)
	{
		ecIds.insert(ec->index);
		for (auto edges : ec->egraph_edges)
		{
			for (auto edge : edges)
			{
				GetAllConnectedEClassId(eg, eg[edge], ecIds);
			}
		}
	}

	return true;
}

boost::uuids::uuid PHELM::GetShapeUID(const TopoDS_Shape& shape)
{
	if (mapShapeToUuid.Seek(shape) != nullptr)
	{
		//std::cerr << mapShapeToUuid(shape) << std::endl;
		return mapShapeToUuid(shape);
	}
	//boost::uuids::random_generator_pure::random_generator_pure(void)
	const auto uid = boost::uuids::random_generator()();
	mapShapeToUuid.Bind(shape, uid);
	return uid;
}

TopoDS_Shape PHELM::GetShapeDS(const boost::uuids::uuid& uid)
{
	/////////////////////////////////
	MapShapeToUuid::Iterator anIter(mapShapeToUuid);
	for (; anIter.More(); anIter.Next())
		if (anIter.Value() == uid)
			return anIter.Key();
	std::cerr << "TopoDS_Shape PHELM::GetShapeDS(const boost::uuids::uuid& uid)" << std::endl;
	system("pause");
	return TopoDS_Shape();
	/////////////////////////////////
}

// shapeA is the newly created shape, and shapeB is the original shape
bool PHELM::SyncTwoShapeUID(const TopoDS_Shape& shapeA, const TopoDS_Shape& shapeB)
{
	auto aUid = mapShapeToUuid.Seek(shapeA);
	auto bUid = mapShapeToUuid.Seek(shapeB);

	if (aUid == nullptr && bUid == nullptr) return false;
	if (bUid == nullptr) mapShapeToUuid.Bind(shapeB, *aUid);
	else mapShapeToUuid.Bind(shapeA, *bUid); 
	return true;
}

void PHELM::WriteProgram(std::list<std::string>& code, const int& cutLine_id,const bool& record)
{
	if (!fastMode)
	{
		if (code.empty() || xmlDoc == nullptr) return;

		//if(record)prog_strs.emplace_back(std::tuple<std::list<std::string>, int>(code, cutLine_id));

		auto equiv = xmlDoc->NewElement("equivalent");

		equiv->SetAttribute("cutLineId", cutLine_id);

		for (auto& c : code)
		{
			auto line = xmlDoc->NewElement("line");
			line->SetText(c.c_str());
			equiv->InsertEndChild(line);
		}
		nProg++;
		prog->InsertEndChild(equiv);

	
	}
	else
	{
		for (auto& c : code) 
		{
			std::cerr << c << std::endl;
		}
	}
}

bool PHELM::GetGeneratedFace(const Vector2d& ns, const Vector2d& ne, const TopoDS_Shape& it, TopoDS_Shape& copiedShape,
	TopoDS_Face& generatedFace)
{
	TopoDS_Shape cutSlice;
	GenerateCuttingSlice(cutSlice, ns, ne, gp_Dir(0, 0, 1));
	GEOMAlgo_Splitter splitterC;
	splitterC.SetFuzzyValue(globalFuzzyVal);
	splitterC.AddArgument(copiedShape);
	splitterC.AddTool(cutSlice);
	splitterC.Perform();

	TopTools_ListIteratorOfListOfShape iterC(splitterC.Modified(copiedShape));
	int ncParts = 0;
	TopoDS_Shape firstComponent;
	for (; iterC.More(); iterC.Next(), ++ncParts)
	{
		if (ncParts == 0) firstComponent = iterC.Value();
	}

	if (ncParts == 0)
	{
		std::cerr << "jump this cutting -- no intersection" << std::endl;

		/*GeomFunc::ExportToIGES(copiedShape, "D:\\Error1.igs");
		std::cerr << ns[0] << " " << ns[1] << " " << ne[0] << " " << ne[1] << std::endl;
		system("pause");*/
		return false;
	}

	for (TopExp_Explorer faceExplorer(firstComponent, TopAbs_FACE); faceExplorer.More(); faceExplorer.Next())
	{
		const TopoDS_Face& face = TopoDS::Face(faceExplorer.Current());
		gp_Pnt p = GeomFunc::GetCentralFacePoint(face);
		// check if point is in base shape
		BRepClass3d_SolidClassifier classifier;
		classifier.Load(it);
		classifier.Perform(p, 10.0 * Precision::Confusion());

		if (classifier.State() == TopAbs_IN)
		{
			generatedFace = face;
			//std::cerr << "found generated face" << std::endl;
		}
	}
	return true;
}

bool PHELM::GetGeneratedFace(const bool& debug_output, const Vector3d& ns, const Vector3d& ne, const gp_Dir& dir, const TopoDS_Shape& it,
	TopoDS_Shape& copiedShape, TopoDS_Face& generatedFace)
{
	TopoDS_Shape cutSlice;
	GenerateCuttingSlice(cutSlice, ns, ne, dir);
	GEOMAlgo_Splitter splitterC;
	splitterC.SetFuzzyValue(globalFuzzyVal);
	splitterC.AddArgument(copiedShape);
	splitterC.AddTool(cutSlice);
	splitterC.Perform();


	if (debug_output)
	{
		const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

		Vector3d3 listMaterial_surfs;
		listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(copiedShape));
		listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(cutSlice,false));
		Math::Functs::OutputObj3d(outputFolder + "\\log\\GetGeneratedFace.obj", listMaterial_surfs, 2);
	}

	TopTools_ListIteratorOfListOfShape iterC(splitterC.Modified(copiedShape));
	int ncParts = 0;
	TopoDS_Shape firstComponent;
	for (; iterC.More(); iterC.Next(), ++ncParts)
	{
		if (ncParts == 0) firstComponent = iterC.Value();
	}

	if (ncParts == 0)
	{
		//std::cerr << "jump this cutting -- no intersection" << std::endl;
		/*GeomFunc::ExportToIGES(copiedShape, "D:\\Error1.igs");
		GeomFunc::ExportToIGES(cutSlice, "D:\\ErrorFace.igs");
		std::cerr << ns[0] << " " << ns[1] << " " << ns[2] << " " << ne[0] << " " << ne[1] << " " << ne[2] << std::endl;
		system("pause");*/
		return false;
	}

	for (TopExp_Explorer faceExplorer(firstComponent, TopAbs_FACE); faceExplorer.More(); faceExplorer.Next())
	{
		const TopoDS_Face& face = TopoDS::Face(faceExplorer.Current());
		gp_Pnt p = GeomFunc::GetCentralFacePoint(face);
		// check if point is in base shape
		BRepClass3d_SolidClassifier classifier;
		classifier.Load(it);
		classifier.Perform(p, globalFuzzyVal);

		if (classifier.State() == TopAbs_IN)
		{
			generatedFace = face;
			//std::cerr << "found generated face" << std::endl;
		}
	}

	return !generatedFace.IsNull();
}