#ifndef CGAL_ONCE
#define CGAL_ONCE
#pragma once

#include <stdio.h>
#include <tchar.h>
#include <stdio.h>
#include <tchar.h>
#include "iostream"
#include <windows.h>
#include <vector>
#include <set>
#include <stdexcept>
#include <cstring>
#include <cstdlib>
#include <ostream>
#include <functional>
#include <queue>
#include <sstream>
#include <math.h>
#include <fstream>
#include "math.hpp"

using namespace std;

//Some necessary CGAL functions encapsulated in the dll (mydll.dll)
typedef void(*CGAL_Vector_Base)(Vector3d n, Vector3d& result);

typedef void(*CGAL_Export_Path_Segment)(std::ofstream &export_file_output, int &export_index,
	std::string s_name, double r, double g, double b, Vector3d start, Vector3d end, double radius);

typedef void(*CGAL_Export_Path_Point)(std::ofstream &export_file_output, int &export_index,
	std::string s_name, double r, double g, double b, Vector3d point, double radius);

typedef bool(*CGAL_2D_Location_Point_Polygon)(Vector2d p, std::vector<Vector2d> py);
typedef bool(*CGAL_2D_Location_Points_Polygon)(const std::vector<Vector2d>& ps, const std::vector<Vector2d>& py);

typedef double(*CGAL_2D_Two_Polygons_Intersection)(const std::vector<Vector2d> &poly_0, const std::vector<Vector2d> poly_1);

typedef double(*CGAL_2D_Two_Polygons_Union)(std::vector<Vector2d> poly_0, std::vector<Vector2d> poly_1,
	std::vector<std::vector<Vector2d>>& inter_polygons);

typedef double (*CGAL_2D_Distance_Point_Line)(Vector2d v, Vector2d l_0, Vector2d l_1);
typedef double (*CGAL_2D_Distance_Point_Segment)(Vector2d v, Vector2d l_0, Vector2d l_1);
typedef double (*CGAL_2D_Distance_Segment_Segment)(Vector2d s_0, Vector2d s_1, Vector2d e_0, Vector2d e_1);
typedef double (*CGAL_2D_Distance_Point_Polygon)(Vector2d p, std::vector<Vector2d> py);

typedef double (*CGAL_2D_Intersection_Line_Line)(const Vector2d& s_0_s, const Vector2d& s_0_e,
	const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);

typedef bool (*CGAL_2D_Intersection_Segment_Polygon)(Vector2d s_s, Vector2d s_e, std::vector<Vector2d>& p);

typedef double (*CGAL_2D_Polygon_One_Offsets)(std::vector<Vector2d>& poly,
	double d, std::vector<std::vector<Vector2d>>& offset_polys);

typedef bool (*CGAL_2D_Polygon_Is_Clockwise_Oriented)(std::vector<Vector2d>& ps);


typedef void (*CGAL_3D_Read_Triangle_Mesh)(std::string path, std::vector<Vector3d>& vecs,
	std::vector<int>& face_id_0, std::vector<int>& face_id_1,
	std::vector<int>& face_id_2);

typedef double(*CGAL_3D_Distance_Point_Segment)(Vector3d p, Vector3d s_s, Vector3d s_e);
typedef void(*CGAL_3D_Plane_Fitting)(std::vector<Vector3d> &points, Vector3d &plane_p, Vector3d &plane_v);
typedef void(*CGAL_3D_Plane_Point_Projection)(Vector3d &plane_p, Vector3d &plane_n, Vector3d &p, Vector3d& result);
typedef void(*CGAL_3D_Plane_Points_Projection)(Vector3d &plane_p, Vector3d &plane_n, std::vector<Vector3d> &points, std::vector<Vector3d> &project_points);

typedef void(*CGAL_3D_Plane_3D_to_2D_Point)(Vector3d &plane_p, Vector3d &plane_n, Vector3d &point_3d, Vector2d& result);
typedef void(*CGAL_3D_Plane_2D_to_3D_Point)(Vector3d &plane_p, Vector3d &plane_n, Vector2d &point_2d, Vector3d& result);

typedef void(*CGAL_3D_Plane_3D_to_2D_Points)(Vector3d &plane_p, Vector3d &plane_n, std::vector<Vector3d> &points_3d, std::vector<Vector2d> &points_2d);
typedef void(*CGAL_3D_Plane_2D_to_3D_Points)(Vector3d &plane_p, Vector3d &plane_n, std::vector<Vector2d> &points_2d, std::vector<Vector3d> &points_3d);

typedef double(*CGAL_2D_Distance_Point_Polygon)(Vector2d p, std::vector<Vector2d> py);
typedef std::vector<int>(*CGAL_Decompose_Polyline)(std::vector<Vector2d>& polyline, double threshold);
typedef bool(*CGAL_Identify_Polycut)(const std::vector<Vector2d>& polygon, const std::vector<Vector2d>& cutLine, std::vector<std::pair<bool, bool>>& result);
typedef bool (*CGAL_Identify_Polycut_Extend)(const std::vector<Vector2d>& polygon, const Vector2d& s, const Vector2d& e, Vector2d& ns, Vector2d& ne);
typedef bool (*CGAL_Identify_Polycut_NotExtend)(const std::vector<Vector2d>& polygon, const Vector2d& s, const Vector2d& e);
typedef double (*CGAL_3D_Distance_Point_Polygon)(const std::vector<Vector3d>& py, const Vector3d& p);
typedef bool (*CGAL_Construct_InOutSide_Polygon)(const std::vector<Vector2d>& py, const Vector2d& p, const Vector2d& q, bool& isPInside, bool& isQInside);
typedef bool (*CGAL_2D_Intersection_Ray_Segment)(const Vector2d& s_0_s, const Vector2d& s_0_e, const Vector2d& s_1_s, const Vector2d& s_1_e, Vector2d& inter);
typedef double (*GetAngleKerfOffsetTan)(const Vector2d& a, const Vector2d& b);
#endif
