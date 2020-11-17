#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

void PCEEnd::PCECUT::Translate(const Vector2d& t)
{
	Vector3d t_3d(t[0], t[1], 0.0);
	origin_s += t;
	origin_e += t;
	origin_c += t;
	offset_s += t;
	offset_e += t;
	offset_c += t;
	u_cutting_line_s += t_3d;
	u_cutting_line_e += t_3d;
	l_cutting_line_s += t_3d;
	l_cutting_line_e += t_3d;
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PCEEnd::PCEEnd(
	const double& kerf, const int& index_, 
	const Vector3d& u_s_, const Vector3d& u_e_, 
	const Vector3d& d_s_, const  Vector3d& d_e_,
	const Vector3d2& con_surfs_,
	const Vector3d1& con_surf_normals_,
	const std::vector<int> &con_surf_indexes_) :
	index(index_), 
	u_s(u_s_), u_e(u_e_),
	l_s(d_s_), l_e(d_e_), 
	con_surfs(con_surfs_), 
	con_surf_normals(con_surf_normals_), 
	con_surf_indexes(con_surf_indexes_)
{
	style = -1;
	GetEdge(kerf);
}

void PCEEnd::GetEdge(const double& kerf)
{
	double angle = Math::Functs::GetAngleBetween(u_e - u_s, l_e - l_s);
	if (angle > Math::Math_PI / 2.0)angle = Math::Math_PI - angle;
	if (!Math::Functs::IsAlmostZero(angle))
	{
		std::cerr << "double angle = Math::Functions::GetAngleBetween(u_e - u_s, l_e - l_s)" << std::endl;
		system("pause");
	}

	auto x = Vector2d((u_e - u_s)[0], (u_e - u_s)[1]);
	auto o = Vector2d(((u_s + u_e) / 2.0)[0], ((u_s + u_e) / 2.0)[1]);
	auto y = Math::Functs::RotationAxis2d(x, -Math::Math_PI / 2.0, Vector2d(0.0, 0.0));
	if (Math::Functs::GetAngleBetween(y, Vector2d(con_surf_normals[0][0], con_surf_normals[0][1])) > Math::Math_PI / 2.0)
		x = -x;

	LocalCS cs(o, x);

	Vector2d l_s_2d = cs.GetLocal(Vector2d(l_s[0], l_s[1]));
	Vector2d l_e_2d = cs.GetLocal(Vector2d(l_e[0], l_e[1]));

	if (Math::Functs::IsAlmostZero(l_s_2d[1]) && Math::Functs::IsAlmostZero(l_e_2d[1]))
		style = VerticalEnd;
	else
		style = TiltingEnd;

	Vector2d2 surfs_2d;
	for (auto& surf : con_surfs)
	{
		Vector2d1 surf_2d;
		for (auto& p : surf)
			surf_2d.emplace_back(cs.GetLocal(Vector2d(p.x, p.y)));
		surfs_2d.emplace_back(surf_2d);
	}

	Vector2d minimal_corner, maximal_corner;
	Math::Functs::GetBoundingBox(surfs_2d, minimal_corner, maximal_corner);
	cut.origin_s = cs.GetGlobalPos(Vector2d(minimal_corner[0], maximal_corner[1]));
	cut.origin_e = cs.GetGlobalPos(Vector2d(maximal_corner[0], maximal_corner[1]));

	double delta_z = u_s[2] - l_s[2];

	if (style == VerticalEnd)
	{
		cut.offset_s = cs.GetGlobalPos(Vector2d(minimal_corner[0], maximal_corner[1] + kerf / 2.0));
		cut.offset_e = cs.GetGlobalPos(Vector2d(maximal_corner[0], maximal_corner[1] + kerf / 2.0));

		cut.cutting_surface_normal = con_surf_normals[0];
		cut.l_cutting_line_s = Vector3d(cut.offset_s[0], cut.offset_s[1], 0.0);
		cut.l_cutting_line_e = Vector3d(cut.offset_e[0], cut.offset_e[1], 0.0);
		cut.u_cutting_line_s = Vector3d(cut.offset_s[0], cut.offset_s[1], delta_z);
		cut.u_cutting_line_e = Vector3d(cut.offset_e[0], cut.offset_e[1], delta_z);
		cut.l_u_vector = Vector3d(0.0, 0.0, delta_z);
	}
	else
	{
		auto tilting_angle = Math::Functs::GetAngleBetween(con_surf_normals[0], Vector3d(0.0, 0.0, 1.0));
		double kerf_distance = 0.5 * kerf / cos(abs(tilting_angle - Math::Math_PI / 2.0));
		cut.offset_s = cs.GetGlobalPos(Vector2d(minimal_corner[0], maximal_corner[1] + kerf_distance));
		cut.offset_e = cs.GetGlobalPos(Vector2d(maximal_corner[0], maximal_corner[1] + kerf_distance));
		cut.cutting_surface_normal = con_surf_normals[0];
		if (tilting_angle < Math::Math_PI / 2.0)
		{
			cut.l_cutting_line_s = Vector3d(cut.offset_s[0], cut.offset_s[1], 0.0);
			cut.l_cutting_line_e = Vector3d(cut.offset_e[0], cut.offset_e[1], 0.0);
			Vector3d cutting_direction = Math::Functs::GetCrossproduct(cut.l_cutting_line_e - cut.l_cutting_line_s, cut.cutting_surface_normal);
			if (Math::Functs::GetAngleBetween(cutting_direction, Vector3d(0.0, 0.0, 1.0)) > Math::Math_PI / 2.0)
				cutting_direction = -cutting_direction;
			Math::Functs::SetVectorLength(cutting_direction, delta_z / cos(Math::Math_PI / 2.0 - tilting_angle));
			cut.u_cutting_line_s = cut.l_cutting_line_s + cutting_direction;
			cut.u_cutting_line_e = cut.l_cutting_line_e + cutting_direction;
			cut.l_u_vector = cutting_direction;
		}
		else
		{
			cut.u_cutting_line_s = Vector3d(cut.offset_s[0], cut.offset_s[1], delta_z);
			cut.u_cutting_line_e = Vector3d(cut.offset_e[0], cut.offset_e[1], delta_z);
			Vector3d cutting_direction = Math::Functs::GetCrossproduct(cut.u_cutting_line_e - cut.u_cutting_line_s, cut.cutting_surface_normal);
			if (Math::Functs::GetAngleBetween(cutting_direction, Vector3d(0.0, 0.0, 1.0)) < Math::Math_PI / 2.0)
				cutting_direction = -cutting_direction;
			Math::Functs::SetVectorLength(cutting_direction, delta_z / cos(tilting_angle - Math::Math_PI / 2.0));
			cut.l_cutting_line_s = cut.u_cutting_line_s + cutting_direction;
			cut.l_cutting_line_e = cut.u_cutting_line_e + cutting_direction;
			cut.l_u_vector = -cutting_direction;
		}
	}
}

void PCEEnd::Clear()
{
	Vector3d2().swap(con_surfs);
}
