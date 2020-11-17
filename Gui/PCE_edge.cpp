#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

PCEEdge::PCEEdge(const Vector3d& l_cutting_line_s_, const Vector3d& l_cutting_line_e_, const Vector3d& cutting_surface_normal_,
	const Vector3d& l_u_vector_)
{
	cut.l_cutting_line_s = l_cutting_line_s_;
	cut.l_cutting_line_e=l_cutting_line_e_;
	cut.cutting_surface_normal = cutting_surface_normal_;
	cut.l_u_vector = l_u_vector_;

	if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(Vector3d(0.0, 0.0, 1.0), cutting_surface_normal_) - Math::Math_PI / 2.0))
	{
		style = VerticalEnd;

		cutting_line_s = cut.l_cutting_line_s;
		cutting_line_e = cut.l_cutting_line_e;

		cut.offset_s = Vector2d(cut.l_cutting_line_s);
		cut.offset_e = Vector2d(cut.l_cutting_line_e);

		cut.u_cutting_line_s = cut.l_cutting_line_s + l_u_vector_;
		cut.u_cutting_line_e = cut.l_cutting_line_e + l_u_vector_;
	}
	else
	{
		style = TiltingEnd;
		cut.u_cutting_line_s = cut.l_cutting_line_s + l_u_vector_;
		cut.u_cutting_line_e = cut.l_cutting_line_e + l_u_vector_;

		if (Math::Functs::GetAngleBetween(cutting_surface_normal_, Vector3d(0.0, 0.0, 1.0)) < Math::Math_PI / 2.0)
		{
			cutting_line_s = cut.l_cutting_line_s;
			cutting_line_e = cut.l_cutting_line_e;
			cut.offset_s = Vector2d(cutting_line_s);
			cut.offset_e = Vector2d(cutting_line_e);
		}
		else
		{
			cutting_line_s = cut.u_cutting_line_s;
			cutting_line_e = cut.u_cutting_line_e;
			cut.offset_s = Vector2d(cutting_line_s);
			cut.offset_e = Vector2d(cutting_line_e);
		}
	}
}

PCEEdge::PCEEdge(const int& index_, const int& shape_index_, const PCECUT& cut_) :
	index(index_), part_index(shape_index_),cut(cut_)
{
	if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(Vector3d(0.0, 0.0, 1.0), cut_.cutting_surface_normal) - Math::Math_PI / 2.0))
		style = VerticalEnd;
	else
		style = TiltingEnd;

	Set(cut_);
}

void PCEEdge::Set(const PCECUT& cut_)
{
	cut = cut_;
	cut.offset_c = (cut.offset_s + cut.offset_e) / 2.0;
	cut.offset_n[0] = cut.offset_s[1] - cut.offset_e[1];//-y
	cut.offset_n[1] = cut.offset_e[0] - cut.offset_s[0];//x
	Math::Functs::SetVectorLength(cut.offset_n, 1.0);

	cut.origin_c = (cut.origin_s + cut.origin_e) / 2.0;
	cut.origin_n[0] = cut.origin_s[1] - cut.origin_e[1];//-y
	cut.origin_n[1] = cut.origin_e[0] - cut.origin_s[0];//x
	Math::Functs::SetVectorLength(cut.origin_n, 1.0);
}
