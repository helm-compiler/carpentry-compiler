#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

LocalCS::LocalCS()
{
	origin = Vector2d(0.0,0.0);
	axis_x = Vector2d(1.0,0.0);
	axis_y = Math::Functs::RotationAxis2d(axis_x, -Math::Math_PI / 2.0, Vector2d(0.0, 0.0));
	Math::Functs::SetVectorLength(axis_x, 1.0);
	Math::Functs::SetVectorLength(axis_y, 1.0);
	angle = Math::Functs::GetAngleBetween(axis_x, Vector2d(1.0, 0.0));
	if (axis_x[1] < 0.0) angle = 2 * M_PI - angle;
	origin_max_y = 0.0;
	origin_min_y = 0.0;
	origin_min_x = 0.0;
	origin_max_x = 0.0;
	offset_max_y = 0.0;
	offset_min_y = 0.0;
	offset_min_x = 0.0;
	offset_max_x = 0.0;
}
LocalCS::LocalCS(const Vector2d& o, const Vector2d& x)
{
	origin = o;
	axis_x = x;
	axis_y = Math::Functs::RotationAxis2d(axis_x, -Math::Math_PI/2.0, Vector2d(0.0, 0.0));
	Math::Functs::SetVectorLength(axis_x, 1.0);
	Math::Functs::SetVectorLength(axis_y, 1.0);
	angle = Math::Functs::GetAngleBetween(axis_x, Vector2d(1.0, 0.0));
	if (axis_x[1] < 0.0) angle = 2 * M_PI - angle;
	origin_max_y = 0.0;
	origin_min_y = 0.0;
	origin_min_x = 0.0;
	origin_max_x = 0.0;
	offset_max_y = 0.0;
	offset_min_y = 0.0;
	offset_min_x = 0.0;
	offset_max_x = 0.0;
}

void LocalCS::SetCS(const Vector2d& o, const Vector2d& x)
{
	origin = o;
	axis_x = x;
	axis_y = Math::Functs::RotationAxis2d(axis_x, -Math::Math_PI / 2.0, Vector2d(0.0, 0.0));
	Math::Functs::SetVectorLength(axis_x, 1.0);
	Math::Functs::SetVectorLength(axis_y, 1.0);
	angle = Math::Functs::GetAngleBetween(axis_x, Vector2d(1.0, 0.0));
	if (axis_x[1] < 0.0) angle = 2 * M_PI - angle;
	origin_max_y = 0.0;
	origin_min_y = 0.0;
	origin_min_x = 0.0;
	origin_max_x = 0.0;
	offset_max_y = 0.0;
	offset_min_y = 0.0;
	offset_min_x = 0.0;
	offset_max_x = 0.0;
}

Vector2d LocalCS::GetLocal(const Vector2d& v) const
{
	Vector2d result = v - origin;
	result = Math::Functs::RotationAxis2d(result, angle, Vector2d(0.0, 0.0));
	return result;
}

Vector2d1 LocalCS::GetLocal(const Vector2d1& points) const
{
	Vector2d1 results;
	for (auto p : points)
		results.emplace_back(GetLocal(p));
	return results;
}

//need to check
Vector2d LocalCS::GetGlobalPos(const Vector2d& v) const
{
	Vector2d result = Math::Functs::RotationAxis2d(v, -angle, Vector2d(0.0, 0.0));
	result = result + origin;
	return result;
}

Vector2d LocalCS::GetGlobalVec(const Vector2d& v) const
{
	Vector2d result = Math::Functs::RotationAxis2d(v, -angle, Vector2d(0.0, 0.0));
	return result;
}

bool LocalCS::GetSegOrigin(const HMODULE& hModule, const double part_match_error, const Vector2d& s, const Vector2d& e, double& min_x, double& max_x) const
{
	auto cgal_2d_inter_line_line = (CGAL_2D_Intersection_Line_Line)GetProcAddress(hModule, "CGAL_2D_Intersection_Line_Line");
	Vector2d low_v = e;
	Vector2d upper_v = s;
	if (s[1] <= e[1])
	{
		low_v = s;
		upper_v = e;
	}

	//outside
	if (low_v[1] > origin_max_y || upper_v[1] < origin_min_y ||
		Math::Functs::IsAlmostZero_Double(low_v[1] - origin_max_y, part_match_error) || Math::Functs::IsAlmostZero_Double(upper_v[1] - origin_min_y, part_match_error))
		return false;

	//all in
	if ((low_v[1] > origin_min_y||Math::Functs::IsAlmostZero_Double(low_v[1]- origin_min_y, part_match_error)) &&
		(upper_v[1] < origin_max_y||Math::Functs::IsAlmostZero_Double(upper_v[1]- origin_max_y, part_match_error)))
	{
		min_x = low_v[0];
		max_x = upper_v[0];
		if (low_v[0] > upper_v[0])
		{
			min_x = upper_v[0];
			max_x = low_v[0];
		}
		return true;
	}

	//cutting lower
	if (low_v[1]< origin_min_y && upper_v[1]>origin_min_y && (upper_v[1] < origin_max_y||
		Math::Functs::IsAlmostZero_Double(upper_v[1]- origin_max_y, part_match_error)))
	{
		Vector2d inter;
		if (cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, origin_min_y), Vector2d(1.0, origin_min_y), inter))
		{
			min_x = inter[0];
			max_x = upper_v[0];
			if (inter[0] > upper_v[0])
			{
				min_x = upper_v[0];
				max_x = inter[0];
			}
			return true;
		}
	}

	//cutting upper
	if (upper_v[1] > origin_max_y && (low_v[1] > origin_min_y||Math::Functs::IsAlmostZero_Double(low_v[1]- origin_min_y, part_match_error)) && low_v[1] < origin_max_y)
	{
		Vector2d inter;
		if (cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, origin_max_y), Vector2d(1.0, origin_max_y), inter))
		{
			min_x = low_v[0];
			max_x = inter[0];
			if (low_v[0] > inter[0])
			{
				min_x = inter[0];
				max_x = low_v[0];
			}
			return true;
		}
	}


	if (upper_v[1] > origin_max_y && low_v[1] < origin_min_y)
	{
		Vector2d inter_0;
		bool b0 = cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, origin_min_y), Vector2d(1.0, origin_min_y), inter_0);
		Vector2d inter_1;
		bool b1 = cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, origin_max_y), Vector2d(1.0, origin_max_y), inter_1);

		if (b0 && b1)
		{
			min_x = inter_0[0];
			max_x = inter_1[0];
			if (inter_0[0] > inter_1[0])
			{
				min_x = inter_1[0];
				max_x = inter_0[0];
			}
			return true;
		}
	}
	return false;
}

bool LocalCS::GetSegOffset(const HMODULE& hModule, const double part_match_error, const Vector2d& s, const Vector2d& e, double& min_x, double& max_x) const
{
	auto cgal_2d_inter_line_line = (CGAL_2D_Intersection_Line_Line)GetProcAddress(hModule, "CGAL_2D_Intersection_Line_Line");

	Vector2d low_v = e;
	Vector2d upper_v = s;
	if (s[1] <= e[1])
	{
		low_v = s;
		upper_v = e;
	}

	//outside
	if (low_v[1] > offset_max_y || upper_v[1] < offset_min_y ||
		Math::Functs::IsAlmostZero_Double(low_v[1] - offset_max_y, part_match_error) ||
		Math::Functs::IsAlmostZero_Double(upper_v[1] - offset_min_y, part_match_error))
		return false;

	//all in
	if (low_v[1] >= offset_min_y && upper_v[1] <= offset_max_y)
	{
		min_x = low_v[0];
		max_x = upper_v[0];
		if (low_v[0] > upper_v[0])
		{
			min_x = upper_v[0];
			max_x = low_v[0];
		}
		return true;
	}

	//cutting lower
	if (low_v[1]< offset_min_y && upper_v[1]>offset_min_y && upper_v[1] <= offset_max_y)
	{
		Vector2d inter;
		if (cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, offset_min_y), Vector2d(1.0, offset_min_y), inter))
		{
			min_x = inter[0];
			max_x = upper_v[0];
			if (inter[0] > upper_v[0])
			{
				min_x = upper_v[0];
				max_x = inter[0];
			}
			return true;
		}
	}

	//cutting upper
	if (upper_v[1] > offset_max_y && low_v[1] >= offset_min_y && low_v[1] < offset_max_y)
	{
		Vector2d inter;
		if (cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, offset_max_y), Vector2d(1.0, offset_max_y), inter))
		{
			min_x = low_v[0];
			max_x = inter[0];
			if (low_v[0] > inter[0])
			{
				min_x = inter[0];
				max_x = low_v[0];
			}
			return true;
		}
	}


	if (upper_v[1] > offset_max_y && low_v[1] < offset_min_y)
	{
		Vector2d inter_0;
		bool b0 = cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, offset_min_y), Vector2d(1.0, offset_min_y), inter_0);
		Vector2d inter_1;
		bool b1 = cgal_2d_inter_line_line(low_v, upper_v, Vector2d(0.0, offset_max_y), Vector2d(1.0, offset_max_y), inter_1);

		if (b0 && b1)
		{
			min_x = inter_0[0];
			max_x = inter_1[0];
			if (inter_0[0] > inter_1[0])
			{
				min_x = inter_1[0];
				max_x = inter_0[0];
			}
			return true;
		}
	}
	return false;
}

Vector2d1 LocalCS::SortX(const Vector2d1& segs) const
{
	Vector2d1 sort_segs;
	if (segs.size() == 0) return sort_segs;
	double minimal_x = segs[0][0];
	double maximal_x = segs[0][1];

	for (int iter = 0; iter < segs.size(); iter++)
	{
		auto& seg = segs[iter];
		if (seg[0] >= minimal_x && seg[0] <= maximal_x)
		{
			if (seg[1] > maximal_x)
				maximal_x = seg[1];
		}
		else
		{
			sort_segs.emplace_back(minimal_x, maximal_x);
			minimal_x = segs[iter][0];
			maximal_x = segs[iter][1];
		}

		if (iter == segs.size() - 1)
		{
			sort_segs.emplace_back(minimal_x, maximal_x);
			break;
		}
	}
	return sort_segs;
}

Vector2d1 LocalCS::EmptyX(const Vector2d1& sort_segs) const
{
	Vector2d1 empty_segs;
	for (int i = 1; i < sort_segs.size(); i++)
		empty_segs.emplace_back(sort_segs[i - 1][1], sort_segs[i][0]);
	return empty_segs;
}