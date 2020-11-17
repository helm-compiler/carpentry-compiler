#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

PCEArrange::
PCEArrange(
	int index_,
	int wlc_index_,
	int wl_index_, 
	const std::vector<PCEPosition> & positions_,
	const std::vector<PCEEdge> & cutting_lines_,
	const double &volume_rate_) :
	index(index_), 
	cutting_number(cutting_lines_.size()),
	wlc_index(wlc_index_),
	wl_index(wl_index_),
	positions(positions_),
	cutting_lines(cutting_lines_),
	volume_rate(volume_rate_), 
	valid(false),
	eclass_index(-1),
	perfect_stock(false)
{
	Encode();
}

void PCEArrange::Encode()
{
	shape_order.clear();

	//encode
	Vector2d1 pos_centers;
	for (auto position : positions)
	{
		pos_centers.emplace_back(position.center);
		shape_order.emplace_back(position.transform->pce_shape->index);
	}

	Vector2d pos_center = Math::Functs::GetCenter(pos_centers);
	for (auto& center : pos_centers)center = center - pos_center;

	//Math::Functions::Vector2d3d();
	auto M = Math::Functs::RotationMatrix(Vector3d(0.0, 0.0, 1.0), Math::Math_PI);
	auto r_pos_centers = Math::Functs::Vector3d2d(Math::Functs::PosApplyM(Math::Functs::Vector2d3d(pos_centers), M));
	auto r_shape_order = shape_order;
	std::reverse(r_shape_order.begin(), r_shape_order.end());

	encode_0 = Vector2dString(pos_centers, shape_order, wl_index);
	encode_1 = Vector2dString(r_pos_centers, r_shape_order, wl_index);
	encode_packing_0 = Vector2dString(pos_centers, shape_order, -1);
	encode_packing_1 = Vector2dString(r_pos_centers, r_shape_order, -1);
}

PCEArrange::PCEArrange():
index(-1), eclass_index(-1), cutting_number(-1), wl_index(-1), wlc_index(-1), valid(false),volume_rate(0), perfect_stock(false)
{
}

void PCEArrange::Clear()
{
	for (auto& position : positions)position.Clear();
	std::vector<PCEPosition>().swap(positions);
	std::vector<PCEEdge>().swap(cutting_lines);
	std::vector<int>().swap(shape_order);
	Vector1i2().swap(cutting_orders);
	std::vector<PCEConstraint>().swap(sharing_constraints);
	std::vector<PCEEdge>().swap(edges);

	std::vector<std::vector<std::tuple<std::list<std::string>,int>>>().swap(prog_strs);
}

Vector3d3 PCEArrange::GetMesh(Vector3d2& surfs) const
{
	Vector3d3 all_points;
	for (auto position : positions)
		all_points.emplace_back(Math::Functs::PosApplyM(position.transform->surfs, position.M));
	return all_points;
}

std::string PCEArrange::Vector2dString(Vector2d1 ps, std::vector<int> ints, int wl_index_)
{
	const auto Comp = [](const Vector2d& a, const Vector2d& b)
	{
		return a[0] < b[0];
	};

	std::sort(ps.begin(), ps.end(), Comp);
	std::string str;
	for (auto p : ps)
	{
		double x = floor(p[0] * 10.0f + 0.5) / 10.0f;
		double y = floor(p[1] * 10.0f + 0.5) / 10.0f;
		str += Math::Functs::DoubleString(x);
		str += Math::Functs::DoubleString(y);
	}
	//sort(ints.begin(), ints.end());
	for (auto i : ints)
		str += Math::Functs::IntString(i);
	if (wl_index_ >= 0)
		str = str + Math::Functs::IntString(wl_index_);
	return str;
}