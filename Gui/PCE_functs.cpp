#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Mod/PartDesign/Gui/CompilerConfig.h>
#include <tchar.h>
#include <iostream>


using namespace PartDesignGui;
using namespace Math;

std::vector<PCEArrange> PCE::LocalPacking(WLCollection& wlc, const PCEArrange& arrange)
{
	WL& wl = wlc.wls[arrange.wl_index];

	std::vector<PCEArrange> opt_arranges;
	std::vector<PCEPosition> fixed_positions;
	std::vector<PCEConstraint> sharing_constraints;
	std::vector<PCEEdge> edges = wl.ExtractEdges();

	for (int i = 0; i < arrange.positions.size(); i++)
		//for (auto& position : arrange.positions)
	{
		auto& position = arrange.positions[i];

		PCEPosition pos;
		pos.transform = position.transform;

		pos.transform->AddNewEdges(edges, pos.edge_indexes);
		ArrangePartbyConstraints(wl, edges, fixed_positions, pos, sharing_constraints);
		if (pos.style != PCEPosition::FreeStyle)
		{
			fixed_positions.emplace_back(pos);
		}
		else
		{
			if (fixed_positions.empty())
			{
				std::cerr << "if (fixed_positions.empty())" << std::endl;
				system("pause");
			}

			PCEArrange packing_arrange = arrange;

			std::vector<PCEEdge> cutting_lines;
			GetCuttingLines(edges, sharing_constraints, fixed_positions, cutting_lines);

			packing_arrange.volume_rate = -1.0;
			packing_arrange.cutting_lines = cutting_lines;
			packing_arrange.positions = fixed_positions;
			packing_arrange.Encode();

			opt_arranges.emplace_back(packing_arrange);


			fixed_positions.clear();
			sharing_constraints.clear();
			edges = wl.ExtractEdges();

			if (i == 0)
			{
				std::cerr << "if (fixed_positions.empty())" << std::endl;
				system("pause");
			}
			i = i - 1;
		}
	}


	PCEArrange packing_arrange = arrange;

	std::vector<PCEEdge> cutting_lines;
	GetCuttingLines(edges, sharing_constraints, fixed_positions, cutting_lines);

	packing_arrange.volume_rate = -1.0;
	packing_arrange.cutting_lines = cutting_lines;
	packing_arrange.positions = fixed_positions;
	packing_arrange.Encode();

	opt_arranges.emplace_back(packing_arrange);


	for (auto& opt_arrange : opt_arranges)
	{
		opt_arrange.Encode();
	}


	//TODO: add computing volume rate
	////volume rate
	//double total_volume = 0.0;
	//for (auto& position : fixed_positions)
	//{
	//	total_volume += position.transform->pce_shape->volume;
	//}
	//total_volume = total_volume / wl.volume;

	return opt_arranges;
}

void PCE::Visualization(
	const std::vector<double>& score,
	const std::vector<WLCollection>& wlcs_, 
	const Vector1i1& order_,
	const std::vector<std::pair<int, int>>& order_stacks_,
	const std::vector<std::pair<int, int>>& order_merges_,
	const Vector1i1& ArrIDs_,
	const Vector1i1& ENodeIDs_,
	const std::vector<std::pair<Vector3d, Vector3d>>& refer_edges_,
	const Vector3d2& refer_faces_,
	const Vector3d1& refer_values_,
	const std::string& path_)
{
	auto GetLUvector = [&](
		const Vector2d& offset_s,
		const Vector2d& offset_e,
		const Vector3d& cutting_surface_normal,
		const double& delta_z,
		Vector3d& l_u_vector,
		Vector3d& l_cutting_line_s,
		Vector3d& l_cutting_line_e)
	{
		auto tilting_angle = Math::Functs::GetAngleBetween(cutting_surface_normal, Vector3d(0.0, 0.0, 1.0));

		if (Math::Functs::IsAlmostZero_Double(tilting_angle - Math::Math_PI / 2.0, angle_match_error))
		{
			l_u_vector = Vector3d(0.0, 0.0, delta_z);
			l_cutting_line_s = Vector3d(offset_s[0], offset_s[1], 0.0);
			l_cutting_line_e = Vector3d(offset_e[0], offset_e[1], 0.0);
		}
		else
		{
			if (tilting_angle < Math::Math_PI / 2.0)
			{
				l_cutting_line_s = Vector3d(offset_s[0], offset_s[1], 0.0);
				l_cutting_line_e = Vector3d(offset_e[0], offset_e[1], 0.0);
				Vector3d cutting_direction = Math::Functs::GetCrossproduct(l_cutting_line_e - l_cutting_line_s, cutting_surface_normal);
				Math::Functs::SetVectorLength(cutting_direction, delta_z / cos(Math::Math_PI / 2.0 - tilting_angle));
				l_u_vector = cutting_direction;
			}
			else
			{
				auto u_ing_line_s = Vector3d(offset_s[0], offset_s[1], delta_z);
				auto u_ing_line_e = Vector3d(offset_e[0], offset_e[1], delta_z);

				Vector3d cutting_direction = Math::Functs::GetCrossproduct(u_ing_line_e - u_ing_line_s, cutting_surface_normal);
				Math::Functs::SetVectorLength(cutting_direction, delta_z / cos(tilting_angle - Math::Math_PI / 2.0));
				if (Math::Functs::GetAngleBetween(cutting_direction, Vector3d(0.0, 0.0, 1.0)) < Math::Math_PI / 2.0)
					l_u_vector = cutting_direction;
				else
					l_u_vector = -cutting_direction;
				l_cutting_line_s = u_ing_line_s - l_u_vector;
				l_cutting_line_e = u_ing_line_e - l_u_vector;
			}
		}
	};

	std::vector<int> all_arranges;
	std::vector<std::vector<PCEEdge>> all_cutting_lines;
	std::vector<std::vector<int>> all_cutting_ids;

	int id_nb = 0;
	for(int iter=0;iter< ArrIDs_.size();iter++)
	{
		int arrange_index= ArrIDs_[iter]; //ArrID
		int order_index = ENodeIDs_[iter]; //ENodeID

		all_arranges.emplace_back(arrange_index);

		auto& ar = ct_egraph.ptr_arranges[arrange_index];
		auto& wlc = wlcs_[ar->wlc_index];
		auto& wl = wlc.wls[ar->wl_index];

		auto cutting_order = ar->cutting_orders[order_index];
		int segment_nb = ar->cutting_number;

		std::vector<PCEEdge> cutting_lines;
		std::vector<int> cutting_line_id;

		for (int i = 0; i < cutting_order.size(); i++)
		{
			auto cutLine = ar->cutting_lines[cutting_order[i]];
			Vector2d s= cutLine.cut.offset_s;
			Vector2d e= cutLine.cut.offset_e;
			Vector3d n= cutLine.cut.cutting_surface_normal;
			Vector3d l_u_vector;
			Vector3d l_cutting_line_s, l_cutting_line_e;
			GetLUvector(s, e, n, wl.size[2], l_u_vector, l_cutting_line_s, l_cutting_line_e);
			cutting_lines.emplace_back(PCEEdge(l_cutting_line_s, l_cutting_line_e, n, l_u_vector));
			cutting_line_id.emplace_back(i + id_nb);
		
		}

		id_nb += segment_nb;
		all_cutting_lines.emplace_back(cutting_lines);
		all_cutting_ids.emplace_back(cutting_line_id);
	}

	std::vector<std::vector<std::pair<Vector3d, Vector3d>>> all_refer_edges;
	Vector3d2 cut_colors;
	Vector3d2 refer_values;
	Vector3d3 refer_faces;

	for (auto& ids : all_cutting_ids)
	{
		all_refer_edges.emplace_back(std::vector<std::pair<Vector3d, Vector3d>>());
		cut_colors.emplace_back(Vector3d1());
		refer_values.emplace_back(Vector3d1());
		refer_faces.emplace_back(Vector3d2());
		for (auto& id : ids)
		{
			id = Math::Functs::VectorIndex(order_, id);
			all_refer_edges.back().emplace_back(refer_edges_[id]);
			double r, g, b;
			Math::Functs::ColorMapping((double)(id) / (double)order_.size(), r, g, b);
			cut_colors.back().emplace_back(r, g, b);
			refer_values.back().emplace_back(refer_values_[id]);
			refer_faces.back().emplace_back(refer_faces_[id]);
		}
	}

	std::ofstream export_fie(path_);
	int export_int = 1;
	Vector3d y_delta(0.0);

	double delta_height = 0.0;
	for (int i = 0; i < all_arranges.size(); i++)
	{
		int arrange_index = all_arranges[i];
		auto& ar = ct_egraph.ptr_arranges[arrange_index];
		int wl_index = ar->wl_index;
		auto& wlc = wlcs_[ar->wlc_index];
		auto& wl = wlc.wls[wl_index];

		OutputFixedPositionsVIS(
			"arr_" + std::to_string(i),
			wl, ar,
			all_cutting_lines[i],
			all_cutting_ids[i],
			cut_colors[i],
			all_refer_edges[i],
			refer_faces[i],
			refer_values[i],
			y_delta,
			export_fie, 
			export_int);

		y_delta[1]+=wl.size[1]*1.5;
		delta_height = wl.size[1];
	}

	//setup
	double own_fab_time = 0.0;
	auto output = (CGAL_Export_Path_Segment)GetProcAddress(*hModule, "CGAL_Export_Path_Segment");
	for (int i = 0; i < order_.size(); i++)
	{
		Vector3d pre_c(delta_height* (i-1), -delta_height * 1.5, 0.0);
		Vector3d cur_c(delta_height* i, -delta_height * 1.5, 0.0);
		Vector3d next_c(delta_height* (i+1), -delta_height * 1.5, 0.0);
		number_label.DrawLabel(NumberLabel::BaseCircle, export_fie, export_int, "cut_" + Math::Functs::IntString(i),i, kerf * 6.0, cur_c);

		pre_c[2] = -kerf * 0.5;
		cur_c[2] = -kerf * 0.5;
		output(export_fie, export_int, "name", 0.0, 0.0, 0.0, pre_c, cur_c, kerf * 0.3);

		if (i == 0)
		{
			auto center = (pre_c + cur_c) / 2.0;
			center[1] += kerf * 3.0;
			number_label.DrawNode(NumberLabel::FullSetup, export_fie, export_int, "setup", kerf * 3.0, center);

			own_fab_time += 1.0;
		}
		else
		{
			bool b0= Math::Functs::CheckContain(order_stacks_, std::pair<int, int>(i - 1, i))|| Math::Functs::CheckContain(order_stacks_, std::pair<int, int>(i, i - 1));
			bool b1 = Math::Functs::CheckContain(order_merges_, std::pair<int, int>(i - 1, i))|| Math::Functs::CheckContain(order_merges_, std::pair<int, int>(i, i - 1));
			auto center = (pre_c + cur_c) / 2.0;
			center[1] += kerf * 3.0;
			int label = -1;
			if (!b0 && !b1) {
				label = NumberLabel::FullSetup; 			
				own_fab_time += 1.0;
			}
			if (!b0 && b1)
			{
				label = Math::Functs::IsAlmostZero(refer_values_[i][2] - refer_values_[i - 1][2]) ? NumberLabel::HalfAngle : NumberLabel::HalfOffset;
				own_fab_time += 0.5;
			}

			if(label>=0)number_label.DrawNode(label, export_fie, export_int, "setup", kerf * 3.0, center);
		}
	}

	//compute fab time


	//score
	number_label.DrawNumber(export_fie, export_int, "score_Mat", score[0], kerf * 8.0, Vector3d(delta_height * 1.5 * 0, -delta_height * 0.75, 0.0));
	number_label.DrawNumber(export_fie, export_int, "score_Pre", score[1], kerf * 8.0, Vector3d(delta_height * 1.5 * 1, -delta_height * 0.75, 0.0));
	number_label.DrawNumber(export_fie, export_int, "score_Time", score[2], kerf * 8.0, Vector3d(delta_height * 1.5 * 2, -delta_height * 0.75, 0.0));

	if(!Math::Functs::IsAlmostZero(own_fab_time- score[2]))
		number_label.DrawNumber(export_fie, export_int, "score_Time", own_fab_time, kerf * 8.0, Vector3d(delta_height * 1.5 * 3, -delta_height * 0.75, 0.0));

	//const std::vector<std::pair<int, int>>& order_stacks_,
		//const std::vector<std::pair<int, int>>& order_merges_,

	export_fie.close();
}

int PCE::GetScore(const std::vector<PCEEdge>& edges, const std::vector<PCEPosition>& fixed_positions,
	const PCEPosition& free_position,
	const std::vector<PCEConstraint>& constraints, std::vector<PCEConstraint>& comb_edges,
	double& total_length, double& boundary_distance)
{
	auto Output = [&]()
	{
		std::ofstream export_fie(packing_o_folder+"\\debug.obj");

		auto output = (CGAL_Export_Path_Segment)GetProcAddress(*hModule, "CGAL_Export_Path_Segment");
		int export_int = 1;
		for (auto iter = 0; iter < edges.size(); iter++)
		{
			Vector3d s = Math::Functs::Vector2d3d(edges[iter].cut.offset_s);
			Vector3d e = Math::Functs::Vector2d3d(edges[iter].cut.offset_e);
			output(export_fie, export_int, "edge_" + Math::Functs::IntString(iter), -1.0, -1.0, -1.0, s, e, kerf * 0.5);
		}

		export_fie.close();
	};
	//Output();


	const auto cgal_2d_distance_segment_segment = 
		(CGAL_2D_Distance_Segment_Segment)GetProcAddress(*hModule, "CGAL_2D_Distance_Segment_Segment");

	auto GetMaximalLength = [](Vector2d e_0_s, Vector2d e_0_e, Vector2d e_1_s, Vector2d e_1_e)
	{
		std::vector<double> lens(1, Math::Functs::GetDistance(e_0_s, e_0_e));
		lens.emplace_back(Math::Functs::GetDistance(e_0_s, e_1_s));
		if (Math::Functs::GetAngleBetween(e_0_e - e_0_s, e_1_s - e_0_s) > Math::Math_PI / 2.0) lens.back() = -lens.back();
		lens.emplace_back(Math::Functs::GetDistance(e_0_s, e_1_e));
		if (Math::Functs::GetAngleBetween(e_0_e - e_0_s, e_1_e - e_0_s) > Math::Math_PI / 2.0) lens.back() = -lens.back();
		std::sort(lens.begin(), lens.end());
		return lens.back() - lens.front();
	};

	auto CollisionDetection = [](const HMODULE& hModule_, const std::vector<PCEPosition>& fixed_positions_, const PCEPosition& free_position_,
		const PCEEdge& fixed_edge, const PCEEdge& free_edge)
	{
		//CGAL_2D_Intersection_Segment_Polygon
		const auto intersection = (CGAL_2D_Intersection_Segment_Polygon)GetProcAddress(hModule_, "CGAL_2D_Intersection_Segment_Polygon");

		Vector2d1 free_u_origin = free_position_.transform->UOrigin(free_position_.t);//collision
		Vector2d1 free_l_origin = free_position_.transform->LOrigin(free_position_.t);//collision

		if (intersection(Vector2d(fixed_edge.cut.u_cutting_line_s), Vector2d(free_edge.cut.u_cutting_line_e) + free_position_.t, free_u_origin))
			return true;
		if (intersection(Vector2d(fixed_edge.cut.l_cutting_line_s), Vector2d(free_edge.cut.l_cutting_line_e) + free_position_.t, free_l_origin))
			return true;

		for (auto& position : fixed_positions_)
		{
			Vector2d1 fixed_u_origin = position.transform->UOrigin(position.t);//collision
			Vector2d1 fixed_l_origin = position.transform->LOrigin(position.t);//collision

			if (intersection(Vector2d(fixed_edge.cut.u_cutting_line_s), Vector2d(free_edge.cut.u_cutting_line_e) + free_position_.t, fixed_u_origin))
				return true;
			if (intersection(Vector2d(fixed_edge.cut.l_cutting_line_s), Vector2d(free_edge.cut.l_cutting_line_e) + free_position_.t, fixed_l_origin))
				return true;
		}
		return false;
	};

	auto Get2DVector = [](Vector3d v) {return Vector2d(v[0], v[1]); };
	total_length = 0.0;
	auto score = 0;
	boundary_distance = 0.0;

	for (auto constraint : constraints)
	{
		const auto& fixed_edge = edges[constraint.edge_0];
		const auto& free_edge = edges[constraint.edge_1];

		double a = Math::Functs::GetAngleBetween(fixed_edge.cut.cutting_surface_normal, free_edge.cut.cutting_surface_normal) - Math::Math_PI;

		if (fixed_edge.style == PCEEdge::TiltingEnd && free_edge.style == PCEEdge::TiltingEnd &&
			Math::Functs::IsAlmostZero_Double(a, angle_match_error))
		{
			auto o1_ = (fixed_edge.cut.u_cutting_line_s + fixed_edge.cut.u_cutting_line_e) / 2.0;
			auto o2_ = (free_edge.cut.u_cutting_line_s + free_edge.cut.u_cutting_line_e) / 2.0;

			auto o1 = Vector2d(o1_[0], o1_[1]);
			auto o2 = Vector2d(o2_[0], o2_[1]);

			const double d = cgal_2d_distance_segment_segment(Get2DVector(fixed_edge.cut.u_cutting_line_s), Get2DVector(fixed_edge.cut.u_cutting_line_e),
				Get2DVector(free_edge.cut.u_cutting_line_s) + free_position.t,
				Get2DVector(free_edge.cut.u_cutting_line_e) + free_position.t);

			auto free_center = (Get2DVector(free_edge.cut.u_cutting_line_s) + Get2DVector(free_edge.cut.u_cutting_line_e)) / 2.0;
			auto fixed_center = (Get2DVector(fixed_edge.cut.u_cutting_line_s) + Get2DVector(fixed_edge.cut.u_cutting_line_e)) / 2.0;

			const double a = Math::Functs::GetDotproduct(free_center + free_position.t - fixed_center, fixed_edge.cut.offset_n);
			if (Math::Functs::IsAlmostZero_Double(a, angle_match_error))
			{
				if (d < maximal_cutting_distance)
				{
					//collision
					if (!CollisionDetection(*hModule, fixed_positions, free_position, fixed_edge, free_edge))
					{
						comb_edges.emplace_back(fixed_edge.index, free_edge.index);
						score++;
					}
				}
			}
		}
		else
		{
			const double d = cgal_2d_distance_segment_segment(fixed_edge.cut.offset_s, fixed_edge.cut.offset_e, free_edge.cut.offset_s + free_position.t,
				free_edge.cut.offset_e + free_position.t);
			const double a = Math::Functs::GetDotproduct(free_edge.cut.offset_c + free_position.t - fixed_edge.cut.offset_c, fixed_edge.cut.offset_n);
			if (Math::Functs::IsAlmostZero_Double(a, angle_match_error))
			{
				if (d < maximal_cutting_distance)
				{
					//collision
					if (!CollisionDetection(*hModule, fixed_positions, free_position, fixed_edge, free_edge))
					{
						if (fixed_edge.style == PCEEdge::VerticalEnd && free_edge.style == PCEEdge::VerticalEnd)
						{
							comb_edges.emplace_back(fixed_edge.index, free_edge.index);
							score++;
						}
					}
				}
			}

		}
	}

	if (false)
	{
		if (!fixed_positions.empty())
		{
			for (auto free_edge_index : free_position.edge_indexes)
			{
				double minimal_dis = Math::MAXDOUBLE;

				for (auto& fixed_position : fixed_positions)
				{
					for (auto& fixed_edge_index : fixed_position.edge_indexes)
					{
						minimal_dis = min(Math::Functs::GetDistance(edges[free_edge_index].cut.offset_c + free_position.t, edges[fixed_edge_index].cut.offset_c), minimal_dis);
					}
				}
				total_length += minimal_dis;
			}
		}
		if (!free_position.edge_indexes.empty())
			total_length = total_length / free_position.edge_indexes.size();
	}
	else
	{
		Vector2d2  all_points;

		all_points.emplace_back(free_position.transform->UOrigin(free_position.t));//bounding box
		if (!fixed_positions.empty())
		{
			for (auto& fixed_position : fixed_positions)
				all_points.emplace_back(fixed_position.transform->UOrigin(fixed_position.t));//bounding box
		}

		Vector2d minimal_corner, maximal_corner;
		Math::Functs::GetBoundingBox(all_points, minimal_corner, maximal_corner);
		Vector2d size_ = maximal_corner - minimal_corner;

		total_length = size_[0] * size_[1];

		for (auto free_edge_index : free_position.edge_indexes)
		{
			auto free_edge_s = edges[free_edge_index].cut.offset_s + free_position.t;
			auto free_edge_e = edges[free_edge_index].cut.offset_e + free_position.t;

			double minimal_dis = Math::MAXDOUBLE;
			for (auto& edge : edges)
			{
				if (edge.part_index < 0)
				{
					minimal_dis = min(cgal_2d_distance_segment_segment
					(free_edge_s, free_edge_e, edge.cut.offset_s, edge.cut.offset_e), minimal_dis);
				}
			}
			boundary_distance += minimal_dis;
		}
		boundary_distance = boundary_distance / free_position.edge_indexes.size();
	}


	return score;
}

//TOTO: the optimal result should be produced...
void PCE::ArrangePartbyConstraints(const WL& wl, std::vector<PCEEdge>& edges, const std::vector<PCEPosition>& fixed_positions,
	PCEPosition& free_position,
	std::vector<PCEConstraint>& sharing_edges)
{
	auto BuildConstraints = [](const WL& wl_,  const std::vector<PCEEdge>& edges_,
		const std::vector<PCEPosition>& fixed_positions_, const PCEPosition& free_position_,
		const std::vector<PCEConstraint>& sharing_edges_, std::vector<PCEConstraint>& constraints_)
	{
		//get saving cutting number
		std::vector<int> edge_labels(edges_.size(), -1);

		for (auto edge : sharing_edges_)
		{
			const auto b0 = edge_labels[edge.edge_0] == -1;
			const auto b1 = edge_labels[edge.edge_1] == -1;

			if (b0 && b1)
			{
				edge_labels[edge.edge_0] = edge.edge_0;
				edge_labels[edge.edge_1] = edge.edge_0;
			}
			if (b0 && !b1)
			{
				edge_labels[edge.edge_0] = edge_labels[edge.edge_1];
			}
			if (!b0 && b1)
			{
				edge_labels[edge.edge_1] = edge_labels[edge.edge_0];
			}
		}



		std::vector<int> fixed_edges;
		for (const auto fixed_position : fixed_positions_)
		{
			//auto& fixed_part = parts[fixed_position.part_index];
			for (auto edge : fixed_position.edge_indexes)fixed_edges.emplace_back(edge);
		}

		if (wl_.index >= 0) for (auto i : wl_.edges) fixed_edges.emplace_back(i);

		for (auto fixed_edge_index : fixed_edges)
		{
			const PCEEdge& fixed_edge = edges_[fixed_edge_index];

			if (edge_labels[fixed_edge.index] == -1 || edge_labels[fixed_edge.index] == fixed_edge.index)
			{
				for (auto free_edge_index : free_position_.edge_indexes)
				{
					const PCEEdge& free_edge = edges_[free_edge_index];
					if (fixed_edge.style == free_edge.style && free_edge.style == PCEEdge::VerticalEnd)
					{
						const double angle = Math::Functs::GetAngleBetween(fixed_edge.cut.offset_n, free_edge.cut.offset_n);
						if (Math::Functs::IsAlmostZero(angle) || Math::Functs::IsAlmostZero(angle - Math::Math_PI))
							constraints_.emplace_back(fixed_edge.index, free_edge.index);
					}
					if (fixed_edge.style == free_edge.style && free_edge.style == PCEEdge::TiltingEnd)
					{
						if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(fixed_edge.cut.cutting_surface_normal, free_edge.cut.cutting_surface_normal) - Math::Math_PI))
							constraints_.emplace_back(fixed_edge.index, free_edge.index);
					}
					/*	if (fixed_edge.style != free_edge.style)
						{
							const double angle = Math::Functions::GetAngleBetween(fixed_edge.cut.offset_n, free_edge.cut.offset_n);
							if (Math::Functions::IsAlmostZero(angle) || Math::Functions::IsAlmostZero(angle - Math::Math_PI))
								constraints_.emplace_back(fixed_edge.index, free_edge.index);
						}*/
				}
			}
		}

		edge_labels.clear();
		std::vector<int>().swap(edge_labels);
	};

	auto ValidDetection = [](const WL& wl_, 
		const std::vector<PCEPosition>& fixed_positions_, const PCEPosition& free_position_,
		HMODULE& hModule_, const double& part_match_error_)
	{
		auto cgal_2d_polygons_intersection = (CGAL_2D_Two_Polygons_Intersection)GetProcAddress(hModule_, "CGAL_2D_Two_Polygons_Intersection");
		auto cgal_2d_location_point_polygon = (CGAL_2D_Location_Point_Polygon)GetProcAddress(hModule_, "CGAL_2D_Location_Point_Polygon");
		auto cgal_2d_distance_point_polygon = (CGAL_2D_Distance_Point_Polygon)GetProcAddress(hModule_, "CGAL_2D_Distance_Point_Polygon");

		Vector2d1 free_u_offset = free_position_.transform->UOffset(free_position_.t);//collision
		Vector2d1 free_l_offset = free_position_.transform->LOffset(free_position_.t);//collision

		//wl
		if (wl_.index >= 0)
		{
			for (const auto p : free_u_offset)
			{
				double dis = cgal_2d_distance_point_polygon(p, wl_.offset_points);

				if (dis > part_match_error_)
				{
					if (!cgal_2d_location_point_polygon(p, wl_.offset_points))
					{
						return false;
					}
				}
			}

			for (const auto p : free_l_offset)
			{
				double dis = cgal_2d_distance_point_polygon(p, wl_.offset_points);

				if (dis > part_match_error_)
				{
					if (!cgal_2d_location_point_polygon(p, wl_.offset_points))
					{
						return false;
					}
				}
			}
		}

		//fixed position
		for (auto position : fixed_positions_)
		{
			Vector2d1 fixed_u_offset = position.transform->UOffset(position.t);//collision
			Vector2d1 fixed_l_offset = position.transform->LOffset(position.t);//collision

			double inter_area = cgal_2d_polygons_intersection(fixed_u_offset, free_u_offset);
			if (inter_area > part_match_error_)
			{
				return false;
			}

			inter_area = cgal_2d_polygons_intersection(fixed_l_offset, free_l_offset);
			if (inter_area > part_match_error_)
			{
				return false;
			}
		}

		return true;
	};

	auto SolveConstraints = [](std::vector<std::pair<PCEEdge, PCEEdge>>& vecPairs, Vector2d& result)
	{
		const int dimension = vecPairs.size();
		Eigen::MatrixXd A(dimension, dimension);
		Eigen::VectorXd x(dimension);
		Eigen::VectorXd b(dimension);

		for (auto i = 0; i < vecPairs.size(); ++i)
		{
			const auto& pair = vecPairs[i];
			auto& ea = pair.first;
			auto& eb = pair.second;

			Vector2d nra, nrb, o1, o2;

			if (ea.style == PCEEdge::TiltingEnd && eb.style == PCEEdge::TiltingEnd)
			{
				if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(ea.cut.cutting_surface_normal, eb.cut.cutting_surface_normal) - Math::Math_PI))
				{
					nra = ea.cut.offset_n;
					nrb = eb.cut.offset_n;

					auto o1_ = (ea.cut.u_cutting_line_s + ea.cut.u_cutting_line_e) / 2.0;
					auto o2_ = (eb.cut.u_cutting_line_s + eb.cut.u_cutting_line_e) / 2.0;

					o1 = Vector2d(o1_[0], o1_[1]);
					o2 = Vector2d(o2_[0], o2_[1]);
				}
				else
				{
					nra = ea.cut.offset_n;
					nrb = eb.cut.offset_n;
					o1 = ea.cut.offset_c;
					o2 = eb.cut.offset_c;
				}
			}
			else
			{
				nra = ea.cut.offset_n;
				nrb = eb.cut.offset_n;
				o1 = ea.cut.offset_c;
				o2 = eb.cut.offset_c;
			}

			Eigen::RowVectorXd tmpRow(dimension);
			tmpRow << nra[0], nra[1];
			A.row(i) = tmpRow;
			double d = (o1.x - o2.x) * nra.x + (o1.y - o2.y) * nra.y;
			b(i) = d;
		}

		Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A);
		if (lu_decomp.rank() != 2) return false;
		//x = lu_decomp.solve(b);
		x = A.colPivHouseholderQr().solve(b);
		result[0] = x[0];
		result[1] = x[1];

		return true;
	};

	auto InsertTs = [&](Vector2d1& ts, Vector2d t)
	{
		for (int i = 0; i < ts.size(); i++)
		{
			double d = Functs::GetLength(ts[i] - t);
			if (Functs::IsAlmostZero_Double(d, part_match_error))
				return false;
		}

		ts.emplace_back(t);
		return true;
	};

	//build constraints
	std::vector<PCEConstraint> constraints;
	BuildConstraints(wl, edges, fixed_positions, free_position, sharing_edges, constraints);

	//std::cerr << constraints.size() << std::endl;
	//output
	if (false)
	{
		std::vector<std::pair<Vector3d, Vector3d>> poly;
		for (auto& edge : edges)
		{
			poly.emplace_back(Math::Functs::Vector2d3d(edge.cut.offset_s), Math::Functs::Vector2d3d(edge.cut.offset_e));
		}

		Output3DSegments_Pairs(*hModule, kerf, packing_o_folder + "\\edges.obj", poly);
	}

	Vector2d1 ts;

	double maximal_score = 0.0;
	PCEPosition maximal_position;
	std::vector<PCEConstraint> optimal_comb_edges;
	for (int i = 0; i < constraints.size(); i++)
	{
		for (int j = i; j < constraints.size(); j++)
		{
			if (i != j)
			{
				std::vector<std::pair<PCEEdge, PCEEdge>> vecPairs;
				vecPairs.emplace_back(edges[constraints[i].edge_0], edges[constraints[i].edge_1]);
				vecPairs.emplace_back(edges[constraints[j].edge_0], edges[constraints[j].edge_1]);
				if (!SolveConstraints(vecPairs, free_position.t))
				{
					vecPairs.clear();
					std::vector<std::pair<PCEEdge, PCEEdge>>().swap(vecPairs);
					continue;
				}


				vecPairs.clear();
				std::vector<std::pair<PCEEdge, PCEEdge>>().swap(vecPairs);

				if (!ValidDetection(wl, fixed_positions, free_position, *hModule, part_match_error)) continue;

				if (!InsertTs(ts, free_position.t))continue;

				std::vector<PCEConstraint> comb_edges;
				double total_score_length = 0.0;
				double boundary_distance = 0.0;

				const int score = GetScore(edges, fixed_positions, free_position, constraints,
					comb_edges, total_score_length, boundary_distance);
				//double d = GetDistance(fixed_positions, free_position);

				double length_bb = (1.0 - total_score_length / (wl.size[0] * wl.size[1] * 2.0)) * 0.5;

				if (length_bb < 0.0 || length_bb>1.0)
				{
					std::cerr << "double length = total_score_length / (wl.size[0] * wl.size[1]*wl.size[2]*5);" << std::endl;
					continue;
					//system("pause");
				}

				double length_bd = (1.0 - boundary_distance / ((wl.size[0] + wl.size[1]) * 2.0)) * 0.5;

				if (length_bd < 0.0 || length_bd>1.0)
				{
					std::cerr << "double length_bd = (1.0 - boundary_distance / ((wl.size[0] + wl.size[1])*2.0)) * 0.5;" << std::endl;
					continue;
					//system("pause");
				}

				auto cur_score = score + length_bb + length_bd;
				if (cur_score > maximal_score)
				{
					maximal_score = cur_score;
					maximal_position = free_position;
					optimal_comb_edges = comb_edges;
				}
				else
				{
					if (Functs::IsAlmostZero_Double(cur_score-maximal_score, part_match_error)&& maximal_position.t[0] > free_position.t[0])
					{
						maximal_score = cur_score;
						maximal_position = free_position;
						optimal_comb_edges = comb_edges;
					}
				}

				comb_edges.clear();
				std::vector<PCEConstraint>().swap(comb_edges);
			}
		}
	}

	if (maximal_score > 0)
	{
		free_position = maximal_position;
		free_position.style = PCEPosition::FixedStyle_2_Constraints;
		free_position.ComputeM();
		free_position.UpdateEdges(edges);
		free_position.score = maximal_score;

		for (auto a : optimal_comb_edges) sharing_edges.emplace_back(a);
	}
	else
	{
		ArrangePartbySlider(wl, edges, constraints, fixed_positions, free_position, sharing_edges);
	}

	constraints.clear();
	optimal_comb_edges.clear();
	std::vector<PCEConstraint>().swap(constraints);
	std::vector<PCEConstraint>().swap(optimal_comb_edges);
}

//this is a bit too complex
void PCE::ArrangePartbySlider(const WL& wl, std::vector<PCEEdge>& edges, const std::vector<PCEConstraint>& constraints,
	const std::vector<PCEPosition>& fixed_positions, PCEPosition& free_position,
	std::vector<PCEConstraint>& sharing_edges)
{
	auto Comp = [](const Vector2d& a, const Vector2d& b) { return a[0] < b[0]; };

	auto OriginOffset = [](const HMODULE& hModule_, const LocalCS& local_cs, const  PCEPosition& free_position_,
		const PCEEdge& fixed_edge_, const PCEEdge& free_edge_,
		Vector2d1& origin_transform_points_, Vector2d1& offset_transform_points_, Vector2d& transform_t)
	{
		auto cgal_2d_distance_point_line = (CGAL_2D_Distance_Point_Line)GetProcAddress(
			hModule_, "CGAL_2D_Distance_Point_Line");
		double distance = cgal_2d_distance_point_line(free_edge_.cut.offset_c, fixed_edge_.cut.offset_s, fixed_edge_.cut.offset_e);
		double angle = Math::Functs::GetAngleBetween(fixed_edge_.cut.offset_n, free_edge_.cut.offset_c - fixed_edge_.cut.offset_c);
		transform_t = angle < M_PI / 2.0 ? -fixed_edge_.cut.origin_n * distance : fixed_edge_.cut.origin_n * distance;
		origin_transform_points_ = free_position_.transform->ULOrigin(transform_t);
		offset_transform_points_ = free_position_.transform->ULOffset(transform_t);

		origin_transform_points_ = local_cs.GetLocal(origin_transform_points_);
		offset_transform_points_ = local_cs.GetLocal(offset_transform_points_);
	};

	double maximal_score = 0.0;
	PCEPosition maximal_position;
	std::vector<PCEConstraint> optimal_comb_edges;

	Vector2d1 ts;

	for (auto constraint : constraints)
	{
		//build local coordinate system
		auto& fixed_edge = edges[constraint.edge_0];
		LocalCS local_cs(fixed_edge.cut.offset_c, fixed_edge.cut.offset_e - fixed_edge.cut.offset_s);

		PCEEdge& free_edge = edges[constraint.edge_1];

		Vector2d1 origin_transform_points;
		Vector2d1 offset_transform_points;
		Vector2d transform_t;
		OriginOffset(*hModule, local_cs, free_position, fixed_edge, free_edge, origin_transform_points, offset_transform_points, transform_t);

		//detect inside / outside
		Vector2d origin_min_bb, origin_max_bb;
		Math::Functs::GetBoundingBox(origin_transform_points, origin_min_bb, origin_max_bb);
		local_cs.origin_max_y = origin_max_bb[1];
		local_cs.origin_min_y = origin_min_bb[1];
		local_cs.origin_min_x = origin_min_bb[0];
		local_cs.origin_max_x = origin_max_bb[0];

		Vector2d offset_min_bb, offset_max_bb;
		Math::Functs::GetBoundingBox(offset_transform_points, offset_min_bb, offset_max_bb);
		local_cs.offset_max_y = offset_max_bb[1];
		local_cs.offset_min_y = offset_min_bb[1];
		local_cs.offset_min_x = offset_min_bb[0];
		local_cs.offset_max_x = offset_max_bb[0];

		if (wl.index >= 0)
		{
			Vector2d wood_min_bb, wood_max_bb;
			Math::Functs::GetBoundingBox(local_cs.GetLocal(wl.origin_points), wood_min_bb, wood_max_bb);
			if (!(local_cs.origin_min_y >= wood_min_bb[1] - part_match_error &&
				local_cs.origin_max_y <= wood_max_bb[1] + part_match_error))
				continue;
		}

		if (wl.index < 0)
		{
			std::cerr << "if (wl.index < 0)" << std::endl;
			system("pause");
		}

		//////////////////

		//projecting
		Vector2d1 proj_segs;
		Vector2d1 proj_kerf_segs;
		if (wl.index >= 0)
		{
			for (auto edge_index : wl.edges)
			{
				PCEEdge& fixed_edge_ = edges[edge_index];
				Vector2d seg_s = local_cs.GetLocal(fixed_edge_.cut.origin_s);
				Vector2d seg_e = local_cs.GetLocal(fixed_edge_.cut.origin_e);
				double min_x;
				double max_x;
				if (local_cs.GetSegOrigin(*hModule,part_match_error, seg_s, seg_e, min_x, max_x))
				{
					proj_segs.emplace_back(min_x, max_x);
					proj_kerf_segs.emplace_back(min_x, max_x);
				}
			}
		}

		std::sort(proj_segs.begin(), proj_segs.end(), Comp);

		Vector2d1 wool_empty_segs = local_cs.EmptyX(local_cs.SortX(proj_segs));
		if (wool_empty_segs.size() != 1)
		{
			continue;
		}

		Vector2d wd_proj_seg;
		if (wool_empty_segs.size() == 1)
			wd_proj_seg = wool_empty_segs[0];

		//origin
		for (auto fixed_position : fixed_positions)
		{
			Vector2d1 ps;
			for (auto edge_index : fixed_position.edge_indexes)
			{
				PCEEdge& fixed_edge_ = edges[edge_index];
				Vector2d seg_s = local_cs.GetLocal(fixed_edge_.cut.origin_s);
				Vector2d seg_e = local_cs.GetLocal(fixed_edge_.cut.origin_e);

				double min_x;
				double max_x;
				if (local_cs.GetSegOrigin(*hModule,part_match_error, seg_s, seg_e, min_x, max_x))
				{
					ps.emplace_back(min_x, max_x);
				}
			}

			if (!ps.empty())
			{
				std::sort(ps.begin(), ps.end(), Comp);
				proj_segs.emplace_back(ps[0][0], ps[ps.size() - 1][1]);
			}
		}

		//offset
		for (auto fixed_position : fixed_positions)
		{
			Vector2d1 ps;
			for (auto edge_index : fixed_position.edge_indexes)
			{
				PCEEdge& fixed_edge = edges[edge_index];
				Vector2d seg_s = local_cs.GetLocal(fixed_edge.cut.offset_s);
				Vector2d seg_e = local_cs.GetLocal(fixed_edge.cut.offset_e);

				double min_x;
				double max_x;
				if (local_cs.GetSegOffset(*hModule,part_match_error, seg_s, seg_e, min_x, max_x))
				{
					ps.emplace_back(min_x, max_x);
				}
			}

			if (!ps.empty())
			{
				std::sort(ps.begin(), ps.end(), Comp);
				proj_kerf_segs.emplace_back(ps[0][0], ps[ps.size() - 1][1]);
			}
		}


		if (proj_segs.empty())continue;
		if (proj_kerf_segs.empty())continue;

		//sort
		std::sort(proj_segs.begin(), proj_segs.end(), Comp);
		std::sort(proj_kerf_segs.begin(), proj_kerf_segs.end(), Comp);

		Vector2d1 empty_segs = local_cs.EmptyX(local_cs.SortX(proj_segs));
		Vector2d1 empty_kerf_segs = local_cs.EmptyX(local_cs.SortX(proj_kerf_segs));

		if (empty_segs.empty())continue;
		if (empty_kerf_segs.empty())continue;

		bool min_end_style = Math::Functs::IsAlmostZero_Double(abs(local_cs.origin_min_x - local_cs.offset_min_x) - kerf, part_match_error); // 0: TiltingEnd; 1:VerticalEnd;
		bool max_end_style = Math::Functs::IsAlmostZero_Double(abs(local_cs.origin_max_x - local_cs.offset_max_x) - kerf, part_match_error); // 0: TiltingEnd; 1:VerticalEnd;

		//get a place to find an optimal one
		for (int i = 0; i < empty_segs.size(); i++)
		{
			auto& empty_seg = empty_segs[i];

			Vector2d empty_kerf_seg;

			bool goon = false;
			for (auto& seg : empty_kerf_segs)
			{
				if (empty_seg[0] <= seg[0] && empty_seg[1] >= seg[1])
				{
					goon = true;
					empty_kerf_seg = seg;
				}
			}
			if (!goon) continue;


			if (!(empty_kerf_seg[0] >= empty_seg[0] && empty_kerf_seg[1] <= empty_seg[1]))
			{
				std::cerr << "if (!(empty_kerf_seg[0] > empty_seg[0]&& empty_kerf_seg[1]< empty_seg[1]))" << std::endl;
				system("pause");
			}

			bool b0 = Math::Functs::IsAlmostZero_Double(empty_seg[0] - wd_proj_seg[0],part_match_error);
			bool b1 = Math::Functs::IsAlmostZero_Double(empty_seg[1] - wd_proj_seg[1],part_match_error);

			if (b0 && b1)
			{
				//empty_seg[0] empty_seg[1]
				//local_cs.origin_min_x local_cs.origin_max_x
				double length = empty_seg[1] - empty_seg[0];

				if (min_end_style)
				{
					if (length > local_cs.origin_max_x - local_cs.origin_min_x || Math::Functs::IsAlmostZero_Double(length - local_cs.origin_max_x + local_cs.origin_min_x,part_match_error))
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[0] - local_cs.origin_min_x, 0.0)));
					}
				}
				else
				{
					//different
					double d0 = local_cs.origin_min_x + empty_seg[0] - local_cs.offset_min_x - empty_seg[0];
					double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;

					double l = (length - d1) - (local_cs.origin_max_x - local_cs.offset_min_x + safe_distance);
					if (l > 0 || Math::Functs::IsAlmostZero_Double(l,part_match_error))
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[0] - local_cs.offset_min_x + d1 + safe_distance, 0.0)));
					}
				}

				if (max_end_style)
				{
					if (length > local_cs.origin_max_x - local_cs.origin_min_x || Math::Functs::IsAlmostZero_Double(length - local_cs.origin_max_x + local_cs.origin_min_x,part_match_error))
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[1] - local_cs.origin_max_x, 0.0)));
					}
				}
				else
				{
					//different
					double d0 = empty_seg[1] - (local_cs.origin_max_x + empty_seg[1] - local_cs.offset_max_x);
					double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;
					double l = (length - d1) - (local_cs.offset_max_x - local_cs.origin_min_x + safe_distance);
					if (l > 0 || Math::Functs::IsAlmostZero_Double(l,part_match_error))
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[1] - local_cs.offset_max_x + d1 + safe_distance, 0.0)));
					}
				}
			}

			if (!b0 && b1)
			{
				double d0 = (local_cs.origin_min_x + empty_kerf_seg[0] - local_cs.offset_min_x) - empty_seg[0];
				double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;

				if (empty_seg[1] - empty_kerf_seg[0] - d1 > local_cs.origin_max_x - local_cs.offset_min_x + safe_distance)
				{
					ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_kerf_seg[0] + d1 - local_cs.offset_min_x + safe_distance, 0.0)));
				}

				if (max_end_style)
				{
					if (empty_seg[1] - empty_kerf_seg[0] - d1 > local_cs.origin_max_x - local_cs.offset_min_x + safe_distance)
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[1] - local_cs.origin_max_x, 0.0)));
					}
				}
				else
				{
					//different
					double d0 = empty_seg[1] - (local_cs.origin_max_x + empty_seg[1] - local_cs.offset_max_x);
					double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;
					double l = (empty_seg[1] - empty_kerf_seg[0] - d1) - (local_cs.offset_max_x - local_cs.origin_min_x + safe_distance);
					if (l > 0 || Math::Functs::IsAlmostZero_Double(l,part_match_error))
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[1] - local_cs.offset_max_x + d1 + safe_distance, 0.0)));
				}
			}

			if (b0 && !b1)
			{
				//empty_seg[0]   empty_kerf_seg[1]
				//local_cs.origin_min_x   local_cs.offset_max_x
				double d2 = empty_seg[1] - (local_cs.origin_max_x + empty_kerf_seg[1] - local_cs.offset_max_x);
				double d3 = ((int)(d2 / kerf) + 1) * kerf - d2;

				if (min_end_style)
				{
					if (empty_kerf_seg[1] - empty_seg[0] - d3 > local_cs.offset_max_x - local_cs.origin_min_x + safe_distance)
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[0] - local_cs.origin_min_x, 0.0)));
					}
				}
				else
				{
					//different
					double d0 = local_cs.origin_min_x + empty_seg[0] - local_cs.offset_min_x - empty_seg[0];
					double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;

					double l = (empty_kerf_seg[1] - empty_seg[0] - d1) - (local_cs.origin_max_x - local_cs.offset_min_x + safe_distance);
					if (l > 0 || Math::Functs::IsAlmostZero_Double(l,part_match_error))
					{
						ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_seg[0] - local_cs.offset_min_x + d1 + safe_distance, 0.0)));
					}
				}


				if (empty_kerf_seg[1] - empty_seg[0] - d3 > local_cs.offset_max_x - local_cs.origin_min_x + safe_distance)
				{
					ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_kerf_seg[1] - d3 - local_cs.offset_max_x - safe_distance, 0.0)));
				}
			}

			if (!b0 && !b1)
			{
				double d0 = (local_cs.origin_min_x + empty_kerf_seg[0] - local_cs.offset_min_x) - empty_seg[0];
				double d1 = ((int)(d0 / kerf) + 1) * kerf - d0;

				if (empty_kerf_seg[1] - empty_kerf_seg[0] - d1 > local_cs.offset_max_x - local_cs.offset_min_x + safe_distance * 2.0)
				{
					ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_kerf_seg[0] + d1 - local_cs.offset_min_x + safe_distance, 0.0)));
				}

				double d2 = empty_seg[1] - (local_cs.origin_max_x + empty_kerf_seg[1] - local_cs.offset_max_x);
				double d3 = ((int)(d2 / kerf) + 1) * kerf - d2;
				if (empty_kerf_seg[1] - empty_kerf_seg[0] - d3 > local_cs.offset_max_x - local_cs.offset_min_x + safe_distance * 2.0)
				{
					ts.emplace_back(transform_t + local_cs.GetGlobalVec(Vector2d(empty_kerf_seg[1] - d3 - local_cs.offset_max_x - safe_distance, 0.0)));
				}
			}
		}

		empty_segs.clear();
		Vector2d1().swap(empty_segs);
	}

	for (int iter = 0; iter < ts.size(); iter++)
		//for (auto& t : ts)
	{
		free_position.t = ts[iter];

		/////////////////////////////////////////////
		std::vector<PCEConstraint> comb_edges;
		double d = 0.0;
		double boundary_distance = 0.0;
		int score = GetScore(edges, fixed_positions, free_position, constraints, comb_edges, d, boundary_distance);
		double length_bb = (1.0 - d / (wl.size[0] * wl.size[1] * 2.0)) * 0.5;

		if (length_bb < 0.0 || length_bb>1.0)
		{
			std::cerr << "double length = total_score_length / (wl.size[0] * wl.size[1]*wl.size[2]*5);" << std::endl;
			//system("pause");
			continue;
			//int score_0 = GetScore(edges, fixed_positions, free_position, constraints, comb_edges, d, boundary_distance);
		}

		double length_bd = (1.0 - boundary_distance / ((wl.size[0] + wl.size[1]) * 2.0)) * 0.5;

		if (length_bd < 0.0 || length_bd>1.0)
		{
			std::cerr << "double length_bd = (1.0 - boundary_distance / ((wl.size[0] + wl.size[1])*2.0)) * 0.5;" << std::endl;
			system("pause");
		}

		if (score + length_bb + length_bd > maximal_score)
		{
			maximal_score = score + length_bb + length_bd;
			maximal_position = free_position;
			optimal_comb_edges = comb_edges;
		}
		comb_edges.clear();
		std::vector<PCEConstraint>().swap(comb_edges);
	}

	if (maximal_score > 0)
	{
		free_position = maximal_position;
		free_position.style = PCEPosition::FixedStyle_1_Constraint;
		free_position.ComputeM();
		free_position.UpdateEdges(edges);
		free_position.score = maximal_score;

		for (auto a : optimal_comb_edges) sharing_edges.emplace_back(a);
	}
	else
	{
		if (constraints.empty())
		{
			//bounding box
			auto edges_ = edges;
			auto sharing_edges_ = sharing_edges;
			auto edge_indexes_ = free_position.edge_indexes;

			edges.erase(edges.begin() + free_position.edge_indexes[0], edges.begin() + edges.size());
			free_position.transform->AddBoundingEdges(edges, free_position.edge_indexes);
			ArrangePartbyConstraints(wl, edges, fixed_positions, free_position, sharing_edges);

			edges = edges_;
			sharing_edges = sharing_edges_;
			free_position.edge_indexes = edge_indexes_;

			free_position.UpdateEdges(edges);
		}
	}

	optimal_comb_edges.clear();
	std::vector<PCEConstraint>().swap(optimal_comb_edges);
}

void PCE::GetCuttingLines(const std::vector<PCEEdge>& edges,
	const std::vector<PCEConstraint>& sharing_edges,
	const std::vector<PCEPosition>& fixed_positions,
	std::vector<PCEEdge>& cutting_lines)
{
	auto Comp = [](const Vector2d& a, const Vector2d& b) { return a[0] < b[0]; };

	auto CollisionDetection = [](const HMODULE& hModule_, const std::vector<PCEPosition>& fixed_positions_,
		const PCEEdge& fixed_edge_, const PCEEdge& free_edge_)
	{
		//CGAL_2D_Intersection_Segment_Polygon
		const auto intersection = (CGAL_2D_Intersection_Segment_Polygon)GetProcAddress(hModule_, "CGAL_2D_Intersection_Segment_Polygon");

		for (auto& position : fixed_positions_)
		{
			Vector2d1 fixed_u_origin = position.transform->UOrigin(position.t);//collision
			Vector2d1 fixed_l_origin = position.transform->LOrigin(position.t);//collision

			if (intersection(Vector2d(fixed_edge_.cut.u_cutting_line_s), Vector2d(free_edge_.cut.u_cutting_line_e), fixed_u_origin))
				return true;
			if (intersection(Vector2d(fixed_edge_.cut.l_cutting_line_s), Vector2d(free_edge_.cut.l_cutting_line_e), fixed_l_origin))
				return true;
		}
		return false;
	};

	auto CombineEdges = [](const std::vector<PCEEdge>& edges_, const std::vector<PCEConstraint>& sharing_edges_,
		std::vector<int>& unique_l_, std::vector<std::vector<int>>& unique_ls_)
	{
		std::vector<int> edge_labels(edges_.size(), -1);

		for (auto edge : sharing_edges_)
		{
			bool b0 = edge_labels[edge.edge_0] == -1;
			bool b1 = edge_labels[edge.edge_1] == -1;

			if (b0 && b1)
			{
				edge_labels[edge.edge_0] = edge.edge_0;
				edge_labels[edge.edge_1] = edge.edge_0;
			}

			if (b0 && !b1)
			{
				edge_labels[edge.edge_0] = edge_labels[edge.edge_1];
			}

			if (!b0 && b1)
			{
				edge_labels[edge.edge_1] = edge_labels[edge.edge_0];
			}
		}

		for (int i = 0; i < edge_labels.size(); i++)
		{
			int index = Math::Functs::VectorIndex(unique_l_, edge_labels[i]);
			if (index < 0)
			{
				unique_l_.emplace_back(edge_labels[i]);
				unique_ls_.emplace_back(1, i);
			}
			else
			{
				unique_ls_[index].emplace_back(i);
			}
		}
	};


	std::vector<int> full_unique_l;
	std::vector<std::vector<int>> full_unique_ls;
	CombineEdges(edges, sharing_edges, full_unique_l, full_unique_ls);

	std::vector<PCEConstraint> sharing_edges_;
	for (int i = 0; i < full_unique_ls.size(); i++)
	{
		if (full_unique_l[i] >= 0)
		{
			if (edges[full_unique_l[i]].part_index < 0)
			{
				if (full_unique_ls[i].size() > 1)
					for (int j = 0; j < full_unique_ls[i].size() - 1; j++)
					{
						sharing_edges_.emplace_back(PCEConstraint(full_unique_ls[i][j], full_unique_ls[i][j + 1]));
					}
			}
			else
			{
				if (full_unique_ls[i].size() > 1)
					for (int j = 0; j < full_unique_ls[i].size() - 1; j++)
						for (int k = j + 1; k < full_unique_ls[i].size(); k++)
						{
							if (!CollisionDetection(*hModule, fixed_positions, edges[full_unique_ls[i][j]], edges[full_unique_ls[i][k]]))
							{
								sharing_edges_.emplace_back(PCEConstraint(full_unique_ls[i][j], full_unique_ls[i][k]));
							}
						}
			}
		}
	}

	std::vector<int> unique_l;
	std::vector<std::vector<int>> unique_ls;
	CombineEdges(edges, sharing_edges_, unique_l, unique_ls);

	//compute cutting lines
	for (int i = 0; i < unique_ls.size(); i++)
	{
		if (unique_l[i] >= 0)
		{
			if (edges[unique_l[i]].part_index < 0)
			{
				for (int j = 0; j < unique_ls[i].size(); j++)
				{
					const PCEEdge& edge = edges[unique_ls[i][j]];
					if (edge.part_index < 0)continue;
					if (edge.style == PCEEdge::VerticalEnd) continue;
					//tilting cutting 
					cutting_lines.emplace_back(PCEEdge(edge.cut.l_cutting_line_s, edge.cut.l_cutting_line_e, edge.cut.cutting_surface_normal, edge.cut.l_u_vector));
					cutting_lines.back().index = cutting_lines.size() - 1;
				}
				continue;
			}
		}

		if (unique_ls[i].size() == 1)
		{
			const PCEEdge& edge = edges[unique_l[i]];
			if (edge.part_index >= 0)
			{
				Vector2d edge_s = edge.cut.offset_s;
				Vector2d edge_e = edge.cut.offset_e;
				cutting_lines.emplace_back(PCEEdge(edge.cut.l_cutting_line_s, edge.cut.l_cutting_line_e, edge.cut.cutting_surface_normal, edge.cut.l_u_vector));
				cutting_lines.back().index = cutting_lines.size() - 1;
			}
		}
		else
		{
			if (unique_l[i] < 0)
			{
				for (int j = 0; j < unique_ls[i].size(); j++)
				{
					if (edges[unique_ls[i][j]].part_index < 0)continue;

					const PCEEdge& edge = edges[unique_ls[i][j]];
					//tilting cutting 
					cutting_lines.emplace_back(PCEEdge(edge.cut.l_cutting_line_s, edge.cut.l_cutting_line_e, edge.cut.cutting_surface_normal, edge.cut.l_u_vector));
					cutting_lines.back().index = cutting_lines.size() - 1;
				}
			}
			else
			{
				const PCEEdge& cs_edge = edges[unique_ls[i][0]];
				LocalCS local_cs(cs_edge.cut.offset_c, cs_edge.cut.offset_e - cs_edge.cut.offset_s);

				//double delta_z=0.0;

				Vector2d1 ps;
				Vector3d1 ps_vectors;
				Vector3d1 lu_vectors;
				for (int j = 0; j < unique_ls[i].size(); j++)
				{
					const PCEEdge& edge = edges[unique_ls[i][j]];

					ps.emplace_back(local_cs.GetLocal(Vector2d(edge.cut.l_cutting_line_s[0], edge.cut.l_cutting_line_s[1])));
					ps.emplace_back(local_cs.GetLocal(Vector2d(edge.cut.l_cutting_line_e[0], edge.cut.l_cutting_line_e[1])));
					ps_vectors.emplace_back(edge.cut.cutting_surface_normal);
					lu_vectors.emplace_back(edge.cut.l_u_vector);
					//delta_z = edge.cutting_line_e[2];
				}

				if (!ps.empty())
				{
					std::sort(ps.begin(), ps.end(), Comp);
					Vector2d s = local_cs.GetGlobalPos(ps[0]);
					Vector2d e = local_cs.GetGlobalPos(ps[ps.size() - 1]);

					cutting_lines.emplace_back(PCEEdge(Vector3d(s[0], s[1], 0.0), Vector3d(e[0], e[1], 0.0), ps_vectors[0], lu_vectors[0]));
					cutting_lines.back().index = cutting_lines.size() - 1;
				}
			}
		}
	}

	unique_ls.clear();
	std::vector<std::vector<int>>().swap(unique_ls);


	auto SortCuttingLine = [&](PCEEdge& cutting_line_0, PCEEdge& cutting_line_1)
	{
		auto s_0 = cutting_line_0.cut.l_cutting_line_s;
		auto s_1 = cutting_line_1.cut.l_cutting_line_s;

		if ( s_0[0]<s_1[0])
		{
			return true;
		}
		else
		{
			if (Math::Functs::IsAlmostZero_Double(s_0[0] - s_1[0], part_match_error))
			{
				if (s_0[1] < s_1[1])
				{
					return true;
				}
			}
		}
		return false;
	};


	std::sort(cutting_lines.begin(), cutting_lines.end(), SortCuttingLine);

	for (int i = 0; i < cutting_lines.size(); i++)
		cutting_lines[i].index = i;
}

void PCE::OutputFixedPositions(
	const WL& wl,
	const std::vector<PCEPosition> &positions, 
	const std::vector<PCEEdge>& cutting_lines,
	const std::string& path)
{
	auto Output = [](std::string path_, const HMODULE& hModule_, const std::vector<std::vector<Vector3d1 > >& points,
		const std::vector<string>& all_names,const Vector1i1& all_shapes,
		const Vector3d1& bound_vecs_, const std::vector<PCEEdge>& cutting_lines_, const double& kerf_, const  WL& wl_,  NumberLabel& number_label)
	{
		std::ofstream export_fie(path_);

		auto output = (CGAL_Export_Path_Segment)GetProcAddress(hModule_, "CGAL_Export_Path_Segment");
		int export_int = 1;
		for (auto iter = 0; iter < bound_vecs_.size(); iter++)
		{
			Vector3d s = bound_vecs_[iter];
			Vector3d e = bound_vecs_[(iter + 1) % bound_vecs_.size()];
			output(export_fie, export_int, "bound_" + Math::Functs::IntString(iter), 0.0, 0.0, 0.0, s, e, kerf_ * 0.5);
		}

		if(false)
		for (auto iter = 0; iter < cutting_lines_.size(); iter++)
		{
			double r, g, b;
			Math::Functs::ColorMapping((double)iter / (double)cutting_lines_.size(), r, g, b);

			Vector3d s = cutting_lines_[iter].cut.u_cutting_line_s;
			Vector3d e = cutting_lines_[iter].cut.u_cutting_line_e;

			output(export_fie, export_int, "cutting_" + Math::Functs::IntString(iter), r, g, b, s, e, kerf_ * 0.9);

			Vector3d normal = cutting_lines_[iter].cut.l_u_vector;
			//Vector3d cross_normal = Math::Functions::GetCrossproduct(normal, e - s);
			output(export_fie, export_int, "cutting_" + Math::Functs::IntString(iter),
				r, g, b, (s + e) / 2.0, (s + e) / 2.0 + Math::Functs::SetVectorLength(normal, 10.0 * kerf_), kerf_ * 0.9);
		}

		for (int i = 0; i < points.size(); i++)
		{
			auto part_center = Math::Functs::GetCenter(points[i]);
			part_center[2] = wl_.size[2] + 2.0 * kerf_;

			number_label.DrawLabel(
				NumberLabel::BaseSquare2,
				export_fie,
				export_int,
				all_names[i],
				all_shapes[i],
				kerf_ * 6.0,
				part_center);
		}


		for (int i = 0; i < points.size(); i++)
		{
			for (int j = 0; j < points[i].size(); j++)
			{
				for (int k = 0; k < points[i][j].size(); k++)
				{
					export_fie << fixed << setprecision(8) << "v " << points[i][j][k][0] << " " << points[i][j][k][1] << " " << points[i][j][k][2] << std::endl;
				}
			}
		}

		int nb = export_int;

		for (int i = 0; i < points.size(); i++)
		{
			export_fie << "g " << all_names[i] << std::endl;

			for (int j = 0; j < points[i].size(); j++)
			{
				export_fie << "f ";
				for (int k = 0; k < points[i][j].size(); k++)
				{
					export_fie << Math::Functs::IntString(nb) << " ";
					nb++;
				}
				export_fie << "" << std::endl;
			}

		}




		export_fie.clear();
		export_fie.close();
	};

	Vector3d3 all_points;

	//fixed positions
	std::vector<string> all_names;
	Vector1i1 all_shapes;
	for (auto position : positions)
	{
		all_points.emplace_back(Math::Functs::PosApplyM(position.transform->surfs, position.M));
		all_names.emplace_back("shape_" + Math::Functs::IntString(position.transform->pce_shape->index));
		all_shapes.emplace_back(position.transform->pce_shape->index);
	}

	//wl
	Vector3d1 bound_vecs;
	if (wl.index >= 0)
	{
		for (auto p : wl.origin_points)
			bound_vecs.emplace_back(p.x, p.y, wl.size[2]);
	}

	Output(path, *hModule, all_points, all_names, all_shapes, bound_vecs, cutting_lines, kerf, wl, number_label);

	all_points.clear();
	std::vector<std::vector<Vector3d1 > >().swap(all_points);
}

void PCE::OutputFixedPositions(const WL& wl, const shared_ptr<PCEArrange> arrange, const std::string& path)
{
	OutputFixedPositions(wl, arrange->positions, arrange->cutting_lines, path);
}

void PCE::OutputFixedPositionsVIS(
	const std::string& pref_str, 
	const WL& wl, 
	const shared_ptr<PCEArrange> arrange,
	const std::vector<PCEEdge>& cutting_lines,
	const Vector1i1& cutting_line_ids,
	const Vector3d1& cutting_line_colors,
	const std::vector<std::pair<Vector3d, Vector3d>>& refer_edges,
	const Vector3d2& refer_faces,
	const Vector3d1& refer_values,
	const Vector3d& delta_y, 
	std::ofstream& export_fie,int& export_int)
{
	auto Output = [](std::ofstream& export_file,
		int &export_int,
		const Vector3d1& refer_values_,
		const Vector3d& delta_y,
		const std::string& pref_str,
		const HMODULE& hModule_,
		const std::vector<std::vector<Vector3d1 > >& points,
		const std::vector<string>& all_names,
		const std::vector<int>& all_shapes_,
		const std::vector<PCEEdge>& cutting_lines_,
		const std::vector<int>& cutting_line_ids_,
		const Vector3d1& cutting_line_colors_,
		const std::vector<std::pair<Vector3d, Vector3d>>& refer_edges_,
		const Vector3d2& refer_faces_,
		const double& kerf_, const  WL& wl_,
		NumberLabel &number_label)
	{
		auto output = (CGAL_Export_Path_Segment)GetProcAddress(hModule_, "CGAL_Export_Path_Segment");
		auto output_point = (CGAL_Export_Path_Point)GetProcAddress(hModule_, "CGAL_Export_Path_Point");
		auto distance3d = (CGAL_3D_Distance_Point_Segment)GetProcAddress(hModule_, "CGAL_3D_Distance_Point_Segment");

		Vector3d1 bound_vecs_;

		//lumber
		if (false)
		{
			bound_vecs_ = wl_.boundary_3d[1];
			for (auto iter = 0; iter < bound_vecs_.size(); iter++)
			{
				Vector3d s = bound_vecs_[iter];
				Vector3d e = bound_vecs_[(iter + 1) % bound_vecs_.size()];
				output(export_file, export_int, pref_str + "_bound_lumber", 0.0, 0.0, 0.0, s + delta_y, e + delta_y, kerf_ * 0.3);
			}
		}

		//wood
		bound_vecs_ = wl_.boundary_3d[0];
		for (auto iter = 0; iter < bound_vecs_.size(); iter++)
		{
			Vector3d s = bound_vecs_[iter];
			Vector3d e = bound_vecs_[(iter + 1) % bound_vecs_.size()];
			output(export_file, export_int, pref_str + "_bound_wood", 0.0, 0.0, 0.0, s + delta_y, e + delta_y, kerf_ * 0.3);
		}

		int upper_refer_int = 0;
		int lower_refer_int = 0;
		for (auto iter = 0; iter < cutting_lines_.size(); iter++)
		{
			auto pref_str_name = pref_str + "_cutting_" + Math::Functs::IntString(cutting_line_ids_[iter]);
			double r, g, b;
			r = cutting_line_colors_[iter][0];
			g = cutting_line_colors_[iter][1];
			b = cutting_line_colors_[iter][2];

			Vector3d cut_s = cutting_lines_[iter].cut.u_cutting_line_s;
			Vector3d cut_e = cutting_lines_[iter].cut.u_cutting_line_e;

			if(cut_s[1] > cut_e[1]) std::swap(cut_s, cut_e);

			//double minimal_distance_s = 10000000.0;
			//double minimal_distance_e = 10000000.0;
			//for (auto iter_ = 0; iter_ < cutting_lines_.size(); iter_++)
			//{
			//	if (iter_ != iter)
			//	{
			//		Vector3d seg_s = cutting_lines_[iter_].cut.u_cutting_line_s;
			//		Vector3d seg_e = cutting_lines_[iter_].cut.u_cutting_line_e;
			//		double d_s = distance3d(s, seg_s, seg_e);
			//		double d_e = distance3d(e, seg_s, seg_e);
			//		minimal_distance_s = min(minimal_distance_s, d_s);
			//		minimal_distance_e = min(minimal_distance_e, d_e);
			//	}
			//}
			//Vector3d s_e = e - s;
			//Math::Functions::SetVectorLength(s_e, 3.5 * kerf_);
			//if (minimal_distance_s < 2.0 * kerf_) s = s + s_e;
			//if (minimal_distance_e < 2.0 * kerf_) e = e - s_e;

			output(export_file, export_int, pref_str_name,r, g, b, cut_s + delta_y, cut_e + delta_y, kerf_ * 0.5);

			//refer edges
			//////////////////////////////////////////////////////////////////////////////////////
			auto refer_length = refer_values_[iter] * 25.4;
			Vector3d refer_s(refer_edges_[iter].first[0], refer_edges_[iter].first[1], refer_edges_[iter].first[2]);
			Vector3d refer_e(refer_edges_[iter].second[0], refer_edges_[iter].second[1], refer_edges_[iter].second[2]);

			//refer_values_[iter];

			auto cut_2d_s = Math::Functs::Vector3d2d(cut_s);
			auto cut_2d_e = Math::Functs::Vector3d2d(cut_e);
			auto refer_2d_s = Math::Functs::Vector3d2d(refer_s);
			auto refer_2d_e = Math::Functs::Vector3d2d(refer_e);

			const auto cgal_2d_distance_point_line = (CGAL_2D_Distance_Point_Line)GetProcAddress(hModule_, "CGAL_2D_Distance_Point_Line");
			auto d_s = cgal_2d_distance_point_line(cut_2d_s, refer_2d_s, refer_2d_e);
			auto d_e = cgal_2d_distance_point_line(cut_2d_e, refer_2d_s, refer_2d_e);

			auto cut_2d_loc = d_s > d_e ? cut_2d_e : cut_2d_s;
			auto cut_3d_loc = d_s > d_e ? cut_e : cut_s;
			auto upper_height = d_s > d_e ? 1 : -1;

			d_s > d_e ? upper_refer_int++ : lower_refer_int++;

			d_s = abs(Math::Functs::GetDistance(cut_2d_loc, refer_2d_s) - refer_values_[iter][2]*25.4);
			d_e = abs(Math::Functs::GetDistance(cut_2d_loc, refer_2d_e) - refer_values_[iter][2] * 25.4);
			auto refer_3d_loc = d_s > d_e ? refer_e : refer_s;

			cut_3d_loc[1] = refer_3d_loc[1];
			cut_3d_loc[2] = cut_s[2];
			refer_3d_loc[2] = cut_s[2];

			Vector3d delta_y_refer(0.0, 0.0,0.0);

			if (upper_height == 1)
				delta_y_refer[1]= kerf_ * 4.0 *upper_refer_int;
			if (upper_height == -1)
				delta_y_refer[1] = -kerf_ * 4.0*lower_refer_int;


			auto refer_edge_center = (cut_3d_loc + refer_3d_loc) / 2.0;

			refer_edge_center[2] = refer_edge_center[2] + 2.0 * kerf_;

			number_label.DrawLabel(
				NumberLabel::BaseSquare1,
				export_file,
				export_int,
				pref_str + "_cutting_refer_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				cutting_line_ids_[iter],
				//3407,
				kerf_ * 4.0,
				refer_edge_center+delta_y + delta_y_refer);

			// cut_3d_loc => refer_3d_loc
			output(export_file,
				export_int,
				pref_str + "_cutting_refer_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				r, g, b,
				cut_3d_loc + delta_y+ delta_y_refer, refer_3d_loc + delta_y+ delta_y_refer, kerf_ * 0.4);

			//cut_3d_loc => delta_y_refer
			output(export_file,
				export_int,
				pref_str + "_cutting_refer_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				r, g, b,
				cut_3d_loc + delta_y, cut_3d_loc + delta_y + delta_y_refer, kerf_ * 0.4);

			//refer_3d_loc=>delta_y_refer
			output(export_file,
				export_int,
				pref_str + "_cutting_refer_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				r, g, b,
				refer_3d_loc + delta_y, refer_3d_loc + delta_y + delta_y_refer, kerf_ * 0.4);
			
			//stop
			output(export_file,
				export_int,
				pref_str + "_cutting_refer_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				r, g, b,
				refer_3d_loc + delta_y + delta_y_refer, refer_3d_loc + delta_y + delta_y_refer*2.0, kerf_ * 0.4);

			if(false)
			output(export_file,
				export_int, 
				pref_str + "_cutting_refer_edge_" + Math::Functs::IntString(cutting_line_ids_[iter]),
				r, g, b, 
				refer_s + delta_y, refer_e + delta_y, kerf_ * 1.2);

			//refer faces
			if(false)
			for (auto iter_ = 0; iter_ < refer_faces_[iter].size(); iter_++)
			{
				auto face_s = refer_faces_[iter][iter_];
				auto face_e = refer_faces_[iter][(iter_+1)% refer_faces_[iter].size()];
				output(export_file,
					export_int,
					pref_str + "_cutting_refer_face_" + Math::Functs::IntString(cutting_line_ids_[iter]),
					r, g, b,
					face_s + delta_y, face_e + delta_y, kerf_ * 1.2);
			}

			///////////////////////////////
			auto label_center = (cut_s + delta_y + cut_e + delta_y) / 2.0;
			label_center[2] = label_center[2] + 2.0 * kerf_;

			number_label.DrawLabel(
				NumberLabel::BaseCircle,
				export_file,
				export_int, 
				pref_str_name,
				cutting_line_ids_[iter], 
				//3407,
				kerf_ * 6.0, 
				label_center);

			///////////////////////////////

			//number_label.DrawNode(NumberLabel::BaseShere, export_fie, export_int, pref_str_name, kerf_ *0.5, cut_s + delta_y);
			//number_label.DrawNode(NumberLabel::BaseShere, export_fie, export_int, pref_str_name, kerf_ * 0.5, cut_e + delta_y);
		}

		for (auto iter = 0; iter < cutting_lines_.size(); iter++)
		{
			double r, g, b;
			Math::Functs::ColorMapping((double)iter / (double)cutting_lines_.size(), r, g, b);

			Vector3d s = cutting_lines_[iter].cut.l_cutting_line_s;
			Vector3d e = cutting_lines_[iter].cut.l_cutting_line_e;

			double minimal_distance_s = 10000000.0;
			double minimal_distance_e = 10000000.0;
			for (auto iter_ = 0; iter_ < cutting_lines_.size(); iter_++)
			{
				if (iter_ != iter)
				{
					Vector3d seg_s = cutting_lines_[iter_].cut.l_cutting_line_s;
					Vector3d seg_e = cutting_lines_[iter_].cut.l_cutting_line_e;

					double d_s = distance3d(s, seg_s, seg_e);
					double d_e = distance3d(e, seg_s, seg_e);

					minimal_distance_s = min(minimal_distance_s, d_s);
					minimal_distance_e = min(minimal_distance_e, d_e);
				}
			}

			Vector3d s_e = e - s;

			Math::Functs::SetVectorLength(s_e, 3.5 * kerf_);
			if (minimal_distance_s < 2.0 * kerf_)
				s = s + s_e;
			if (minimal_distance_e < 2.0 * kerf_)
				e = e - s_e;

			//output(export_fie, export_int, pref_str + "_cutting_edge", r, g, b, s + delta_y, e + delta_y, kerf_ * 1.2);
			//number_label.DrawNode(NumberLabel::BaseShere, export_fie, export_int, pref_str + "_cutting_edge", kerf_ * 4.0, s + delta_y);
			//number_label.DrawNode(NumberLabel::BaseShere, export_fie, export_int, pref_str + "_cutting_edge", kerf_ * 4.0, e + delta_y);
		}

		for (auto iter = 0; iter < cutting_lines_.size(); iter++)
		{
			double r, g, b;
			Math::Functs::ColorMapping((double)iter / (double)cutting_lines_.size(), r, g, b);

			Vector3d s = cutting_lines_[iter].cut.l_cutting_line_s;
			Vector3d e = cutting_lines_[iter].cut.l_cutting_line_e;

			double minimal_distance_s = 10000000.0;
			double minimal_distance_e = 10000000.0;
			for (auto iter_ = 0; iter_ < cutting_lines_.size(); iter_++)
			{
				if (iter_ != iter)
				{
					Vector3d seg_s = cutting_lines_[iter_].cut.l_cutting_line_s;
					Vector3d seg_e = cutting_lines_[iter_].cut.l_cutting_line_e;

					double d_s = distance3d(s, seg_s, seg_e);
					double d_e = distance3d(e, seg_s, seg_e);

					minimal_distance_s = min(minimal_distance_s, d_s);
					minimal_distance_e = min(minimal_distance_e, d_e);
				}
			}

			Vector3d s_e = e - s;

			Math::Functs::SetVectorLength(s_e, 3.5 * kerf_);
			if (minimal_distance_s < 2.0 * kerf_)
			{
				s = s + s_e;
			}
			if (minimal_distance_e < 2.0 * kerf_)
			{
				e = e - s_e;
			}

			Vector3d sss = s + cutting_lines_[iter].cut.l_u_vector;
			Vector3d eee = e + cutting_lines_[iter].cut.l_u_vector;

			//output(export_fie, export_int, pref_str+"_cutting_edge", r, g, b, s + delta_y, sss + delta_y, kerf_ * 1.2);
			//output(export_fie, export_int, pref_str + "_cutting_edge", r, g, b, e + delta_y, eee + delta_y, kerf_ * 1.2);
		}

		//part shape label


		for (int i = 0; i < points.size(); i++)
		{
			auto part_center = Math::Functs::GetCenter(points[i]);
			part_center[2] = wl_.size[2] + 2.0 * kerf_;
			part_center += delta_y;

			number_label.DrawLabel(
				NumberLabel::BaseSquare2,
				export_file,
				export_int,
				"name",
				all_shapes_[i],
				//3407,
				kerf_ * 6.0,
				part_center);
		}

		for (int i = 0; i < points.size(); i++)
		{
			for (int j = 0; j < points[i].size(); j++)
				for (int k = 0; k < points[i][j].size(); k++)
					export_file << fixed << setprecision(8) << "v " << points[i][j][k][0] + delta_y[0] << " " << points[i][j][k][1] + delta_y[1] << " " << points[i][j][k][2] + delta_y[2] << std::endl;
		}
		

		for (int i = 0; i < points.size(); i++)
		{
			export_file << "g " << pref_str + "_" + all_names[i] << std::endl;
			for (int j = 0; j < points[i].size(); j++)
			{
				export_file << "f ";
				for (int k = 0; k < points[i][j].size(); k++)
				{
					export_file << Math::Functs::IntString(export_int) << " ";
					export_int++;
				}
				export_file << "" << std::endl;
			}
		}

	};

	std::vector<std::vector<Vector3d1 > > all_points;

	//fixed positions
	std::vector<string> all_names;
	std::vector<int> all_shapes;
	for (auto position : arrange->positions)
	{
		Vector3d2 points = position.transform->surfs;

		for (auto& points_ : points)
		{
			for (auto& p : points_)
			{
				glm::vec4 v(p, 1.0);
				Vector3d t_v = Vector3d(position.M * v);
				p = Vector3d(position.M * v);
			}
		}
		all_points.emplace_back(points);
		all_names.emplace_back("shape_" + Math::Functs::IntString(position.transform->pce_shape->index));
		all_shapes.emplace_back(position.transform->pce_shape->index);
	}

	//wl
	Output(export_fie, export_int, refer_values, delta_y, pref_str,*hModule, all_points, all_names, all_shapes, cutting_lines, cutting_line_ids, cutting_line_colors, refer_edges, refer_faces, kerf, wl, number_label);

	all_points.clear();
	std::vector<std::vector<Vector3d1 > >().swap(all_points);
}

void PCE::PartsRegistration(GLMATS& glmats,const std::vector<shared_ptr<PCEShape>> &ptr_shapes_) 
{
	std::vector<int> parts_encodes_;

	//init lpce_shapes
	std::vector<PCEShape> lpce_shapes;
	for (auto& glmat : glmats)
	{
		lpce_shapes.emplace_back(PCEShape());
		auto& local_shape = lpce_shapes.back();
		local_shape.ashape = std::get<1>(glmat);
		local_shape.index = lpce_shapes.size() - 1;
		GeomFunc::LoadTopoDSGeometry(std::get<1>(glmat), local_shape.surfs, local_shape.surf_normals, local_shape.volume, local_shape.holes, local_shape.center);
	}


	Vector3d1 rotations = Math::Functs::EmumerateRotations();

	auto shapes_ = ptr_shapes_;
	std::map<int, glm::dmat4> part_LUMs;
	for (auto& local_shape : lpce_shapes)
	{
		int existing_shape_index = -1;
		if (!shapes_.empty())
		{
			for (const auto& rotation : rotations)
			{
				glm::dmat4 PM =
					Math::Functs::TranslationMatrix(local_shape.center) *
					Math::Functs::RotationMatrix(Vector3d(0, 0, 1), -rotation[2]) *
					Math::Functs::RotationMatrix(Vector3d(0, 1, 0), -rotation[1]) *
					Math::Functs::RotationMatrix(Vector3d(1, 0, 0), -rotation[0]) *
					Math::Functs::TranslationMatrix(-local_shape.center);

				auto surfs = Math::Functs::PosApplyM(local_shape.surfs, PM);
				Vector3d maximal_corner, minimal_corner;
				Math::Functs::GetBoundingBox(surfs, minimal_corner, maximal_corner);
				for (auto& surf : surfs) for (auto& p : surf) p = p - minimal_corner;

				for (auto& shape : shapes_)
				{
					//dis_match_error
					double dis = Math::Functs::GetDistance(Math::Functs::PosApplyM(shape->surfs, shape->init_LUM), surfs);

					if (dis < CompilerConfig::Instance().GetPartMatchError()) {
						existing_shape_index = shape->index;
						auto M = Math::Functs::TranslationMatrix(-minimal_corner) * PM;
						part_LUMs.insert(std::pair<int, glm::dmat4>(local_shape.index, M));
						break;
					}
				}
				if (existing_shape_index >= 0) break;
			}
		}

		if (existing_shape_index < 0)
		{
			Vector3d maximal_corner, minimal_corner;
			Math::Functs::GetBoundingBox(local_shape.surfs, minimal_corner, maximal_corner);
			shapes_.emplace_back(shared_ptr<PCEShape>(new PCEShape(local_shape)));
			shapes_.back()->index = shapes_.size() - 1;
			auto M = Math::Functs::TranslationMatrix(-minimal_corner);
			part_LUMs.insert(std::pair<int, glm::dmat4>(local_shape.index, M));
			shapes_.back()->init_LUM = M;
			existing_shape_index = shapes_.back()->index;
		}

		parts_encodes_.emplace_back(existing_shape_index);
	}

	for (int i = 0; i < glmats.size(); i++)
	{
		auto mat = std::get<2>(glmats[i]) * glm::inverse(part_LUMs.at(i));
		const auto trsf = GeomFunc::ConvertTogpTrsf(part_LUMs.at(i));
		auto l_shape = BRepBuilderAPI_GTransform(std::get<1>(glmats[i]), trsf).Shape();
		std::get<2>(glmats[i]) = mat;
		std::get<1>(glmats[i]) = l_shape;
	}
}


void PCE::BuildShapes(
	std::vector<WLCollection>& wlcs_,
	const int& ct_egraph_heads_size_,
	const GLMATS& glmats,
	const std::string& packing_o_folder_,
	std::vector<shared_ptr<PCEShape>>& ptr_shapes_,
	std::vector<int>& parts_encodes_,
	const bool WLCFree)
{
	const auto Comp = [](const PCETransform::WLRATE& wl_rate_0, const PCETransform::WLRATE& wl_rate_1)
	{
		if (wl_rate_0.perfect_matching > wl_rate_1.perfect_matching)
			return true;
		else
		{
			if (wl_rate_0.perfect_matching == wl_rate_1.perfect_matching)
				return wl_rate_0.volume_rate > wl_rate_1.volume_rate;
			else
				return false;
		}
	};

	auto ReadPCEPart = [&](
		const int& part_index,
		const TopoDS_Shape& local_aShape,
		const TopoDS_Shape& global_aShape,
		const glm::dmat4& M)
	{
		shared_ptr<PCEShape> ptr_cur_shape(new PCEShape(local_aShape));

		int existing_shape_index = -1;
		if (!ptr_shapes_.empty())
		{
			auto surfs = ptr_cur_shape->surfs;
			for (auto& ptr_shape_ : ptr_shapes_)
			{
				//dis_match_error
				double dis = Math::Functs::GetDistance(ptr_shape_->surfs, surfs);
				if (dis < CompilerConfig::Instance().GetPartMatchError())
				{
					existing_shape_index = ptr_shape_->index;
					break;
				}
			}
		}

		if (existing_shape_index < 0)
		{
			ptr_shapes_.emplace_back(ptr_cur_shape);
			ptr_shapes_.back()->index = ptr_shapes_.size() - 1;
			existing_shape_index = ptr_shapes_.back()->index;

			//assignment
			auto ptr_shape_ = ptr_shapes_.back();

			for (auto& wlc : wlcs_)
			{
				std::vector<PCETransform> transforms = wlc.DetectSize(
					CompilerConfig::Instance().GetHModule(), 
					CompilerConfig::Instance().GetKerf(),
					CompilerConfig::Instance().GetPartMatchError(),
					ptr_shape_, WLCFree);

				if (!transforms.empty() && (ptr_shape_->transforms.empty() ||
					Comp(transforms.front().wl_rates.front(), ptr_shape_->transforms.front().wl_rates.front())))
				{
					ptr_shape_->transforms = transforms;
					ptr_shape_->WLC_index = wlc.index;
				}
			}

			if (ptr_shape_->WLC_index < 0)
			{
				std::cerr << "Cannot find any library to place this shape: " << ptr_shape_->index << std::endl;
				system("pause");
				return;
			}

			wlcs_[ptr_shape_->WLC_index].input_shapes.emplace_back(ptr_shape_->index);

			if (!packing_o_folder_.empty() && CompilerConfig::Instance().GetObjOutput())
			{
				for (auto& transform : ptr_shape_->transforms)
				{
					std::string path = packing_o_folder_ + "\\transform\\shape_transform_" + Math::Functs::IntString(ptr_shape_->index) + "_" + Math::Functs::IntString(transform.index) + ".obj";
					Math::Functs::OutputObj3d(path, transform.surfs,1);
				}
			}
		}

		parts_encodes_.emplace_back(existing_shape_index);
	};

	
	//std::cerr << "BuildShapes******BuildShapes********BuildShapes******BuildShapes********" << std::endl;
	//if (!packing_o_folder_.empty()) Math::Functs::ClearFolder(packing_o_folder_+"\\transform");
	if(!Functs::DetectExisting(packing_o_folder_ + "\\transform"))Math::Functs::ClearFolder(packing_o_folder_ + "\\transform");

	Vector3d3 g_vec3d3, l_vec3d3;

	for(int i=0;i<glmats.size();i++)
	//for (int i = 0; i < g_shapes.size(); i++)
	{
		int nb = ptr_shapes_.size();
		ReadPCEPart(i, std::get<1>(glmats[i]), std::get<0>(glmats[i]), std::get<2>(glmats[i]));

		if (ptr_shapes_.size() > nb)
		{
			auto ptr_shape_ = ptr_shapes_.back();
			if (!packing_o_folder_.empty() && CompilerConfig::Instance().GetObjOutput())
			{
				std::string path = packing_o_folder_ + "\\shape_" + std::to_string(ptr_shape_->index) + ".obj";
				Math::Functs::OutputObj3d(path, ptr_shape_->surfs, false);
			}

			std::cerr << "Shape_Index: " << ptr_shape_->index << " wlc: " << ptr_shape_->WLC_index;
			std::cerr << " #transforms: " << ptr_shape_->transforms.size() << " transforms: wl_index(volume_rate) " << std::endl;
			for (auto& transform : ptr_shape_->transforms)
			{
				std::cerr << "            Rotation: " << transform.rotation[0] << " " << transform.rotation[1] << " " << transform.rotation[2] << " : Perfect_Matching&Volume_Rate";
				for (auto wl_rate : transform.wl_rates) std::cerr << wl_rate.wl->index << "(" << wl_rate.perfect_matching << "/" << wl_rate.volume_rate << ") ";
				std::cerr << "\n";
			}
			std::cerr << "" << std::endl;
		}

		if (!packing_o_folder_.empty())
		{
			g_vec3d3.emplace_back(GeomFunc::LoadTopoDSGeometry(std::get<0>(glmats[i])));
			l_vec3d3.emplace_back(GeomFunc::LoadTopoDSGeometry(std::get<1>(glmats[i])));

#ifndef MAINOBJ
			std::string global_iges_path = packing_o_folder_ + "\\part_global_" + Math::Functs::IntString(i) + ".iges";
			std::string local_iges_path = packing_o_folder_ + "\\part_local_" + Math::Functs::IntString(i) + ".iges";
			GeomFunc::ExportToIGES(std::get<0>(glmats[i]), global_iges_path.c_str());
			GeomFunc::ExportToIGES(std::get<1>(glmats[i]), local_iges_path.c_str());
#endif
		}
	}

	if (!packing_o_folder_.empty()&& CompilerConfig::Instance().GetObjOutput())
	{
		std::string global_path = packing_o_folder_ + "\\" + std::to_string(ct_egraph_heads_size_) + "_input_global.obj";
		std::string local_path = packing_o_folder_ + "\\" + std::to_string(ct_egraph_heads_size_) + "_input_local.obj";
		Math::Functs::OutputObj3d(global_path, g_vec3d3, 2);
		Math::Functs::OutputObj3d(local_path, l_vec3d3, 2);
	}

	//output mat
#ifndef MAINOBJ
	if (!packing_o_folder_.empty())
	{
		std::string mat_path = packing_o_folder_ + "\\mat.txt";
		std::ofstream file(mat_path);
		file << glmats.size() << std::endl;

		for (auto& glmat : glmats)
		{
			auto& m = std::get<2>(glmat);
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					if (j == 3)
						file << m[i][j];
					else
						file << m[i][j] << " ";
				}
				file << "\n";
			}
		}

		file.close();
	}
#endif
}

VectorPI1 PCE::GetCTArrangeEgraphIndex(const ArrNode& arr_node)
{
	if(!arr_node.order_levels.empty()&& !arr_node.order_levels.back().empty())
		return  GetCTArrangeEgraphIndex(arr_node.order_levels, arr_node.order_levels.back().front());
	return VectorPI1();
}

VectorPI1 PCE::GetCTArrangeEgraphIndex(const std::vector<std::vector<OrderNode>>& order_levels, const OrderNode& order_node)
{
	VectorPI1 arrange_eclass;
	OrderNode node = order_node;
	while (true) {
		if (node.level < 0)break;

		if (node.ct_arrange_eclass_index.first >= 0 && node.ct_arrange_eclass_index.second >= 0)
			arrange_eclass.emplace_back(node.ct_arrange_eclass_index);

		if (node.level == 0) break;
		node = order_levels[node.level - 1][node.parent];
	}
	std::reverse(arrange_eclass.begin(), arrange_eclass.end());
	return arrange_eclass;
}

void PCE::GetCTArrangeEgraphIndex(
	const std::vector<std::vector<OrderNode>>& order_levels,
	const OrderNode& order_node,
	std::vector<int>& ct_arrange_indexes,
	std::vector<int>& ct_eclass_indexes)
{
	OrderNode node = order_node;
	while (true) {
		if (node.level < 0)break;

		if (node.ct_arrange_eclass_index.first >= 0 && node.ct_arrange_eclass_index.second >= 0)
		{
			ct_arrange_indexes.emplace_back(node.ct_arrange_eclass_index.first);
			ct_eclass_indexes.emplace_back(node.ct_arrange_eclass_index.second);
		}

		if (node.level == 0) break;
		node = order_levels[node.level - 1][node.parent];
	}
	std::reverse(ct_arrange_indexes.begin(), ct_arrange_indexes.end());
	std::reverse(ct_eclass_indexes.begin(), ct_eclass_indexes.end());
};


int PCEEncode::Insert(std::map<int, int>& encode, const int& shape)
{
	if (encode.find(shape) == encode.end())
		encode[shape] = 0;
	encode[shape]++;

	return encode[shape];
}

std::map<int, int> PCEEncode::GetEncode(const std::vector<int>& shapes)
{
	std::map<int, int> encode;
	for (auto& shape : shapes) Insert(encode, shape);
	return encode;
}

int PCEEncode::GetEncodeSize(const std::map<int, int>& encode)
{
	int nb = 0;
	for (auto iter = encode.rbegin(); iter != encode.rend(); iter++)
		nb += iter->second;
	return nb;
}

int PCEEncode::GetEncodeSize(const Vector1i1& shapes)
{
	return shapes.size();
}

std::vector<int> PCEEncode::GetShapes(const std::map<int, int>& encode)
{
	std::vector<int> shapes;
	for (auto iter = encode.rbegin(); iter != encode.rend(); iter++)
		for (int i = 0; i < iter->second; i++)shapes.emplace_back(iter->first);
	return shapes;
}
bool PCEEncode::CheckInclude(const Vector1i1& shapes_0, const Vector1i1& shapes_1)
{
	return CheckInclude(GetEncode(shapes_0), GetEncode(shapes_1));
}

bool PCEEncode::CheckInclude(const std::map<int, int>& big_encodes, const std::map<int, int>& small_encodes)
{
	for (auto iter = small_encodes.rbegin(); iter != small_encodes.rend(); iter++)
	{
		if (big_encodes.find(iter->first) == big_encodes.end())
		{
			return false;
		}
		else
		{
			if (small_encodes.at(iter->first) > big_encodes.at(iter->first))
				return false;
		}
	}
	return true;
}

std::map<int, int> PCEEncode::Subtract(const Vector1i1& shapes_0, const Vector1i1& shapes_1)
{
	return Subtract(GetEncode(shapes_0), GetEncode(shapes_1));
}

std::map<int, int> PCEEncode::Subtract(const std::map<int, int>& big_encodes, const std::map<int, int>& small_encodes)
{
	std::map<int, int> sub_map;
	for (auto iter = big_encodes.rbegin(); iter != big_encodes.rend(); iter++)
	{
		if (small_encodes.find(iter->first) == small_encodes.end())
		{
			sub_map[iter->first] = iter->second;
		}
		else
		{
			if (big_encodes.at(iter->first) > small_encodes.at(iter->first))
			{
				sub_map[iter->first] = big_encodes.at(iter->first) - small_encodes.at(iter->first);
			}
		}
	}
	return sub_map;
}

bool PCEEncode::CheckIdentify(const Vector1i1& shapes_0, const Vector1i1& shapes_1)
{
	return CheckIdentify(GetEncode(shapes_0), GetEncode(shapes_1));
}

bool PCEEncode::CheckIdentify(const std::map<int, int>& encode_0, const std::map<int, int>& encode_1)
{
	if (encode_0.size() != encode_1.size()) return false;
	return CheckInclude(encode_0, encode_1) && CheckInclude(encode_1, encode_0);
}

std::map<int, int> PCEEncode::Union(const Vector1i1& shapes_0, const Vector1i1& shapes_1)
{
	return Union(GetEncode(shapes_0), GetEncode(shapes_1));
}

std::map<int, int> PCEEncode::Union(const std::map<int, int>& encode_0, const std::map<int, int>& encode_1)
{
	std::vector<std::map<int, int>> encodes{encode_0,encode_1};
	return Union(encodes);
}

std::map<int, int> PCEEncode::Union(const std::vector<std::map<int, int>>& encodes)
{
	std::map<int, int> union_shapes;
	for (auto& encode : encodes)
	{
		for (auto iter = encode.rbegin(); iter != encode.rend(); iter++)
		{
			if (union_shapes.find(iter->first) == union_shapes.end())
				union_shapes[iter->first] = iter->second;
			else
				union_shapes[iter->first] += iter->second;
		}
	}
	return union_shapes;
}

Vector1i2 PCE::ComputeGroupCombOrderAdd(const Vector1i2& shape_groups_, const Vector1i1& add_shapes, const int max_nb)
{
	Vector1i2 shape_groups(shape_groups_.size(), add_shapes);

	auto selections = ComputeGroupCombOrderDel(shape_groups, add_shapes, max_nb);

	Vector1i1 all;
	for (int i = 0; i < add_shapes.size(); i++) all.emplace_back(i);

	Vector1i2 outputs;
	for (int i = 0; i < selections.size(); i++)
	{
		outputs.emplace_back(Vector1i1(add_shapes.size(), -1));
		for (int j = 0; j < selections[i].size(); j++)
		{
			auto a = PCEEncode::GetShapes(PCEEncode::Subtract(all, selections[i][j]));
			for (int k = 0; k < a.size(); k++) outputs.back()[a[k]] = j;
		}
	}
	
	//std::cerr << Math::Functs::IntString(outputs,false,",","\n") << std::endl;
	random_shuffle(outputs.begin(), outputs.end());

	return outputs;
}

Vector1i3 PCE::ComputeGroupCombOrderDel(const Vector1i2& shape_groups, const Vector1i1& delete_shapes, const int max_nb)
{
	//check
	Vector1i2 group_comb_order;
	Vector1i1 shapes;
	std::map<int, int> temp;
	for (auto& shape_group : shape_groups)
	{
		group_comb_order.emplace_back(Vector1i1());
		for (auto& shape : shape_group)
		{
			shapes.emplace_back(shape);
			group_comb_order.back().emplace_back(PCEEncode::Insert(temp, shape));
		}
	}

	//get delete_shapes_combs
	auto shapes_encode = PCEEncode::GetEncode(shapes);
	auto delete_shapes_encode = PCEEncode::GetEncode(delete_shapes);
	Vector1i3 combs;
	std::map<int, int> delete_shapes_combs;
	for (auto iter = delete_shapes_encode.rbegin(); iter != delete_shapes_encode.rend(); iter++)
	{
		delete_shapes_combs.insert(std::pair<int, int>(iter->first, delete_shapes_combs.size()));
		combs.emplace_back(Math::Functs::CombNonRepeat(shapes_encode.at(iter->first), iter->second));
	}


	Vector1i3 all_selections;
	Math::Functs::Combination(combs, 0, Vector1i2(), max_nb,all_selections);

	Vector1i3 all_outputs;
	for (int selection_id = 0; selection_id < all_selections.size(); selection_id++)
	{
		const auto select = all_selections[selection_id];

		all_outputs.emplace_back(Vector1i2());
		for (int group_id = 0; group_id < shape_groups.size(); group_id++)
		{
			all_outputs.back().emplace_back(Vector1i1());
			for (int shape_id = 0; shape_id < shape_groups[group_id].size(); shape_id++)
			{
				auto shape = shape_groups[group_id][shape_id];

				if (delete_shapes_combs.find(shape) == delete_shapes_combs.end())
				{
					all_outputs.back().back().emplace_back(shape_id);
				}
				else
				{
					auto t = select[delete_shapes_combs.at(shape)];
					auto c = group_comb_order[group_id][shape_id];
					if (std::find(t.begin(), t.end(), c) == t.end())
					{
						all_outputs.back().back().emplace_back(shape_id);
					}
				}
			}
		}
		//std::cerr << Math::Functs::IntString(vs, false, ",", "/") << std::endl;
	}
	return all_outputs;
}

ArrNode PCE::AddOperation(const ArrNode& old_arr_node, const Vector1i1& add_shape_ids)
{
	//add shape
	auto& vec_shapes_ = CompilerConfig::Instance().GetVecShape();

	//wlc wl
	auto& wlc = CompilerConfig::Instance().GetVecWLC()[old_arr_node.wlc_index];
	auto& wl = wlc.wls[old_arr_node.wl_index];

	ArrNode add_arr_node;
	add_arr_node.wlc_index = old_arr_node.wlc_index;
	add_arr_node.wl_index = old_arr_node.wl_index;

	//order levels
	add_arr_node.order_levels = old_arr_node.order_levels;
	if (!add_arr_node.order_levels.empty() && !add_arr_node.order_levels.back().empty())
		add_arr_node.order_levels.back().front().ct_arrange_eclass_index = std::pair<int, int>(-1, -1);

	//seq_shapes_transforms
	std::vector<std::pair<int, Vector1i1>> seq_shapes_transforms(old_arr_node.shape_order.size(), std::pair<int, Vector1i1>());
	Vector1i1 seq_shapes = old_arr_node.shape_order;
	for (int i = 0; i < add_shape_ids.size(); i++)
	{
		add_arr_node.order_levels.emplace_back(std::vector<OrderNode>());

		seq_shapes_transforms.emplace_back(std::pair<int, Vector1i1>());
		seq_shapes_transforms.back().first = add_shape_ids[i];
		for (auto& transform : vec_shapes_[add_shape_ids[i]]->transforms)
			seq_shapes_transforms.back().second.emplace_back(transform.index);

		seq_shapes.emplace_back(add_shape_ids[i]);
	}

	//packing according to the shape_transform
	OrderNode* order_node;
	if (add_arr_node.order_levels.size() >= 2)
		order_node = &add_arr_node.order_levels[add_arr_node.order_levels.size() - (1+ add_shape_ids.size())].front();
	else
		order_node = new OrderNode(wl.ExtractEdges());


	Order2Arrange(std::vector<std::string>(), wlc, wl, *order_node,
		1, PCEEncode::GetEncode(seq_shapes), 
		std::map<int, Vector1i1>(),
		seq_shapes_transforms, 
		add_arr_node.order_levels, false, -1);



	if (!add_arr_node.order_levels.empty() && !add_arr_node.order_levels.back().empty())
	{
		//add_arr_node.arrange_id = add_arr_node.order_levels.back().front().ct_arrange_eclass_index.first;
		add_arr_node.shape_order = seq_shapes;
	}

	auto add_arrange_eclss = GetCTArrangeEgraphIndex(add_arr_node);
	if (add_arrange_eclss.size() != 1)
	{
		int d = 0;
	}

	return add_arr_node;
};

//sequence packing
ArrNode PCE::SequencePacking(const int& wlc_index_, const int& wl_index_, const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms)
{
	//wlc wl
	auto& wlc = CompilerConfig::Instance().GetVecWLC()[wlc_index_];
	auto& wl = wlc.wls[wl_index_];

	ArrNode arr_node;
	arr_node.wlc_index = wlc_index_;
	arr_node.wl_index = wl_index_;

	//Vector1i1 seq_shapes;
	for (auto& shape_transform : seq_shapes_transforms)
		arr_node.shape_order.emplace_back(shape_transform.first);

	//packing according to the shape_transform
	arr_node.order_levels = std::vector<std::vector<OrderNode>>(arr_node.shape_order.size(), std::vector<OrderNode>());//order levels
	Order2Arrange(
		std::vector<std::string>(), 
		wlc, wl, OrderNode(wl.ExtractEdges()),
		1, PCEEncode::GetEncode(arr_node.shape_order), 
		std::map<int, Vector1i1>(), 
		seq_shapes_transforms, 
		arr_node.order_levels, 
		false, -1);

	//shape transforms
	if (!arr_node.order_levels.empty()&& !arr_node.order_levels.back().empty())
	{
		//arr_node.arrange_id = arr_node.order_levels.back().front().ct_arrange_eclass_index.first;
	}

	return arr_node;
};

Vector1i1 PCE::ArrangeDAI(
	const Vector1i1& arrange_ids,
	const Vector1i1& delete_shapes, 
	const Vector1i1& add_shapes)
{
	//check
	Vector1i2 shape_groups;
	for (auto& arrange_id : arrange_ids)
		shape_groups.emplace_back(ct_egraph.ptr_arranges[arrange_id]->shape_order);

	Vector1i3 all_selections = ComputeGroupCombOrderDel(shape_groups, delete_shapes,1);

	auto used_selections = all_selections.front();


	//delete operations
	std::vector<ArrNode> delete_nodes;
	Vector1i1 output_old_arrange_ids;
	for (int arrange_iter = 0; arrange_iter < used_selections.size(); arrange_iter++)
	{
		auto& arrange = ct_egraph.ptr_arranges[arrange_ids[arrange_iter]];

		std::vector<std::pair<int, Vector1i1>> seq_shapes_transforms;
		for (int selection_iter = 0; selection_iter < used_selections[arrange_iter].size(); selection_iter++)
		{
			auto& position = arrange->positions[used_selections[arrange_iter][selection_iter]];
			std::pair<int, Vector1i1> shape_transform(position.transform->pce_shape->index, Vector1i1(1, position.transform->index));
			seq_shapes_transforms.emplace_back(shape_transform);
		}

		if (seq_shapes_transforms.size() < arrange->shape_order.size())
			delete_nodes.emplace_back(SequencePacking(arrange->wlc_index, arrange->wl_index, seq_shapes_transforms));
		else
			output_old_arrange_ids.emplace_back(arrange_ids[arrange_iter]);
	}

	//add operation
	auto add_nodes = delete_nodes;
	for (auto& add_shape : add_shapes)
	{
		bool goon = false;
		for (auto& arr_node : add_nodes)
		{
			auto add_arr_node = AddOperation(arr_node, Vector1i1(1,add_shape));
			auto arrange_eclss = GetCTArrangeEgraphIndex(arr_node);
			auto add_arrange_eclss = GetCTArrangeEgraphIndex(add_arr_node);
			if (add_arrange_eclss.size() == 1)
			{
				std::cerr <<"Add(arrange,eclass): "<< Math::Functs::IntString(arrange_eclss, false, ",", "/") << " => " << Math::Functs::IntString(add_arrange_eclss, false, ",", "/") << std::endl;
				arr_node = add_arr_node;
				goon = true;
				break;
			}
		}

		if (!goon)
		{
			for (auto& arr_node : add_nodes)
			{
				auto add_arr_node = AddOperation(arr_node, Vector1i1(1, add_shape));
				auto arrange_eclss = GetCTArrangeEgraphIndex(arr_node);
				auto add_arrange_eclss = GetCTArrangeEgraphIndex(add_arr_node);
				if (add_arrange_eclss.size()>0)
				{
					std::cerr << "Add(arrange,eclass): "<< Math::Functs::IntString(arrange_eclss, false, ",", "/") << " => " << Math::Functs::IntString(add_arrange_eclss, false, ",", "/") << std::endl;
					arr_node = add_arr_node;
					goon = true;
					break;
				}
			}
		}

		if (!goon)
		{
			std::cerr << "for (auto& add_shape : add_shapes)" << std::endl;
			system("pause");
		}
	}

	Vector1i1 update_arrange_ids= output_old_arrange_ids;
	for (auto& add_node : add_nodes)
	{
		auto add_arrange_eclss = GetCTArrangeEgraphIndex(add_node);
		for (auto& arrange_eclass : add_arrange_eclss)update_arrange_ids.emplace_back(arrange_eclass.first);
	}
	
	return update_arrange_ids;
}
