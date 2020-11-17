#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <iostream>
#include <Mod/PartDesign/Gui/CompilerConfig.h>
#include <omp.h>
#include <Eigen/SVD>  
#include <Eigen/Dense>

using namespace PartDesignGui;

PCE::PCE(bool b, const std::string& packing_o_folder_)
{
	//load the dll built from carpentry-geom-lib
	hModule = &CompilerConfig::Instance().GetHModule();

	//parameter initialization
	kerf = 3.175;
	part_match_error = 0.01;
	angle_match_error = 0.08;
	maximal_cutting_distance = 200.0;
	total_order_number = 100;
	prunning_limitation = 10;
	random_pruning = false;
	safe_distance = 10.0;
	tree_searching_branch = 3;
	output_all_log = false;

	Init(b,packing_o_folder_);
}
PCE::PCE(const std::string& pref_path)
{
	//load the dll built from carpentry-geom-lib
	hModule = &CompilerConfig::Instance().GetHModule();

	//parameter initialization
	kerf = 3.175;
	part_match_error = 0.01;
	angle_match_error = 0.08;
	maximal_cutting_distance = 200.0;
	total_order_number = 100;
	prunning_limitation = 10;
	random_pruning = false;
	safe_distance = 10.0;
	tree_searching_branch = 3;
	output_all_log = false;

	Init(pref_path);
}

void PCE::Init(bool b, const std::string& packing_o_folder_)
{
	ClearAll();

	auto& config = CompilerConfig::Instance();
	if (hModule && config.HasLoaded())
	{
		kerf = config.GetKerf();
		part_match_error = config.GetPartMatchError();
		maximal_cutting_distance = config.GetMaximalCuttingDistance();
		total_order_number = config.GetTotalOrderNumber();
		prunning_limitation = config.GetPruningLimitation();
		safe_distance = config.GetSafeDistance();
		random_pruning = config.IsRandomPruning();
		packing_o_folder = packing_o_folder_;
		safe_distance = ((int)(safe_distance / kerf) + 1) * kerf;
	}
	else
	{
		std::cerr << "Load configuration file error..." << std::endl;
		system("pause");
		return;
	}

	Math::Functs::ClearFolder(packing_o_folder);
	Math::Functs::ClearFolder(packing_o_folder + "\\temp");
}

void PCE::Init(const std::string& pref_path)
{
	ClearAll();

	auto& config = CompilerConfig::Instance();
	if (hModule && config.HasLoaded())
	{
		kerf = config.GetKerf();
		part_match_error = config.GetPartMatchError();
		maximal_cutting_distance = config.GetMaximalCuttingDistance();
		total_order_number = config.GetTotalOrderNumber();
		prunning_limitation = config.GetPruningLimitation();
		safe_distance = config.GetSafeDistance();
		random_pruning = config.IsRandomPruning();
		packing_o_folder = config.GetPackingOutputFolder();
		if (pref_path.size() != 0) packing_o_folder = packing_o_folder + "\\" + pref_path;
		safe_distance = ((int)(safe_distance / kerf) + 1) * kerf;
	}
	else
	{
		std::cerr << "Load configuration file error..." << std::endl;
		system("pause");
		return;
	}

	Math::Functs::ClearFolder(packing_o_folder);
	Math::Functs::ClearFolder(packing_o_folder + "\\temp");
}

//2. egraph arrange
void PCE::EnumerateOrders(std::vector<WLCollection>& wlcs_, const std::vector<int>& shape_encodes_)
{
	auto GetStartNodes = [](const WLCollection& wlc_, const int& total_order_number_, std::vector<int>& start_nodes)
	{
		auto shape_encode = PCEEncode::GetEncode(wlc_.shapes);
		auto iter = shape_encode.rbegin();
		for (int i = 0; i < total_order_number_ && i < shape_encode.size(); i++)
		{
			start_nodes.emplace_back(iter->first);
			iter++;
		}
		return total_order_number_ <= shape_encode.size() ?
			1 : (int)((double)total_order_number_ / (double)shape_encode.size());
	};

	auto OrderLevels2Edges = [](std::vector<std::vector<OrderNode>>& order_levels, const std::string& path) {
		int node_nb = 0;
		std::vector<int> edges;
		std::vector<string> labels;
		for (int i = 0; i < order_levels.size(); i++)
		{
			auto& order_level = order_levels[i];
			for (auto& order_node : order_level) {
				order_node.index = node_nb;
				labels.emplace_back(std::to_string(node_nb) + "_" + std::to_string(order_node.position.transform->pce_shape->index));
				node_nb++;
				edges.emplace_back(order_node.index);
				if (i == 0)
					edges.emplace_back(-1);
				else
					edges.emplace_back(order_levels[i - 1][order_node.parent].index);
			}
		}

		Math::Functs::Output_tree(node_nb, edges, path, labels);
	};

	//build head
	if(ct_egraph.heads.empty())
		ct_egraph.heads.emplace_back(shared_ptr<EGraphClass>(new EGraphClass(-1, shape_encodes_)));
	else
		ct_egraph.heads.emplace_back(shared_ptr<EGraphClass>(new EGraphClass(ct_egraph.heads.back()->index - 1, shape_encodes_)));

	int head_index = ct_egraph.heads.back()->index;

	std::vector<int> head_edges;
	int ec_nb = ct_egraph.ptr_egraph.size();

	//shape transforms
	auto unique_shapes_ = Math::Functs::RemoveDuplicate(shape_encodes_);
	auto& ptr_shapes_ = CompilerConfig::Instance().GetVecShape();
	std::map<int, Vector1i1> shape_transforms;
	for (auto& shape_ : unique_shapes_)
	{
		Vector1i1 transforms_;
		for (auto& transform_ : ptr_shapes_[shape_]->transforms)
			transforms_.emplace_back(transform_.index);
		shape_transforms.insert(std::pair<int, Vector1i1>(shape_, transforms_));
	}

	//wlc
	for (auto& wlc : wlcs_)
	{
		for (auto& ec_shape : shape_encodes_)
			if (Math::Functs::CheckContain(wlc.input_shapes, ec_shape))
				wlc.shapes.emplace_back(ec_shape);

		//initialization
		if (wlc.shapes.empty()) continue;

		std::cerr << "WLC: " << wlc.index << " Shapes: " << wlc.shapes.size() << " Total orders: " << total_order_number << std::endl;

		auto ec_index = PCEEGraphContainer::CTNewEClass(ct_egraph.ptr_egraph, wlc.shapes);
		if(ec_index< ec_nb) ct_egraph.ptr_egraph[ec_index]->head_index = head_index;

		head_edges.emplace_back(ec_index);

		//select the starting nodes
		std::cerr << "////////////////////////////////////////////////////////" << std::endl;
		std::cerr << "Select the starting nodes..." << std::endl;
		std::vector<int> start_shapes;
		int each_nb = GetStartNodes(wlc, total_order_number, start_shapes);
		std::cerr << Math::Functs::IntString(start_shapes, false, " ") << std::endl;

		std::cerr << "Each start node should be following with " << each_nb << "'s arranges.\n";

		for (int start_shape : start_shapes)
		{
			std::cerr << "////////////////////////////////////////////////////////" << std::endl;
			std::cerr << "start Shape: " << start_shape <<" maximal_arrange_nb: "<< each_nb << std::endl;
			std::cerr << "eclass number: " << ct_egraph.ptr_egraph.size() << " arrange number: " << ct_egraph.ptr_arranges.size() << std::endl;
			std::cerr << "Full_part_number: " <<ct_egraph.ptr_egraph[0]->arranges.size() + ct_egraph.ptr_egraph[0]->egraph_edges.size() << std::endl;

			//each_nb
			auto& wl = wlc.wls[0];
			std::map<int, int> all_targeting_shapes = PCEEncode::GetEncode(wlc.shapes);
			OrderNode order_node(wl.ExtractEdges());
			std::vector<std::vector<OrderNode>> order_levels(wlc.shapes.size(), std::vector<OrderNode>());//order levels
			Order2Arrange(enumerate_order_strs_,wlc, wl, order_node, each_nb, all_targeting_shapes, shape_transforms, std::vector<std::pair<int, Vector1i1>>(), order_levels,true, start_shape);
			BreadthFirstOrder(enumerate_order_strs_,wlc, wl, each_nb, all_targeting_shapes, shape_transforms, std::vector<std::pair<int, Vector1i1>>(), order_levels,true);
			wl.all_enumerating_nb = 0;

			if (packing_o_folder.size() != 0)
				OrderLevels2Edges(order_levels, packing_o_folder + "\\order_levels_" + Math::Functs::IntString(wlc.index) + ".gml");

			std::cerr << "////////////////////////////////////////////////////////" << std::endl;
			std::cerr <<"eclass number: "<< ct_egraph.ptr_egraph.size() << " arrange number: " << ct_egraph.ptr_arranges.size() << std::endl;
		}

		wlc.shapes.clear();
	}

	ct_egraph.heads.back()->PushComp(head_edges);

	Vector1i1 rel_ecs;
	PCEEGraphContainer::GetRelatedEclass(ct_egraph.ptr_egraph, ct_egraph.GetHead(head_index), rel_ecs);
	
	for (auto& ptr_ec : ct_egraph.ptr_egraph)
		ptr_ec->head_index = Functs::CheckContain(rel_ecs, ptr_ec->index)? head_index : 0;

	//ct_egraph.ComputeMoreSubSets();
	ct_egraph.SetEclassValid(ct_egraph.heads.back());//only for the current head
	ct_egraph.Pruning(ct_egraph.heads.back()->index,prunning_limitation, random_pruning);//only for the current head
	ct_egraph.SetArrangeValid(ct_egraph.heads.back()->index);//only for the current head

	PostProcess(wlcs_);
}

void PCE::BreadthFirstOrder(
	std::vector<std::string>& enumerate_order_strs,
	WLCollection& wlc, WL& wl,
	const int& maximal_arrange_nb,
	const std::map<int, int>& all_targeting_shapes,
	const std::map<int, Vector1i1>& shape_transforms,
	//const VectorPI1& seq_shapes_transforms,
	const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms,
	std::vector<std::vector<OrderNode>>& order_levels,
	const bool depth_searching_bool)
{
	auto GetFixedPositions = [&](const OrderNode& order_node)
	{
		OrderNode node = order_node;
		std::vector<PCEPosition> fixed_positions;
		while (true) {
			fixed_positions.emplace_back(node.position);
			if (node.level == 0) break;
			node = order_levels[node.level - 1][node.parent];
		}
		std::reverse(fixed_positions.begin(), fixed_positions.end());
		return fixed_positions;
	};

	//check all_enumerating_nb
	if (wl.all_enumerating_nb >= maximal_arrange_nb)return;

	for (int level = 0; level < order_levels.size() - 1; level++)
	{
		if (depth_searching_bool)
		{
			if (wl.index == 0 || output_all_log)
				std::cerr << "/////////////////////_Level: " << level << "////////////////////////" << std::endl;
		}

		auto& order_level = order_levels[level];
		for (auto& order_node : order_level)
		{
			if (order_node.visited)continue;
			Order2Arrange(enumerate_order_strs,wlc, wl, order_node, maximal_arrange_nb, all_targeting_shapes, shape_transforms, seq_shapes_transforms, order_levels, depth_searching_bool);
			if (wl.all_enumerating_nb >= maximal_arrange_nb)return;
		}
	}
}



void PCE::Order2Arrange(
	std::vector<std::string>& enumerate_order_strs,
	WLCollection& wlc, WL& wl,
	OrderNode& parent_order_node,
	const int& maximal_arrange_nb,
	const std::map<int, int>& all_targeting_shapes,
	const std::map<int, Vector1i1>& shape_transforms,
	//const VectorPI1& seq_shapes_transforms,
	const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms,
	std::vector<std::vector<OrderNode>>& order_levels,
	const bool depth_searching_bool,
	const int start_shape_index)
{
	auto GetFixedGroupPositions = [&](const OrderNode& order_node)
	{
		OrderNode node = order_node;
		std::vector<PCEPosition> fixed_positions;
		while (true) {
			if (node.position.style != PCEPosition::FreeStyle)
				fixed_positions.emplace_back(node.position);
			else
				break;
			if (node.level == 0) break;
			node = order_levels[node.level - 1][node.parent];
			if (node.ct_arrange_eclass_index.first >= 0 || node.ct_arrange_eclass_index.second >= 0)break;
		}
		std::reverse(fixed_positions.begin(), fixed_positions.end());
		return fixed_positions;
	};

	auto GetVisitedParts = [&](const OrderNode& order_node, std::vector<int>& fixed_shapes)
	{
		OrderNode node = order_node;
		while (true) {
			if (node.level < 0)break;
			fixed_shapes.emplace_back(node.position.transform->pce_shape->index);
			if (node.level == 0) break;
			node = order_levels[node.level - 1][node.parent];
		}
		std::reverse(fixed_shapes.begin(), fixed_shapes.end());

		return PCEEncode::GetEncode(fixed_shapes);
	};



	const auto GetCTEclassIndex = [](
		const std::vector<int>& arrange_,
		std::vector<shared_ptr<EGraphClass>>& ptr_egraph_)
	{
		string str = Math::Functs::IntString(arrange_, true,",");
		for (auto& eclass : ptr_egraph_)
		{
			if (eclass->str == str)
				return eclass->index;
		}
		return -1;
	};

	auto GetNextPositions = [&](
		const std::vector<PCEPosition>& fixed_group_positions,
		const std::map<int, int>& visited_shapes,
		const std::vector<PCEConstraint>& fixed_group_sharing_constraints,
		const std::vector<PCEEdge>& fixed_group_edges)
	{
		//fixed_volume
		double fixed_volume = 0.0;
		for (auto& position : fixed_group_positions)
			fixed_volume += position.transform->pce_shape->volume;

		auto& ptr_shapes_ = CompilerConfig::Instance().GetVecShape();

		//next positions
		std::vector<PCEPosition> next_positions;

		if (depth_searching_bool)
		{
			if (start_shape_index < 0 || !visited_shapes.empty()) {

				auto sub_shapes = PCEEncode::GetShapes(PCEEncode::Subtract(all_targeting_shapes, visited_shapes));

				sort(sub_shapes.begin(), sub_shapes.end());
				sub_shapes.erase(unique(sub_shapes.begin(), sub_shapes.end()), sub_shapes.end());

				for (const auto& shape : sub_shapes)
				{
					if (ptr_shapes_[shape]->volume < wl.volume - fixed_volume)//volume
					{
						for (auto& transform : ptr_shapes_[shape]->transforms)
						{
							if (Math::Functs::CheckContain(shape_transforms.at(shape), transform.index))
							{
								//maximal volume
								if (transform.CheckWL(wl.index)) {
									next_positions.emplace_back(PCEPosition());
									next_positions.back().transform = &transform;
								}
							}
						}
					}
				}
			}
			else
			{
				if (ptr_shapes_[start_shape_index]->volume < wl.volume - fixed_volume)//volume
					for (auto& transform : ptr_shapes_[start_shape_index]->transforms)
					{
						if (Math::Functs::CheckContain(shape_transforms.at(start_shape_index), transform.index))
						{
							//maximal volume
							if (transform.CheckWL(wl.index)) {
								next_positions.emplace_back(PCEPosition());
								next_positions.back().transform = &transform;
							}
						}
					}
			}
		}
		else
		{
			int visited_nb = PCEEncode::GetShapes(visited_shapes).size();
			//fixed_group_positions
			//int visited_nb = fixed_group_positions.size();

			auto shape = seq_shapes_transforms[visited_nb].first;

			if (ptr_shapes_[shape]->volume < wl.volume - fixed_volume)//volume
			{
				for (auto& transform : ptr_shapes_[shape]->transforms)
				{
					auto t = seq_shapes_transforms[visited_nb].second;
					if (Math::Functs::CheckContain(t, transform.index))
					{
						if (transform.CheckWL(wl.index)) {
							next_positions.emplace_back(PCEPosition());
							next_positions.back().transform = &transform;
						}
					}
				}
			}


			int dsad = 10;
		}


		

		std::vector<OrderNode> next_order_nodes;
		for (auto& next_position : next_positions) {
			auto sharing_constraints_ = fixed_group_sharing_constraints;
			auto edges_ = fixed_group_edges;
			
			next_position.transform->pce_shape->index;
			next_position.transform->index;

			next_position.transform->AddNewEdges(edges_, next_position.edge_indexes);
			ArrangePartbyConstraints(wl, edges_, fixed_group_positions, next_position, sharing_constraints_);
			if (next_position.style != PCEPosition::FreeStyle)
				next_order_nodes.emplace_back(OrderNode(next_position, sharing_constraints_, edges_));
		}
		auto Comp = [](OrderNode& next_0, OrderNode& next_1)
		{
			//if (Math::Functs::IsAlmostZero_Double(next_0.position.score - next_1.position.score, CompilerConfig::Instance().GetPartMatchError()))
			//{
			//	return next_0.position.t[0] < next_1.position.t[0];
			//}
			//else return next_0.position.score > next_1.position.score;
			return next_0.position.score > next_1.position.score;
		};
		std::sort(next_order_nodes.begin(), next_order_nodes.end(), Comp);

		return next_order_nodes;
	};

	auto DepthSearching = [&](const OrderNode& selected_node)
	{
		//subset

		std::vector<int> ct_arrange_indexes, ct_eclass_indexes;
		GetCTArrangeEgraphIndex(order_levels, selected_node, ct_arrange_indexes, ct_eclass_indexes);
		
		std::vector<int> encode_shapes = PCEEncode::GetShapes(all_targeting_shapes);

		const int ct_father_eclass_index = GetCTEclassIndex(encode_shapes,ct_egraph.ptr_egraph);
		if (ct_father_eclass_index >= 0)
			ct_egraph.ptr_egraph[ct_father_eclass_index]->PushComp(ct_eclass_indexes);

		//recursive
		for (auto ct_arrange_index : ct_arrange_indexes)
		{
			auto& arrange = ct_egraph.ptr_arranges[ct_arrange_index];
			auto backup_shape_order = arrange->shape_order;

			bool goon = true;
			if (arrange->shape_order.size() == wlc.shapes.size())
			{
				string str = Math::Functs::IntString(arrange->shape_order, true,",");
				for (auto arrange_index_ :ct_egraph.full_part_arranges)
				{
					string str_ = Math::Functs::IntString(ct_egraph.ptr_arranges[arrange_index_]->shape_order, true,",");
					if (str == str_)
					{
						goon = false;
						break;
					}
				}
				if (goon) ct_egraph.full_part_arranges.emplace_back(ct_arrange_index);
			}

			if (goon)
			{
				for (auto& wl_ : wlc.wls)
				{
					if (wl_.index > wl.index)
					{
						//order levels
						std::vector<std::vector<OrderNode>> depth_order_levels(backup_shape_order.size(), std::vector<OrderNode>());

						//each_nb
						//std::cerr << "\n**************" << std::endl;
						OrderNode order_node(wl_.ExtractEdges());

						std::map<int, int> all_targeting_shapes = PCEEncode::GetEncode(backup_shape_order);
						auto& depth_all_targeting_shapes = PCEEncode::GetEncode(backup_shape_order);

						Order2Arrange(enumerate_order_strs,wlc, wl_, order_node, maximal_arrange_nb, depth_all_targeting_shapes, shape_transforms, seq_shapes_transforms, depth_order_levels, depth_searching_bool);
						BreadthFirstOrder(enumerate_order_strs,wlc, wl_, maximal_arrange_nb, depth_all_targeting_shapes, shape_transforms, seq_shapes_transforms, depth_order_levels, depth_searching_bool);
						wl_.all_enumerating_nb = 0;
						//std::cerr << "**************\n" << std::endl;
					}
				}
			}
		}
	};

	auto OutputArrange = [&](shared_ptr<PCEArrange>& arrange)
	{
		std::string o_path = packing_o_folder +
			"\\temp\\wlc_" + Math::Functs::IntString(wlc.index) +
			"_wl_"+ Math::Functs::IntString(wl.index) +
			"_arrange_" + Math::Functs::IntString(arrange->index) +
			"_valid_" + Math::Functs::IntString(arrange->valid) +
			"_shapes_" + Math::Functs::IntString(arrange->shape_order.size()) +
			"(" + Math::Functs::IntString(arrange->shape_order,false,",") + ")" +
			"_cutting_" + Math::Functs::IntString(arrange->cutting_number) + ".obj";

		OutputFixedPositions(wlc.wls[arrange->wl_index], arrange, o_path);
	};

	//init
	//get fixed positions and visited shapes
	std::vector<PCEPosition> fixed_group_positions = GetFixedGroupPositions(parent_order_node);
	std::vector<int> fixed_shapes;
	std::map<int, int> visited_shapes = GetVisitedParts(parent_order_node, fixed_shapes);
	
	if (wl.index == 0 && start_shape_index < 0&&fixed_group_positions.empty())
	{
		auto sub_shapes = PCEEncode::GetShapes(PCEEncode::Subtract(all_targeting_shapes, visited_shapes));
		auto sub_encode = Math::Functs::IntString(wlc.index)+ Math::Functs::IntString(wl.index)+ Math::Functs::IntString(sub_shapes, true,",");

		if (!Functs::VectorInsertNoDuplicate(enumerate_order_strs, sub_encode)) return;
	}


	//log
	//output visited shapes and all targeting shapes 

	if (depth_searching_bool)
	{
		if (wl.index == 0 || output_all_log)
		{
			std::cerr << parent_order_node.level << "," << parent_order_node.level_index << " : ";
			std::cerr << "(" << fixed_shapes.size() << ")" << Math::Functs::IntString(fixed_shapes, false, ",") << " : ";
			std::cerr << "(" << PCEEncode::GetShapes(all_targeting_shapes).size() << ")" << Math::Functs::IntString(PCEEncode::GetShapes(all_targeting_shapes), true, ",") << " : ";
		}
	}



	//check terminal condition
	if (visited_shapes.size() == PCEEncode::GetEncodeSize(all_targeting_shapes))return;

	std::vector<PCEConstraint> fixed_group_sharing_constraints = parent_order_node.sharing_constraints;
	std::vector<PCEEdge> fixed_group_edges = parent_order_node.edges;

	//OrderNode
	std::vector<OrderNode> next_order_nodes = GetNextPositions(fixed_group_positions, visited_shapes, fixed_group_sharing_constraints, fixed_group_edges);


	if (next_order_nodes.empty()&& !fixed_group_positions.empty())
	{
		if (depth_searching_bool)
		{
			if (wl.index == 0 || output_all_log)std::cerr << " :switch to new wood: ";
		}

		//save push
		std::vector<PCEEdge> cutting_lines;
		GetCuttingLines(parent_order_node.edges, parent_order_node.sharing_constraints, fixed_group_positions, cutting_lines);

		int nb = ct_egraph.ptr_arranges.size();

		parent_order_node.ct_arrange_eclass_index =  
			PCEEGraphContainer::CTNewArrange(ct_egraph.ptr_arranges,ct_egraph.ptr_egraph,
			wlc, wl, fixed_group_positions, cutting_lines,parent_order_node.sharing_constraints, parent_order_node.edges);
		
		if (ct_egraph.ptr_arranges.size() != nb && packing_o_folder.size() != 0 && CompilerConfig::Instance().GetObjOutput())
			OutputArrange(ct_egraph.ptr_arranges.back());

		fixed_group_positions.clear();
		fixed_group_sharing_constraints.clear();
		fixed_group_edges = wl.ExtractEdges(); //wl edges
		next_order_nodes = GetNextPositions(fixed_group_positions, visited_shapes, fixed_group_sharing_constraints, fixed_group_edges);

		//log
		if (depth_searching_bool)
		{
			if (wl.index == 0 || output_all_log)
			{
				std::cerr << "(" << next_order_nodes.size() << ")";
				for (auto& order_node : next_order_nodes)
					if (order_node.position.style != PCEPosition::FreeStyle)
						std::cerr << order_node.position.transform->pce_shape->index << "," << order_node.position.transform->index << "/";
				std::cerr << "" << std::endl;
			}

			if (ct_egraph.ptr_arranges.size() != nb)
				std::cerr << ":::::::::::New arrange: " << ct_egraph.ptr_arranges.back()->index << " shape: " << Math::Functs::IntString(ct_egraph.ptr_arranges.back()->shape_order, true, ",") << std::endl;
		}

	}
	else
	{
		//log
		if (depth_searching_bool)
		{
			if (wl.index == 0 || output_all_log)
			{
				std::cerr << "(" << next_order_nodes.size() << ")";
				for (auto& order_node : next_order_nodes)
					if (order_node.position.style != PCEPosition::FreeStyle)
						std::cerr << order_node.position.transform->pce_shape->index << "," << order_node.position.transform->index << "/";
				std::cerr << "" << std::endl;
			}
		}
	}

	//build order level
	int next_parent_level_index = -1;
	auto& order_level = order_levels[parent_order_node.level + 1];
	for (int i = 0; i < tree_searching_branch && i < next_order_nodes.size(); i++)
	{
		auto& next_position = next_order_nodes[i];

		if (next_parent_level_index < 0)
			next_parent_level_index = order_level.size();

		order_level.emplace_back(OrderNode(
			-1,
			PCEEncode::GetEncodeSize(visited_shapes),
			order_level.size(),
			parent_order_node.level_index,
			!((bool)i),
			next_position.position,
			next_position.sharing_constraints,
			next_position.edges));

		if (!depth_searching_bool)break;
	}

	//next
	auto visited_shapes_ = visited_shapes;
	if (!next_order_nodes.empty())
	{
		//select the first one
		auto& selected_node = order_level[next_parent_level_index];
		if (selected_node.position.style != PCEPosition::FreeStyle)
		{
			auto fixed_positions_ = fixed_group_positions;
			fixed_positions_.emplace_back(selected_node.position);
			PCEEncode::Insert(visited_shapes_, selected_node.position.transform->pce_shape->index);
			if (PCEEncode::GetEncodeSize(visited_shapes_) == PCEEncode::GetEncodeSize(all_targeting_shapes))
			{
				//save push
				std::vector<PCEEdge> cutting_lines;
				GetCuttingLines(selected_node.edges, selected_node.sharing_constraints, fixed_positions_, cutting_lines);
				//int nb = wlc.wlc_arranges.size();
				int nb = ct_egraph.ptr_arranges.size();

				selected_node.ct_arrange_eclass_index = 
					PCEEGraphContainer::CTNewArrange(ct_egraph.ptr_arranges, ct_egraph.ptr_egraph, 
						wlc, wl, fixed_positions_, cutting_lines, selected_node.sharing_constraints, selected_node.edges);

				if (depth_searching_bool)
				{
					if (ct_egraph.ptr_arranges.size() != nb)
						std::cerr << ":::::::::::New arrange: " << ct_egraph.ptr_arranges.back()->index << " shape: " << Math::Functs::IntString(ct_egraph.ptr_arranges.back()->shape_order, true, ",") << std::endl;
				}
			
				if (ct_egraph.ptr_arranges.size() != nb && packing_o_folder.size() != 0 && CompilerConfig::Instance().GetObjOutput())
					OutputArrange(ct_egraph.ptr_arranges.back());

				wl.all_enumerating_nb++;

				if (depth_searching_bool)
				{
					if (wl.index == 0 || output_all_log)
						std::cerr << "-----All_enumerating_nb: " << wl.all_enumerating_nb << " ------max: " <<maximal_arrange_nb<< std::endl;
					//else
						//std::cerr << wl.all_enumerating_nb << "/" << maximal_arrange_nb << std::endl;
				}

				//depth searching
				if(depth_searching_bool) DepthSearching(selected_node);
			}
			else
			{
				Order2Arrange(enumerate_order_strs,wlc, wl, selected_node, maximal_arrange_nb, all_targeting_shapes, shape_transforms, seq_shapes_transforms, order_levels, depth_searching_bool);
			}
		}
	}
}

//3. post process
void PartDesignGui::PCE::PostProcess(const std::vector<WLCollection>& wlcs_)
{
	auto OutputArrange = [&](const WLCollection& wlc, shared_ptr<PCEArrange> arrange)
	{
		std::string o_path = packing_o_folder +
			"\\wlc_" + Math::Functs::IntString(wlc.index) +
			"_wl_" + Math::Functs::IntString(arrange->wl_index) +
			"_arrange_" + Math::Functs::IntString(arrange->index) +
			"_valid_" + Math::Functs::IntString(arrange->valid) +
			"_shapes_" + Math::Functs::IntString(arrange->shape_order.size()) +
			"(" + Math::Functs::IntString(arrange->shape_order,false,",") + ")" +
			"_cutting_" + Math::Functs::IntString(arrange->cutting_number) + ".obj";

		OutputFixedPositions(wlc.wls[arrange-> wl_index], arrange, o_path);
	};

	auto OutputLog = [&]()
	{
		if (output_all_log)
		{
			std::cerr << "" << std::endl;
			for (auto& eclass : ct_egraph.ptr_egraph)
			{
				std::cerr << "***********" << "E-class: " << eclass->index << "***********" << std::endl;
				std::cerr << "Valid: " << eclass->valid << std::endl;
				std::cerr << "shapes: " << Math::Functs::IntString(eclass->ec_shapes, false, " ") << std::endl;
				std::cerr << "Comp sets: " << eclass->egraph_edges.size() << " : ";
				for (auto& edges : eclass->egraph_edges)
				{
					std::cerr << "(";
					for (auto& edge : edges) std::cerr << edge << " ";
					std::cerr << ")";
				}
				std::cerr << std::endl;

				std::cerr << "Arranges: " << eclass->arranges.size() << " index(wl_index):  ";
				for (auto& arrange : eclass->arranges)
					std::cerr << arrange << "(" << ct_egraph.ptr_arranges[arrange]->wl_index << ")";
				std::cerr << "\n";
			}
		}

		std::ofstream file(packing_o_folder + "\\egraph.txt");

		for (auto& head : ct_egraph.heads)
		{
			file << "Head index: " << head->index <<" head_index:"<<head->head_index<< " shapes: " << Math::Functs::IntString(head->ec_shapes, true, ",");
			file << " edges: " << head->egraph_edges.size() << " : ";
			for (auto& edges : head->egraph_edges)
				file << "(" << Math::Functs::IntString(edges, true, ",") << ")";

			file << " arranges: " << head->arranges.size() << " : " << Math::Functs::IntString(head->arranges, false, "/") << std::endl;
		}

		file << "\n\n";
		for (auto& eclass : ct_egraph.ptr_egraph)
		{
			file << "E-class: " << eclass->index 
				<< " head_index "<< eclass->head_index 
				<< " valid: " << eclass->valid
				<<" pre_eclass: "<<eclass->pre_eclass
				<< " shapes: " 
				<< Math::Functs::IntString(eclass->ec_shapes, true, ",");

			file << " edges: " << eclass->egraph_edges.size() << " : ";
			for (auto& edges : eclass->egraph_edges)
				file << "("<<Math::Functs::IntString(edges,true,",") << ")";

			file << " arranges: " << eclass->arranges.size() << " : " << Math::Functs::IntString(eclass->arranges, false, "/");
			
			file << std::endl;
		}

		file << "\n\n";

		for (auto& arrange :ct_egraph.ptr_arranges)
		{
			file << "Arrange: " << arrange->index<<" valid: "<<arrange->valid<<" eclass_index: "<<arrange->eclass_index <<" wlc: "<<arrange->wlc_index<<" wl:"<<arrange->wl_index;
			file << " shapes: " << Math::Functs::IntString(arrange->shape_order, true, ",") << std::endl;
		}

		file.close();
	};

	//wlc
		//wlc.arranges
	for (auto& arrange : ct_egraph.ptr_arranges)
		if (arrange->valid && CompilerConfig::Instance().GetObjOutput())
			OutputArrange(wlcs_[arrange->wlc_index], arrange);

	OutputLog();
}

void PCE::ClearAll()
{
	std::map<std::pair<int, int>, double>().swap(ptr_shapes_packing);

	number_label.Clear();
	ct_egraph.Clear();

	std::vector<std::string>().swap(enumerate_order_strs_);
}

int PCEEGraphContainer::CTNewEClass(std::vector<shared_ptr<EGraphClass>>& ptr_egraph_, Vector1i1& shape_order)
{
	string str = Math::Functs::IntString(shape_order, true,",");
	for (auto& eclass : ptr_egraph_)
		if (eclass->str == str)
			return eclass->index;

	ptr_egraph_.emplace_back(shared_ptr<EGraphClass>(new EGraphClass(ptr_egraph_.size(), shape_order)));
	
	return ptr_egraph_.back()->index;
}

std::pair<int, int> PCEEGraphContainer::CTNewArrange(
	std::vector<shared_ptr<PCEArrange>>& ptr_arranges_, 
	std::vector<shared_ptr<EGraphClass>>& ptr_egraph_, 
	const WLCollection& wlc_, const WL& wl_,
	const std::vector<PCEPosition>& fixed_positions_,
	const std::vector<PCEEdge>& cutting_lines_,
	const std::vector<PCEConstraint>& sharing_constraints_,
	const std::vector<PCEEdge>& edges_)
{
	double total_volume = 0.0;
	for (const auto& position : fixed_positions_)
		total_volume += position.transform->pce_shape->volume;
	double volume_rate = total_volume / wl_.volume;

	//input arrange
	shared_ptr<PCEArrange> arr(new PCEArrange(ptr_arranges_.size(), wlc_.index, wl_.index, fixed_positions_, cutting_lines_, volume_rate));

	int arrange_index = -1;
	for (auto& ptr_arrange_ : ptr_arranges_)
	{
		//if (ptr_arrange_->encode_0 == arr->encode_0 || ptr_arrange_->encode_1 == arr->encode_1 ||
		//	ptr_arrange_->encode_0 == arr->encode_1 || ptr_arrange_->encode_1 == arr->encode_0)
		//{
		//	arrange_index = ptr_arrange_->index;
		//	break;
		//}
		if (ptr_arrange_->encode_0 == arr->encode_0)
		{
			arrange_index = ptr_arrange_->index;
			break;
		}

	}

	if (arrange_index < 0)
	{
		ptr_arranges_.emplace_back(arr);
		arrange_index = ptr_arranges_.back()->index;

		ptr_arranges_.back()->edges = edges_;
		ptr_arranges_.back()->sharing_constraints = sharing_constraints_;
	}

	//input egraph
	auto ptr_arrange_ = ptr_arranges_[arrange_index];

	if (ptr_arrange_->eclass_index < 0)
	{
		ptr_arrange_->eclass_index = CTNewEClass(ptr_egraph_, ptr_arrange_->shape_order);
		ptr_egraph_[ptr_arrange_->eclass_index]->arranges.emplace_back(ptr_arrange_->index);

		if (!ptr_egraph_[ptr_arrange_->eclass_index]->xmlDoc_path.empty()) ptr_egraph_[ptr_arrange_->eclass_index]->xmlDoc_path = "";
	}
	return std::pair<int, int>(ptr_arrange_->index, ptr_arrange_->eclass_index);
}

void PCEEGraphContainer::Clear()
{
	for (auto& ptr_arrange : ptr_arranges)ptr_arrange->Clear();
	std::vector<shared_ptr<PCEArrange>>().swap(ptr_arranges);
	for (auto& ptr_ec : ptr_egraph)ptr_ec->Clear();
	std::vector<shared_ptr<EGraphClass>>().swap(ptr_egraph);
	std::vector<int>().swap(full_part_arranges);
}

void PartDesignGui::PCEEGraphContainer::OutputArrange(
	const std::vector<WLCollection>& wlcs_, 
	const shared_ptr<PCEArrange>& ptr_arrange_, 
	std::vector<TopoDS_Shape>& transfShapes, 
	std::vector<PCEEdge>& cutting_lines, 
	std::vector<PCEHole>& holes, 
	std::vector<std::vector<Vector2d>>& offset_polygons, 
	Vector3d& woods, boost::uuids::uuid& uid)
{
	cutting_lines = ptr_arrange_->cutting_lines;

	woods = wlcs_[ptr_arrange_->wlc_index].wls[ptr_arrange_->wl_index].size;
	uid = wlcs_[ptr_arrange_->wlc_index].wls[ptr_arrange_->wl_index].uid;

	for (auto position : ptr_arrange_->positions)
	{
		//auto& part = parts[position.part_index];
		auto& shape = position.transform->pce_shape;

		glm::dmat4 M = Math::Functs::TranslationMatrix(Vector3d(position.t.x - position.transform->tranform_corner[0], position.t.y - position.transform->tranform_corner[1], -position.transform->tranform_corner[2]))
			* Math::Functs::TranslationMatrix(Vector3d(position.transform->center[0], position.transform->center[1], position.transform->center[2]))
			* Math::Functs::RotationMatrix(Vector3d(0.0, 0.0, 1.0), -position.transform->rotation[2])
			* Math::Functs::RotationMatrix(Vector3d(0.0, 1.0, 0.0), -position.transform->rotation[1])
			* Math::Functs::RotationMatrix(Vector3d(1.0, 0.0, 0.0), -position.transform->rotation[0])
			* Math::Functs::TranslationMatrix(Vector3d(-position.transform->center[0], -position.transform->center[1], -position.transform->center[2]));

		const auto trsf = GeomFunc::ConvertTogpTrsf(M);
		// Transform shape
		TopoDS_Shape curShape = shape->ashape;
		curShape = BRepBuilderAPI_GTransform(curShape, trsf).Shape();
		transfShapes.emplace_back(curShape);

		//offset
		std::vector<Vector2d> offset_polygon;
		for (auto p : position.transform->cut_ends.ul_offset_2d)
		{
			glm::vec4 v(Vector3d(p.x, p.y, 0.0), 1.0);
			Vector3d t_v = Vector3d(position.M * v);
			offset_polygon.emplace_back(Vector2d(t_v.x, t_v.y));
		}
		offset_polygons.emplace_back(offset_polygon);

		//holes
		if (!shape->holes.empty())
		{
			for (auto& hole : shape->holes)
			{
				Vector3d center_0 = Vector3d(M * glm::vec4(hole.first[0], 1.0));
				Vector3d center_1 = Vector3d(M * glm::vec4(hole.first[1], 1.0));
				holes.emplace_back(std::vector<Vector3d>{center_0, center_1}, hole.second, transfShapes.size() - 1);
			}
		}
	}
}

shared_ptr<PartDesignGui::EGraphClass> PCEEGraphContainer::GetHead(const int& head_index)
{
	for (auto& head : heads)
	{
		if (head->index == head_index)
			return head;
	}
	return heads.front();
}

void PCEEGraphContainer::Pruning(const int head_index, const int& prunning_limitation, const bool& random_pruning)
{
	for (auto& eclass : ptr_egraph)
	{
		if (eclass->head_index == head_index)
		{
			eclass->Pruning(ptr_arranges, prunning_limitation, random_pruning);
		}
	}
}

//without check valid
Vector1i1 PCEEGraphContainer::GetHeadArranges(const int& head_index)
{
	Vector1i1 arranges;
	for (const auto& eclass : ptr_egraph)
	{
		if (eclass->head_index == head_index)
			for (auto& arrange : eclass->arranges)
				Math::Functs::VectorInsertNoDuplicate(arranges, ptr_arranges[arrange]->index);
	}
	return arranges;
}

//without check valid
Vector1i1 PCEEGraphContainer::GetHeadEClasses(const int& head_index)
{
	Vector1i1 eclasses;
	for (const auto& eclass : ptr_egraph)
		if (eclass->head_index == head_index)
			Math::Functs::VectorInsertNoDuplicate(eclasses, eclass->index);
	return eclasses;
}

void PCEEGraphContainer::SetArrangeValid(const int& head_index)
{
	for (const auto& eclass : ptr_egraph)
	{
		if (eclass->valid && eclass->head_index == head_index)
			for (auto& arrange : eclass->arranges)
				ptr_arranges[arrange]->valid = true;
	}
}

void PCEEGraphContainer::SetEclassValid(shared_ptr<EGraphClass> eclass)
{
	eclass->valid = true;
	for (auto edges : eclass->egraph_edges)
		for (auto edge : edges)
			if (!ptr_egraph[edge]->valid)
				SetEclassValid(ptr_egraph[edge]);
}

void PCEEGraphContainer::GetRelatedEclass(const std::vector<shared_ptr<EGraphClass>>& ptr_egraph, const shared_ptr<EGraphClass> eclass, Vector1i1& rel_ecs)
{
	if(std::find(rel_ecs.begin(),rel_ecs.end(),eclass->index)==rel_ecs.end())
		rel_ecs.emplace_back(eclass->index);
	for (auto edges : eclass->egraph_edges)
		for (auto edge : edges)
			GetRelatedEclass(ptr_egraph,ptr_egraph[edge], rel_ecs);
}

void PCEEGraphContainer::ComputeMoreSubSets()
{
	auto EClassUnion = [](shared_ptr<EGraphClass> target_ec, std::vector<shared_ptr<EGraphClass>> union_ecs)
	{
		auto& target_shapes = target_ec->encode_shapes;

		std::vector<std::map<int, int>> encodes;
		for (auto& union_ec : union_ecs) encodes.emplace_back(union_ec->encode_shapes);
		std::map<int, int> union_shapes = PCEEncode::Union(encodes);

		return PCEEncode::CheckIdentify(union_shapes, target_shapes);
	};

	std::vector<std::vector<int>> sub_sets(ptr_egraph.size(), std::vector<int>());

	for (auto i = 1; i < ptr_egraph.size(); i++)
	{
		auto& eclass_0 = ptr_egraph[i];
		for (auto j = i + 1; j < ptr_egraph.size(); j++)
		{
			auto& eclass_1 = ptr_egraph[j];

			if (PCEEncode::CheckInclude(eclass_0->encode_shapes, eclass_1->encode_shapes))
			{
				Math::Functs::VectorInsertNoDuplicate(sub_sets[eclass_0->index], eclass_1->index);
			}

			if (PCEEncode::CheckInclude(eclass_1->encode_shapes, eclass_0->encode_shapes))
			{
				Math::Functs::VectorInsertNoDuplicate(sub_sets[eclass_1->index], eclass_0->index);
			}
		}
	}

	for (auto& eclass : ptr_egraph)
	{
		if (eclass->index == 0)continue;

		//2
		for (int i = 0; i < sub_sets[eclass->index].size(); i++)
		{
			int sub_set_0 = sub_sets[eclass->index][i];
			for (int j = 0; j < sub_sets[eclass->index].size(); j++)
			{
				auto sub_set_1 = sub_sets[eclass->index][j];

				std::vector<shared_ptr<EGraphClass>> union_2_ecs = { ptr_egraph[sub_set_0], ptr_egraph[sub_set_1] };
				if (EClassUnion(eclass, union_2_ecs))
				{
					std::vector<int> edges;
					edges.emplace_back(sub_set_0);
					edges.emplace_back(sub_set_1);
					eclass->PushComp(edges);
				}

				//three sub sets
				for (int k = 0; k < sub_sets[eclass->index].size(); k++)
				{
					int sub_set_2 = sub_sets[eclass->index][k];

					std::vector<shared_ptr<EGraphClass>> union_3_ecs = { ptr_egraph[sub_set_0], ptr_egraph[sub_set_1], ptr_egraph[sub_set_2] };

					if (EClassUnion(eclass, union_3_ecs))
					{
						std::vector<int> edges;
						edges.emplace_back(sub_set_0);
						edges.emplace_back(sub_set_1);
						edges.emplace_back(sub_set_2);
						eclass->PushComp(edges);
					}
				}
			}
		}
	}
}

void PCE::ComputeShapeScore(
	std::vector<WLCollection>& wlcs_,
	const std::vector<shared_ptr<PCEShape>>& ptr_shapes_,
	std::map<std::pair<int, int>, double>& ptr_shapes_packing_)
{
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

	for (auto shape_0 : ptr_shapes_)
	{
		for (auto shape_1 : ptr_shapes_)
		{
			if (shape_0->WLC_index != shape_1->WLC_index)continue;

			std::pair<int, int> shape_pair(shape_0->index, shape_1->index);
			std::pair<int, int> shape_pair_inverse(shape_1->index, shape_0->index);

			if (ptr_shapes_packing_.find(shape_pair) != ptr_shapes_packing_.end()) continue;
			if (ptr_shapes_packing_.find(shape_pair_inverse) != ptr_shapes_packing_.end()) continue;

			auto& wl = wlcs_[shape_0->WLC_index].wls[0];

			double maximal_score = 0;
			int score_nb = 0;
			std::cerr << "Shape index: " << shape_0->index << " " << shape_1->index << std::endl;

			for (auto& transform_0 : shape_0->transforms)
			{
				for (auto& transform_1 : shape_1->transforms)
				{
					auto a = Math::Functs::GetCenter(wl.offset_points);
					auto b = Math::Functs::GetCenter(transform_0.cut_ends.ul_offset_2d);

					PCEPosition position;
					position.transform = &transform_0;
					std::vector<PCEConstraint> constraints;
					position.t = a - b;;
					position.style = PCEPosition::FixedStyle_2_Constraints;
					position.ComputeM();
					auto edges = wl.ExtractEdges();
					position.transform->AddNewEdges(edges, position.edge_indexes);
					position.UpdateEdges(edges);
					position.score = 0.0;

					std::vector<PCEPosition> positions(1,position);
					PCEPosition second_position;
					second_position.transform = &transform_1;
					second_position.transform->AddNewEdges(edges, second_position.edge_indexes);
					ArrangePartbyConstraints(wl, edges, positions, second_position, constraints);

					std::vector<int> full_unique_l;
					std::vector<std::vector<int>> full_unique_ls;
					CombineEdges(edges, constraints, full_unique_l, full_unique_ls);

					//only consider the inner edges 
					int shared_edges=0;
					for (int i = 0; i < full_unique_l.size(); i++)
					{
						if (full_unique_l[i] == -1)continue;

						bool goon = true;
						for (int j = 0; j < full_unique_ls[i].size(); j++)
						{
							if (full_unique_ls[i][j] < 4)
							{
								goon = false;
								break;
							}
						}
						if (goon)
						{
							shared_edges = shared_edges + full_unique_ls[i].size()-1;
						}
					}

					//std::cerr << "Transform index: " << transform_0.index << " " << transform_1.index << " : " << shared_edges << std::endl;;
					if (second_position.style == PCEPosition::FreeStyle)
					{
						std::cerr << "if (second_position.style == PCEPosition::FreeStyle)" << std::endl;
						system("pause");
					}

					maximal_score += shared_edges;
					score_nb++;
					//if (shared_edges > maximal_score)
					{
					
						//output arrange
						positions.emplace_back(second_position);
						std::vector<PCEEdge> cutting_lines;
						GetCuttingLines(edges, constraints, positions, cutting_lines);
						auto path = packing_o_folder + "\\score_"
							+ std::to_string(shape_0->index) + "_" 
							+ std::to_string(shape_1->index) + "_"
							+std::to_string(transform_0.index)+"_"
							+std::to_string(transform_1.index)+"_score_"
							+std::to_string(shared_edges)+".obj";

						OutputFixedPositions(wl, positions, cutting_lines, path);
					}
				}
			}

			//void PCE::OutputFixedPositions(const WL& wl, const shared_ptr<PCEArrange> arrange, const std::string& path)
			ptr_shapes_packing_[shape_pair] = maximal_score/ score_nb;

			int dsad = 0;

		}
	}
}

void NumberLabel::DrawNode(const int style, std::ofstream& export_fie, int& export_int, const std::string& name, const double radius, const Vector3d& label_center)
{
	if (!load) Init();
	//enum { BaseCircle, BaseSquare, FullSetup, HalfAngle, HalfOffset};

	Vector3d1 style_vecs;
	Vector1i2 style_faces;

	switch (style)
	{
	case BaseCircle: {
		style_vecs = base_circle_vecs;
		style_faces = base_circle_faces;
		break;
	}
	case BaseShere: {
		style_vecs = sphere_vecs;
		style_faces = sphere_faces;
		break;
	}
	case BaseSquare1: {
		style_vecs = base_square_1_vecs;
		style_faces = base_square_1_faces;
		break;
	}
	case BaseSquare2: {
		style_vecs = base_square_2_vecs;
		style_faces = base_square_2_faces;
		break;
	}
	case FullSetup: {
		style_vecs = full_setup_vecs;
		style_faces = full_setup_faces;
		break;
	}
	case HalfAngle: {
		style_vecs = half_angle_setup_vecs;
		style_faces = half_angle_setup_faces;
		break;
	}
	case HalfOffset: {
		style_vecs = half_offset_setup_vecs;
		style_faces = half_offset_setup_faces;
		break;
	} 
	default:
		break;
	}

	for (auto& vec : style_vecs)
	{
		auto v = vec * radius + label_center;
		export_fie << fixed << setprecision(8) << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}

	export_fie << "g " << name << std::endl;
	for (int iter_ = 0; iter_ < style_faces[0].size(); iter_++)
		export_fie << "f " << style_faces[0][iter_] + export_int << " " << style_faces[1][iter_] + export_int << " " << style_faces[2][iter_] + export_int << std::endl;
	export_int = export_int + style_vecs.size();
}

void NumberLabel::DrawNumber1(std::ofstream& export_fie, int& export_int, const std::string& name, const double value, const double radius, const Vector3d& label_center)
{
	if (!load) Init();

	auto GetNumb = [](double number_double)
	{
		int number = number_double;
		Vector1i1 num;
		int unit_place = 1;
		int all = 0;
		while (true)
		{
			int unitPlace1 = number / unit_place % 10;
			if (unitPlace1 == 0 && all == number)break;
			all += unit_place * unitPlace1;
			num.emplace_back(unitPlace1);
			unit_place = unit_place * 10;
		}

		if (num.size() == 0) num.emplace_back(0);
		
		return num;
	};

	double integer_double = std::floor(value);
	double decimal_double = value - integer_double;

	auto integer_num = GetNumb(integer_double);
	auto decimal_num = GetNumb(decimal_double*1000.0);
	int dot_index = decimal_num.size();
	for (auto& n : integer_num) decimal_num.emplace_back(n);

	auto num = decimal_num;

	double d = 1.0;
	bool goon = true;
	for (int i = 0; i < num.size(); i++)
	{
		double delta_x = -i * d + (num.size() * d) / 2.0 - 0.5 * d;

		if (i == dot_index - 1)
		{
			auto dot_center = label_center;
			dot_center[0] += (delta_x + 0.5 * d) * radius;
			dot_center[1] += -0.5 * d * radius;
			DrawNode(BaseCircle, export_fie, export_int, name, radius * 0.07, dot_center);
		}

		if (i == num.size() - 2)goon = false;

		auto number_ = num[i];

		if (number_ == 0 && goon)
			continue;
		else
			goon = false;

		for (auto base_vec : numbers_vecs[number_])
		{
			base_vec[0] = base_vec[0] + delta_x;
			auto v = (base_vec)*radius + label_center;
			export_fie << fixed << setprecision(8) << "v " << v[0] << " " << v[1] << " " << v[2] << " 1.0 0.0 0.0" << std::endl;
		}

		export_fie << "g " << name << std::endl;

		for (int iter_ = 0; iter_ < numbers_faces[number_][0].size(); iter_++)
		{
			export_fie << "f " << numbers_faces[number_][0][iter_] + export_int << " " << numbers_faces[number_][1][iter_] + export_int << " " << numbers_faces[number_][2][iter_] + export_int << std::endl;
		}
		export_int = export_int + numbers_vecs[number_].size();

	}


}

void NumberLabel::DrawNumber(std::ofstream& export_fie, int& export_int, const std::string& name, const double value, const double radius, const Vector3d& label_center)
{
	//if (!load) Init();

	//int number = value * 10000;

	////number
	//Vector1i1 num;
	//int unit_place = 1;
	//int all = 0;
	//while (true)
	//{
	//	int unitPlace1 = number / unit_place % 10;
	//	if (unitPlace1 == 0 && all == number)break;
	//	all += unit_place * unitPlace1;
	//	num.emplace_back(unitPlace1);
	//	unit_place = unit_place * 10;
	//}

	//if (number == 0)
	//{ 
	//	num.emplace_back(0);
	//	num.emplace_back(0);
	//}

	//double d = 1.0;
	//bool goon = true;
	//for (int i = 0; i < num.size(); i++)
	//{
	//	double delta_x = -i * d + (num.size() * d) / 2.0 - 0.5 * d;

	//	if (i == num.size()-1)
	//	{
	//		auto dot_center = label_center;
	//		dot_center[0] += (delta_x + 0.5 * d) * radius;
	//		dot_center[1] += -0.5 * d * radius;
	//		DrawNode(BaseCircle, export_fie, export_int, name, radius*0.07, dot_center);
	//	}

	//	if (i == num.size() - 2)goon = false;

	//	auto number_ = num[i];

	//	if (number_ == 0&&goon) 
	//		continue; 
	//	else
	//		goon=false;

	//	for (auto base_vec : numbers_vecs[number_])
	//	{
	//		base_vec[0] = base_vec[0] + delta_x;
	//		auto v = (base_vec)*radius + label_center;
	//		export_fie << fixed << setprecision(8) << "v " << v[0] << " " << v[1] << " " << v[2] << " 1.0 0.0 0.0" << std::endl;
	//	}

	//	export_fie << "g " << name << std::endl;

	//	for (int iter_ = 0; iter_ < numbers_faces[number_][0].size(); iter_++)
	//	{
	//		export_fie << "f " << numbers_faces[number_][0][iter_] + export_int << " " << numbers_faces[number_][1][iter_] + export_int << " " << numbers_faces[number_][2][iter_] + export_int << std::endl;
	//	}
	//	export_int = export_int + numbers_vecs[number_].size();

	//}


	if (!load) Init();

	auto GetNumb = [](double number_double)
	{
		int number = number_double;
		Vector1i1 num;
		int unit_place = 1;
		int all = 0;
		while (true)
		{
			int unitPlace1 = number / unit_place % 10;
			if (unitPlace1 == 0 && all == number)break;
			all += unit_place * unitPlace1;
			num.emplace_back(unitPlace1);
			unit_place = unit_place * 10;
		}

		if (num.size() == 0) num.emplace_back(0);

		return num;
	};

	double integer_double = std::floor(value);
	double decimal_double = value - integer_double;

	auto integer_num = GetNumb(integer_double);
	auto decimal_num = GetNumb(decimal_double * 1000.0);
	int dot_index = decimal_num.size();
	for (auto& n : integer_num) decimal_num.emplace_back(n);

	auto num = decimal_num;

	double d = 1.0;
	bool goon = true;
	for (int i = 0; i < num.size(); i++)
	{
		double delta_x = -i * d + (num.size() * d) / 2.0 - 0.5 * d;

		if (i == dot_index)
		{
			auto dot_center = label_center;
			dot_center[0] += (delta_x + 0.5 * d) * radius;
			dot_center[1] += -0.5 * d * radius;
			DrawNode(BaseCircle, export_fie, export_int, name, radius * 0.07, dot_center);
		}

		if (i == num.size() - 2)goon = false;

		auto number_ = num[i];

		if (number_ == 0 && goon)
			continue;
		else
			goon = false;

		for (auto base_vec : numbers_vecs[number_])
		{
			base_vec[0] = base_vec[0] + delta_x;
			auto v = (base_vec)*radius + label_center;
			export_fie << fixed << setprecision(8) << "v " << v[0] << " " << v[1] << " " << v[2] << " 1.0 0.0 0.0" << std::endl;
		}

		export_fie << "g " << name << std::endl;

		for (int iter_ = 0; iter_ < numbers_faces[number_][0].size(); iter_++)
		{
			export_fie << "f " << numbers_faces[number_][0][iter_] + export_int << " " << numbers_faces[number_][1][iter_] + export_int << " " << numbers_faces[number_][2][iter_] + export_int << std::endl;
		}
		export_int = export_int + numbers_vecs[number_].size();

	}
}

void NumberLabel::DrawLabel(const int base_style, std::ofstream& export_fie, int& export_int, const std::string& name, const int number, const double radius, const Vector3d& label_center)
{
	if (!load) Init();

	//base
	DrawNode(base_style, export_fie, export_int, name, radius, label_center);

	//number

	Vector1i1 num;
	int unit_place = 1;
	int all = 0;
	while (true)
	{
		int unitPlace1 = number / unit_place % 10;
		if (unitPlace1 == 0 && all == number)break;
		all += unit_place * unitPlace1;
		num.emplace_back(unitPlace1);
		unit_place = unit_place * 10;
	}

	if (number == 0)num.emplace_back(0);

	double d = 1.0;
	for (int i = 0; i < num.size(); i++)
	{
		double delta_x = -i * d + (num.size() * d) / 2.0-0.5*d;

		auto number_ = num[i];
		for (auto base_vec : numbers_vecs[number_])
		{
			base_vec[0] = base_vec[0] + delta_x;
			auto v = (base_vec)*radius + label_center;
			export_fie << fixed << setprecision(8) << "v " << v[0] << " " << v[1] << " " << v[2] <<" 1.0 0.0 0.0"<< std::endl;
		}

		export_fie << "g " << name << std::endl;

		for (int iter_ = 0; iter_ < numbers_faces[number_][0].size(); iter_++)
		{
			export_fie << "f " << numbers_faces[number_][0][iter_] + export_int << " " << numbers_faces[number_][1][iter_] + export_int << " " << numbers_faces[number_][2][iter_] + export_int << std::endl;
		}
		export_int = export_int + numbers_vecs[number_].size();

	}
}

void NumberLabel::Init()
{
	std::cerr<<"NumberLabel Initialization..."<<std::endl;

	base_circle_faces = Vector1i2(3, std::vector<int>());
	sphere_faces = Vector1i2(3, std::vector<int>());
	base_square_1_faces = Vector1i2(3, std::vector<int>());
	base_square_2_faces = Vector1i2(3, std::vector<int>());
	full_setup_faces = Vector1i2(3, std::vector<int>());
	half_angle_setup_faces = Vector1i2(3, std::vector<int>());
	half_offset_setup_faces = Vector1i2(3, std::vector<int>());
	numbers_vecs = Vector3d2(10, Vector3d1());
	numbers_faces = Vector1i3(10, Vector1i2(3, std::vector<int>()));

	std::string pref_path = "../../../../../src/Mod/PartDesign/Ext/obj/";


	std::string cur_path = _pgmptr;
	cur_path = cur_path.substr(0, cur_path.find_last_of('\\'));

	if (!Functs::DetectExisting(cur_path + "\\Ext\\obj\\"))
	{
		std::cerr << "std::string cur_path = _pgmptr" << std::endl;
		system("pause");
	}

	pref_path = cur_path + "\\Ext\\obj\\";

	auto read_mesh = (CGAL_3D_Read_Triangle_Mesh)GetProcAddress(CompilerConfig::Instance().GetHModule(), "CGAL_3D_Read_Triangle_Mesh");
	read_mesh(pref_path + "base_circle.obj", base_circle_vecs, base_circle_faces[0], base_circle_faces[1], base_circle_faces[2]);
	read_mesh(pref_path + "sphere.obj", sphere_vecs, sphere_faces[0], sphere_faces[1], sphere_faces[2]);
	read_mesh(pref_path + "base_square_1.obj", base_square_1_vecs, base_square_1_faces[0], base_square_1_faces[1], base_square_1_faces[2]);
	read_mesh(pref_path + "base_square_2.obj", base_square_2_vecs, base_square_2_faces[0], base_square_2_faces[1], base_square_2_faces[2]);
	read_mesh(pref_path + "full_setup.obj", full_setup_vecs, full_setup_faces[0], full_setup_faces[1], full_setup_faces[2]);
	read_mesh(pref_path + "half_angle_setup.obj", half_angle_setup_vecs, half_angle_setup_faces[0], half_angle_setup_faces[1], half_angle_setup_faces[2]);
	read_mesh(pref_path + "half_offset_setup.obj", half_offset_setup_vecs, half_offset_setup_faces[0], half_offset_setup_faces[1], half_offset_setup_faces[2]);

	for (int i = 0; i < 10; i++)
	{
		read_mesh(pref_path + std::to_string(i) + ".obj", numbers_vecs[i], numbers_faces[i][0], numbers_faces[i][1], numbers_faces[i][2]);

		Vector3d min_corner,max_corner;
		Math::Functs::GetBoundingBox(numbers_vecs[i], min_corner, max_corner);

		auto t = (min_corner + max_corner) / 2.0;
		t[2] = 0.0;
		for (auto& vec : numbers_vecs[i])
		{ vec = vec - t; }
	}
	load = true;
}
