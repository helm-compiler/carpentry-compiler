#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>
#include <CompilerConfig.h>

using namespace PartDesignGui;

PCETransform::PCETransform(
	Vector3d2 surfs_,
	Vector3d1 surf_normals_,
	Vector3d center_,
	Vector3d rotation_,
	glm::dmat4 PM_,
	glm::dmat4 RM_,
	Vector3d tranform_corner_) :
	surfs(surfs_),
	surf_normals(surf_normals_),
	center(center_),
	rotation(rotation_),
	PM(PM_),
	RM(RM_),
	tranform_corner(tranform_corner_) {

};

PCETransform::PCETransform()
{
	index = -1;
}

bool PCETransform::CheckWL(const int& wl_index_) const
{
	for (auto wl_rate : wl_rates)
	{
		if (wl_rate.wl->index == wl_index_)
			return true;
	}
	return false;
}

void PCETransform::Clear()
{
	std::vector<WLRATE>().swap(wl_rates);
	cut_ends.Clear();
	Vector3d2().swap(surfs);
	Vector3d1().swap(surf_normals);
}

void PCETransform::AddNewEdges(std::vector<PCEEdge>& edges_, std::vector<int>& edge_indexes)
{
	edge_indexes.clear();
	if (cut_ends.valid)
	{
		for (const auto& end : cut_ends.ends)
		{
			edges_.emplace_back(PCEEdge(edges_.size(), pce_shape->index, end.cut));
			edge_indexes.emplace_back(edges_.back().index);
		}
	}
	else
		AddBoundingEdges(edges_, edge_indexes);
}

void PCETransform::AddBoundingEdges(std::vector<PCEEdge>& edges_, std::vector<int>& edge_indexes)
{
	edge_indexes.clear();

	auto& offset_points = bounding_box.offset_2d;
	auto& origin_points = bounding_box.origin_2d;

	double height = wl_rates[0].wl->size[2];
	for (int i = 0; i < offset_points.size(); i++)
	{
		const int index_0 = i;
		const int index_1 = (i + 1) % offset_points.size();

		const Vector2d origin_s = origin_points[index_0];
		const Vector2d origin_e = origin_points[index_1];
		const Vector2d offset_s = offset_points[index_0];
		const Vector2d offset_e = offset_points[index_1];

		PCECUT cut;
		cut.origin_s = origin_s;
		cut.origin_e = origin_e;
		cut.offset_s = offset_s;
		cut.offset_e = offset_e;
		cut.cutting_surface_normal = bounding_box.normals[i];
		cut.u_cutting_line_s = Vector3d(offset_s[0], offset_s[1], height);
		cut.u_cutting_line_e = Vector3d(offset_e[0], offset_e[1], height);
		cut.l_cutting_line_s = Vector3d(offset_s[0], offset_s[1], 0.0);
		cut.l_cutting_line_e = Vector3d(offset_e[0], offset_e[1], 0.0);

		cut.l_u_vector = Vector3d(0.0, 0.0, height);

		edges_.emplace_back(PCEEdge(edges_.size(), pce_shape->index, cut));
		edge_indexes.emplace_back(edges_.back().index);
	}
}

//only for collision detection
Vector2d1 PCETransform::UOffset(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.u_offset_2d : bounding_box.offset_2d;
	for (auto& p : points)p += t;
	return points;
}

//only for collision detection
Vector2d1 PCETransform::LOffset(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.l_offset_2d : bounding_box.offset_2d;
	for (auto& p : points)p += t;
	return points;
}

//only for collision detection
Vector2d1 PCETransform::UOrigin(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.u_origin_2d : bounding_box.origin_2d;
	for (auto& p : points) p += t;
	return points;
}

//only for collision detection
Vector2d1 PCETransform::LOrigin(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.l_origin_2d : bounding_box.origin_2d;
	for (auto& p : points) p += t;
	return points;
}

Vector2d1 PCETransform::ULOffset(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.ul_offset_2d : bounding_box.offset_2d;
	for (auto& p : points)p += t;
	return points;
}

Vector2d1 PCETransform::ULOrigin(const Vector2d t) const
{
	auto points = cut_ends.valid ? cut_ends.ul_origin_2d : bounding_box.origin_2d;
	for (auto& p : points)p += t;
	return points;
}

PCETransform::CutEnds PCETransform::BuildEnd(const double& kerf, const HMODULE& hModule,
	const Vector3d2& surfs_, Vector3d upper_normal, Vector3d lower_normal)
{
	//build data structure
	struct Node
	{
	public:
		enum { EdgeNode, SurfNode };
		enum { UpperEdge, LowerEdge };

		int index, style;
		std::vector<int> neighbors;
		std::vector<int> ends;

		int surface_index;

		bool used;

		//edge
		Vector3d edge_s, edge_e;
		int edge_style;
		int friend_node;
		int edge_index;

		//surface
		Vector3d1 surf_points;
		Vector3d surf_normal;

		//EdgeNode
		Node(int index_, Vector3d s_, Vector3d e_, int edge_style_, int edge_index_, int surface_index_) :
			index(index_), style(EdgeNode), edge_s(s_), edge_e(e_), edge_style(edge_style_), edge_index(edge_index_), surface_index(surface_index_)
		{
			friend_node = -1;
			used = false;
		};

		//SurfNode
		Node(int index_, const Vector3d1& points_, const Vector3d& surf_normal_, int surface_index_) :
			index(index_), style(SurfNode), surf_points(points_), surf_normal(surf_normal_), surface_index(surface_index_)
		{
			edge_style = -1;
			friend_node = -1;
			used = false;
		};

		void PushEnd(int end)
		{
			if (Math::Functs::VectorIndex(ends, end) == -1) ends.emplace_back(end);
		}

		void PushEnd(std::vector<int> ints)
		{
			for (auto i : ints) PushEnd(i);
		}
	};
	struct Edge
	{
	public:
		int index;
		int node_0, node_1;

		Edge(int index_, int n_0, int n_1) : index(index_), node_0(n_0), node_1(n_1)
		{
		};
	};
	auto CheckConnecting0 = [&](const Vector3d& s, const Vector3d& e, const Vector3d1& contour)
	{
		for (int i = 0; i < contour.size(); i++)
		{
			Vector3d s_ = contour[i];
			Vector3d e_ = contour[(i + 1) % contour.size()];
			bool ss = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(s, s_));
			bool se = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(s, e_));
			bool es = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(e, s_));
			bool ee = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(e, e_));
			if ((ss && ee) || (se && es)) return true;
		}
		return false;
	};
	auto CheckConnecting1 = [&](Vector3d1& contour0, Vector3d1& contour1)
	{
		for (int i = 0; i < contour0.size(); i++)
		{
			Vector3d s = contour0[i];
			Vector3d e = contour0[(i + 1) % contour0.size()];

			for (int j = 0; j < contour1.size(); j++)
			{
				Vector3d s_ = contour1[j];
				Vector3d e_ = contour1[(j + 1) % contour1.size()];
				bool ss = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(s, s_));
				bool se = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(s, e_));
				bool es = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(e, s_));
				bool ee = Math::Functs::IsAlmostZero(Math::Functs::GetDistance(e, e_));
				if ((ss && ee) || (se && es)) return true;
			}
		}
		return false;
	};

	//searching for the upper surface and the lower surface
	auto BuildNodes = [&](std::vector<Node>& nodes, Vector3d1& surf_normals_, PCETransform::CutEnds& cut_ends_) {

		bool ub(false), lb(false);
		for (int i = 0; i < surf_normals_.size(); i++)
		{
			if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(upper_normal, surf_normals_[i])))
			{
				if (ub)
				{
					std::cerr << "Find not single surfaces which are upper surfaces..." << std::endl;
					system("pause");
					return false;
				}

				cut_ends_.u_surf = surfs_[i];
				cut_ends_.u_surf_index = i;

				ub = true;
				for (int j = 0; j < cut_ends_.u_surf.size(); j++)
					nodes.emplace_back(Node(nodes.size(), cut_ends_.u_surf[j], cut_ends_.u_surf[(j + 1) % cut_ends_.u_surf.size()], Node::UpperEdge, j, i));
			}
			else
			{
				if (Math::Functs::IsAlmostZero(Math::Functs::GetAngleBetween(lower_normal, surf_normals_[i])))
				{
					if (lb)
					{
						std::cerr << "Find not single surfaces which are lower surfaces..." << std::endl;
						system("pause");
						return false;
					}

					cut_ends_.l_surf = surfs_[i];
					cut_ends_.l_surf_index = i;

					lb = true;
					for (int j = 0; j < cut_ends_.l_surf.size(); j++)
						nodes.emplace_back(Node(nodes.size(), cut_ends_.l_surf[j], cut_ends_.l_surf[(j + 1) % cut_ends_.l_surf.size()], Node::LowerEdge, j, i));
				}
				else
				{
					nodes.emplace_back(Node(nodes.size(), surfs_[i], surf_normals_[i], i));
				}
			}
		}
		if (!ub || !lb) return false;
		return true;
	};

	//build edges
	auto BuildEdges = [&](std::vector<Node>& nodes)
	{
		std::vector<Edge> edges;
		for (auto& n : nodes)
		{
			for (auto& nn : nodes)
			{
				if (n.index < nn.index && n.style == Node::EdgeNode && nn.style == Node::SurfNode)
				{
					if (CheckConnecting0(n.edge_s, n.edge_e, nn.surf_points))
					{
						edges.emplace_back(edges.size(), n.index, nn.index);
						n.neighbors.emplace_back(nn.index);
						nn.neighbors.emplace_back(n.index);
					}
				}

				if (n.index < nn.index && n.style == Node::SurfNode && nn.style == Node::EdgeNode)
				{
					if (CheckConnecting0(nn.edge_s, nn.edge_e, n.surf_points))
					{
						edges.emplace_back(edges.size(), n.index, nn.index);
						n.neighbors.emplace_back(nn.index);
						nn.neighbors.emplace_back(n.index);
					}
				}

				if (n.index < nn.index && n.style == Node::SurfNode && nn.style == Node::SurfNode)
				{
					if (CheckConnecting1(n.surf_points, nn.surf_points))
					{
						edges.emplace_back(edges.size(), n.index, nn.index);
						n.neighbors.emplace_back(nn.index);
						nn.neighbors.emplace_back(n.index);
					}
				}
			}
		}
		return edges;
	};

	//Note: there are challenged case I can not handle
	//This algorithm cannot guarantee to decompose all of surface into valid ends
	//multi merge
	auto RegionGrow = [&](std::vector<Node>& nodes, PCETransform::CutEnds& cut_ends_)
	{
		for (auto& n : nodes)
		{
			if (n.style == Node::EdgeNode)
			{
				nodes[n.neighbors[0]].PushEnd(n.index);
				n.used = true;
			}
		}

		while (true)
		{
			for (auto& n : nodes)
			{
				if (n.style == Node::SurfNode)
				{
					if (n.ends.size() > 2)
					{
						std::cerr << "Nodes' ends > 2..." << std::endl;
						system("pause");
						return false;
					}
					if (n.ends.size() == 2)
					{
						Node& n0 = nodes[n.ends[0]];
						Node& n1 = nodes[n.ends[1]];
						if (n0.friend_node < 0)n0.friend_node = n1.index;
						if (n1.friend_node < 0)n1.friend_node = n0.index;
						n0.PushEnd(n.index);
						n1.PushEnd(n.index);
						n.used = true;
					}
					if (n.ends.size() == 1)
					{
						Node& n_ = nodes[n.ends[0]];
						n_.PushEnd(n.index);
					}
				}
			}

			bool goon = false;
			for (auto& n : nodes)
			{
				if (n.style == Node::SurfNode && n.ends.size() == 0)
				{
					goon = true;
					break;
				}
				if (n.style == Node::EdgeNode && n.friend_node < 0)
				{
					goon = true;
					break;
				}
			}

			if (!goon) break;

			for (auto& n : nodes)
			{
				if (n.style == Node::SurfNode && !n.used && n.ends.size() > 0)
				{
					for (auto i : n.neighbors)
					{
						if (!nodes[i].used && nodes[i].ends.size() < 2)
							nodes[i].PushEnd(n.ends);
					}
					n.used = true;
				}
			}
		}

		//valid detection
		for (auto& n : nodes)
		{
			if (n.style == Node::EdgeNode && n.friend_node == -1)
			{
				std::cerr << "End computing error...: n.friend_node==-1" << std::endl;
				system("pause");
				return false;
			}
			if (n.style == Node::SurfNode && n.ends.size() == 0)
			{
				std::cerr << "End computing error...:  n.ends.size()==0" << std::endl;
				system("pause");
				return false;
			}
		}

		//build ends
		for (auto& n : nodes)
		{
			if (n.style == Node::EdgeNode && n.edge_style == Node::UpperEdge)
			{
				Vector3d2 con_surfs;
				std::vector<int> con_surf_indexes;
				Vector3d1 con_surf_normals;
				n.PushEnd(nodes[n.friend_node].ends);
				for (auto i : n.ends)
				{
					con_surfs.emplace_back(nodes[i].surf_points);
					con_surf_normals.emplace_back(nodes[i].surf_normal);
					con_surf_indexes.emplace_back(nodes[i].surface_index);
				}

				if (con_surfs.size() != 1)
				{
					std::cerr << "Two ends cutting...." << std::endl;
					return false;
				}
				cut_ends_.ends.emplace_back(PCEEnd(kerf, n.edge_index, n.edge_s, n.edge_e,
					nodes[n.friend_node].edge_s, nodes[n.friend_node].edge_e, con_surfs, con_surf_normals, con_surf_indexes));
			}
		}
	};

	//Note: this function can only handle very simple case
	auto ConnectContours = [](const HMODULE& hModule_, const std::vector<std::pair<Vector2d, Vector2d>>& contours)
	{
		Vector2d1 results;
		auto line_itersection = (CGAL_2D_Intersection_Line_Line)GetProcAddress(hModule_, "CGAL_2D_Intersection_Line_Line");
		for (int i = 0; i < contours.size(); i++)
		{
			int index_0 = i;
			int index_1 = (i + 1) % contours.size();
			Vector2d inter;
			line_itersection(contours[index_0].first, contours[index_0].second, contours[index_1].first, contours[index_1].second, inter);
			results.emplace_back(inter);
		}
		return results;
	};

	PCETransform::CutEnds cut_ends_;

	cut_ends_.valid = false;
	Vector3d1 surf_normals_;
	for (const auto& surf_ : surfs_)
	{
		auto n = Math::Functs::ComputeNormalFromPolyline(surf_);
		Math::Functs::SetVectorLength(n, 1.0);
		surf_normals_.emplace_back(-n[0], -n[1], -n[2]);
	}

	//input valid detection
	for (auto& s : surfs_)
		if (s.size() < 3)
		{
			std::cerr << "Input surface is not valid..." << std::endl;
			system("pause");
			return cut_ends_;
		}

	if (surfs_.size() != surf_normals_.size())
	{
		std::cerr << "Input surface and normals are not valid..." << std::endl;
		system("pause");
		return cut_ends_;
	}

	//build nodes and edges
	std::vector<Node> nodes;
	if (!BuildNodes(nodes, surf_normals_, cut_ends_)) return cut_ends_;
	std::vector<Edge> edges = BuildEdges(nodes);

	//valid detection
	for (auto& n : nodes)
	{
		if (n.style == Node::EdgeNode && n.neighbors.size() != 1)
		{
			std::cerr << "Edge computing error..." << std::endl;
			system("pause");
			return cut_ends_;
		}
		if (n.style == Node::SurfNode && n.neighbors.size() != n.surf_points.size())
		{
			std::cerr << "Edge computing error..." << std::endl;
			system("pause");
			return cut_ends_;
		}
	}

	//Note: there are challenged case I can not handle
	//This algorithm cannot guarantee to decompose all of surface into valid ends
	//multi merge
	RegionGrow(nodes, cut_ends_);

	//sort ends
	auto SortEnd = [](PCEEnd& end0, PCEEnd& end1) { return  end0.index < end1.index; };
	std::sort(cut_ends_.ends.begin(), cut_ends_.ends.end(), SortEnd);

	//kerf offset
	auto check_orientation = (CGAL_2D_Polygon_Is_Clockwise_Oriented)GetProcAddress(hModule, "CGAL_2D_Polygon_Is_Clockwise_Oriented");
	auto PolyOffset = (CGAL_2D_Polygon_One_Offsets)GetProcAddress(hModule, "CGAL_2D_Polygon_One_Offsets");
	auto PolyUnion = (CGAL_2D_Two_Polygons_Union)GetProcAddress(hModule, "CGAL_2D_Two_Polygons_Union");

	//u_origin_2d  l_origin_2d
	for (auto p : cut_ends_.u_surf) cut_ends_.u_origin_2d.emplace_back(p.x, p.y);
	for (auto p : cut_ends_.l_surf) cut_ends_.l_origin_2d.emplace_back(p.x, p.y);

	//ul_origin_2d
	Vector2d2 polygons;
	polygons.clear();
	PolyUnion(cut_ends_.u_origin_2d, cut_ends_.l_origin_2d, polygons);
	cut_ends_.ul_origin_2d = polygons[0];

	//u_offset_2d
	std::vector<std::pair<Vector2d, Vector2d>> contours;
	for (auto& end : cut_ends_.ends)
	{
		contours.emplace_back(Vector2d(end.cut.u_cutting_line_s[0], end.cut.u_cutting_line_s[1]),
			Vector2d(end.cut.u_cutting_line_e[0], end.cut.u_cutting_line_e[1]));
	}
	cut_ends_.u_offset_2d = ConnectContours(hModule, contours);

	//l_offset_2d
	contours.clear();
	for (auto& end : cut_ends_.ends)
	{
		contours.emplace_back(Vector2d(end.cut.l_cutting_line_s[0], end.cut.l_cutting_line_s[1]),
			Vector2d(end.cut.l_cutting_line_e[0], end.cut.l_cutting_line_e[1]));
	}
	cut_ends_.l_offset_2d = ConnectContours(hModule, contours);

	//valid detection
	if (cut_ends_.u_offset_2d.size() != cut_ends_.l_offset_2d.size())
	{
		std::cerr << "if (u_offset_2d.size() != l_offset_2d.size())" << std::endl;
		system("pause");
		return cut_ends_;
	}

	//ul_offset_2d
	polygons.clear();
	PolyUnion(cut_ends_.u_offset_2d, cut_ends_.l_offset_2d, polygons);
	cut_ends_.ul_offset_2d = polygons[0];

	cut_ends_.ul_origin_2d=Math::Functs::Polygon_Clear(cut_ends_.ul_origin_2d, CompilerConfig::Instance().GetAngleMatchError(), CompilerConfig::Instance().GetPartMatchError());
	cut_ends_.ul_offset_2d=Math::Functs::Polygon_Clear(cut_ends_.ul_offset_2d, CompilerConfig::Instance().GetAngleMatchError(), CompilerConfig::Instance().GetPartMatchError());

	cut_ends_.valid = true;
	return cut_ends_;
}


bool PCETransform::BuildEnd(const double& kerf, const HMODULE& hModule, Vector3d upper_normal, Vector3d lower_normal)
{
	cut_ends = BuildEnd(kerf, hModule, surfs, upper_normal, lower_normal);
	return cut_ends.valid;
}


bool PCETransform::BuildBoundingBox(const double& kerf, const HMODULE& hModule)
{
	auto GetBoundingBox_2D = [](const Vector2d2& points, double kerf_ = 0.0)
	{
		Vector2d corner_0, corner_2;
		Math::Functs::GetBoundingBox(points, corner_0, corner_2);
		Vector2d1 box;
		box.emplace_back(Vector2d(corner_0[0] - kerf_ / 2.0, corner_0[1] - kerf_ / 2.0));
		box.emplace_back(Vector2d(corner_2[0] + kerf_ / 2.0, corner_0[1] - kerf_ / 2.0));
		box.emplace_back(Vector2d(corner_2[0] + kerf_ / 2.0, corner_2[1] + kerf_ / 2.0));
		box.emplace_back(Vector2d(corner_0[0] - kerf_ / 2.0, corner_2[1] + kerf_ / 2.0));
		return box;
	};

	Vector2d2 surfs_2d = Math::Functs::Vector3d2d(surfs);

	bounding_box.origin_2d = GetBoundingBox_2D(surfs_2d);
	bounding_box.offset_2d = GetBoundingBox_2D(surfs_2d, kerf);
	bounding_box.center_2d = (bounding_box.origin_2d[0] + bounding_box.origin_2d[2]) / 2.0;

	bounding_box.normals.emplace_back(Vector3d(0.0, -1.0, 0.0));
	bounding_box.normals.emplace_back(Vector3d(1.0, 0.0, 0.0));
	bounding_box.normals.emplace_back(Vector3d(0.0, 1.0, 0.0));
	bounding_box.normals.emplace_back(Vector3d(-1.0, 0.0, 0.0));

	return true;
}

