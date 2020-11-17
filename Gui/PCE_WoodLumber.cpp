#include "PCE.h"
#include "GeomCommonFunctions.h"
#include <Gui/Application.h>
#include <BRepBuilderAPI_GTransform.hxx>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tchar.h>
#include <iostream>

using namespace PartDesignGui;

//WL
////////////////////////////////////////////////////////////////////////////////////////////
WLCollection::WL::WL(int index_, string style_, double x_, double y_, double z_, string u_, double kerf) :
	index(index_), style(style_), volume(x_* y_* z_)
{
	std::vector<double> xyz{ x_,y_,z_ };
	std::sort(xyz.begin(), xyz.end());

	size = Vector3d(xyz[2], xyz[1], xyz[0]);

	const boost::uuids::string_generator gen;
	uid = gen(u_);

	offset_points.emplace_back(-kerf * 0.5, -kerf * 0.5);
	offset_points.emplace_back(size[0] + kerf * 0.5, -kerf * 0.5);
	offset_points.emplace_back(size[0] + kerf * 0.5, size[1] + kerf * 0.5);
	offset_points.emplace_back(-kerf * 0.5, size[1] + kerf * 0.5);

	origin_points.emplace_back(0.0, 0.0);
	origin_points.emplace_back(size[0], 0.0);
	origin_points.emplace_back(size[0], size[1]);
	origin_points.emplace_back(0.0, size[1]);

	surf_normals.emplace_back(Vector3d(0.0, -1.0, 0.0));
	surf_normals.emplace_back(Vector3d(1.0, 0.0, 0.0));
	surf_normals.emplace_back(Vector3d(0.0, 1.0, 0.0));
	surf_normals.emplace_back(Vector3d(-1.0, 0.0, 0.0));

	Vector3d1 bound_wood{
	Vector3d(-kerf,-kerf,0.0),
	Vector3d(size[0] + kerf,-kerf,0.0),
	Vector3d(size[0] + kerf,size[1] + kerf,0.0),
	Vector3d(-kerf,size[1] + kerf,0.0)
	};

	boundary_3d.emplace_back(bound_wood);

	Vector3d1 bound_lumber{
	Vector3d(-kerf,0.0,-kerf),
	Vector3d(size[0] + kerf,0.0,-kerf),
	Vector3d(size[0] + kerf,0.0,size[2] + kerf),
	Vector3d(-kerf,0.0,size[2] + kerf)
	};
	boundary_3d.emplace_back(bound_lumber);
}

std::vector<PCEEdge> WLCollection::WL::ExtractEdges()
{
	std::vector<PCEEdge> edges_;

	edges.clear();
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
		cut.cutting_surface_normal = surf_normals[i];
		cut.u_cutting_line_s = Vector3d(offset_s[0], offset_s[1], size[2]);
		cut.u_cutting_line_e = Vector3d(offset_e[0], offset_e[1], size[2]);
		cut.l_cutting_line_s = Vector3d(offset_s[0], offset_s[1], 0.0);
		cut.l_cutting_line_e = Vector3d(offset_e[0], offset_e[1], 0.0);

		cut.l_u_vector = Vector3d(0.0, 0.0, size[2]);

		edges_.emplace_back(PCEEdge(edges_.size(), -1, cut));
		edges.emplace_back(edges_.back().index);
	}
	return edges_;
}

////////////////////////////////////////////////////////////////////////////////////////////
WLCollection::WLCollection(int index_) : index(index_)
{
}

void WLCollection::SortWLS()
{
	const auto Comp = [](WL & wl0, WL & wl1)
	{
		return wl0.volume > wl1.volume;
	};
	std::sort(wls.begin(), wls.end(), Comp);
	for (auto i = 0; i < wls.size(); i++) wls[i].index = i;
}

void WLCollection::Clear()
{
	std::vector<WL>().swap(wls);
	//std::vector<int>().swap(shapes);
	//std::vector<int>().swap(input_shapes);
}

std::vector<PCETransform> WLCollection::DetectSize
(const HMODULE& hModule,
	const double& kerf,
	const double part_match_error,
	const shared_ptr<PCEShape>& ptr_shape,
	const bool WLCFree)
{
	const auto CompTransform = [](const PCETransform & transform_0, const PCETransform & transform_1)
	{
		return transform_0.wl_rates.front().volume_rate > transform_1.wl_rates.front().volume_rate;
	};

	const auto CompWlrate = [](const PCETransform::WLRATE & wl_rate_0, const PCETransform::WLRATE & wl_rate_1)
	{
		if (wl_rate_0.perfect_matching > wl_rate_1.perfect_matching)
		{
			return true;
		}
		else
		{
			if (wl_rate_0.perfect_matching == wl_rate_1.perfect_matching)
				return wl_rate_0.volume_rate > wl_rate_1.volume_rate;
			else
				return false;
		}
	};

	auto GetRP = [](const glm::dmat4 &PM, const Vector3d2& in_surfs)
	{
		Vector3d2 r_ps;
		for (const auto surf : in_surfs)
		{
			Vector3d1 ps;
			for (const auto p : surf)
				ps.emplace_back(Vector3d(PM * glm::vec4(p,1.0)));
			r_ps.emplace_back(ps);
		}
		return r_ps;
	};

	auto GetRN = [](const glm::dmat4& RM, const Vector3d1 & normals)
	{
		Vector3d1 r_normals;
		for (const auto normal : normals)
			r_normals.emplace_back(Vector3d(RM * glm::vec4(normal, 1.0)));
		return r_normals;
	};

	auto GetDistance = [](const HMODULE & hModule_, const Vector3d2 & ps_0, const Vector3d2 & ps_1)
	{
		auto cgal_3d_distance_point_polygon = (CGAL_3D_Distance_Point_Polygon)GetProcAddress(hModule_, "CGAL_3D_Distance_Point_Polygon");

		Vector3d2 ps_0_ = ps_0;
		Vector3d2 ps_1_ = ps_1;

		double total_distance = 0.0;
		int nb = 0;
		for (auto& ps : ps_0_)
			for (auto& p : ps)
			{
				double min_dis = 100000000000.0;
				for (auto& poly : ps_1_)
					min_dis = min(min_dis, cgal_3d_distance_point_polygon(poly, p));
				total_distance += min_dis;
				nb++;
			}
		total_distance = total_distance / (double)nb;

		return total_distance;
	};

	auto OutputRectangle = [](std::string path, const PCETransform & transform, const HMODULE & hModule, const double& kerf)
	{
		auto output = (CGAL_Export_Path_Segment)GetProcAddress(hModule, "CGAL_Export_Path_Segment");
		const auto& points = transform.surfs;

		std::ofstream file(path);
		for (int i = 0; i < points.size(); i++)
			for (int j = 0; j < points[i].size(); j++)
				file << "v " << points[i][j][0] << " " << points[i][j][1] << " " << points[i][j][2] << std::endl;

		int nb = 1;
		for (int i = 0; i < points.size(); i++)
		{
			file << "f ";
			for (int j = 0; j < points[i].size(); j++)
			{
				file << Math::Functs::IntString(nb) << " ";
				nb++;
			}
			file << "" << std::endl;
		}
		file.clear();
		file.close();
	};

	

	std::vector<PCETransform> transforms;
	if (part_match_error < 0.0) return transforms;
	Vector3d1 rotations= Math::Functs::EmumerateRotations();

	//angle 0.0
	for (const auto& rotation : rotations)
	{
		glm::dmat4 PM = 
			Math::Functs::TranslationMatrix(ptr_shape->center) *
			Math::Functs::RotationMatrix(Vector3d(0, 0, 1), -rotation[2]) *
			Math::Functs::RotationMatrix(Vector3d(0, 1, 0), -rotation[1]) *
			Math::Functs::RotationMatrix(Vector3d(1, 0, 0), -rotation[0]) *
			Math::Functs::TranslationMatrix(-ptr_shape->center);

		glm::dmat4 RM =
			Math::Functs::RotationMatrix(Vector3d(0, 0, 1), -rotation[2]) *
			Math::Functs::RotationMatrix(Vector3d(0, 1, 0), -rotation[1]) *
			Math::Functs::RotationMatrix(Vector3d(1, 0, 0), -rotation[0]);

		auto surfs = GetRP(PM, ptr_shape->surfs);
		auto surf_normals = GetRN(RM, ptr_shape->surf_normals);

		Vector3d maximal_corner, minimal_corner;
		Math::Functs::GetBoundingBox(surfs, minimal_corner, maximal_corner);
		Vector3d shape_size = maximal_corner - minimal_corner;

		for (auto& surf : surfs)
			for (auto& p : surf)
				p = p - minimal_corner;

		PM = Math::Functs::TranslationMatrix(-minimal_corner) * PM;

		std::vector<PCETransform::WLRATE> wl_rates;

		for (auto& wl : wls)
		{
			double ratio = -1.0;
			int perfect_matching = 0;


			if (WLCFree)
			{
				ratio = shape_size[0] * shape_size[1] * shape_size[2] / wl.volume;

				if (Math::Functs::IsAlmostZero_Double(wl.size[0] - shape_size[0], part_match_error)) perfect_matching++;
				if (Math::Functs::IsAlmostZero_Double(wl.size[1] - shape_size[1], part_match_error)) perfect_matching++;
				if (Math::Functs::IsAlmostZero_Double(wl.size[2] - shape_size[2], part_match_error)) perfect_matching++;
			}
			else
			{
				if (wl.size[0] + part_match_error >= shape_size[0] && wl.size[1] + part_match_error >= shape_size[1] && wl.size[2] + part_match_error >= shape_size[2])
				{
					ratio = shape_size[0] * shape_size[1] * shape_size[2] / wl.volume;

					if (Math::Functs::IsAlmostZero_Double(wl.size[0] - shape_size[0], part_match_error)) perfect_matching++;
					if (Math::Functs::IsAlmostZero_Double(wl.size[1] - shape_size[1], part_match_error)) perfect_matching++;
					if (Math::Functs::IsAlmostZero_Double(wl.size[2] - shape_size[2], part_match_error)) perfect_matching++;
					else
						ratio = -1.0;
				}
			}

			if (ratio >= 0)
			{
				bool b = false;
				
				for (auto& trf : transforms)
				{
					if (GetDistance(hModule, surfs, trf.surfs) < part_match_error)
					{
						b = true;
						if (!trf.CheckWL(wl.index))
							trf.wl_rates.emplace_back(PCETransform::WLRATE(wl,ratio, perfect_matching));
						break;
					}
				}

				if (!b)
				{
					transforms.emplace_back(PCETransform(surfs, surf_normals, ptr_shape->center, rotation, PM, RM, minimal_corner));
					transforms.back().wl_rates.emplace_back(PCETransform::WLRATE(wl,ratio, perfect_matching));
				}

			}
		}
	}

	for (auto& transform : transforms) std::sort(transform.wl_rates.begin(), transform.wl_rates.end(), CompWlrate);

	std::sort(transforms.begin(), transforms.end(), CompTransform);

	for (int iter = 0; iter < transforms.size(); iter++)
		transforms[iter].index = iter;

	std::vector<PCETransform> return_transforms;
	for (auto& transform : transforms)
	{
		transform.BuildBoundingBox(kerf, hModule);
		if (transform.BuildEnd(kerf, hModule))
			return_transforms.emplace_back(transform);
	}

	for (int i = 0; i < return_transforms.size(); i++)
	{
		return_transforms[i].index = i;
		return_transforms[i].pce_shape = ptr_shape;
	}
		
	return return_transforms;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//EGraphClass
//////////////////////////////////////////////////////////////////////////////////////////////////////////

EGraphClass::EGraphClass(int index_,const std::vector<int>& ec_shapes_) :
	index(index_), str(Math::Functs::IntString(ec_shapes_,true,",")), ec_shapes(ec_shapes_), head_index(0)
{
	if (ec_shapes_.empty())
	{
		std::cerr << "if (ec_shapes_.empty())" << std::endl;
		system("pause");
	}

	valid = false;
	compiledNodesB=false,
	ShapeEncode();
}

void EGraphClass::ShapeEncode()
{
	if (ec_shapes.empty())
	{
		std::cerr << "if (ec_shapes.empty())" << std::endl;
		system("pause");
	}
	encode_shapes=PCEEncode::GetEncode(ec_shapes);
};

void EGraphClass::PushComp(std::vector<int> comp_)
{
	if (!(comp_.size() == 1 && comp_[0] == index))
	{
		string str_ = Math::Functs::IntString(comp_,true,",");

		for (auto edges : egraph_edges)
		{
			string str = Math::Functs::IntString(edges,true,",");
			if (str_ == str)return;
		}
		egraph_edges.emplace_back(comp_);
	}
}

void EGraphClass::Pruning(std::vector<shared_ptr<PCEArrange>>& egraph_arranges,
	const int& maximal_number, const bool& random_pruning)
{
	if (arranges.empty())return;

	std::vector<std::tuple<int, int, double>> cns;
	for (int i = 0; i < arranges.size(); i++)
	{
		const auto& arrange = egraph_arranges[arranges[i]];
		cns.emplace_back(arranges[i], arrange->cutting_number, arrange->volume_rate);
	}

	if (random_pruning)
		std::random_shuffle(cns.begin(), cns.end());
	else
	{
		std::sort(cns.begin(), cns.end(), [](const std::tuple<int, int, double>& cn0, const std::tuple<int, int, double>& cn1)
			{return std::get<1>(cn0) < std::get<1>(cn1) || (std::get<1>(cn0) == std::get<1>(cn1) && std::get<2>(cn0) > std::get<2>(cn1));});
	}

	std::vector<int>().swap(arranges);

	for (int i = 0; i < maximal_number && i < cns.size() && std::get<1>(cns[i]) == std::get<1>(cns.front()); i++)
	{
		int arrange_index = std::get<0>(cns[i]);

		const auto& arrange = egraph_arranges[arrange_index];
		bool b = true;
		for (auto& index_ : arranges)
		{
			const auto& arrange_ = egraph_arranges[index_];
			if (arrange->encode_packing_0 == arrange_->encode_packing_0 ||
				arrange->encode_packing_0 == arrange_->encode_packing_1 ||
				arrange->encode_packing_1 == arrange_->encode_packing_0 ||
				arrange->encode_packing_1 == arrange_->encode_packing_1)
			{
				b = false;
				break;
			}
		}
		if (b) arranges.emplace_back(arrange_index);
	}
}

void EGraphClass::Pruning(std::vector<PCEArrange> & egraph_arranges, 
	const int& maximal_number, const bool& random_pruning)
{
	struct CN
	{
	public:
		int index;
		int cutting_number;
		double volume_rate;
		CN()
		{
			index = -1;
			cutting_number = -1;
		}
		CN(int index_, int cutting_number_, double volume_rate_)
		{
			index = index_;
			cutting_number = cutting_number_;
			volume_rate = volume_rate_;
		}
	};

	auto Comp = [](const CN & cn0, const CN & cn1)
	{
		if (cn0.cutting_number < cn1.cutting_number)
			return true;
		else
			if (cn0.cutting_number == cn1.cutting_number && cn0.volume_rate > cn1.volume_rate)
				return true;
		return false;
	};

	std::vector<CN> cns;
	for (int i = 0; i < arranges.size(); i++)
	{
		const auto& arrange = egraph_arranges[arranges[i]];
		cns.emplace_back(CN(arranges[i], arrange.cutting_number, arrange.volume_rate));
	}

	if (random_pruning)
		std::random_shuffle(cns.begin(), cns.end());
	else
		sort(cns.begin(), cns.end(), Comp);

	std::vector<int>().swap(arranges);

	for (int i = 0; i < maximal_number && i < cns.size(); i++)
	{
		const auto& arrange = egraph_arranges[cns[i].index];
		bool b = true;
		for (auto& index_ : arranges)
		{
			const auto& arrange_ = egraph_arranges[index_];
			if (arrange.encode_packing_0 == arrange_.encode_packing_0||
				arrange.encode_packing_0 == arrange_.encode_packing_1||
				arrange.encode_packing_1 == arrange_.encode_packing_0 || 
				arrange.encode_packing_1 == arrange_.encode_packing_1)
			{
				b = false;
				break;
			}
		}
		if (b) arranges.emplace_back(cns[i].index);
	}
}

void EGraphClass::Clear()
{
	arranges.clear();
	egraph_edges.clear();
}
