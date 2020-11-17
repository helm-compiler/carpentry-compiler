#ifndef _PCE_COMBINE_H
#define _PCE_COMBINE_H
#include "PreCompiled.h"

#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <GeomLProp_SLProps.hxx>
#include <gp_Dir.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepExtrema_ExtFF.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TCollection_AsciiString.hxx>
#include <Geom_CartesianPoint.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepTools_WireExplorer.hxx>
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <ShapeFix_Wire.hxx>
#include <vector>
#include <map>
#include <unordered_map>
#include <boost/uuid/uuid.hpp>
#include "Math.hpp"
#include "CGAL.h"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/string_generator.hpp>
#include "slvs.h"
#include "tinyxml2.h"

namespace PartDesignGui {

	class PCE;
	struct PCEEdge;
	struct WLCollection;
	struct PCEArrange;
	struct PCEPosition;
	struct PCETransform;
	struct PCEEnd;
	struct PCEShape;
	struct ArrNode;

	//std::vector<TopoDS_Shape>& g_shapes, std::vector<TopoDS_Shape>& l_shapes, std::vector<glm::dmat4>& LGMs
	typedef std::tuple<TopoDS_Shape, TopoDS_Shape, glm::dmat4> GLMAT;
	typedef std::vector<std::tuple<TopoDS_Shape, TopoDS_Shape, glm::dmat4>> GLMATS;

	static void Output3DSegments_Pairs(
		const HMODULE& hModule,
		const double kerf,
		const std::string& path, 
		const std::vector<std::pair<Vector3d, Vector3d>>& poly)
	{
		auto output = (CGAL_Export_Path_Segment)GetProcAddress(hModule, "CGAL_Export_Path_Segment");
		std::ofstream export_fie(path);
		int export_int = 1;
		for (int i = 0; i < poly.size(); i++)
			output(export_fie,
				export_int,
				"edge_" + Math::Functs::IntString(i),
				0.0, 0.0, 0.0,
				poly[i].first, poly[i].second, kerf * 0.5);

		export_fie.clear();
		export_fie.close();
	};

	static void Output2DSegments(
		const HMODULE& hModule,
		const double kerf,
		const std::string& path, const Vector2d1& poly)
	{
		auto output = (CGAL_Export_Path_Segment)GetProcAddress(hModule, "CGAL_Export_Path_Segment");
		std::ofstream export_fie(path);
		int export_int = 1;
		for (int i = 0; i < poly.size(); i++)
			output(export_fie,
				export_int,
				"edge_" + Math::Functs::IntString(i),
				0.0, 0.0, 0.0,
				Math::Functs::Vector2d3d(poly[i]), Math::Functs::Vector2d3d(poly[(i + 1) % poly.size()]), kerf * 0.5);

		export_fie.clear();
		export_fie.close();
	};

	struct NumberLabel
	{
		enum {BaseCircle,BaseShere, BaseSquare1, BaseSquare2, FullSetup, HalfAngle, HalfOffset};

		Vector3d1 base_circle_vecs;
		Vector1i2 base_circle_faces;

		Vector3d1 sphere_vecs;
		Vector1i2 sphere_faces;

		Vector3d1 base_square_1_vecs;
		Vector1i2 base_square_1_faces;

		Vector3d1 base_square_2_vecs;
		Vector1i2 base_square_2_faces;

		Vector3d2 numbers_vecs;
		Vector1i3 numbers_faces;

		Vector3d1 full_setup_vecs;
		Vector1i2 full_setup_faces;

		Vector3d1 half_angle_setup_vecs;
		Vector1i2 half_angle_setup_faces;

		Vector3d1 half_offset_setup_vecs;
		Vector1i2 half_offset_setup_faces;

		bool load = false;
		void Clear()
		{
			Vector3d1().swap(base_circle_vecs);
			Vector1i2().swap(base_circle_faces);
			Vector3d1().swap(sphere_vecs);
			Vector1i2().swap(sphere_faces);
			Vector3d1().swap(base_square_1_vecs);
			Vector1i2().swap(base_square_1_faces);
			Vector3d1().swap(base_square_2_vecs);
			Vector1i2().swap(base_square_2_faces);
			Vector3d1().swap(full_setup_vecs);
			Vector1i2().swap(full_setup_faces);
			Vector3d1().swap(half_angle_setup_vecs);
			Vector1i2().swap(half_angle_setup_faces);
			Vector3d1().swap(half_offset_setup_vecs);
			Vector1i2().swap(half_offset_setup_faces);
			Vector3d2().swap(numbers_vecs);
			Vector1i3().swap(numbers_faces);
		}

		NumberLabel()
		{
			load = false;
		}

		void Init();

		void DrawNumber1(std::ofstream& export_fie, int& export_int, const std::string& name, const double number, const double radius, const Vector3d& label_center);
		void DrawNumber(std::ofstream& export_fie, int& export_int, const std::string& name, const double number, const double radius, const Vector3d& label_center);
		void DrawLabel(const int base_style, std::ofstream& export_fie, int& export_int, const std::string& name, const int number, const double radius, const Vector3d& label_center);
		void DrawNode(const int style, std::ofstream& export_fie, int& export_int, const std::string& name, const double radius, const Vector3d& label_center);
	};

	struct PCEHole
	{
	public:
		PCEHole(const std::vector<Vector3d>& _rp, const double _r, const int _p) : refPoints(_rp), radius(_r), partId(_p) {}
		std::vector<Vector3d> refPoints;
		double radius;
		int partId;
	};

	struct LocalCS
	{
	public:
		Vector2d origin;
		Vector2d axis_x;
		Vector2d axis_y;

		double origin_max_y;
		double origin_min_y;
		double origin_min_x;
		double origin_max_x;

		double offset_max_y;
		double offset_min_y;
		double offset_min_x;
		double offset_max_x;

		double angle;
		LocalCS();
		LocalCS(const Vector2d& o, const Vector2d& x);
		void SetCS(const Vector2d& o, const Vector2d& x);
		Vector2d GetLocal(const Vector2d& v) const;
		Vector2d1 GetLocal(const Vector2d1& vecs) const;

		//need to check
		Vector2d GetGlobalPos(const Vector2d& v) const;
		Vector2d GetGlobalVec(const Vector2d& v) const;

		bool GetSegOrigin(
			const HMODULE& hModule,
			const double part_match_error,
			const Vector2d& s, 
			const Vector2d& e,
			double& min_x, 
			double& max_x) const;

		bool GetSegOffset(
			const HMODULE& hModule, 
			const double part_match_error,
			const Vector2d& s, 
			const Vector2d& e, 
			double& min_x, 
			double& max_x) const;

		Vector2d1 SortX(const Vector2d1& segs) const;
		Vector2d1 EmptyX(const Vector2d1& sort_segs) const;
	};

	struct PCEEnd
	{
		enum { VerticalEnd, TiltingEnd };

		struct PCECUT
		{
			Vector2d origin_s;
			Vector2d origin_e;
			Vector2d origin_c; // center point
			Vector2d origin_n; // normal

			Vector2d offset_s;
			Vector2d offset_e;
			Vector2d offset_c; // center point
			Vector2d offset_n; // normal
			Vector3d cutting_surface_normal;

			Vector3d u_cutting_line_s;
			Vector3d u_cutting_line_e;
			Vector3d l_cutting_line_s;
			Vector3d l_cutting_line_e;

			Vector3d l_u_vector;

			void Translate(const Vector2d& t);
		};

		PCEEnd(const double& kerf,
			const int& index_, 
			const Vector3d& u_s_, 
			const Vector3d& u_e_, 
			const Vector3d& d_s_, 
			const Vector3d& d_e_,
			const Vector3d2& con_surfs_,
			const Vector3d1& con_surf_normals_,
			const std::vector<int>& con_surf_indexes_);

		void GetEdge(const double& kerf);
		void Clear();
		
		int index;
		Vector3d u_s, u_e, l_s, l_e;
		Vector3d2 con_surfs;
		Vector3d1 con_surf_normals;
		std::vector<int> con_surf_indexes;
		int style;
		PCECUT cut;
	};

	typedef PCEEnd::PCECUT PCECUT;

	struct PCEConstraint
	{
		int edge_0;
		int edge_1;
		PCEConstraint(const int& edge_0_, const int& edge_1_) :edge_0(edge_0_), edge_1(edge_1_) {};
	};

	struct PCEEdge
	{
	public:
		enum { VerticalEnd, TiltingEnd};

		int index;
		int part_index;

		Vector3d cutting_line_s;
		Vector3d cutting_line_e;

		PCECUT cut;
		int style=-1;

		PCEEdge(const Vector3d& l_cutting_line_s_,
			const Vector3d& l_cutting_line_e_, 
			const Vector3d& cutting_surf_normal_, 
			const Vector3d& l_u_vector_);

		PCEEdge(const int &index_, const int &part_index_, const PCECUT &cut_);

		void Set(const PCECUT& cut_);
	};

	struct EGraphClass
	{
	public:
		int index;
		std::vector<int> ec_shapes;
		std::vector<int> arranges;
		std::vector<std::vector<int>> egraph_edges;
		std::map<int, int> encode_shapes;
		int head_index;

		int pre_eclass=-1;
		void ShapeEncode();

		string str;
		bool valid;

		bool compiledNodesB;
		std::string xmlDoc_path;

		EGraphClass() :index(0), str(""), valid(false), compiledNodesB(false), head_index(0) {}
		EGraphClass(int index_, const std::vector<int>& ec_shapes_);
		void PushComp(std::vector<int> comp_);
		void Pruning(std::vector<PCEArrange>& egraph_arranges, const int& maximal_number, const bool& random_pruning);
		void Pruning(std::vector<shared_ptr<PCEArrange>>& egraph_arranges, const int& maximal_number, const bool& random_pruning);
		void Clear();
	};

	struct WLCollection
	{
	public:
		struct WL
		{
		public:
			WL() :index(-1) {};
			WL(int index_, string style_, double x_, double y_, double z_, string u_, double kerf);
			std::vector<PCEEdge> ExtractEdges();
			int index;
			Vector3d size;
			string style;
			boost::uuids::uuid uid;
			double volume=0.0;
			std::vector<int> edges;
			Vector2d1 offset_points;
			Vector2d1 origin_points;
			Vector3d1 surf_normals;
			Vector3d2 boundary_3d;
			int all_enumerating_nb = 0;
		};

		

		WLCollection(int index_);
		void Clear();

		std::vector<PCETransform> DetectSize
			(const HMODULE& hModule,
			const double& kerf,
			const double part_match_error,
			const shared_ptr<PCEShape>& ptr_shape,
			const bool WLCFree=false);

		void SortWLS();

		int index;
		std::vector<WL> wls;
		std::vector<int> shapes;
		std::vector<int> input_shapes;//only include the shape type
	};

	typedef WLCollection::WL WL;
	//typedef WLCollection::EGraphClass EGraphClass;

	struct PCETransform
	{
		struct WLRATE
		{
			WL* wl;
			double volume_rate;
			int perfect_matching;
			WLRATE(WL& wl_, double volume_rate_, int perfect_matching_) :
				wl(&wl_), volume_rate(volume_rate_), perfect_matching(perfect_matching_) {};
		};
		
		struct Boundingbox
		{
			Vector2d center_2d;
			Vector2d1 offset_2d;
			Vector2d1 origin_2d;
			Vector3d1 normals;
		};

		struct CutEnds
		{
			bool valid=false;
			std::vector<PCEEnd> ends;
			Vector3d1 u_surf;
			Vector3d1 l_surf;
			Vector2d1 u_offset_2d;
			Vector2d1 u_origin_2d;
			Vector2d1 l_offset_2d;
			Vector2d1 l_origin_2d;
			Vector2d1 ul_offset_2d;
			Vector2d1 ul_origin_2d;
			int u_surf_index;
			int l_surf_index;
			void Clear()
			{
				for (auto& e : ends) e.Clear();
				std::vector<PCEEnd>().swap(ends);
				Vector3d1().swap(u_surf);
				Vector3d1().swap(l_surf);
				Vector2d1().swap(ul_offset_2d);
				Vector2d1().swap(ul_origin_2d);
				Vector2d1().swap(u_offset_2d);
				Vector2d1().swap(u_origin_2d);
				Vector2d1().swap(l_offset_2d);
				Vector2d1().swap(l_origin_2d);
			}

			Vector3d2 Mesh()
			{
				Vector3d2 mesh;
				mesh.emplace_back(u_surf);

				for (auto& end : ends)
					mesh.emplace_back(Vector3d1{ end.u_e, end.u_s, end.l_e, end.l_s });

				mesh.emplace_back(l_surf);
				return mesh;
			};

			Vector3d1 Normals()
			{
				Vector3d1 normals;
				normals.emplace_back(Vector3d(0.0,0.0,1.0));
				for (auto& end : ends)
					normals.emplace_back(end.con_surf_normals[0]);
				normals.emplace_back(Vector3d(0.0, 0.0, -1.0));
				return normals;
			};
		};

		int index;
		
		//transform
		Vector3d2 surfs;
		Vector3d1 surf_normals;

		Vector3d center;
		Vector3d rotation;
		Vector3d tranform_corner;
		glm::dmat4 PM;
		glm::dmat4 RM;

		//int shape_index;

		shared_ptr<PCEShape> pce_shape;

		std::vector<WLRATE> wl_rates;

		//process
		CutEnds cut_ends;
		Boundingbox bounding_box;

		PCETransform(Vector3d2 surfs_, 
			Vector3d1 surf_normals_,
			Vector3d center_, 
			Vector3d rotation_,
			glm::dmat4 PM_,
			glm::dmat4 RM_,
			Vector3d tranform_corner_);

		PCETransform();
		bool CheckWL(const int& wl_index_) const;
		void Clear();
		void AddNewEdges(std::vector<PCEEdge>& edges_, std::vector<int>& edge_indexes);
		void AddBoundingEdges(std::vector<PCEEdge>& edges_, std::vector<int>& edge_indexes);

		Vector2d1 UOffset(const Vector2d t) const;
		Vector2d1 LOffset(const Vector2d t) const;
		Vector2d1 UOrigin(const Vector2d t) const;
		Vector2d1 LOrigin(const Vector2d t) const;
		Vector2d1 ULOffset(const Vector2d t) const;
		Vector2d1 ULOrigin(const Vector2d t) const;

		bool BuildEnd(const double& kerf, const HMODULE& hModule, 
			Vector3d upper_normal=Vector3d(0.0, 0.0, 1.0), Vector3d lower_normal = Vector3d(0.0, 0.0, -1.0));

		static PCETransform::CutEnds BuildEnd(const double& kerf, const HMODULE& hModule,
			const Vector3d2& surfs_, 
			Vector3d upper_normal = Vector3d(0.0, 0.0, 1.0), Vector3d lower_normal = Vector3d(0.0, 0.0, -1.0));

		bool BuildBoundingBox(const double& kerf, const HMODULE& hModule);
	};


	struct PCEEncode
	{
		static int Insert(std::map<int, int>& encode, const int& shape);

		static std::map<int, int> GetEncode(const std::vector<int> &shapes);

		static int GetEncodeSize(const std::map<int, int>& encode);

		static int GetEncodeSize(const Vector1i1& shapes);

		static std::vector<int> GetShapes(const std::map<int, int>& encode);

		static bool CheckInclude(const std::map<int, int>& big_encodes, const std::map<int, int>& small_encodes);

		static bool CheckInclude(const Vector1i1& shapes_0, const Vector1i1& shapes_1);

		static std::map<int,int> Subtract(const std::map<int, int>& big_encodes, const std::map<int, int>& small_encodes);

		static std::map<int, int> Subtract(const Vector1i1& shapes_0, const Vector1i1& shapes_1);
		static bool CheckIdentify(const std::map<int, int>& encode_0, const std::map<int, int>& encode_1);
		static bool CheckIdentify(const Vector1i1& shapes_0, const Vector1i1& shapes_1);

		static std::map<int, int> Union(const std::vector<std::map<int, int>>& encodes);
		static std::map<int, int> Union(const std::map<int, int>& encode_0, const std::map<int, int>& encode_1);

		static std::map<int, int> Union(const Vector1i1& shapes_0, const Vector1i1& shapes_1);
	};

	struct PCEShape
	{
	public:

		PCEShape();
		PCEShape(const TopoDS_Shape& ashape_);

		PCEShape(const PCEShape& pce_shape);

		void Clear();

		int index=-1;

		TopoDS_Shape ashape;
		Vector3d2 surfs;
		std::vector<std::pair<Vector3d1, double>> holes;
		Vector3d1 surf_normals;
		Vector3d center;
		glm::dmat4 LGM;//design space => assembling space in freecad

		double volume;

		glm::dmat4 init_LUM;
		int WLC_index;
		std::vector<PCETransform> transforms;
	};



	struct PCEPosition
	{
	public:

		enum { FreeStyle, FixedStyle_1_Constraint, FixedStyle_2_Constraints };

		Vector2d t;
		int style;
		glm::dmat4 M;
		Vector2d center;
		double score;

		PCETransform* transform;
		std::vector<int> edge_indexes;

		void Clear();
		void ComputeM();
		void UpdateEdges(std::vector<PCEEdge>& edges_);
		PCEPosition();
		//PCEPosition(int part_index_);
	};

	struct PCEArrange
	{
	public:
		int index;
		int eclass_index;
		int cutting_number;
		int wl_index;
		int wlc_index;
		bool valid;
		string encode_0;
		string encode_1;
		string encode_packing_0;
		string encode_packing_1;
		double volume_rate;
		std::vector<PCEPosition> positions;
		std::vector<PCEEdge> cutting_lines;
		std::vector<int> shape_order;

		std::string Vector2dString(Vector2d1 ps, std::vector<int> ints, int wl_index_);

		Vector1i2 cutting_orders;

		bool perfect_stock=false;

		std::vector<PCEConstraint> sharing_constraints;
		std::vector<PCEEdge> edges;

		std::vector<std::vector<std::tuple<std::list<std::string>, int>>> prog_strs;

		PCEArrange(
			int index_, 
			int wlc_index_,
			int wl_index_, 
			const std::vector<PCEPosition>& positions_, 
			const std::vector<PCEEdge>& cutting_lines_, 
			const double& volume_rate_);
		PCEArrange();

		void Clear();
		Vector3d3 GetMesh(Vector3d2 &surfs = Vector3d2()) const;
		void Encode();
	};

	struct OrderNode
	{
	public:
		int index=-1;
		int level=-1;
		int level_index=-1;
		int parent=-1;
		bool visited=false;
		PCEPosition position;
		std::vector<PCEConstraint> sharing_constraints;
		std::vector<PCEEdge> edges;
		
		//std::pair<int, int> arrange_eclass_index;
		std::pair<int, int> ct_arrange_eclass_index;

		OrderNode(const int &index_,
			const int &level_,
			const int &level_index_,
			const int &parent_,
			const bool& visited_,
			const PCEPosition& position_, 
			const std::vector<PCEConstraint>& sharing_constraints_,
			const std::vector<PCEEdge>& edges_) :
			index(index_), 
			level(level_),
			level_index(level_index_),
			parent(parent_),
			visited(visited_),
			position(position_), 
			sharing_constraints(sharing_constraints_),
			edges(edges_),
			//arrange_eclass_index(std::pair<int, int>(-1,-1)),
			ct_arrange_eclass_index(std::pair<int, int>(-1, -1)){};

		OrderNode(
			const PCEPosition& position_,
			const std::vector<PCEConstraint>& sharing_constraints_,
			const std::vector<PCEEdge>& edges_) :
			index(-1),
			level(-1),
			level_index(-1),
			parent(-1),
			visited(false),
			position(position_),
			sharing_constraints(sharing_constraints_),
			edges(edges_),
			//arrange_eclass_index(std::pair<int, int>(-1, -1)),
			ct_arrange_eclass_index(std::pair<int, int>(-1, -1)) {};

		OrderNode(const std::vector<PCEEdge>& edges_) :
			index(-1),
			level(-1),
			level_index(-1),
			parent(-1),
			visited(false),
			edges(edges_),
			//arrange_eclass_index(std::pair<int, int>(-1, -1)),
			ct_arrange_eclass_index(std::pair<int, int>(-1, -1)) {};

		OrderNode():
			index(-1),
			level(-1),
			level_index(-1),
			parent(-1),
			visited(false),
			//arrange_eclass_index(std::pair<int, int>(-1, -1)),
			ct_arrange_eclass_index(std::pair<int, int>(-1, -1)) {}; 
		
		void Clear(){
			std::vector<PCEConstraint>().swap(sharing_constraints);
			std::vector<PCEEdge>().swap(edges);
			position.Clear();
		}


	};

	struct PCEEGraphContainer
	{
	public:
		std::vector<shared_ptr<EGraphClass>> heads;

		std::vector<shared_ptr<PCEArrange>> ptr_arranges;
		std::vector<shared_ptr<EGraphClass>> ptr_egraph;

		std::vector<int> full_part_arranges;

		void Clear();

		static std::pair<int, int> CTNewArrange(
			std::vector<shared_ptr<PCEArrange>>& ptr_arranges_,
			std::vector<shared_ptr<EGraphClass>>& ptr_egraph_,
			const WLCollection& wlc_, const WL& wl_,
			const std::vector<PCEPosition>& fixed_positions_,
			const std::vector<PCEEdge>& cutting_lines_,
			const std::vector<PCEConstraint>& sharing_constraints_,
			const std::vector<PCEEdge>& edges_
			);

		static int CTNewEClass(std::vector<shared_ptr<EGraphClass>>& ptr_egraph_, Vector1i1& shapes);

		//old version
		static void OutputArrange(
			const std::vector<WLCollection>& wlcs_, 
			const shared_ptr<PCEArrange> &ptr_arrange_,
			std::vector<TopoDS_Shape>& transfShapes, 
			std::vector<PCEEdge>& cutting_lines, 
			std::vector<PCEHole>& holes,
			std::vector<std::vector<Vector2d>>& offset_polygons, 
			Vector3d& woods, boost::uuids::uuid& uid);

		void SetEclassValid(shared_ptr<EGraphClass> eclass);

		static void GetRelatedEclass(const std::vector<shared_ptr<EGraphClass>>&ptr_egraph,const shared_ptr<EGraphClass> eclass, Vector1i1 &rel_ecs);
	
		void SetArrangeValid(const int& head_index);

		Vector1i1 GetHeadArranges(const int& head_index);
		Vector1i1 GetHeadEClasses(const int& head_index);

		void Pruning(const int head_index, const int& prunning_limitation, const bool& random_pruning);

		void ComputeMoreSubSets();

		shared_ptr<EGraphClass> GetHead(const int& head_index);
	};

	struct ArrNode
	{
		//int arrange_id = -1;
		std::vector<std::vector<OrderNode>> order_levels;
		int wlc_index = -1;
		int wl_index = -1;
		Vector1i1 shape_order;

		void Clear()
		{
			Vector1i1().swap(shape_order);
			std::vector<std::vector<OrderNode>>().swap(order_levels);
		}
	};

	class PCE
	{
	public:
		PCE(const std::string& pref_path = "");

		PCE(bool b, const std::string& packing_o_folder_);


		void Visualization(
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
			const std::string& path_);

		PCEEGraphContainer& GetEGraphCT() { return ct_egraph; };
		std::string& GetPackingFolder() { return packing_o_folder; };
		std::map<std::pair<int, int>, double>& GetPtrShapesPacking() { return ptr_shapes_packing; };

		static void PartsRegistration(GLMATS& glmats, const std::vector<shared_ptr<PCEShape>>& ptr_shapes_);


		void ComputeShapeScore(
			std::vector<WLCollection>& wlcs_,
			const std::vector<shared_ptr<PCEShape>>& ptr_shapes_,
			std::map<std::pair<int, int>, double>& ptr_shapes_packing_);

		static void BuildShapes(
			std::vector<WLCollection>& wlcs_,
			const int &ct_egraph_heads_size_,
			const GLMATS& glmats,
			const std::string& packing_o_folder_,
			std::vector<shared_ptr<PCEShape>>& ptr_shapes_,
			std::vector<int>& parts_encodes_,
			const bool WLCFree = false);

		std::vector<PCEArrange> LocalPacking(WLCollection& wlc, const PCEArrange& arrange);

		void EnumerateOrders(std::vector<WLCollection>& wlcs_, const std::vector<int>& parts_encodes_);//2. egraph arrange

		void ClearAll();

		//PHELM* langCompiler;

		//void SetPHELM(PHELM& langCompile_)
		//{
		//	langCompile = &langCompile_;
		//};

		Vector1i3 ComputeGroupCombOrderDel(const Vector1i2& shape_groups, const Vector1i1& delete_shapes, const int max_nb);
		Vector1i2 ComputeGroupCombOrderAdd(const Vector1i2& shape_groups, const Vector1i1& delete_shapes, const int max_nb);

		Vector1i1 ArrangeDAI(const Vector1i1& arrange_ids, const Vector1i1& delete_shapes, const Vector1i1& add_shapes);


		void PostProcess(const std::vector<WLCollection>& wlcs_);//3. post process

		ArrNode AddOperation(const ArrNode& old_arr_node, const Vector1i1& add_shape_ids);
		ArrNode SequencePacking(const int& wlc_index_, const int& wl_index_, const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms);

		void GetCTArrangeEgraphIndex(const std::vector<std::vector<OrderNode>>& order_levels,const OrderNode& order_node,
			std::vector<int>& ct_arrange_indexes,std::vector<int>& ct_eclass_indexes);

		VectorPI1 GetCTArrangeEgraphIndex(const std::vector<std::vector<OrderNode>>& order_levels, const OrderNode& order_node);
		VectorPI1 GetCTArrangeEgraphIndex(const ArrNode& arr_node);

	private:
		
		void Init(const std::string& pref_path = "");
		void Init(bool b, const std::string& packing_o_folder_);




		void OutputFixedPositions(const WL& wl, const shared_ptr<PCEArrange> arrange, const std::string& path);
		void OutputFixedPositions(const WL& wl, const std::vector<PCEPosition>& positions, const std::vector<PCEEdge>& cutting_lines, const std::string& path);
		void OutputFixedPositionsVIS(
			const std::string& pref_str,
			const WL& wl,
			const shared_ptr<PCEArrange> arrange,
			const std::vector<PCEEdge>& cutting_lines,
			const Vector1i1& cutting_line_ids,
			const Vector3d1& cutting_line_colors,
			const std::vector<std::pair<Vector3d, Vector3d>>& refer_edges,
			const Vector3d2& refer_faces,
			const Vector3d1 &refer_values,
			const Vector3d &delta_y, 
			std::ofstream &export_fie, 
			int& export_int);
		
		void ArrangePartbyConstraints(const WL& wl, std::vector<PCEEdge>& edges, const std::vector<PCEPosition>& fixed_positions, PCEPosition& free_position,
			std::vector<PCEConstraint>& sharing_edges);

		void ArrangePartbySlider(const WL& wl, std::vector<PCEEdge>& edges, const std::vector<PCEConstraint>& constraints,
			const std::vector<PCEPosition>& fixed_positions, PCEPosition& free_position, std::vector<PCEConstraint>& sharing_edges);

		void GetCuttingLines(const std::vector<PCEEdge>& edges, const std::vector<PCEConstraint>& sharing_edges,
			const std::vector<PCEPosition>& fixed_positions, std::vector<PCEEdge>& cutting_lines);

		int GetScore(const std::vector<PCEEdge>& edges, const std::vector<PCEPosition>& fixed_positions, const PCEPosition& free_position,
			const std::vector<PCEConstraint>& constraints, std::vector<PCEConstraint>& comb_edges, 
			double& total_length, double& boundary_distance);

		void Order2Arrange(
			std::vector<std::string> &enumerate_order_strs,
			WLCollection& wlc, WL& wl,
			OrderNode& parent_order_node,
			const int& maximal_arrange_nb,
			const std::map<int, int>& all_targeting_shapes,
			const std::map<int, Vector1i1>& shape_transforms,
			//const VectorPI1& seq_shapes_transforms,
			const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms,
			std::vector<std::vector<OrderNode>>& order_levels,
			const bool depth_searching_bool,
			const int start_shape_index=-1);

		void BreadthFirstOrder(
			std::vector<std::string>& enumerate_order_strs,
			WLCollection& wlc, WL& wl,
			const int& maximal_arrange_nb,
			const std::map<int, int>& all_targeting_shapes,
			const std::map<int, Vector1i1>& shape_transforms,
			//const VectorPI1& seq_shapes_transforms,
			const std::vector<std::pair<int, Vector1i1>>& seq_shapes_transforms,
			std::vector<std::vector<OrderNode>>& order_levels,
			const bool depth_searching_bool);

		HMODULE *hModule;
		//PHELM* langCompiler;
		//std::vector<WLCollection> wlcs;
		
		std::map<std::pair<int, int>, double> ptr_shapes_packing;
		
		PCEEGraphContainer ct_egraph;

		NumberLabel number_label;

		//parameters
		double angle_match_error = 0.08;
		double part_match_error;
		double kerf;
		double maximal_cutting_distance;
		int    total_order_number;
		int    prunning_limitation;
		bool   random_pruning;
		double safe_distance;
		bool   output_all_log;

		std::string packing_o_folder;

		//delete later
		int tree_searching_branch;

		std::vector<std::string> enumerate_order_strs_;
	};
}


#endif // _PCE_COMBINE_H
