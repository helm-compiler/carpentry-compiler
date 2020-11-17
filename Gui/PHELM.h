#ifndef _PHELM_H
#define _PHELM_H
#include "PreCompiled.h"

#include <vector>
#include <unordered_map>

#include <Mod/PartDesign/App/FeaturePrimitive.h>
#include <Mod/PartDesign/App/FeatureHole.h>
#include <Mod/Sketcher/App/SketchObject.h>
#include <Mod/PartDesign/Gui/PCE.h>

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <boost/uuid/random_generator.hpp>      // random_generator etc.

#include "tinyxml2.h"
#include "GeomCommonFunctions.h"
#include "Geom2d_Edge.h"
#include "../App/GEOMAlgo_Splitter.hxx"

const double TOOL_EPS = 1e-4;
const double CHOPSAW_WIDTH_MIN = 19.9;
const double CHOPSAW_WIDTH_MAX = 40.1;
const double CHOPSAW_HEIGHT = 50.1;

typedef std::list<std::tuple<TopoDS_Shape, std::vector<Vector2d>, gp_Pnt> > ListPairShapePoly;
typedef NCollection_DataMap<TopoDS_Shape, boost::uuids::uuid, TopTools_ShapeMapHasher> MapShapeToUuid;

class PHELM
{
	class RefPnt
	{
	public:
		RefPnt(const gp_Pnt2d& p, const double d, const Handle(Geom2d_Edge) e, int id) : pnt(p), dist(d), edge(e), edgeId(id) {};

		gp_Pnt2d pnt;
		double dist;
		Handle(Geom2d_Edge) edge;
		int edgeId;
		bool operator<(const RefPnt& other) const
		{
			return dist < other.dist;
		}
	};

	enum RefType
	{
		REFChopsaw,
		REFTracksaw,
		REFJigsaw,
		REFBandsaw
	};

	struct CSRef
	{
		CSRef(double a1, double a2, double s, double h, RefType t) : angle1(a1), angle2(a2), stop(s), type(t)
		{
			pangle = std::max(GeomFunc::GetAnglePrecision(angle1), GeomFunc::GetAnglePrecision(angle2));

			// consider stop precision 1/16 inch
			pstop = GeomFunc::GetDimensionPrecision(stop);
			iangle1 = Round(angle1 * 10000.);
			iangle2 = Round(angle2 * 10000.);
			istop = Round(stop * 10000.);
			iheight = Round(h * 10000.);
			hash = std::to_string(iangle1) + "," + std::to_string(iangle2) + "," + std::to_string(istop) + "," + std::
				to_string(iheight) + "," + std::to_string(type);
		};

		bool operator==(const CSRef& rhs) const
		{
			return iangle1 == rhs.iangle1 && iangle2 == rhs.iangle2 && istop == rhs.istop && iheight == rhs.iheight;
		}

		void SetType(RefType t)
		{
			type = t;
		}

		RefType GetType()
		{
			return type;
		}

		bool IsType(RefType t)
		{
			return type == t;
		}

		double angle1;
		double angle2;
		double stop;
		int iangle1;
		int iangle2;
		int istop;
		int iheight;
		int pangle;
		int pstop;
		RefType type;
		std::string hash;
	};

public:
	typedef PartDesignGui::EGraphClass PCEEClass;
	typedef std::vector<PCEEClass> PCEEGraph;
	PHELM();
	
	static bool CheckENodeEffectiveness(std::vector<shared_ptr<PCEEClass>> eg, shared_ptr<PCEEClass> ec);
	bool GetAllConnectedEClassId(std::vector<shared_ptr<PCEEClass>> eg, shared_ptr<PCEEClass> ec, std::set<int>& ecIds);
	void TestPost(PartDesignGui::PCE* pce, shared_ptr<PCEEClass> ec, std::vector<int>& nodes);

	void PostprocessEGraph(PartDesignGui::PCE* pce, const int& head_index, std::unordered_set<int>& compiledENodes, const std::string& egraph_folder);
	void ExportResultEGraphEnhanced(PartDesignGui::PCE* pce, const int& head_index, std::unordered_set<int>& compiledNodes, const std::string& egraph_folder);
	QString CompileEgraphEnhanced(PartDesignGui::PCE* pce, const int& head_index, std::unordered_set<int>& compiledNodes, std::unordered_set<int>& nonCompiledNodes, 
		std::unordered_set<int>& nonCompiledArranges, const std::string& egraph_folder);
	void SetKerfDistance(double k) { kerfDistance = k; kerfDistance2 = k / 2.0; }
	bool CompileToChopSaw(const bool& debug_output, TopoDS_Shape& input, const TopoDS_Face& surface, std::list<TopoDS_Shape>& listShape,
		std::list<std::string>& compiledCode, bool performCut = true);
	bool CompileToBandSaw(const Vector3d& cutTowards, const std::string& refS, const std::string& refE, double height,
		const boost::uuids::uuid& uid, const std::vector<Vector2d>& poly, const Vector2d& s,
		const Vector2d& e, std::list<std::string>& compiledCode);
	bool CompileToJigSaw(const Vector3d& cutTowards, const std::string& refS, const std::string& refE, double height,
		const boost::uuids::uuid& uid, const std::vector<Vector2d>& poly, const Vector2d& s,
		const Vector2d& e, std::list<std::string>& compiledCode);
	bool CompileToDrill(TopoDS_Shape& shape, PartDesignGui::PCEHole& hole, std::list<std::string>& compiledCode);
	bool CompileHoles(std::list<TopoDS_Shape>& pceParts, std::vector<PartDesignGui::PCEHole>& holes, ListPairShapePoly& listShapes);

	bool CompileSuccessiveCutting(const bool& debug_output, std::list<TopoDS_Shape>& listMaterial, std::vector<PartDesignGui::PCEEdge>& cuttingLines,
		ListPairShapePoly& listShapes, int &i, 
		std::vector<string>& output_infos,
		std::vector<Vector2d>& output_ns, std::vector<Vector2d>& output_ne, std::vector<Vector3d>& output_normals);
	
	bool SimpleCompileSuccCutting(std::list<TopoDS_Shape>& listMaterial, std::vector<PartDesignGui::PCEEdge>& cuttingLines,
		ListPairShapePoly& listShapes);

	void SetFastMode(bool _f);

	static void GetUpDownFaces(const TopoDS_Shape& aShape, const bool& up_down, std::vector<Vector2d>& points, TopoDS_Face& face);
	static void GetUpDownFaces(const TopoDS_Shape& aShape, const bool& up_down, std::vector<Vector2d>& points);

	TopoDS_Shape GetShapeDS(const boost::uuids::uuid& uid);

	MapShapeToUuid& GetMapShapeToUuid() { return mapShapeToUuid; };

	void ClearAll();
protected:
	std::string GetEdgeReference(TopoDS_Face& face, const std::vector<Vector2d>& facePoly, const Vector2d& s);
	void GetPolygonOfFace(TopoDS_Face& face, std::vector<Vector2d>& points);
	
	void FaceToVertexList(const TopoDS_Face& face, const gp_Dir& normal, TopoDS_Wire& wire,
		std::vector<gp_Pnt>& vecPolygon, std::vector<Vector2d>& vecPolygon2D);
	void GetOffsetUpFaces(const TopoDS_Shape& aShape, std::vector<Vector2d>& points);
	void GenerateCuttingPlane(TopoDS_Shape& prism, const Vector2d& s, const Vector2d& e, const gp_Dir& dir);
	void GenerateCuttingPlane(TopoDS_Shape& prism, const std::vector<Vector3d>& poly, const double& offset);
	void GenerateCuttingPlane(TopoDS_Shape& prism, const std::vector<Vector3d>& poly, TopoDS_Face& face,
		const double& offset);
	void GenerateCuttingPlane(TopoDS_Shape& prism, const Vector3d& s, const Vector3d& e, const gp_Dir& dir);
	void GenerateCuttingSlice(TopoDS_Shape& prism, const Vector2d& s, const Vector2d& e, const gp_Dir& dir);
	void GenerateCuttingSlice(TopoDS_Shape& prism, const Vector3d& s, const Vector3d& e, const gp_Dir& dir);
	bool ConvertVecUidToBuf(std::vector<boost::uuids::uuid>& vecUids);
	bool GetGeneratedFace(const Vector2d& ns, const Vector2d& ne, const TopoDS_Shape& it, TopoDS_Shape& copiedShape,
		TopoDS_Face& generatedFace);
	bool GetGeneratedFace(const bool& debug_output, const Vector3d& ns, const Vector3d& ne, const gp_Dir& dir, const TopoDS_Shape& it,
		TopoDS_Shape& copiedShape,
		TopoDS_Face& generatedFace);
	bool CGALIdentifyPolycutExtend(const std::vector<Vector2d>& polygon, const Vector2d& s, const Vector2d& e,
		Vector2d& ns, Vector2d& ne, Vector2d& refNs, Vector2d& refNe);

private:
	boost::uuids::uuid GetShapeUID(const TopoDS_Shape& shape);

	bool SyncTwoShapeUID(const TopoDS_Shape& shapeA, const TopoDS_Shape& shapeB);
	void WriteProgram(std::list<std::string>& code, const int & cutLine_id, const bool& record=true);

	std::vector<std::string> CADCode;
	std::vector<std::string> HELMCode;
	std::vector<std::pair<double, double> > vecDrillTools;
	std::string retBuffer;
	int curVarIndex;
	bool materialLoaded;
	bool toolLoaded;
	MapShapeToUuid mapShapeToUuid;
	std::shared_ptr<tinyxml2::XMLDocument> xmlDoc;
	tinyxml2::XMLElement* prog{};
	std::vector<std::tuple<std::list<std::string>,int>> prog_strs;

	int nProg{};

	HMODULE *hModule;

	// configurations
	int maxPerm;
	bool randomReferences;
	bool PLY_WOOD;

	double CHOPSAW_MAX_PERP_ANGLE = 89.01;//Haisen: default: 60.01 
	double CHOPSAW_MIN_PERP_ANGLE = -89.01;//Haisen: default: -50.01 
	double CHOPSAW_MAX_BEVE_ANGLE = 89.01;//Haisen: default: 45.1 
	double CHOPSAW_MAX_HEIGHT = 101.61;
	double CHOPSAW_MAX_Y = 152.41;
	double CHOPSAW_MAX_X = 2438.41;
	double CHOPSAW_MAX_X_PRIME = 2500.0; //Haisen: default:914.41
	//double CHOPSAW_MIN_X_PRIME = 148;
	double CHOPSAW_MIN_X_PRIME = -1e-6;
	//double CHOPSAW_MIN_Z = 50.8;
	double CHOPSAW_MIN_Z = 0.0;
	double JIGSAW_MAX_HEIGHT = 25.41;
	double BANDSAW_MAX_HEIGHT = 152.41;
	double BANDSAW_MAX_X = 660.41;
	double BANDSAW_MAX_X_PRIME = 330.21;
	double BANDSAW_MAX_Y = 609.61;
	double TRACKSAW_MAX_HEIGHT = 25.41;
	double TRACKSAW_MAX_X = 2438.41;
	double TRACKSAW_MAX_X_PRIME = 914.41;
	double TRACKSAW_MAX_BEVE_ANGLE = 45.01;
	double TRACKSAW_MAX_Y = 1219.21;
	double SPLIT_ANGLE_THRES = M_PI / 6.0;
	double SPLIT_ERROR_THRES = 1.0;
	double kerfDistance = 3.175;
	double kerfDistance2 = 1.5875;
	double globalFuzzyVal = 1e-2;

	bool isChopSawEnabled;
	bool isBandSawEnabled;
	bool isJigSawEnabled;
	bool isTrackSawEnabled;

	bool visFlag = false;
	bool fastMode = false;
};

#endif // _PHELM_H
