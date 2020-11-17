// Processing-specific HELM

#include "PHELM.h"
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRepTools.hxx>
#include <gp_Dir.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TCollection_AsciiString.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRepOffset_Tool.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <TopOpeBRepBuild_Tools.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <QString>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <boost/uuid/random_generator.hpp>      // random_generator etc.

#include "GeomCommonFunctions.h"
#include "CompilerConfig.h"

#include <random>
#include <BRepPrimAPI_MakeCylinder.hxx>

#define sq QString::fromUtf8
#define u(x) QString::fromUtf8(x)

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;

inline int frac(int n)
{
	if (n > 14) n = 14;
	int res = 1;
	for (auto i = 2; i <= n+1; ++i)
	{
		res *= i;
	}
	return res;
}

Vector3d ComputeFaceCenter(const std::vector<Vector3d>& poly)
{
	Vector3d c(0, 0, 0);
	for (const auto& p : poly)
	{
		c = p + c;
	}
	c = c / static_cast<double>(poly.size());
	return c;
}

void OffsetCutFace(std::vector<Vector3d>& poly, const Vector3d& c)
{
	double eps = 1e-4;
	for (auto& p : poly)
	{
		const auto dir = eps * normalize(p - c);
		p += dir;
	}
}

TopoDS_Face MergeFaces(const std::vector<TopoDS_Face>& mergCand)
{
	BRep_Builder builder;
	TopoDS_Compound Comp;

	for (auto& sh : mergCand)
		builder.Add(Comp, sh);

	return TopoDS::Face(Comp);
}

bool MergeTopoDSFaceFromComb(
	const std::vector<TopoDS_Face>& faceSet,
	const std::vector<int>& vecInt,
	std::vector<std::vector<TopoDS_Face> >& resultFace)
{
	std::vector<TopoDS_Face> mergeCand;
	for (auto i : vecInt)
	{
		if (i == 1)
			mergeCand.push_back(faceSet[i]);
		else
			resultFace.back().push_back(faceSet[i]);
	}

	resultFace.back().emplace_back(MergeFaces(mergeCand));

	return !resultFace.back().back().IsNull();
}

void GetMFromN(int M, int N, const std::vector<TopoDS_Face>& faceSet, std::vector<std::vector<TopoDS_Face> >& resultFace)
{
	std::vector<int> vecInt(N, 0);
	for (auto i = 0; i < M; ++i) vecInt[i] = 1;

	for (auto i = 0; i < M - 1; ++i)
	{
		if (vecInt[i] == 1 && vecInt[i + 1] == 0)
		{
			//1. first exchange 1 and 0 to 0 1
			swap(vecInt[i], vecInt[i + 1]);

			//2.move all 1 before vecInt[i] to left
			sort(vecInt.begin(), vecInt.begin() + i, [](const int a, const int b)
			{
				return a > b;
			});

			MergeTopoDSFaceFromComb(faceSet, vecInt, resultFace);

			//try do step 1 and 2 from front
			i = -1;
		}
	}
}

// Two faces can be merged: if their normals are in the same plane and the have the exact same edge
QString CompilePackedParts(
	TopoDS_Shape& originalPart,
	std::vector<TopoDS_Shape>& resultParts,
	std::vector<TopoDS_Face>& cutFaces)
{
	// 0. Construct and initialize a graph
	Graph G;
	QString resultPCEHELM;
	// 1. Project the `cutFaces' to some faces with the same normals
	for (auto i = 0; i < cutFaces.size(); ++i)
	{
		for (auto j = 1 + 1; j < cutFaces.size(); ++j)
		{
			// Exist a common edge -> Wrong
			// 			TopTools_ListOfShape commonEdges;
			// 			BRepOffset_Tool::FindCommonShapes(cutFaces[i], cutFaces[j], TopAbs_EDGE, commonEdges);
			// 			if (commonEdges.Extent() >= 1)
			// 			{
			// 				boost::add_edge(i, j, G);
			// 			}

			// Should only consider co-planar case
			Handle(Geom_Surface) surfacei = BRep_Tool::Surface(cutFaces[i]);
			Handle(Geom_Surface) surfacej = BRep_Tool::Surface(cutFaces[j]);
			if (surfacei->DynamicType() == STANDARD_TYPE(Geom_Plane) &&
				surfacej->DynamicType() == STANDARD_TYPE(Geom_Plane))
			{
				Handle(Geom_Plane) planei = Handle(Geom_Plane)::DownCast(surfacei);
				Handle(Geom_Plane) planej = Handle(Geom_Plane)::DownCast(surfacej);
				if (planei->Position().IsCoplanar(planej->Position(), Precision::Confusion(), Precision::Confusion()))
					boost::add_edge(i, j, G);
			}
		}
	}

	std::vector<int> component(num_vertices(G));
	int num = connected_components(G, &component[0]);
	std::vector<int>::size_type i;
	std::cout << "Total number of components: " << num << std::endl;
	for (i = 0; i != component.size(); ++i)
		std::cout << "Vertex " << i << " is in component " << component[i] << std::endl;
	std::cout << std::endl;

	return resultPCEHELM;
}

PHELM::PHELM()
{
	hModule = &CompilerConfig::Instance().GetHModule();

	curVarIndex = 0;
	materialLoaded = false;
	toolLoaded = false;
	isChopSawEnabled = false;
	isJigSawEnabled = false;
	isBandSawEnabled = false;
	isTrackSawEnabled = false;
	randomReferences = false;
	fastMode = false;
	maxPerm = 10;

	auto& config = CompilerConfig::Instance();
	if (config.HasLoaded())
	{
		maxPerm = config.GetMaximalPerm();
		randomReferences = config.IsRandomPruning();
		PLY_WOOD = true;
		kerfDistance = config.GetKerf();
		kerfDistance2 = config.GetKerf2();
		globalFuzzyVal = 1e-2;

		if (config.HasProcess("Chopsaw"))
		{
			auto tcs = dynamic_pointer_cast<ChopsawConfig>(config.GetProcess("Chopsaw"));
			isChopSawEnabled = true;
			CHOPSAW_MAX_PERP_ANGLE = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_PERP_ANGLE);
			CHOPSAW_MIN_PERP_ANGLE = tcs->GetValue(ChopsawConfig::CHOPSAW_MIN_PERP_ANGLE);
			CHOPSAW_MAX_BEVE_ANGLE = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_BEVE_ANGLE);
			CHOPSAW_MAX_HEIGHT = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_HEIGHT);
			CHOPSAW_MAX_Y = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_Y);
			CHOPSAW_MAX_X = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_X);
			CHOPSAW_MAX_X_PRIME = tcs->GetValue(ChopsawConfig::CHOPSAW_MAX_X_PRIME);
			CHOPSAW_MIN_X_PRIME = tcs->GetValue(ChopsawConfig::CHOPSAW_MIN_X_PRIME);
			CHOPSAW_MIN_Z = tcs->GetValue(ChopsawConfig::CHOPSAW_MIN_Z);
		}

		if (config.HasProcess("Bandsaw"))
		{
			auto tbs = dynamic_pointer_cast<BandsawConfig>(config.GetProcess("Bandsaw"));
			isBandSawEnabled = true;
		}

		if (config.HasProcess("Jigsaw"))
		{
			auto tjs = dynamic_pointer_cast<JigsawConfig>(config.GetProcess("Jigsaw"));
			JIGSAW_MAX_HEIGHT = tjs->GetValue(JigsawConfig::JIGSAW_MAX_HEIGHT);
			isJigSawEnabled = true;
		}

		if (config.HasProcess("Tracksaw"))
		{
			auto tts = dynamic_pointer_cast<TracksawConfig>(config.GetProcess("Tracksaw"));
			isTrackSawEnabled = true;
		}
	}
};

bool PHELM::CompileToJigSaw(const Vector3d& cutTowards, const std::string& refS, const std::string& refE,
	double height, const boost::uuids::uuid& uid,
	const std::vector<Vector2d>& poly, const Vector2d& s,
	const Vector2d& e, std::list<std::string>& compiledCode)
{
	if (glm::abs(dot(cutTowards, Vector3d(0, 0, 1))) < 0.9999)
	{
		return false;
	}

	if (refS.empty() || refE.empty()) return false;

	if (height < JIGSAW_MAX_HEIGHT)
	{
		// Generate output
		std::stringstream ss;
		ss.setf(ios::fixed);

		ss << std::setprecision(4) << retBuffer << "Jigsaw(" << uid << ", Path(" << s.x << "/" << s.y << ", " << e.x <<
			"/" << e.y <<
			"), Ref(" << refS << ", " << refE << "), " << height / 25.4 << ")" << std::endl;
		compiledCode.push_back(ss.str());
		return true;
	}
	return false;
}

bool PHELM::CompileToBandSaw(const Vector3d& cutTowards, const std::string& refS, const std::string& refE,
	double height, const boost::uuids::uuid& uid,
	const std::vector<Vector2d>& poly, const Vector2d& s,
	const Vector2d& e, std::list<std::string>& compiledCode)
{
	if (glm::abs(dot(cutTowards, Vector3d(0, 0, 1))) < 0.9999)
	{
		return false;
	}

	if (refS.empty() || refE.empty()) return false;

	bool canCompile = false;

	const Vector2d cutDir = normalize(e - s);
	const Vector2d verDir(cutDir.y, -cutDir.x);

	double minY = DBL_MAX, maxY = -DBL_MAX;
	auto minX = DBL_MAX, maxX = -DBL_MAX;
	const double curX = dot(verDir, s);
	double curY = dot(cutDir, s);

	for (auto& v : poly)
	{
		const double tx = dot(verDir, v);
		const double ty = dot(cutDir, v);

		if (tx < minX) minX = tx;
		if (tx > maxX) maxX = tx;
		if (ty < minY) minY = ty;
		if (ty > maxY) maxY = ty;
	}

	// TODO: Use a more comprehensive way to reference points
	const double xprime = std::min(curX - minX, maxX - curX);
	const double x = maxX - minX;
	const double y = maxY - minY;
	const double l = length(e - s);
	if (x < BANDSAW_MAX_X && y < BANDSAW_MAX_Y && xprime < BANDSAW_MAX_X_PRIME && height < BANDSAW_MAX_HEIGHT)
	{
		// Generate output
		std::stringstream ss;
		ss.setf(ios::fixed);

		ss << std::setprecision(4) << retBuffer << "Bandsaw(" << uid << ", Path(" << s.x << "/" << s.y << ", " << e.x <<
			"/" << e.y <<
			"), Ref(" << refS << ", " << refE << "), " << height / 25.4 << ")" << std::endl;
		compiledCode.push_back(ss.str());

		canCompile = true;
	}

	if (isJigSawEnabled && height < JIGSAW_MAX_HEIGHT)
	{
		// Generate output
		std::stringstream ss;
		ss.setf(ios::fixed);

		ss << std::setprecision(4) << retBuffer << "Jigsaw(" << uid << ", Path(" << s.x << "/" << s.y << ", " << e.x <<
			"/" << e.y <<
			"), Ref(" << refS << ", " << refE << "), " << height / 25.4 << ")" << std::endl;
		compiledCode.push_back(ss.str());

		canCompile = true;
	}

	return canCompile;
}

bool PHELM::CompileToDrill(TopoDS_Shape& shape, PartDesignGui::PCEHole& hole, std::list<std::string>& compiledCode)
{
	if (hole.refPoints.size() != 2) return false;
	bool canCompile = false;
	// Search for the location of faces
	Vector3d cutTowards = hole.refPoints[1] - hole.refPoints[0];
	TopoDS_Face refFace;
	Vector3d faceNormal;
	for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
	{
		gp_Dir dir = GeomFunc::GetFaceNormal(TopoDS::Face(exp.Current()));
		faceNormal =  Vector3d(dir.X(), dir.Y(), dir.Z());

		double angle = Math::Functs::GetAngleBetween(cutTowards, faceNormal);
		if (Math::Functs::IsAlmostZero(angle) || Math::Functs::IsAlmostZero(angle - Math_PI))
		{
			refFace = TopoDS::Face(exp.Current());
			break;
		}
	}

	std::vector<Vector2d> polyFace;
	GetPolygonOfFace(refFace, polyFace);
	const auto base1 = GeomFunc::GetAnyPerpendicularVector(faceNormal);
	const auto base2 = glm::cross(faceNormal, base1);

	// Compute references
	std::string refN = "";
	if (glm::abs(dot(cutTowards, Vector3d(0, 0, 1))) > 0.9999)
	{
		auto px = glm::dot(hole.refPoints[0], base1);
		auto py = glm::dot(hole.refPoints[0], base2);
		Vector2d ns(px, py);
		refN = GetEdgeReference(refFace, polyFace, ns);
		std::stringstream ss;
		ss.setf(ios::fixed);
		auto sourceUid = GetShapeUID(shape);
		auto faceUid = GetShapeUID(refFace);

		ss << std::setprecision(4) << retBuffer << "Drill(" << sourceUid << ", " << faceUid << ", Ref(" << refN  << "), " << 2 * hole.radius / 25.4 << ")" << std::endl;
		compiledCode.push_back(ss.str());

		canCompile = true;
	}
	return canCompile;
}

bool PHELM::CompileToChopSaw(const bool& debug_output, TopoDS_Shape& base, const TopoDS_Face& cut_surface, std::list<TopoDS_Shape>& listShape,
	std::list<std::string>& compiledCode, bool performCut)
{
	const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

	if (debug_output)
	{

		Vector3d3 listShapes_surfs;
		listShapes_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(base));
		listShapes_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(cut_surface));

		Math::Functs::OutputObj3d(outputFolder + "\\log\\CompileToChopSaw.obj", listShapes_surfs);
	}

	//Vector3d3 debug_surfs;
	//std::vector<std::pair<Vector3d, Vector3d>> debug_edges;

	try
	{
		// neighboring face
		TopTools_IndexedDataMapOfShapeListOfShape edgeFaceMap;
		TopExp::MapShapesAndAncestors(base, TopAbs_EDGE, TopAbs_FACE, edgeFaceMap);

		// start to look at every edge
		NCollection_DataMap<TopoDS_Edge, std::vector<TopoDS_Face>, TopTools_ShapeMapHasher> mapEdgeToFace;//cut_edges and base_surfaces
		TopTools_IndexedMapOfShape cut_originMap;
		TopExp::MapShapes(cut_surface, TopAbs_EDGE, cut_originMap);

		const auto projection = (CGAL_3D_Plane_Points_Projection)GetProcAddress(CompilerConfig::Instance().GetHModule(), "CGAL_3D_Plane_Points_Projection");


		for (int i = 1; i <= cut_originMap.Extent(); ++i)
		{
			const TopoDS_Edge& cut_edge = TopoDS::Edge(cut_originMap(i));
			std::vector<TopoDS_Face> tmpFaces;
			auto pnt = GeomFunc::GetCentralEdgePoint(cut_edge);
			auto mkVertex = BRepBuilderAPI_MakeVertex(pnt);
			auto vert = mkVertex.Shape();



			Vector3d3 tmpFaces_surfs;
			for (TopExp_Explorer faceExplorer(base, TopAbs_FACE); faceExplorer.More(); faceExplorer.Next())
			{
				const TopoDS_Face& face = TopoDS::Face(faceExplorer.Current());
				//BRepExtrema_DistShapeShape distShapeVertex(face, vert);
				//distShapeVertex.Perform();
				//const double distVal = distShapeVertex.Value();

				Vector3d face_normal;
				Vector3d1 face_points;
				GeomFunc::GetOutsideContour(face, face_normal,face_points);
				Vector3d face_center = Math::Functs::GetCenter(face_points);
				Vector3d1 ps0(1, Vector3d(pnt.X(), pnt.Y(), pnt.Z()));
				Vector3d1 ps1;
				projection(face_center, face_normal, ps0, ps1);
				const double distVal_1 = Math::Functs::GetLength(ps0.back(),ps1.back());

				if (distVal_1 < CompilerConfig::Instance().GetPartMatchError())
				{
					tmpFaces.push_back(face);
					tmpFaces_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(face));
				}

				////1e-5
				//if (distVal < CompilerConfig::Instance().GetPartMatchError())
				////if (distVal < 1e-5)
				//{
				//	tmpFaces.push_back(face);
				//	tmpFaces_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(face));
				//}


				/*
				BRepClass3d_SolidClassifier classifier;
				classifier.Load(face);
				classifier.Perform(pnt, 1e-4);

				if (classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON)
				{
					mapEdgeToFace.Bind(edge, face);
					break;
				}*/
			}

			if (debug_output) Math::Functs::OutputObj3d(outputFolder + "\\log\\tmpFaces_surfs_"+std::to_string(i)+".obj", tmpFaces_surfs);

			mapEdgeToFace.Bind(cut_edge, tmpFaces);
		}

		const TopoDS_Face& cut_face = TopoDS::Face(cut_surface);
		gp_Dir cut_normal = GeomFunc::GetFaceNormal(cut_face);
		// start to look at every edge
		cut_originMap.Clear();
		TopExp::MapShapes(cut_face, TopAbs_EDGE, cut_originMap);

		bool canCompiled = false;
		std::vector<CSRef> vecCSRefs;
		for (int i = 1; i <= cut_originMap.Extent(); ++i)
		{
			const TopoDS_Edge& curEdge = TopoDS::Edge(cut_originMap(i));
			gp_Pnt st, nd;
			GeomFunc::GetEdgeEndPnt(curEdge, st, nd);

			// select the face on which the edge lies
			if (mapEdgeToFace.Seek(curEdge) != nullptr)
			{
				auto& tmpFaces = mapEdgeToFace(curEdge);
				for (const TopoDS_Face& liedFace : tmpFaces)
				{
					gp_Dir liedFaceNormal = GeomFunc::GetFaceNormal(liedFace);
					// std::cerr << "for " << i << " lied face" << std::endl;

					// see if there is a long edge that can be used to align
					for (TopExp_Explorer edgeExplorer(liedFace, TopAbs_EDGE); edgeExplorer.More(); edgeExplorer.Next())
					{
						const TopoDS_Edge& alignEdge = TopoDS::Edge(edgeExplorer.Current());

						if (GeomFunc::GetEdgeLength(alignEdge) <= 1.0)continue;

						gp_Dir alignDir[2];
						alignDir[0] = GeomFunc::GetEdgeDirection(alignEdge);
						alignDir[1] = -alignDir[0];

						for (auto k = 0; k < 2; ++k)
						{
							//debug_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(liedFace));
							//debug_edges.emplace_back(Vector3d(st.X(), st.Y(), st.Z()), Vector3d(nd.X(), nd.Y(), nd.Z()));
							//Math::Functions::OutputObj3d(outputFolder + "\\debug_surfs.obj", debug_surfs);
							//PartDesignGui::Output3DSegments_Pairs(CompilerConfig::Instance().GetHModule(), CompilerConfig::Instance().GetKerf(), outputFolder + "\\debug_edges.obj", debug_edges);

							gp_Dir zAxis = alignDir[k].Crossed(liedFaceNormal);
							gp_Pnt midPntEdge = GeomFunc::GetCentralEdgePoint(alignEdge);
							Standard_Real rmin = zAxis.XYZ().Dot(midPntEdge.XYZ());
							bool isMinimal = true;

							// make sure the aligned edge has is minimal along the z-axis direction
							for (TopExp_Explorer vExp(base, TopAbs_VERTEX); vExp.More(); vExp.Next())
							{
								const TopoDS_Vertex& mv = TopoDS::Vertex(vExp.Current());
								auto mp = BRep_Tool::Pnt(mv);
								const double tempDist = zAxis.XYZ().Dot(mp.XYZ());
								if (tempDist + CompilerConfig::Instance().GetPartMatchError() < rmin) // rounding-error
								{
									isMinimal = false;
									break;
								}
							}

							if (!isMinimal)
							{
								/*std::cerr << "not the minimal point" << std::endl; */
								continue;
							}

							// plane doesn't have direction
							double dotWithLiedFace = liedFaceNormal.Dot(cut_normal);

							/*
							if (dotWithLiedFace < 0)
							{
								cutNormal = -cutNormal;
								dotWithLiedFace = -dotWithLiedFace;
							}
							*/

							const double dotWithAlignDir = alignDir[k].Dot(cut_normal);
							// which means cut direction is the same as z-axis
							if (Math::Functs::IsAlmostZero(dotWithLiedFace) && Math::Functs::IsAlmostZero(dotWithAlignDir))
								continue;

							/////////////////////////////////////////////////////////////
							// Construct a orthogonal relationship
							TopoDS_Face targetFace;
							bool faceFound = TopOpeBRepBuild_Tools::GetAdjacentFace(liedFace, alignEdge, edgeFaceMap, targetFace);

							if (!faceFound)
							{
								continue;
							}

							auto targetFaceNormal = GeomFunc::GetFaceNormal(targetFace);
							if (!Math::Functs::IsAlmostZero(liedFaceNormal.Dot(targetFaceNormal)))
							{
								continue;
							}
							// Construction ends
							/////////////////////////////////////////////////////////////

							double perpendicularRad = std::acos(std::abs(dotWithAlignDir));
							if ((dotWithAlignDir < 0 && dotWithLiedFace > 0) ||
								(dotWithAlignDir > 0 && dotWithLiedFace < 0))
							{
								perpendicularRad = -perpendicularRad;
							}
							//double perpendicularAngle = 90.0 - std::acos(dotWithLiedFace) / M_PI * 180.0;

							/*
							// old version: 06-30-2019
							const gp_Dir projectedDir(dotWithAlignDir * alignDir[k].XYZ() + dotWithLiedFace * liedFaceNormal.XYZ());
							if (alignDir[k].Dot(cutNormal) < 0)
							{
								perpendicularRad = -perpendicularRad;
							}
							double perpendicularRad = M_PI_2 - std::acos(projectedDir.Dot(liedFaceNormal));
							*/
							double perpendicularAngle = 0.0, beveledAngle = 0.0;
							const double dotWithZAxis = zAxis.Dot(cut_normal);
							double beveledRad = M_PI_2 - std::acos(std::abs(dotWithZAxis));
							if ((dotWithAlignDir < 0 && dotWithZAxis > 0) ||
								(dotWithAlignDir > 0 && dotWithZAxis < 0))
							{
								beveledRad = -beveledRad;
							}
							if (Math::Functs::IsAlmostZero(dotWithLiedFace)) // beveled cut
							{
								beveledAngle = beveledRad / M_PI * 180.0;
								perpendicularAngle = 0.0;
								perpendicularRad = 0.0;
							}
							else
							{
								perpendicularAngle = perpendicularRad / M_PI * 180.0;
								beveledAngle = 0.0;
								beveledRad = 0.0;
							}

							{
								double depth = liedFaceNormal.XYZ().Dot(midPntEdge.XYZ());
								const double minHeight = zAxis.XYZ().Dot(midPntEdge.XYZ());
								double maxHeight = minHeight;
								double cutMinDepth = DBL_MAX;
								bool isDepth = true;
								std::vector<TopoDS_Vertex> vertices;
								for (TopExp_Explorer vExp(base, TopAbs_VERTEX); vExp.More(); vExp.Next())
								{
									const TopoDS_Vertex& mv = TopoDS::Vertex(vExp.Current());
									auto mp = BRep_Tool::Pnt(mv);

									// The height constraint makes sure that stackability can be verified in e-graph
									const double curHeight = zAxis.XYZ().Dot(mp.XYZ());
									if (curHeight > maxHeight) maxHeight = curHeight;

									// The depth constraint makes sure that no other vertices are beyond the lied face
									const double curDepth = liedFaceNormal.XYZ().Dot(mp.XYZ());
									if (curDepth > depth + 0.1)
									{
										isDepth = false;
										//break;
									}

									if (curDepth < cutMinDepth)
									{
										cutMinDepth = curDepth;
									}
								}

								// compute the right-most point, maximal x and maximal x'
								double cutMaxX = -DBL_MAX;
								double cutMinX = DBL_MAX;
								double stopY = DBL_MAX;
								for (TopExp_Explorer vExp(base, TopAbs_VERTEX); vExp.More(); vExp.Next())
								{
									const TopoDS_Vertex& mv = TopoDS::Vertex(vExp.Current());
									auto mp = BRep_Tool::Pnt(mv);
									const double tempDist = alignDir[k].XYZ().Dot(mp.XYZ());
									if (tempDist < cutMinX) cutMinX = tempDist;

									// compute the depth and update it
									if (std::fabs(tempDist - cutMaxX) < Precision::Confusion())
									{
										// same point, select a closer one
										stopY = std::min(stopY, depth - liedFaceNormal.XYZ().Dot(mp.XYZ()));
									}
									else if (tempDist > cutMaxX)
									{
										cutMaxX = tempDist;
										stopY = depth - liedFaceNormal.XYZ().Dot(mp.XYZ());
									}
								}

								//std::cerr << "Candidate x1 = " << alignDir[i].XYZ().Dot(st.XYZ()) << " " << alignDir[i].XYZ().Dot(nd.XYZ())  << " " << cutMinX << " " << cutMaxX << std::endl;

								// here might have problem, need to identify st/nd height along axis_z
								double stopX = std::min(alignDir[k].XYZ().Dot(st.XYZ()) - cutMinX,
									alignDir[k].XYZ().Dot(nd.XYZ()) - cutMinX);

								//std::cerr << "original stopx" << stopX  << " inch " << stopX /25.4 << std::endl;
								//std::cerr << "before " << stopX / 25.4 <<  stopX;
								if (!GeomFunc::IsAlmostZero(beveledRad))
								{
									stopX -= kerfDistance2 / std::cos(beveledRad);
									// std::cerr << "error " << kerfDistance2 / std::cos(beveledRad) << std::endl;
								}
								else
								{
									stopX -= kerfDistance2 / std::cos(perpendicularRad);
								}
								//std::cerr << " after" << stopX / 25.4 << std::endl;

								if (stopX < 0) stopX = 0.0;
								const double height = maxHeight - minHeight;



								/*
								std::cerr << retBuffer
									<< std::setprecision(4) << "Chopsaw(" << GetShapeUID(base)
									<< ", " << faceUID
									<< ", " << edgeUID << ", Ref("
									<< perpendicularAngle << ", "	// Ref -> Perpendicular Angle
									<< beveledAngle << ", "		// Ref-> Beveled Angle
									<< Round(stopX / 0.00254) / 10000.0 << "), " 	// Ref -> StopY
									<< Round((maxHeight - minHeight) / 0.00254) / 10000.0 << ")" << std::endl;	// Ref -> Height
								*/

								bool tb0 = perpendicularAngle < CHOPSAW_MAX_PERP_ANGLE;
								bool tb1 = perpendicularAngle > CHOPSAW_MIN_PERP_ANGLE;
								bool tb2 = std::abs(beveledAngle) < CHOPSAW_MAX_BEVE_ANGLE;
								bool tb3 = isDepth;
								bool tb4 = depth - cutMinDepth < CHOPSAW_MAX_Y;
								bool tb5 = cutMaxX - cutMinX < CHOPSAW_MAX_X;
								bool tb6 = stopX < CHOPSAW_MAX_X_PRIME;
								bool tb7 = stopX >= CHOPSAW_MIN_X_PRIME;
								bool tb8 = height < CHOPSAW_MAX_HEIGHT;
								bool tb9 = height > CHOPSAW_MIN_Z;

								if (perpendicularAngle < CHOPSAW_MAX_PERP_ANGLE &&
									perpendicularAngle > CHOPSAW_MIN_PERP_ANGLE &&
									std::abs(beveledAngle) < CHOPSAW_MAX_BEVE_ANGLE &&
									isDepth && depth - cutMinDepth < CHOPSAW_MAX_Y &&
									cutMaxX - cutMinX < CHOPSAW_MAX_X &&
									stopX < CHOPSAW_MAX_X_PRIME &&
									stopX >= CHOPSAW_MIN_X_PRIME &&
									height < CHOPSAW_MAX_HEIGHT &&
									height > CHOPSAW_MIN_Z)
								{
									canCompiled = true;
									std::stringstream ss;
									ss.setf(ios::fixed);
									auto roundPerpendicularAngle = Round(perpendicularAngle * 100.0) / 100.0;
									auto roundBeveledAngle = Round(beveledAngle * 100.0) / 100.0;

									if (stopX < 0) stopX = 0.0; // a workaround for some extreme cases

									auto faceUID = GetShapeUID(liedFace);
									auto edgeUID = GetShapeUID(alignEdge);

									ss << retBuffer
										<< std::setprecision(4) << "Chopsaw(" << GetShapeUID(base)
										<< ", " << faceUID
										<< ", " << edgeUID << ", Ref("
										<< (Math::Functs::IsAlmostZero(roundPerpendicularAngle) ? 0 : roundPerpendicularAngle) << ", " // Ref -> Perpendicular Angle
										<< (Math::Functs::IsAlmostZero(roundBeveledAngle) ? 0 : roundBeveledAngle) << ", "
										// Ref-> Beveled Angle
										<< Round(stopX / 0.00254) / 10000.0 << "), " // Ref -> StopY
										<< Round((maxHeight - minHeight) / 0.00254) / 10000.0 << ")" << std::endl;// Ref -> Height

									compiledCode.push_back(ss.str());
									vecCSRefs.emplace_back(roundPerpendicularAngle, roundBeveledAngle, stopX, height, REFChopsaw);
									//std::cerr << ss.str();
								}

								// This is only for track saw!
								if (isTrackSawEnabled && Math::Functs::IsAlmostZero(perpendicularAngle))
								{
									//if (Math::Functions::IsAlmostZero(beveledAngle))
									if (std::abs(beveledAngle) < TRACKSAW_MAX_BEVE_ANGLE)
									{
										if (stopX < TRACKSAW_MAX_X_PRIME &&
											isDepth &&
											((depth - cutMinDepth < TRACKSAW_MAX_Y && cutMaxX - cutMinX <
												TRACKSAW_MAX_X) ||
												(depth - cutMinDepth < TRACKSAW_MAX_X && cutMaxX - cutMinX <
													TRACKSAW_MAX_Y)) &&
											height < TRACKSAW_MAX_HEIGHT)
										{
											canCompiled = true;
											std::stringstream ss;
											ss.setf(ios::fixed);
											auto roundPerpendicularAngle = Round(perpendicularAngle * 100.0) /
												100.0;
											auto roundBeveledAngle = Round(beveledAngle * 100.0) / 100.0;

											if (stopX < 0) stopX = 0.0; // a workaround

											auto faceUID = GetShapeUID(liedFace);
											auto edgeUID = GetShapeUID(alignEdge);

											ss << retBuffer
												<< std::setprecision(4) << "Tracksaw(" << GetShapeUID(base)
												<< ", " << faceUID
												<< ", " << edgeUID << ", Ref("
												<< (Math::Functs::IsAlmostZero(roundPerpendicularAngle)
													? 0
													: roundPerpendicularAngle) << ", "
												// Ref -> Perpendicular Angle
												<< (Math::Functs::IsAlmostZero(roundBeveledAngle) ? 0 : roundBeveledAngle)
												<< ", " // Ref-> Beveled Angle
												<< Round(stopX / 0.00254) / 10000.0 << "), " // Ref -> StopY
												<< Round((maxHeight - minHeight) / 0.00254) / 10000.0 << ")" << std
												::endl; // Ref -> Height

											compiledCode.push_back(ss.str());
											vecCSRefs.emplace_back(roundPerpendicularAngle, roundBeveledAngle,
												stopX, height, REFTracksaw);
										}
									}
								}
							}
						}
					}
				}
			}
		}


		if (canCompiled)
		{
			// heuristic pruning, erase non-dominated references
			std::vector<std::list<std::string>::iterator> eraseRefs;
			auto cit = compiledCode.begin();
			auto p_a = 4, p_s = 4;

			// Get maximal precision of a set of references
			for(int i=0;i< vecCSRefs.size();i++)
			{
				auto ref = vecCSRefs[i];
				const auto tp_a = ref.pangle;
				const auto& tp_s = ref.pstop;
				//std::cerr << ref.angle1 <<" " << ref.angle2 << " " << tp_a << " " << ref.stop << " " << tp_s << std::endl;;

				if (tp_a > p_a || tp_s > p_s)
					continue;
				p_a = tp_a, p_s = tp_s;
			}

			// identify the references that either are duplicated or non-dominated
			std::unordered_map<std::string, int> mapDuplicated;
			for (auto it = vecCSRefs.begin(); it != vecCSRefs.end(); ++it, ++cit)
			{
				if (!GeomFunc::IsAlmostZero(it->angle1) && !GeomFunc::IsAlmostZero(it->angle2))
				{
					eraseRefs.push_back(cit);
					continue;
				}
				const auto tp_a = it->pangle;
				const auto& tp_s = it->pstop;

				if (tp_a > p_a || tp_s > p_s)
				{
					eraseRefs.push_back(cit);
					continue;
				}

				if (mapDuplicated.find(it->hash) != mapDuplicated.end())
				{
					eraseRefs.push_back(cit);
					continue;
				}

				mapDuplicated[it->hash] = 1;
			}

			if (randomReferences)
			{
				const int nbErased = eraseRefs.size();
				const int nbAllCodes = compiledCode.size();

				if (nbAllCodes != 0)
				{
					std::random_device rd;
					std::mt19937 rng(rd());
					std::uniform_int_distribution<> dis(0, nbAllCodes - 1);
					std::set<unsigned> duplicateAvoider;
					eraseRefs.clear();

					while (eraseRefs.size() < nbErased)
					{
						int nToErase = dis(rng);
						if (duplicateAvoider.find(nToErase) == duplicateAvoider.end())
						{
							eraseRefs.push_back(std::next(compiledCode.begin(), nToErase));
							duplicateAvoider.insert(nToErase);
						}
					}
				}
			}

			// erase unnecessary refs
			for (auto eit : eraseRefs)
			{
				compiledCode.erase(eit);
			}


			if (!performCut)
			{
				for (auto& c : compiledCode)
					std::cerr << c;
			}
		}
		else
		{
			//std::cerr << "Cannot compile, return" << std::endl;
			return false;
		}
	}
	catch (Standard_Failure & e)
	{
		std::cerr << e.GetMessageString() << std::endl;
		return false;
	}
	catch (Base::Exception & e)
	{
		std::cerr << e.getMessage() << std::endl;
		return false;
	}

	return true;
}

bool PHELM::SimpleCompileSuccCutting(std::list<TopoDS_Shape>& listMaterial, std::vector<PartDesignGui::PCEEdge>& cuttingLines,
	ListPairShapePoly& listShapes)
{
	std::vector<Vector2d> output_ns;
	std::vector<Vector2d> output_ne;
	std::vector<Vector3d> output_normals;
	std::vector<string> output_infos;
	int i = 0;
	return CompileSuccessiveCutting(false,listMaterial, cuttingLines, listShapes, i, output_infos, output_ns, output_ne, output_normals);
}

void PHELM::SetFastMode(bool _f)
{
	fastMode = _f;

	if (fastMode)
	{
		kerfDistance = 0.;
		kerfDistance2 = 0.;
	}
	else
	{
		auto& config = CompilerConfig::Instance();
		kerfDistance = config.GetKerf();
		kerfDistance2 = config.GetKerf2();
	}
}

bool PHELM::CompileSuccessiveCutting(
	const bool& debug_output,
	std::list<TopoDS_Shape>& listMaterial,
	std::vector<PartDesignGui::PCEEdge>& cuttingLines,
	ListPairShapePoly& listShapes,
	int &i,	// for recursive calls
	std::vector<string>& output_infos,
	std::vector<Vector2d>& output_ns, std::vector<Vector2d>& output_ne, std::vector<Vector3d>& output_normals) // for visualization
{
	// see if any part in listShapes has been finished, if any, erase it and push a returned value
	auto CheckFinish = [&]()
	{
		if (debug_output)
		{
			const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

			Vector3d3 listMaterial_surfs;
			Vector3d3 listShapes_surfs;
			for (auto newIter = listMaterial.begin(); newIter != listMaterial.end(); ++newIter)
				listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(*newIter));
			for (auto it = listShapes.begin(); it != listShapes.end(); ++it)
				listShapes_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(std::get<0>(*it)));
			Math::Functs::OutputObj3d(outputFolder + "\\log\\listMaterial_surfs_" + std::to_string(i) + ".obj", listMaterial_surfs, 2);
			Math::Functs::OutputObj3d(outputFolder + "\\log\\listShapes_surfs_" + std::to_string(i) + ".obj", listShapes_surfs, 2);
		}

		std::vector<std::list<TopoDS_Shape>::iterator> eraseListMaterial_;
		std::vector<ListPairShapePoly::iterator> eraseListPairSP_;
		for (auto newIter = listMaterial.begin(); newIter != listMaterial.end(); ++newIter)
		{
			for (auto it = listShapes.begin(); it != listShapes.end(); ++it)
			{
				bool b = GeomFunc::AreTwoShapesSame_Mesh(*newIter, std::get<0>(*it), CompilerConfig::Instance().GetPartMatchError());
				if (b)
				{
					eraseListMaterial_.push_back(newIter);
					eraseListPairSP_.push_back(it);
				}
			}
		}

		for (auto& it : eraseListMaterial_)
			listMaterial.erase(it);
		for (auto& it : eraseListPairSP_)
			listShapes.erase(it);
	};

	if (cuttingLines.empty())
	{
		if (listMaterial.size() != 1 || listShapes.size() != 1)
		{
			std::cerr << "if (listMaterial.size() != 1 || listShapes.size() != 1)" << std::endl;
			system("pause");
		}

		auto d = GetShapeUID(listMaterial.back());
		auto c = GetShapeUID(std::get<0>(listShapes.back()));

		ConvertVecUidToBuf(std::vector<boost::uuids::uuid>(1,c));

		std::stringstream ss;
		ss.setf(ios::fixed);
		ss << retBuffer
			<< std::setprecision(4) << "(" << d<< ")" << std::endl;// Ref -> Height

		WriteProgram(std::list<std::string>(1, ss.str()), -1);

		CheckFinish();

		output_infos.emplace_back("There are no cutting lines in this arrange.");
		return true;
	}

	if (i == cuttingLines.size())
	{
		output_infos.emplace_back("Have finished carving the final cutting line.");
		return true;
	}
	// first, extract all up faces on listMaterial
	std::vector<std::vector<Vector2d> > upPolylineListMaterial(listMaterial.size());
	const auto& cutLine = cuttingLines[i];
	const auto& cutLine_id = cutLine.index;
	const Vector2d cutLineDir = normalize(cutLine.cut.offset_e - cutLine.cut.offset_s);
	Vector3d ne3d = cutLine.cutting_line_e;
	Vector3d ns3d = cutLine.cutting_line_s;


	// given a cutting line and its normal, compute its cross-product as its direction
	const auto cutTowards = cross(cutLine.cut.cutting_surface_normal, Vector3d(cutLineDir[0], cutLineDir[1], 0));
	const gp_Dir cutTowardsGpDir(cutTowards[0], cutTowards[1], cutTowards[2]);
	//std::cerr << "Current Cutting Line = " << ne3d[0] << " " << ne3d[1] << " " << ne3d[2] << " " << ns3d[0] << " " << ns3d[1] << " " << ns3d[2] << std::endl;
	//std::cout << "Cut towards = " << cutTowards[0] << " " << cutTowards[1] << " " << cutTowards[2] << std::endl;


	int cnt = 0;
	auto cgal_identify_polycut_extend = (CGAL_Identify_Polycut_Extend)GetProcAddress(
		*hModule, "CGAL_Identify_Polycut_Extend");
	auto cgal_identify_polycut_notextend = (CGAL_Identify_Polycut_NotExtend)GetProcAddress(
		*hModule, "CGAL_Identify_Polycut_NotExtend");
	auto cgal_2d_poly_union = (CGAL_2D_Two_Polygons_Union)GetProcAddress(*hModule, "CGAL_2D_Two_Polygons_Union");

	// for each up faces on listMaterial
	std::vector<std::list<TopoDS_Shape>::iterator> eraseListMaterial;

	std::list<TopoDS_Shape> createdShapes;


	if (debug_output)
	{
		const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

		Vector3d3 listMaterial_surfs;
		for (auto newIter = listMaterial.begin(); newIter != listMaterial.end(); ++newIter)
			listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(*newIter));
		Math::Functs::OutputObj3d(outputFolder + "\\log\\listMaterial_input_" + std::to_string(i) + ".obj", listMaterial_surfs, 2);
	}
	for (auto it = listMaterial.begin(); it != listMaterial.end(); ++it, ++cnt)
	{
		// -in for loop- second, shorten this line if necessary
		std::vector<Vector2d> polyUpFace, polyDownFace;
		TopoDS_Face refFace;
		GetUpDownFaces(*it,true, polyUpFace, refFace);
		GetUpDownFaces(*it,false, polyDownFace);

		if (polyUpFace.size() < 3 || polyDownFace.size() < 3) continue;

		std::vector<std::vector<Vector2d> > unionedPolys;

		cgal_2d_poly_union(polyUpFace, polyDownFace, unionedPolys);

		if (unionedPolys.empty()) continue;

		upPolylineListMaterial[cnt] = unionedPolys[0];

		//std::cerr << "number of points in the up face " << upPolylineListMaterial[cnt].size() << std::endl;
		Vector2d ns, ne, nsr, ner;
		double vs = dot(cutLine.cut.offset_s, cutLineDir), ve = dot(cutLine.cut.offset_e, cutLineDir);

		//std::cerr << std::setprecision(10) << "Debug" << cutLine.s.x << " " << cutLine.s.y << " " << cutLine.e.x << " " << cutLine.e.y << std::endl;
		// -in for loop- third, identify if two ends are on edges

		if (!CGALIdentifyPolycutExtend(upPolylineListMaterial[cnt], cutLine.cut.offset_s, cutLine.cut.offset_e, ns, ne, nsr, ner))
		{
			//std::cerr << "pass a cut because of false value returned by cgal_extend" << std::endl;
			continue;
		}

		// add kerf here
		//ns = ns - kerfDistance2 * cutLineDir;
		//ne = ne + kerfDistance2 * cutLineDir;

		double nvs = dot(ns, cutLineDir), nve = dot(ne, cutLineDir);

		//std::cerr << std::setprecision(10) << ns.x << " " << ns.y << " " << ne.x << " " << ne.y << std::endl;
		//std::cerr << std::setprecision(10)  << vs << " " << ve << " " << nvs << " " << nve << std::endl;

		CompilerConfig::Instance().GetPartMatchError();

		bool b1 = true, b2 = true;
		for (auto& psp : listShapes)
		{
			const auto& poly = std::get<1>(psp);
			// -in for loop- if left end not on, cast a ray!
			if (!Math::Functs::IsAlmostZero_Double(nvs - vs, CompilerConfig::Instance().GetPartMatchError()))
			{
				if (nvs < vs)
				{
					if (b1 && !cgal_identify_polycut_notextend(poly, cutLine.cut.offset_s, ns))
					{
						/*std::vector<std::vector<Vector2d> > polyT = { poly };
						OutputRectangle("D:\\texttext.obj", polyT);
						system("pause");*/
						b1 = false;
					}
				}
			}

			// -in for loop- if right end not on, cast a ray!
			if (!Math::Functs::IsAlmostZero_Double(nve - ve, CompilerConfig::Instance().GetPartMatchError()))
			{
				if (nve > ve)
				{
					if (b2 && !cgal_identify_polycut_notextend(poly, cutLine.cut.offset_e, ne))
					{
						b2 = false;
					}
				}
			}
		}

		if (!b1)
		{
			ns = cutLine.cut.offset_s;
			nsr = ns;
		}
		if (!b2)
		{
			ne = cutLine.cut.offset_e;
			ner = ne;
		}

		if (visFlag)
		{
			output_ns.emplace_back(ns);
			output_ne.emplace_back(ne);
			output_normals.emplace_back(cutLine.cut.cutting_surface_normal);
		}

		// Compute references
		std::string refNs = "";
		std::string refNe = "";
		if (glm::abs(dot(cutTowards, Vector3d(0, 0, 1))) > 0.9999)
		{
			refNs = GetEdgeReference(refFace, polyUpFace, nsr);
			refNe = GetEdgeReference(refFace, polyUpFace, ner);
		}

		double height = std::abs(cutLine.cut.u_cutting_line_s[2] - cutLine.cut.l_cutting_line_s[2]);

		// replace with 3d vectors
		ns3d = Vector3d(ns[0], ns[1], ns3d[2]);
		ne3d = Vector3d(ne[0], ne[1], ne3d[2]);
		//std::cerr << std::setprecision(10) << "Output = " << ns.x << " " << ns.y << " " << ne.x << " " << ne.y << std::endl;

		if (GeomFunc::IsAlmostZero(length2(ne - ns))) continue;

		// std::cerr << "b1/b2 = " << b1 << " " << b2 << std::endl;

		auto sourceUid = GetShapeUID(*it);

		BRepBuilderAPI_Copy copier;
		copier.Perform(*it);
		TopoDS_Shape copiedShape = copier.Shape();
		std::list<std::string> compiledCode;
		TopoDS_Face generatedFace;
		//std::cerr << "CT: " << cutTowards[0] << " " << cutTowards[1] << " " << cutTowards[2] << std::endl;
		if (!GetGeneratedFace(debug_output, ns3d, ne3d, cutTowardsGpDir, *it, copiedShape, generatedFace))
		{
			Vector3d1 cut_poly{ cutLine.cut.u_cutting_line_s,cutLine.cut.u_cutting_line_e, cutLine.cut.l_cutting_line_e, cutLine.cut.l_cutting_line_s};
			generatedFace = GeomFunc::Poly2TopoDSFace(CompilerConfig::Instance().GetHModule(), cut_poly);
			//continue;
		}


		//static TopoDS_Face Poly2TopoDSFace(const HMODULE& hModule, const Vector3d1& poly)


		if (false)
		{
			GeomFunc::ExportToIGES(generatedFace, std::string("D:\\CutFace").append(std::to_string(i)).append(".igs").c_str());
			GeomFunc::ExportToIGES(*it, std::string("D:\\CutPart").append(std::to_string(i)).append(".igs").c_str());
		}


		bool successed = false;


		auto start = clock();
		// Execute a virtual cutting
		std::vector<boost::uuids::uuid> retUids;
		TopoDS_Shape cutPlane;
		GenerateCuttingPlane(cutPlane, ns3d, ne3d, cutTowardsGpDir);
		GEOMAlgo_Splitter splitter;
		splitter.SetFuzzyValue(globalFuzzyVal);
		splitter.AddArgument(*it);
		splitter.AddTool(cutPlane);
		splitter.Perform();

		if (debug_output)
		{
			const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

			Vector3d3 listMaterial_surfs;
			listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(*it));
			listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(cutPlane));
			Math::Functs::OutputObj3d(outputFolder + "\\log\\cut_input.obj", listMaterial_surfs,2);
		}


		TopTools_ListIteratorOfListOfShape iter(splitter.Modified(*it));
		TopoDS_Shape firstComponent;
		int ncParts = 0;

		const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();
		
		Vector3d3 cut_output_surfs;

		for (; iter.More(); iter.Next(), ++ncParts)
		{
			cut_output_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(iter.Value()));

			if (ncParts == 0) firstComponent = iter.Value();
			bool hasAnyPartIncluded = false;
			BRepClass3d_SolidClassifier classifier;
			classifier.Load(iter.Value());

			for (auto& listShape : listShapes)
			{
				classifier.Perform(std::get<2>(listShape), globalFuzzyVal);
				if (classifier.State() == TopAbs_IN || classifier.State() == TopAbs_ON)
				{
					hasAnyPartIncluded = true;
					break;
				}
			}

			if (hasAnyPartIncluded)
			{
				createdShapes.push_back(iter.Value());

				for (auto it = listShapes.begin(); it != listShapes.end(); ++it)
				{
					if (GeomFunc::AreTwoShapesSame_Mesh(createdShapes.back(), std::get<0>(*it), CompilerConfig::Instance().GetPartMatchError()))
					{
						SyncTwoShapeUID(createdShapes.back(), std::get<0>(*it));
					}
				}
				retUids.push_back(GetShapeUID(createdShapes.back()));
			}
		}

		if (ncParts == 0)continue;

		if(debug_output)
		Math::Functs::OutputObj3d(outputFolder + "\\log\\cut_output.obj", cut_output_surfs,2);

		ConvertVecUidToBuf(retUids);

		// -in for loop-  compile to Bandsaw cutting or Chopsaw cutting
		if (!b1 && !b2)
		{
			// only Jigsaw
			if (isJigSawEnabled)
				successed |= CompileToJigSaw(cutTowards, refNs, refNe, height, sourceUid, upPolylineListMaterial[cnt],
					ns, ne, compiledCode);

			if (!successed)
			{
				output_infos.emplace_back("Fail to compile for jig saw");
				return false;
			}
		}
		else if (b1 && !b2)
		{
			if (isBandSawEnabled)
				successed |= CompileToBandSaw(cutTowards, refNs, refNe, height, sourceUid,
					upPolylineListMaterial[cnt], ns, ne, compiledCode);

			if (!successed) 
			{
				output_infos.emplace_back("Fail to compile for band saw");
				return false;
			}
		}
		else if (!b1 && b2)
		{
			if (isBandSawEnabled)
				successed |= CompileToBandSaw(cutTowards, refNs, refNe, height, sourceUid,
					upPolylineListMaterial[cnt], ne, ns, compiledCode);

			if (!successed) 
			{
				output_infos.emplace_back("Fail to compile for band saw");
				return false;
			}
		}
		else // b1 and b2 are both true
		{
			bool tryChopTrack = (ncParts != 1);
			if (fastMode) 
			{
				tryChopTrack = (ncParts != 1);
			}

			if (tryChopTrack)
			{
				std::list<TopoDS_Shape> placeHolder;

				if (isChopSawEnabled)
				{
					successed |= CompileToChopSaw(debug_output, *it, generatedFace, placeHolder, compiledCode, true);
				}
					
				if (isBandSawEnabled)
				{
					successed |= CompileToBandSaw(cutTowards, refNs, refNe, height, sourceUid,
						upPolylineListMaterial[cnt], ns, ne,
						compiledCode);
				}
					
				if (!successed) 
				{
					output_infos.emplace_back("Fail to compile for chop saw and band saw");
					return false;
				}
			}
			else
			{
				// TODO: what is the case?
				output_infos.emplace_back("// ncParts == 2");
				output_infos.emplace_back("// TODO: what is the case?");
				return false;
			}
		}

		WriteProgram(compiledCode,cutLine_id);//CompileSuccessiveCutting
		eraseListMaterial.push_back(it);
	}

	for (auto& it : eraseListMaterial)
	{
		listMaterial.erase(it);
	}

	listMaterial.splice(listMaterial.end(), createdShapes);



	

	CheckFinish();

	i++;
	return CompileSuccessiveCutting(debug_output,listMaterial, cuttingLines, listShapes, i, output_infos,output_ns, output_ne, output_normals);
}

bool PHELM::CompileHoles(std::list<TopoDS_Shape>& listMaterial, std::vector<PartDesignGui::PCEHole>& holes, ListPairShapePoly& listShapes)
{


	std::list<std::string> compiledCode;
	// for each up faces on listMaterial
	std::vector<std::list<TopoDS_Shape>::iterator> eraseListMaterial;
	std::vector<ListPairShapePoly::iterator> eraseListPairSP;
	std::list<TopoDS_Shape> createdShapes;

	int hole_index = 0;
	for (auto& hole : holes)
	{
		int cnt = 0;
		auto& holePnts = hole.refPoints;
		auto pntA = gp_Pnt(holePnts[0].x, holePnts[0].y, holePnts[0].z);
		auto pntB = gp_Pnt(holePnts[1].x, holePnts[1].y, holePnts[1].z);
		auto dirHole = gp_Dir(gp_Vec(pntA, pntB));
		auto mkVertexA = BRepBuilderAPI_MakeVertex(pntA);
		auto vertA = mkVertexA.Shape();
		auto mkVertexB = BRepBuilderAPI_MakeVertex(pntB);
		auto vertB = mkVertexB.Shape();
		const TopoDS_Face* faceA = nullptr, *faceB = nullptr;

		for (auto it = listMaterial.begin(); it != listMaterial.end(); ++it, ++cnt)
		{
			for (TopExp_Explorer faceExplorer(*it, TopAbs_FACE); faceExplorer.More(); faceExplorer.Next())
			{
				const TopoDS_Face& face = TopoDS::Face(faceExplorer.Current());
				if (faceA == nullptr)
				{
					BRepExtrema_DistShapeShape distShapeVertex(face, vertA);
					distShapeVertex.Perform();
					const double distVal = distShapeVertex.Value();
					if (distVal < CompilerConfig::Instance().GetPartMatchError()) faceA = &face;
				}
				if (faceB == nullptr)
				{
					BRepExtrema_DistShapeShape distShapeVertex(face, vertB);
					distShapeVertex.Perform();
					const double distVal = distShapeVertex.Value();
					if (distVal < CompilerConfig::Instance().GetPartMatchError()) faceB = &face;
				}
			}
			
			if (faceA != nullptr && faceB != nullptr)
			{
				gp_Ax1 ax1(pntA, dirHole);
				gp_Ax2 ax2(pntA, dirHole);
				ax2.Rotate(ax1, Math_PI / 2.0);

				BRepPrimAPI_MakeCylinder mkCylinder(ax2, hole.radius, pntA.Distance(pntB));
				std::vector<boost::uuids::uuid> retUids;
				GEOMAlgo_Splitter splitter;
				splitter.SetFuzzyValue(globalFuzzyVal);
				splitter.AddArgument(*it);
				splitter.AddTool(mkCylinder.Shape());
				splitter.Perform();
				TopTools_ListIteratorOfListOfShape iter(splitter.Modified(*it));
				createdShapes.push_back(iter.Value());
				
				// TODO: Identify if two shapes are same
				// TODO: Sync uids here

				retUids.push_back(GetShapeUID(createdShapes.back()));
				ConvertVecUidToBuf(retUids);

				if (!CompileToDrill(*it, hole, compiledCode)) return false;
				WriteProgram(compiledCode, hole_index);//CompileHoles
				eraseListMaterial.push_back(it);
				break;
			}
		}


		hole_index++;
	}

	for (auto& it : eraseListMaterial)
	{
		listMaterial.erase(it);
	}

	if (!listMaterial.empty()) return false;
	listShapes.clear();
	return true;
}

int debug_debug_nb = 0;
QString PHELM::CompileEgraphEnhanced(
	PartDesignGui::PCE* pce, 
	const int& head_index,
	std::unordered_set<int>& compiledNodes, 
	std::unordered_set<int> &nonCompiledNodes,
	std::unordered_set<int>& nonCompiledArranges,
	const std::string& egraph_folder)
{
	bool debug_output = false;
	int debug_output_eclass = 52;
	int debug_output_arrange = 88;


	//clear folder
	Math::Functs::ClearFolder(egraph_folder);

	std::cerr << "\n\n";
	std::cerr << "Start to working on CompileEgraphEnhanced ... " << std::endl;
	std::cerr << debug_debug_nb << std::endl;

	debug_debug_nb++;
	QString resultPCEHELM;
	auto& wlcs = CompilerConfig::Instance().GetVecWLC();

	auto seed = 1;
	std::random_device rd;
	std::mt19937 rng(seed);

	// Blank file for the root node
	auto rtXmlFileName = std::string(egraph_folder + "\\-1.xml");
	auto rtDeclaration = R"(<?xml version="1.0" encoding="UTF-8"?>)";
	auto rtXmlDoc = std::make_shared<tinyxml2::XMLDocument>(rtXmlFileName.c_str());
	rtXmlDoc->Parse(rtDeclaration);
	auto rtRoot = rtXmlDoc->NewElement("root");
	rtXmlDoc->InsertEndChild(rtRoot);
	rtXmlDoc->SaveFile(rtXmlFileName.c_str());

	//for (auto& wlc : wlcs)
	{
		auto& allEClasses = pce->GetEGraphCT().ptr_egraph;
		auto& allArranges = pce->GetEGraphCT().ptr_arranges;
		auto& allHeads = pce->GetEGraphCT().heads;

		//allHeads.back()->egraph_edges;
		Vector1i1 rel_ecs;
		PartDesignGui::PCEEGraphContainer::GetRelatedEclass(allEClasses, pce->GetEGraphCT().GetHead(head_index), rel_ecs);

		std::cerr <<" eclass_nb:  " << allEClasses.size() << std::endl;

		for (auto i = 0; i < allEClasses.size(); ++i)
		{
			auto& ec = allEClasses[i];
			
			if (!ec->valid) continue;
			if (debug_output&&i != debug_output_eclass)continue;
			//if (ec->head_index !=head_index)continue;
			if (std::find(rel_ecs.begin(), rel_ecs.end(), ec->index) == rel_ecs.end()) continue;
			if (!ec->xmlDoc_path.empty())continue;

			const auto enIndex = ec->index;
			auto curECiter = 0;
			// Create xml file for each eclass
			auto xmlFileName = std::string(egraph_folder + "\\").append(std::to_string(enIndex)).append(".xml");
			auto declaration = R"(<?xml version="1.0" encoding="UTF-8"?>)";
			xmlDoc = std::make_shared<tinyxml2::XMLDocument>(xmlFileName.c_str());
			xmlDoc->Parse(declaration);
			tinyxml2::XMLElement* root = xmlDoc->NewElement("root");
			xmlDoc->InsertEndChild(root);

			std::cerr <<i <<" E-Class Index: " << enIndex << " Arrange_nb: " << ec->arranges.size() << std::endl;

			auto allNProgs = 0;
			for (auto& arrange : ec->arranges)
			{
				auto ptr_arrange= allArranges[arrange];

				if (debug_output&&arrange != debug_output_arrange)continue;

				std::cerr << "E-Class Index: " << enIndex << " , Arrange Index = " << arrange << std::endl;

				std::vector<TopoDS_Shape> pceParts;
				std::vector<std::vector<Vector2d> > pceOffsetFaces;
				std::vector<PartDesignGui::PCEEdge> cuttingLines;

				std::vector<double> pceDsFaces;
				Vector3d stockStyle;
				boost::uuids::uuid stockUid{};
				std::vector<PartDesignGui::PCEHole> holes;

				pce->GetEGraphCT().OutputArrange(wlcs, ptr_arrange, pceParts, cuttingLines, holes, pceOffsetFaces,stockStyle, stockUid);

				// output arrange
				/*
				Vector3d3 o_surfs;
				for (auto& pcepart : pceParts)
				{
					o_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(pcepart));
				}
				Math::Functions::OutputObj3d(egraph_folder + "\\" + std::to_string(enIndex) + "_" + std::to_string(arrange) + ".obj", o_surfs,2);*/

				
				// TODO: solve no need to cut problem
				//if (cuttingLines.empty()) continue;

				std::cerr << "Nb of cutting lines is " << cuttingLines.size() << std::endl;
				const auto maxOrder = std::min(frac(pceParts.size()), maxPerm);
				std::cerr << "max order = " << maxOrder << std::endl;

				auto countPerm = 0;
				auto cnt = 0;

				NCollection_DataMap<TopoDS_Shape, boost::uuids::uuid, TopTools_ShapeMapHasher> goalShapeMap;
				for (auto it = pceParts.begin(); it != pceParts.end(); ++it, ++cnt)
				{
					auto uid = GetShapeUID(*it);
					goalShapeMap.Bind(*it, uid);
				}

				//remove duplicated cutting order
				std::vector<int> order;
				for (auto& cuttingLine : cuttingLines)
					order.emplace_back(cuttingLine.index);

				auto r_order = order;
				std::reverse(r_order.begin(), r_order.end());

				Vector1i2 orders{order, r_order};
				while (next_permutation(order.begin(), order.end()))
				{
					auto str = Math::Functs::IntString(order,false,",");
					bool goon = true;
					for (auto& o : orders)
					{
						if (str == Math::Functs::IntString(o, false, ","))
						{
							goon = false;
							break;
						}
					}
					if (goon)
					{ 
						orders.emplace_back(order);
						if (orders.size() >= maxOrder)break;
					}
				}

				auto original_cuttingLines = cuttingLines;
				bool compile_arrange_success = false;

				bool goon_while = true;
				if(false)
				if (!ptr_arrange->prog_strs.empty())
				{
					goon_while = false;
					/*			
								{
									WriteProgram(std::get<0>(prog_str), std::get<1>(prog_str));
								}*/

					for (int str_index = 0; str_index < ptr_arrange->prog_strs.size(); str_index++)
					{
						prog = xmlDoc->NewElement("Program");
						prog->SetAttribute("id", str_index);
						prog->SetAttribute("wlc", ptr_arrange->wlc_index);
						prog->SetAttribute("arr", ptr_arrange->index);

						nProg = 0, cnt = 0;
						for (auto& prog_str : ptr_arrange->prog_strs[str_index])
						{
							WriteProgram(std::get<0>(prog_str), std::get<1>(prog_str), false);
						}
						root->InsertEndChild(prog);

						std::cerr << nProg << " programs and at node " << allNProgs << std::endl;
						++allNProgs;
						compile_arrange_success = true;
						std::cerr << "---------------------------------" << std::endl;
					}
				}
				
				Vector1i1 failed_orders_nb;
				std::vector<std::string> failed_orders_strs;
				while(goon_while)
				{
					++curECiter;
					
					if (countPerm >= maxOrder) break;
					if (!(countPerm < orders.size()))break;
					order = orders[countPerm];

					std::vector<PartDesignGui::PCEEdge> cuttingLines_;
					for (auto order_ : order) cuttingLines_.emplace_back(original_cuttingLines[order_]);
					cuttingLines = cuttingLines_;

					std::cerr << "---------------------------------" << std::endl;
					std::cerr << countPerm << " : " << Math::Functs::IntString(order, false, ",") << std::endl;

					//check in failed orders
					bool goon = true;
					for (int index=0;index<failed_orders_nb.size()&& goon;index++)
					{
						Vector1i1 temp(order.begin(), order.begin() + failed_orders_nb[index]);
						auto str = Math::Functs::IntString(temp,false,",");
						if (str == failed_orders_strs[index])
						{
							std::cerr << "This order has been failed before..." << std::endl;
							goon = false;
						}
					}

					if (goon)
					{
						std::vector<std::tuple<std::list<std::string>, int>>().swap(prog_strs);
						prog = xmlDoc->NewElement("Program");
						prog->SetAttribute("id", (int)ptr_arrange->cutting_orders.size());
						prog->SetAttribute("wlc", ptr_arrange->wlc_index);
						prog->SetAttribute("arr", ptr_arrange->index);

						//mapShapeToUuid.Clear(false);
						BRepPrimAPI_MakeBox mkBox(stockStyle.x, stockStyle.y, stockStyle.z);
						TopoDS_Shape inputShape = mkBox.Shape();
						std::list<TopoDS_Shape> pceMaterial = { inputShape };
						ListPairShapePoly listParts(pceParts.size());
						mapShapeToUuid.Bind(pceMaterial.front(), stockUid);

						nProg = 0, cnt = 0;
						for (auto it = listParts.begin(); it != listParts.end(); ++it, ++cnt)
						{
							std::get<0>(*it) = pceParts[cnt];
							std::get<1>(*it) = pceOffsetFaces[cnt];
							Base::Vector3d cp;
							GeomFunc::GetCenterOfGravity(std::get<0>(*it), cp);
							std::get<2>(*it) = gp_Pnt(cp.x, cp.y, cp.z);
							mapShapeToUuid.Bind(std::get<0>(*it), *goalShapeMap.Seek(pceParts[cnt]));
						}

						if (debug_output)
						{
							const auto& outputFolder = CompilerConfig::Instance().GetEGraphOutputFolder();

							Vector3d3 listMaterial_surfs;
							for (auto it = pceMaterial.begin(); it != pceMaterial.end(); ++it, ++cnt)
							{
								listMaterial_surfs.emplace_back(GeomFunc::LoadTopoDSGeometry(*it));
							}
							Math::Functs::OutputObj3d(outputFolder + "\\pceMaterial.obj", listMaterial_surfs);

							Math::Functs::ClearFolder(outputFolder + "\\log");
						}

						std::vector<Vector2d> output_ns;
						std::vector<Vector2d> output_ne;
						std::vector<Vector3d> output_normals;
						std::vector<string> output_infos;
						int cuttingLine_index = 0;
						auto tempResult = CompileSuccessiveCutting(debug_output, pceMaterial, cuttingLines, listParts, cuttingLine_index,
							output_infos, output_ns, output_ne, output_normals);


						if (!tempResult)
						{
							failed_orders_nb.emplace_back(cuttingLine_index + 1);
							Vector1i1 temp(order.begin(), order.begin() + failed_orders_nb.back());
							failed_orders_strs.emplace_back(Math::Functs::IntString(temp,false,","));

							std::cerr << "Failed order: " << failed_orders_strs.back() << std::endl;
							std::cerr << "!tempResult" << std::endl;
							std::cerr << "---------------------------------" << std::endl;
							for (auto& info : output_infos) std::cerr << info << std::endl;
							goto JumpOrder;
						}

						if (listParts.size() != holes.size())
						{
							failed_orders_nb.emplace_back(cuttingLine_index + 1);
							Vector1i1 temp(order.begin(), order.begin() + failed_orders_nb.back());
							failed_orders_strs.emplace_back(Math::Functs::IntString(temp,false,","));

							std::cerr << "Failed order: " << failed_orders_strs.back() << std::endl;
							std::cerr << "Fail: listParts.size() != holes.size()" << std::endl;
							std::cerr << "---------------------------------" << std::endl;
							for (auto& info : output_infos) std::cerr << info << std::endl;
							goto JumpOrder;
						}

						auto tempResultHole = CompileHoles(pceMaterial, holes, listParts);

						if (!listParts.empty())
						{
							std::cerr << "Fail: !listParts.empty()" << std::endl;
							std::cerr << "---------------------------------" << std::endl;
							for (auto& info : output_infos) std::cerr << info << std::endl;
							goto JumpOrder;
						}

						/*
						if (!listParts.empty())
						{
							std::cerr << "error: remaining " << listParts.size() << std::endl;
							auto ncParts = 0;
							for (const auto& kkk : pceMaterial)
							{
								GeomFunc::ExportToIGES(
									kkk, std::string("D:\\REMAINED-")
										 .append(std::to_string(ncParts++).append(".igs")).c_str());
							}
							for (const auto& kkk : listParts)
							{
								GeomFunc::ExportToIGES(std::get<0>(kkk),
													   std::string("D:\\GOAL-")
													   .append(std::to_string(ncParts++).append(".igs")).c_str());
							}
							system("pause");
							continue;
						}*/

						ptr_arrange->cutting_orders.emplace_back(order);
						root->InsertEndChild(prog);
						//ptr_arrange->prog_strs.emplace_back(prog_strs);

						std::cerr << nProg << " programs and at node " << allNProgs << std::endl;
						++allNProgs;
						compile_arrange_success = true;
						std::cerr << "---------------------------------" << std::endl;
					}
					
				JumpOrder:
					++countPerm;
					/*			if (next_permutation(order.begin(), order.end()))
								{
									std::vector<PartDesignGui::PCEEdge> cuttingLines_;
									for (auto order_ : order) cuttingLines_.emplace_back(original_cuttingLines[order_]);
									cuttingLines = cuttingLines_;
								}
								else
								{
									break;
								}*/
					//////

		

				};

				if (!compile_arrange_success)
				{
					nonCompiledArranges.insert(arrange);
				}
			}


			if (allNProgs != 0)
			{
				compiledNodes.insert(enIndex);
				ec->compiledNodesB = true;
			}
			else
			{
				if (!ec->arranges.empty())
					nonCompiledNodes.insert(enIndex);
				ec->compiledNodesB = false;
			}
			xmlDoc->SaveFile(xmlFileName.c_str());
			ec->xmlDoc_path = xmlFileName.c_str();
		}
	}

	return resultPCEHELM;
}

void PHELM::PostprocessEGraph(PartDesignGui::PCE* pce, const int& head_index, std::unordered_set<int>& compiledENodes,const std::string& egraph_folder)
{
	auto& wlcs = CompilerConfig::Instance().GetVecWLC();

	auto& eg = pce->GetEGraphCT().ptr_egraph;
	auto& heads = pce->GetEGraphCT().heads;
	auto& head = pce->GetEGraphCT().GetHead(head_index);
	
	Vector1i1 rel_ecs;
	PartDesignGui::PCEEGraphContainer::GetRelatedEclass(eg, head,rel_ecs);
	std::sort(rel_ecs.begin(),rel_ecs.end());

	for (auto& ec : eg)
	{
		//if (ec->head_index != head_index)
		{
			if (std::find(rel_ecs.begin(), rel_ecs.end(), ec->index) != rel_ecs.end())
			{
				if (!Math::Functs::DetectExisting(ec->xmlDoc_path))
				{
					
					auto str = ec->xmlDoc_path.substr(0, ec->xmlDoc_path.find_last_of("\\") + 1) + "programs\\" + ec->xmlDoc_path.substr(ec->xmlDoc_path.find_last_of("\\") + 1);
					auto xmlFileName = std::string(egraph_folder + "\\").append(std::to_string(ec->index)).append(".xml");
					std::string cmd_str = "copy " + str + " " + xmlFileName;
					system(cmd_str.c_str());
				}
			}
		}
	}

	for (auto& ec : eg) ec->valid = ec->compiledNodesB;

	head->valid = false;
	head->valid = CheckENodeEffectiveness(eg, head);

	//erase invalid edges

	for (int i = 0; i < eg.size(); ++i)
	{
		auto& ec = eg[i];

		if (ec->egraph_edges.empty())continue;

		std::vector<std::vector<int> > egraph_edges_;
		for (auto edges : ec->egraph_edges)
		{
			bool b = true;
			for (auto edge : edges)
			{
				auto& ec_ = eg[edge];
				b &= ec_->valid;
			}

			if (b)
			{
				egraph_edges_.emplace_back(edges);
			}
		}
		ec->egraph_edges = egraph_edges_;
	}

	std::set<int> allValidECIds;
	
	GetAllConnectedEClassId(eg, head, allValidECIds);

	for(auto ec:eg)
	{
		if (!ec->valid) continue;
		if (allValidECIds.find(ec->index) == allValidECIds.end())
			ec->valid = false;
	}
}

void PHELM::ExportResultEGraphEnhanced(PartDesignGui::PCE* pce, const int& head_index, std::unordered_set<int>& compiledNodes, const std::string& egraph_folder)
{
	auto& allEClasses = pce->GetEGraphCT().ptr_egraph;
	auto& heads = pce->GetEGraphCT().heads;
	auto& wlcs = CompilerConfig::Instance().GetVecWLC();
	auto& head = pce->GetEGraphCT().GetHead(head_index);

	Vector1i1 rel_ecs;
	PartDesignGui::PCEEGraphContainer::GetRelatedEclass(allEClasses, head, rel_ecs);

	const auto xmlFileName = std::string(egraph_folder + "\\egraph.xml");
	const auto declaration = R"(<?xml version="1.0" encoding="UTF-8"?>)";
	xmlDoc = std::make_shared<tinyxml2::XMLDocument>(xmlFileName.c_str());
	xmlDoc->Parse(declaration);
	auto root = xmlDoc->NewElement("root");

	// Combined node, id = -1
	auto cbXmlEC = xmlDoc->NewElement("EClass");
	cbXmlEC->SetAttribute("ID", -1);

	auto cbXmlEN = xmlDoc->NewElement("ENode");

	for (auto edges : head->egraph_edges)
	{
		for (auto& edge : edges)
		{
			auto cbXmlCombo = xmlDoc->NewElement("C");
			cbXmlCombo->SetText(edge);
			cbXmlEN->InsertEndChild(cbXmlCombo);
		}
	}

	cbXmlEC->InsertEndChild(cbXmlEN);
	root->InsertEndChild(cbXmlEC);

	// Combine all WLC together to get a top-layer E-Graph
	
	for (auto i = 0; i < allEClasses.size(); ++i)
	{
		const auto& ec = allEClasses[i];
		if (!ec->valid) continue;

		if (std::find(rel_ecs.begin(), rel_ecs.end(), ec->index) == rel_ecs.end())continue;

		auto xmlEC = xmlDoc->NewElement("EClass");
		//if (uncompiledNodes.find(ec.index) == uncompiledNodes.end())
		{
			xmlEC->SetAttribute("ID", ec->index);

			for (auto& en : ec->egraph_edges)
			{
				auto xmlEN = xmlDoc->NewElement("ENode");
				for (auto& combo : en)
				{
					auto xmlCombo = xmlDoc->NewElement("C");
					xmlCombo->SetText(combo);
					xmlEN->InsertEndChild(xmlCombo);
				}

				xmlEC->InsertEndChild(xmlEN);
			}
			root->InsertEndChild(xmlEC);
		}
	}

	xmlDoc->InsertEndChild(root);
	xmlDoc->SaveFile(xmlFileName.c_str());

	if (!Math::Functs::DetectExisting(xmlFileName))
	{
		std::cerr << "if (!Math::Functs::DetectExisting(xmlFileName)): " <<xmlFileName<< std::endl;
		system("pause");
	}
}