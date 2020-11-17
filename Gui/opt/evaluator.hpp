#pragma once
#include <iostream>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/edge_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/bimap.hpp>
#include <parsers.hpp>
#include <egraph.hpp>
#include <pagmo/types.hpp>

namespace OptCompiler {

	class Evaluator
	{
		typedef boost::adjacency_list<boost::hash_setS, boost::hash_setS, boost::directedS,
			boost::property<boost::vertex_index_t, std::size_t, boost::property<
			boost::vertex_name_t, int>>,
			boost::no_property> Graph;
		typedef boost::graph_traits<Graph>::vertex_descriptor GVertex;
		typedef boost::graph_traits<Graph>::edge_descriptor GEdge;
		typedef boost::bimap<int, GVertex>::value_type position;
		typedef std::vector<GVertex> container;

		typedef std::map<int, std::shared_ptr<std::vector<int>>> MapIntToUnion;

		struct cycle_detector : public boost::dfs_visitor<>
		{
			cycle_detector(bool& has_cycle)
				: _has_cycle(has_cycle)
			{
			}

			template <class Edge, class Graph>
			void back_edge(Edge, Graph&)
			{
				_has_cycle = true;
			}

		protected:
			bool& _has_cycle;
		};

		class bfs_in_edges_vis : public boost::dfs_visitor<>
		{
		public:
			bfs_in_edges_vis(GVertex& v, std::vector<GVertex>& ie) : targetV(v), inEdges(ie)
			{
			}

			void examine_edge(GEdge e, const Graph& g)
			{
				auto targetVertex = target(e, g);

				if (targetVertex == targetV)
				{
					inEdges.push_back(source(e, g));
				}
			}

			GVertex& targetV;
			std::vector<GVertex>& inEdges;
		};

	public:
		Evaluator();

		virtual ~Evaluator();

		static Evaluator& GetEvaluator(void)
		{
			return *_pcSingleton;
		}

		static void InitEvaluator(void);

		static Evaluator* _pcSingleton;

		// non-static members
		static void ConstructDependencyGraph(
			std::unordered_map<std::string, int>& stockQuery,
			std::vector<std::string> inputMaterials,
			int N,
			Graph& g);

		static void ConstructPartialDependencyGraph(std::unordered_map<std::string, int>& stockQuery,
			std::vector<std::pair<std::string, int>>& inputMaterials,
			Graph& g);

		static double EvaluateCostForSetups(Graph& g,
			Vector2<int>& sharedSetups,
			Vector2<double>& weights,
			std::vector<std::pair<int, int>>& mergedResults,
			std::vector<LLHelm*> references,
			std::unordered_map<int, int>& mergedDict,
			MapIntToUnion& mapIntToUnion,
			bool isStack = false);

		static double EvaluateCostForSetups(Graph& g, std::vector<std::pair<int, int>>& setups);

		static bool ExportProgramsToFile(std::vector<ENode*>& progs, std::ostream& file);

		static std::vector<int> PartialTopologicalSort(const std::vector<int>& elemSet,
			const std::vector<std::string>& inputMaterials,
			int curIndex,
			const boost::property_map<Graph, boost::vertex_name_t>::type&
			depGraphIndex,
			std::unordered_map<std::string, int>& stockQuery);

		static std::vector<int> PartialTopologicalSort(const std::vector<int>& elemSet,
			const boost::property_map<Graph, boost::vertex_name_t>::
			type& depGraphIndex,
			container& oriSortResult);

		static std::vector<double> EvaluatePrograms(std::vector<ENode*>& progs);

		static std::vector<double> EvaluatePrograms(const EncodedVec& vec);

		static void CollapseEdges(GVertex& b, GVertex& a, Graph& g);

		static double GetDimensionPrecision(double num);

		static double GetAnglePrecision(double num);

		static double PrecisionToDouble(int p);
	};
}