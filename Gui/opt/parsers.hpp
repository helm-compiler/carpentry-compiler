#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/edge_list.hpp>
#include <llhelm.hpp>
#include <egraph.hpp>

#ifndef _PARSERS_HEADER_
#define _PARSERS_HEADER_

namespace OptCompiler
{
	namespace SexpParser
	{
		typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, LLHelm*, boost::property<
			boost::vertex_index_t, std::size_t>,
			boost::no_property> Graph;
		typedef boost::graph_traits<Graph>::vertex_descriptor GVertex;
		typedef boost::graph_traits<Graph>::edge_descriptor GEdge;
		typedef std::unordered_map<GVertex, std::string> VTS;
		// For color map
		typedef std::vector<boost::default_color_type>::iterator RAIter;
		typedef boost::vec_adj_list_vertex_id_map<LLHelm*, size_t> ID;
		typedef boost::iterator_property_map<RAIter, ID, std::iterator_traits<RAIter>::value_type, std::iterator_traits<
			RAIter>::reference> CM;

		class Found_vertex_exception : public std::exception
		{
		};

		class bfs_process_visitor : public boost::default_dfs_visitor
		{
		public:
			bfs_process_visitor(GVertex& process, std::vector<std::string>& sq, std::string& im, std::string& rf, std::string& re)
				: vertx(process), outputParts(sq), inputMaterial(im),referFace(rf),referEdge(re)
			{
			}

			void discover_vertex(GVertex u, const Graph& g)
			{
				if (g[u]->get_type() == TYPE_CHOPSAW || g[u]->get_type() == TYPE_BANDSAW || g[u]->get_type() ==
					TYPE_TRACKSAW || g[u]->get_type() == TYPE_JIGSAW)
				{
					vertx = u;
				}
				else if (g[u]->get_type() == TYPE_LUMBER)
				{
					const auto uid = target(*out_edges(u, g).first, g);
					if (g[uid]->get_type() == TYPE_LLHELM)
					{
						inputMaterial = g[uid]->ToString();
						//if (startLumber.empty()) startLumber = inputMaterials.back();
					}
				}
				else if (g[u]->get_type() == TYPE_ASSIGN)
				{
					auto vars = out_edges(u, g);
					vars = out_edges(target(*vars.first, g), g);
					for (auto var = vars.first; var != vars.second; ++var)
					{
						const auto uid = target(*out_edges(target(*var, g), g).first, g);
						outputParts.push_back(g[uid]->ToString());
						//stockQuery[g[uid]->ToString()] = inputMaterials.size();
					}
				}
				else if(g[u]->get_type() == TYPE_FACE)
				{
					const auto uid = target(*out_edges(u, g).first, g);
					if (g[uid]->get_type() == TYPE_LLHELM)
					{
						referFace = g[uid]->ToString();
					}
				}
				else if (g[u]->get_type() == TYPE_EDGE)
				{
					const auto uid = target(*out_edges(u, g).first, g);
					if (g[uid]->get_type() == TYPE_LLHELM)
					{
						referEdge = g[uid]->ToString();
					}
				}
			}

			// std::unordered_map<std::string, int>& stockQuery;
			GVertex& vertx;
			std::vector<std::string>& outputParts;
			std::string& inputMaterial;

			std::string& referFace;
			std::string& referEdge;
		};

		class bfs_ref_visitor : public boost::default_bfs_visitor
		{
		public:
			bfs_ref_visitor(CM& cm, std::pair<bool, Ref*>& r)
				: colorMap(cm), ref(r)
			{
				ref.first = false;
			}

			void discover_vertex(GVertex u, const Graph& g)
			{
				const auto& curType = g[u]->get_type();
				if (curType == TYPE_REFC || curType == TYPE_REFB || curType == TYPE_REFK || curType == TYPE_REFJ || curType == TYPE_REFD)
				{
					ref.first = true;
					ref.second = static_cast<Ref*>(g[u]);
					if (ref.second->IsComplete())
						return;

					typeRef = curType;
				}

				if (curType == TYPE_HEIGHT)
				{
					auto h = target(*out_edges(u, g).first, g);
					h = target(*out_edges(h, g).first, g);
					const auto ht = static_cast<Float*>(g[h]);
					ref.second->height = ht->num;
					ref.second->hasHeight = true;
				}

				if (ref.first)
				{
					if (typeRef == TYPE_REFC)
					{
						if (g[u]->get_type() == TYPE_ANGLE || g[u]->get_type() == TYPE_LENGTH)
						{
							// TODO: more than one?
							const auto var = target(*out_edges(target(*out_edges(u, g).first, g), g).first, g);
							if (g[var]->get_type() == TYPE_FLOAT)
							{
								auto f = static_cast<Float*>(g[var]);
								auto refCont = dynamic_cast<Refc*>(ref.second);
								refCont->values.push_back(f->num);
								refCont->hasValues = true;
							}
						}
					}
					else if (typeRef == TYPE_REFB)
					{
						if (g[u]->get_type() == TYPE_TUPLE)
						{
							boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
							tie(vi, vi_end) = out_edges(u, g);
							for (next = vi; vi != vi_end; vi = next)
							{
								++next;
								auto var = target(*out_edges(target(*vi, g), g).first, g);
								if (g[var]->get_type() == TYPE_FLOAT)
								{
									auto f = dynamic_cast<Float*>(g[var]);
									auto refCont = dynamic_cast<Refb*>(ref.second);
									refCont->values.push_back(f->num);
									refCont->hasValues = true;
								}
							}
						}
					}
					else if (typeRef == TYPE_REFK)
					{
						if (g[u]->get_type() == TYPE_ANGLE || g[u]->get_type() == TYPE_LENGTH)
						{
							// TODO: more than one?
							auto var = target(*out_edges(target(*out_edges(u, g).first, g), g).first, g);
							if (g[var]->get_type() == TYPE_FLOAT)
							{
								auto f = static_cast<Float*>(g[var]);
								auto refCont = dynamic_cast<Refk*>(ref.second);
								refCont->values.push_back(f->num);
								refCont->hasValues = true;
							}
						}
					}
					else if (typeRef == TYPE_REFJ)
					{
						if (g[u]->get_type() == TYPE_TUPLE)
						{
							boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
							tie(vi, vi_end) = out_edges(u, g);
							for (next = vi; vi != vi_end; vi = next)
							{
								++next;
								auto var = target(*out_edges(target(*vi, g), g).first, g);
								if (g[var]->get_type() == TYPE_FLOAT)
								{
									auto f = dynamic_cast<Float*>(g[var]);
									auto refCont = dynamic_cast<Refj*>(ref.second);
									refCont->values.push_back(f->num);
									refCont->hasValues = true;
								}
							}
						}
					}
					else if (typeRef == TYPE_REFD)
					{
						if (g[u]->get_type() == TYPE_TUPLE)
						{
							boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
							tie(vi, vi_end) = out_edges(u, g);
							for (next = vi; vi != vi_end; vi = next)
							{
								++next;
								auto var = target(*out_edges(target(*vi, g), g).first, g);
								if (g[var]->get_type() == TYPE_FLOAT)
								{
									auto f = dynamic_cast<Float*>(g[var]);
									auto refCont = dynamic_cast<Refd*>(ref.second);
									refCont->values.push_back(f->num);
									refCont->hasValues = true;
								}
							}
						}
					}
				}
			}

			void examine_edge(GEdge e, const Graph& g)
			{
				auto targetVertex = target(e, g);

				if (ref.second != nullptr && ref.second->IsComplete())
					put(colorMap, vertex(targetVertex, g), boost::color_traits<Graph>::black());

				if (g[targetVertex]->get_type() == TYPE_CHOPSAW)
				{
					//std::cerr << "we found chop saw" << std::endl;
					//put(colorMap, vertex(targetVertex, g), boost::color_traits<SexpParser::Graph>::black());
				}

				if (g[targetVertex]->get_type() == TYPE_REFC)
				{
					//std::cerr << "we found ref" << std::endl;
				}
			}

			CM& colorMap;
			std::pair<bool, Ref*>& ref;
			Atom typeRef;
		};

		class AST
		{
		public:
			AST() = default;

			void CreateVertex(const GVertex& v, const std::string& str)
			{
				content[v] = str;
			}

			void PreProcess();

			LLHelm* GetProcess() const { return process_; }

			Ref* GetReference() const { return ref_; }

			std::string& GetInputMaterial() { return inputMaterial_; }

			std::vector<std::string>& GetOutputParts() { return outputParts_; }

			std::string& GetReferFace() { return referFace_; }
			std::string& GetReferEdge() { return referEdge_; }

			Graph g;
			VTS content;

		private:
			LLHelm* process_;
			Ref* ref_;
			std::string inputMaterial_;
			std::vector<std::string> outputParts_;
			std::string referFace_;
			std::string referEdge_;
		};

		AST Parse(const std::string& str);

		AST ParseFile(const std::string& filename);
	}

	class EClass;
	class EGraph;

	class EGraphParser
	{
	public:
		typedef std::pair<std::string, std::string> PairStrs;
		typedef std::vector<std::vector<std::vector<PairStrs>>> SubModule;

		static void ParseXML(const std::string& path, std::unordered_map<int, SubModule>& all_progs,
			std::unordered_set<int>& emptyProgs);

		static EClass* ParseEGraph(const std::string& path, EGraph* eg);
	};
}
#endif