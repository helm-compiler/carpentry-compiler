#include "evaluator.hpp"
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/copy.hpp>
#include <graphmap.hpp>

namespace OptCompiler
{
	/*static std::unordered_map<std::string, double> PRICE_DICT = {
		{"a032d738-d738-4880-8f50-36f4250db666", 0.45}, // 2x4x48
		{"652f5cd5-9e6e-48fa-b3f8-176f10fcdb56", 0.2}, // 2x4x24
		{"b024b206-900c-4e9a-a934-00293b4bd186", 1.0}, // 2x4x96
		{"b024b206-900c-4e9a-a934-00293b4bd188", 1.0},
		{"d7e0368b-4bcd-43a1-8037-1ddb6cecd4e2", 2.2}
	};*/

	static std::unordered_map<std::string, double> PRICE_DICT = {
		{"652f5cd5-9e6e-48fa-b3f8-176f10fcdb56", 0.75}, // 2x4x24
		{"a032d738-d738-4880-8f50-36f4250db666", 1.375}, // 2x4x48
		{"b024b206-900c-4e9a-a934-00293b4bd186", 2.5}, // 2x4x96
		{"764917a9-4b5a-4e52-a402-5cdfc7b0d055", 0.75}, // 2x4x24
		{"906eb3e2-1df2-4764-a150-e736a3a37483", 1.375}, // 2x4x48
		{"7eaa6d88-3d87-4afa-a5f5-2d0cae939480", 2.5}, // 2x4x96
		{"b054fd41-d023-4cde-b5ca-cf6356ec6b14", 0.3}, // 2x4x24
		{"32b733c8-9bb2-4c2b-b3af-4844992d5807", 0.55}, // 2x4x48
		{"805da52d-c8e4-4ec6-8cfe-8a304f52ff24", 1.0}, // 2x4x96
		{"15411ad4-b2b1-4dbe-b5f3-f6e5df128d31", 0.3}, // 2x4x24
		{"d01ecf45-d768-404c-8518-2dae6444cdff", 0.55}, // 2x4x48
		{"1df1a8cb-d30e-4c7f-93ab-7de601a610ad", 1.0}, // 2x4x96

		{"e49f6c42-16fe-4b9c-8198-0469f33182da", 0.55}, // 2x4x48
		{"aac77b14-f1cf-4950-b48c-c4d9e012fddb", 1.0}, // 2x4x96

		{"af330a8f-db6b-46ee-9427-2af4ac227167", 0.55}, // 2x4x48
		{"44335d53-2916-4f71-81b3-2a5d76d71fff", 1.0}, // 2x4x96
		{"9d786f94-62a6-4e98-9c9c-0a2a2c9c27dd", 6.5}, // 2x4x96

		{"8538f8db-adb6-4b5d-8ec0-6d6bad325ac4", 0.7}, // 2x4x48
		{"786c393f-b7fa-4910-9fc7-1492021f1a77", 1.2}, // 2x4x96
		{"69a58769-b39e-4267-a28e-aa8a09a93669", 8.5}, // 2x4x96
	};

	constexpr auto PRECISION_CHOPSAW = 1.0;
	constexpr auto PRECISION_BANDSAW = 2.5;
	constexpr auto PRECISION_JIGSAW = 6.25;
	constexpr auto PRECISION_TRACKSAW = 1.0;
	constexpr auto PRECISION_DRILL = 1.0;
	constexpr auto TIME_CHOPSAW = 1.0;
	constexpr auto TIME_BANDSAW = 4.0;
	constexpr auto TIME_JIGSAW = 10.0;
	constexpr auto TIME_TRACKSAW = 5.0;
	constexpr auto TIME_DRILL = 2.0;

	Evaluator* Evaluator::_pcSingleton = nullptr;

	Evaluator::Evaluator()
	{
	}

	Evaluator::~Evaluator()
	{
	}

	void Evaluator::InitEvaluator(void)
	{
		_pcSingleton = new Evaluator();
	}

	void Evaluator::CollapseEdges(GVertex& b, GVertex& a, Graph& g)
	{
		// First step: get all in edges
		std::vector<GVertex> inEdges;
		bfs_in_edges_vis vis(b, inEdges);
		boost::depth_first_search(g, visitor(vis));
		// std::cerr << "We found overall " << inEdges.size() << " in vertices" << std::endl;

		for (auto v : inEdges)
		{
			boost::graph_traits<Graph>::out_edge_iterator vi, vi_end, next;
			tie(vi, vi_end) = out_edges(v, g);

			for (next = vi; vi != vi_end; vi = next)
			{
				++next;
				if (target(*vi, g) == b)
				{
					remove_edge(*vi, g);
				}
			}

			boost::add_edge(v, a, g);
		}

		const auto outEdges = out_edges(b, g);
		for (auto beit = outEdges.first; beit != outEdges.second; ++beit)
		{
			if (a != beit->m_target)
			{
				boost::add_edge(a, beit->m_target, g);
			}
		}

		clear_out_edges(b, g);
		remove_vertex(b, g);
	}

	void Evaluator::ConstructDependencyGraph(
		std::unordered_map<std::string, int>& stockQuery,
		std::vector<std::string> inputMaterials,
		int N,
		Graph& g)
	{
		boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, g);
		boost::property_map<Graph, boost::vertex_name_t>::type names = get(boost::vertex_name, g);
		boost::graph_traits<Graph>::vertex_iterator vi, vend;
		boost::graph_traits<Graph>::vertices_size_type cnt = 0;

		std::unordered_map<int, GVertex> mapIntToVert;
		for (int i = 0; i < N; ++i)
		{
			auto v = boost::add_vertex(g);
			put(index, v, cnt);
			mapIntToVert[i] = v;
			put(names, v, cnt++);
		}

		for (int i = 0; i < inputMaterials.size(); ++i)
		{
			auto& im = inputMaterials[i];
			if (stockQuery.find(im) != stockQuery.end())
			{
				boost::add_edge(mapIntToVert[stockQuery[im]], mapIntToVert[i], g);
				// std::cerr << "Connect " << stockQuery[im] << " -> " << i << std::endl;
			}
		}
	}

	void Evaluator::ConstructPartialDependencyGraph(
		std::unordered_map<std::string, int>& stockQuery,
		std::vector<std::pair<std::string, int>>& inputMaterials,
		Graph& g)
	{
		boost::property_map<Graph, boost::vertex_index_t>::type index = get(boost::vertex_index, g);
		boost::property_map<Graph, boost::vertex_name_t>::type names = get(boost::vertex_name, g);
		boost::graph_traits<Graph>::vertex_iterator vi, vend;
		boost::graph_traits<Graph>::vertices_size_type cnt = 0;

		std::unordered_map<int, GVertex> mapIntToVert;
		for (auto& im : inputMaterials)
		{
			auto v = boost::add_vertex(g);
			put(index, v, cnt++);
			mapIntToVert[im.second] = v;
			put(names, v, im.second);
		}

		for (auto& im : inputMaterials)
		{
			if (stockQuery.find(im.first) != stockQuery.end())
			{
				const auto& a = stockQuery[im.first];
				const auto& b = im.second;

				if (mapIntToVert.find(a) != mapIntToVert.end() &&
					mapIntToVert.find(b) != mapIntToVert.end())
				{
					boost::add_edge(mapIntToVert[a], mapIntToVert[b], g);
					std::cerr << "Connect " << a << " -> " << b << std::endl;
				}
				else
				{
					std::cerr << "big error" << std::endl;
				}
			}
		}
	}

	double Evaluator::EvaluateCostForSetups(Graph& g,
		Vector2<int>& sharedSetups,
		Vector2<double>& weights,
		std::vector<std::pair<int, int>>& mergedResults,
		std::vector<LLHelm*> references,
		std::unordered_map<int, int>& mergedDict,
		MapIntToUnion& mapIntToUnion,
		bool isStack)
	{
		double sumOfSharedSetups = 0.0;

		// Why 2D vector here? Because A-> (b1, b2, b3)
		for (auto i = 0; i < sharedSetups.size(); ++i)
		{
			for (auto j = 0; j < sharedSetups[i].size(); ++j)
			{
				Graph tempG;
				copy_graph(g, tempG);
				boost::property_map<Graph, boost::vertex_name_t>::type copyIndex = get(boost::vertex_name, tempG);
				boost::graph_traits<Graph>::vertex_iterator vi, vend;

				std::unordered_map<int, GVertex> mapIndexVert;
				for (boost::tie(vi, vend) = vertices(tempG); vi != vend; ++vi)
				{
					// std::cerr << get(copyIndex, *vi) << std::endl;
					mapIndexVert[get(copyIndex, *vi)] = *vi;
				}
				int a = i, b = sharedSetups[i][j];
				// std::cerr << "looking for " << a << " -> " << b << std::endl;

				if (isStack)
				{
					while (mapIndexVert.find(a) == mapIndexVert.end())
					{
						if (mergedDict.find(a) != mergedDict.end())
						{
							a = mergedDict[a];
						}
					}

					while (mapIndexVert.find(b) == mapIndexVert.end())
					{
						if (mergedDict.find(b) != mergedDict.end())
						{
							b = mergedDict[b];
						}
					}

					if (a == b) continue;
				}
				else
				{
					// check if the source and target vertices are both existing
					// todo: if it is collapsed, we also need to find it
					if (mapIndexVert.find(a) == mapIndexVert.end() ||
						mapIndexVert.find(b) == mapIndexVert.end())
					{
						continue;
					}
				}

				CollapseEdges(mapIndexVert[a], mapIndexVert[b], tempG);
				// todo: if a and b are close to each other, then we must keep the order
				// otherwise, we can use arbitrary order
				// 07-05-2019 1am

				// obtain a property map for the vertex_index property
				boost::property_map<Graph, boost::vertex_index_t>::type
					index = get(boost::vertex_index, tempG);
				// initialize the vertex_index property values
				// boost::graph_traits<Graph>::vertex_iterator vi, vend;
				boost::graph_traits<Graph>::vertices_size_type cnt = 0;

				//boost::remove_vertex(Vertex(A), g);

				for (boost::tie(vi, vend) = vertices(tempG); vi != vend; ++vi)
				{
					put(index, *vi, cnt++);
					// std::cerr << copyIndex[*vi] << std::endl;
				}

				std::vector<boost::default_color_type> vertex_color(num_vertices(tempG));
				auto idmap = get(boost::vertex_index, tempG);
				auto vcmap = make_iterator_property_map(vertex_color.begin(), idmap);

				bool has_cycle = false;
				cycle_detector vis(has_cycle);
				boost::depth_first_search(tempG, visitor(vis).color_map(vcmap));
				// std::cerr << "The graph has a cycle? " << has_cycle << std::endl;

				if (!has_cycle)
				{
					// Results
					mergedResults.emplace_back(a, b);

					if (GraphMap::get_instance().GetLogMerged())
					{
						// always have a small subset of union
						if (mapIntToUnion.find(b) == mapIntToUnion.end())
						{
							mapIntToUnion[b] = std::make_shared<std::vector<int>>();
						}

						mapIntToUnion[b]->push_back(a);

						if (mapIntToUnion.find(a) != mapIntToUnion.end())
						{
							mapIntToUnion[b]->insert(mapIntToUnion[b]->end(),
								mapIntToUnion[a]->begin(),
								mapIntToUnion[a]->end());

							mapIntToUnion.erase(a);
						}
					}

					g.clear();
					// std::cerr << "Successfully merge " << a << " and " << b << std::endl;
					copy_graph(tempG, g);

					double timeSaved = 0;
					switch (references[a]->get_type())
					{
					case TYPE_REFB:
						timeSaved += TIME_BANDSAW;
						break;
					case TYPE_REFC:
						timeSaved += TIME_CHOPSAW;
						break;
					case TYPE_REFJ:
						timeSaved += TIME_JIGSAW;
						break;
					case TYPE_REFK:
						timeSaved += TIME_TRACKSAW;
						break;
					case TYPE_REFD:
						timeSaved += TIME_DRILL;
						break;
					default:
						break;
					}

					sumOfSharedSetups += timeSaved * weights[i][j];
					mergedDict[a] = b;
				}
			}
		}

		return sumOfSharedSetups;
	}

	// deprecated
	double Evaluator::EvaluateCostForSetups(Graph& g, std::vector<std::pair<int, int>>& setups)
	{
		auto sumOfSetups = 0;
		for (int i = 0; i < setups.size(); ++i)
		{
			auto& s = setups[i];
			Graph tempG;
			copy_graph(g, tempG);

			// printf("test");

			boost::property_map<Graph, boost::vertex_name_t>::type copyIndex = get(boost::vertex_name, tempG);
			boost::graph_traits<Graph>::vertex_iterator vi, vend;

			std::unordered_map<int, GVertex> mapIndexVert;
			for (boost::tie(vi, vend) = vertices(tempG); vi != vend; ++vi)
			{
				// std::cerr << get(copyIndex, *vi) << std::endl;
				mapIndexVert[get(copyIndex, *vi)] = *vi;
			}

			// check if the source and target vertices are both existing
			if (mapIndexVert.find(s.first) == mapIndexVert.end() ||
				mapIndexVert.find(s.second) == mapIndexVert.end())
			{
				continue;
			}

			CollapseEdges(mapIndexVert[s.first], mapIndexVert[s.second], tempG);

			// obtain a property map for the vertex_index property
			boost::property_map<Graph, boost::vertex_index_t>::type
				index = get(boost::vertex_index, tempG);
			// initialize the vertex_index property values
			// boost::graph_traits<Graph>::vertex_iterator vi, vend;
			boost::graph_traits<Graph>::vertices_size_type cnt = 0;

			//boost::remove_vertex(Vertex(A), g);

			for (boost::tie(vi, vend) = vertices(tempG); vi != vend; ++vi)
			{
				put(index, *vi, cnt++);
				// std::cerr << copyIndex[*vi] << std::endl;
			}

			std::vector<boost::default_color_type> vertex_color(num_vertices(tempG));
			auto idmap = get(boost::vertex_index, tempG);
			auto vcmap = make_iterator_property_map(vertex_color.begin(), idmap);

			bool has_cycle = false;
			cycle_detector vis(has_cycle);
			boost::depth_first_search(tempG, visitor(vis).color_map(vcmap));
			// std::cerr << "The graph has a cycle? " << has_cycle << std::endl;

			if (!has_cycle)
			{
				g.clear();
				//	std::cerr << "Successfully merge " << s.first << " and " << s.second << std::endl;
				copy_graph(tempG, g);
				++sumOfSetups;
			}
		}

		boost::property_map<Graph, boost::vertex_name_t>::type
			index = get(boost::vertex_name, g);
		boost::graph_traits<Graph>::vertex_iterator vi, vend;
		container c;
		topological_sort(g, std::back_inserter(c));

		return sumOfSetups;
	}

	double Evaluator::GetAnglePrecision(double num)
	{
		return std::abs(num - std::round(num));
	}

	double Evaluator::GetDimensionPrecision(double num)
	{
		double num_0 = num;
		double tmp_0 = num_0 - std::floor(num_0);
		tmp_0 = tmp_0 / 0.0625;
		auto d0 =std::abs(tmp_0 - std::round(tmp_0));

		return d0;
	}

	double Evaluator::PrecisionToDouble(int p)
	{
		return (1.0 - 1.0 / (static_cast<double>(p) + 1.0));
	}

	bool Evaluator::ExportProgramsToFile(std::vector<ENode*>& progs, std::ostream& file)
	{
		for (auto i = 0; i < progs.size(); ++i)
		{
			auto str = GraphMap::get_instance().cont[&(progs[i]->ast.g)];
			file << str.moduleID << " " << str.programID << std::endl;
		}
		return true;
	}

	std::vector<int> Evaluator::PartialTopologicalSort(const std::vector<int>& elemSet,
		const boost::property_map<Graph, boost::vertex_name_t>::type&
		depGraphIndex,
		container& oriSortResult)
	{
		std::set<int> waitForOutput;
		std::vector<int> result;

		for (auto& v : elemSet) waitForOutput.insert(v);

		for (auto rc = oriSortResult.rbegin(); rc != oriSortResult.rend(); ++rc)
		{
			int p = depGraphIndex[*rc];
			if (waitForOutput.find(p) != waitForOutput.end())
			{
				result.push_back(p);
				waitForOutput.erase(p);
			}
		}
		return result;
	}

	std::vector<int> Evaluator::PartialTopologicalSort(const std::vector<int>& elemSet,
		const std::vector<std::string>& inputMaterials, const int curIndex,
		const boost::property_map<Graph, boost::vertex_name_t>::type&
		depGraphIndex,
		std::unordered_map<std::string, int>& stockQuery)
	{
		std::vector<int> result;
		std::vector<std::pair<std::string, int>> localInputMaterial;
		localInputMaterial.reserve(elemSet.size() + 1);

		localInputMaterial.emplace_back(inputMaterials[curIndex], curIndex);
		for (auto& v : elemSet)
		{
			localInputMaterial.emplace_back(inputMaterials[v], v);
		}

		Graph localG_;
		ConstructPartialDependencyGraph(stockQuery, localInputMaterial, localG_);
		container cont_;
		topological_sort(localG_, std::back_inserter(cont_));

		for (auto rc = cont_.rbegin(); rc != cont_.rend(); ++rc)
		{
			//std::cerr << "(T3): " << depGraphIndex[*rc] << " ";
			result.push_back(depGraphIndex[*rc]);
		}

		return result;
	}
	int debug_nb = 0;
	// Thread-safe evaluator
	std::vector<double> Evaluator::EvaluatePrograms(std::vector<ENode*>& progs)
	{
		double result = 0.0;
		double materialCost = 0.0;


		std::vector<ENode*> progs_;
		for (auto i = 0; i < progs.size(); ++i)
		{
			auto& prog = progs[i];
			if (prog->ast.GetInputMaterial().empty())
			{
				auto _str = prog->_str;
				std::string output_part = _str.substr(1, _str.find_first_of("=") - 3);
				std::string inputMat = _str.substr(_str.find_first_of("=") + 3, _str.size() - _str.find_first_of("=") - 5);
				materialCost += PRICE_DICT[inputMat];
				continue;
			}
			progs_.emplace_back(progs[i]);
		}


		const auto nProgs = progs_.size();
		std::vector<LLHelm*> references;
		std::vector<int> eachProgRefs(nProgs);

		std::vector<std::string> inputMaterials;
		std::unordered_map<std::string, int> stockQuery;
		for (auto i = 0; i < nProgs; ++i)
		{
			auto& prog = progs_[i];
			auto process = prog->ast.GetProcess();
			auto reference = prog->ast.GetReference();

			references.push_back(reference);

			// construct stockquery by checking its output
			for (auto& op : prog->ast.GetOutputParts())
			{
				stockQuery[op] = inputMaterials.size();
			}

			// take input and push it back to input materials
			const auto& inputMat = prog->ast.GetInputMaterial();
			inputMaterials.push_back(inputMat);

			if (PRICE_DICT.find(inputMat) != PRICE_DICT.end())
			{
				materialCost += PRICE_DICT[inputMat];
			}

			eachProgRefs[i] = 1;
		}

		// Problem 1: every enode represents a line of code, then we have to get its process information

		Graph dependencyG;

		ConstructDependencyGraph(stockQuery, inputMaterials, nProgs, dependencyG);

		Graph dupDepenGraph;
		copy_graph(dependencyG, dupDepenGraph);
		container oriSortResult;
		topological_sort(dupDepenGraph, std::back_inserter(oriSortResult));
		auto oriDepGraphIndex = get(boost::vertex_name, dupDepenGraph);

		Vector2<int> stackSetups(references.size());
		Vector2<int> sharedSetups(references.size());
		Vector2<double> stackWeights(references.size());
		Vector2<double> sharedWeights(references.size());
		std::vector<std::pair<int, int>> resultShared, resultStack;

		auto soFarRef = 0;
		for (auto m : eachProgRefs)
		{
			for (auto i = soFarRef; i < soFarRef + m; ++i)
			{
				for (auto j = 0; j < references.size(); ++j)
				{
					// skip the references in the region it self
					if (j >= soFarRef && j < soFarRef + m)
						continue;

					auto& r1 = references[i];
					auto& r2 = references[j];
					if (r1->get_type() == r2->get_type())
					{
						// Hints for chop saw references
						if (r1->get_type() == TYPE_REFC)
						{
							auto cr1 = dynamic_cast<Refc*>(r1);
							auto cr2 = dynamic_cast<Refc*>(r2);

							auto tmpWeight = 1.0;
							//if (!IsAlmostZero(cr1->height - cr2->height)) tmpWeight = 0.4;

							if (cr1->intValues[0] - cr2->intValues[0] == 0 &&
								cr1->intValues[1] - cr2->intValues[1] == 0 &&
								cr1->intValues[2] - cr2->intValues[2] == 0)
							{
								stackSetups[i].push_back(j);
								stackWeights[i].push_back(tmpWeight);
							}
							else if ((cr1->intValues[0] - cr2->intValues[0] == 0 && cr1->intValues[1] - cr2->intValues[1] ==
								0)
								|| cr1->intValues[2] - cr2->intValues[2] == 0)
							{
								sharedSetups[i].push_back(j);
								sharedWeights[i].push_back(tmpWeight);
							}
						}
						// Hints for band saw references
						else if (r1->get_type() == TYPE_REFK)
						{
							auto cr1 = dynamic_cast<Refk*>(r1);
							auto cr2 = dynamic_cast<Refk*>(r2);

							auto tmpWeight = 1.0;
							//if (!IsAlmostZero(cr1->height - cr2->height)) tmpWeight = 0.4;

							if (cr1->intValues[0] - cr2->intValues[0] == 0 &&
								cr1->intValues[1] - cr2->intValues[1] == 0 &&
								cr1->intValues[2] - cr2->intValues[2] == 0)
							{
								stackSetups[i].push_back(j);
								stackWeights[i].push_back(tmpWeight);
							}
							else if ((cr1->intValues[0] - cr2->intValues[0] == 0 && cr1->intValues[1] - cr2->intValues[1] ==
								0)
								|| cr1->intValues[2] - cr2->intValues[2] == 0)
							{
								sharedSetups[i].push_back(j);
								sharedWeights[i].push_back(tmpWeight);
							}
						}
						// Hints for band saw references
						else if (r1->get_type() == TYPE_REFB)
						{
							auto cr1 = dynamic_cast<Refb*>(r1);
							auto cr2 = dynamic_cast<Refb*>(r2);

							auto tmpWeight = 1.0;

							const auto& nRefCr1 = cr1->intValues.size();
							const auto& nRefCr2 = cr2->intValues.size();
							auto nbLarger = std::max(nRefCr1, nRefCr2);

							auto tmpIntValue = cr1->intValues;
							tmpIntValue.insert(tmpIntValue.end(), cr2->intValues.begin(), cr2->intValues.end());
							std::sort(tmpIntValue.begin(), tmpIntValue.end());
							auto lastUnqItr = std::unique(tmpIntValue.begin(), tmpIntValue.end());
							auto nUnique = std::distance(tmpIntValue.begin(), lastUnqItr);

							auto repRatio = static_cast<double>(nUnique) / static_cast<double>(nbLarger);

							if (repRatio > 0.7)
							{
								//							sharedSetups[i].push_back(j);
								//							sharedWeights[i].push_back(repRatio);
							}
						}
						// Hints for jig saw references
						else if (r1->get_type() == TYPE_REFJ)
						{
							auto cr1 = dynamic_cast<Refj*>(r1);
							auto cr2 = dynamic_cast<Refj*>(r2);

							auto tmpWeight = 1.0;

							const auto& nRefCr1 = cr1->intValues.size();
							const auto& nRefCr2 = cr2->intValues.size();
							auto nbLarger = std::max(nRefCr1, nRefCr2);

							auto tmpIntValue = cr1->intValues;
							tmpIntValue.insert(tmpIntValue.end(), cr2->intValues.begin(), cr2->intValues.end());
							std::sort(tmpIntValue.begin(), tmpIntValue.end());
							auto lastUnqItr = std::unique(tmpIntValue.begin(), tmpIntValue.end());
							auto nUnique = std::distance(tmpIntValue.begin(), lastUnqItr);

							auto repRatio = static_cast<double>(nUnique) / static_cast<double>(nbLarger);

							if (repRatio > 0.7)
							{
								//							sharedSetups[i].push_back(j);
								//							sharedWeights[i].push_back(repRatio);
							}
						}
					}
				}
			}
			soFarRef += m;
		}

		/****************************************************************************/
		// Evaluate costs for setups
		std::unordered_map<int, int> mergedDict, stackDict;
		std::vector<std::pair<int, int>> mergedDict_pair, stackDict_pair;
		MapIntToUnion mapIntToUnion, mapIntToUnionShared;
		result += EvaluateCostForSetups(dependencyG, stackSetups, stackWeights, resultStack, references, mergedDict,
			mapIntToUnion, true);
		stackDict = mergedDict; // Replace dictionary
		result += 0.5 * EvaluateCostForSetups(dependencyG, sharedSetups, sharedWeights, resultShared, references,
			mergedDict, mapIntToUnionShared);

		/****************************************************************************/
		// Output results
		//if(true)
		if (GraphMap::get_instance().GetLogMerged())
		{
			if (true)
			{
				std::cerr << "------------------------------------" << std::endl;
				std::cerr << "The following is stacked instructions" << std::endl;
				for (auto& m : stackDict)
				{
					std::cerr << "(" << m.first << "->" << m.second << ") ";
					stackDict_pair.emplace_back(m.first, m.second);
				}

				std::cerr << std::endl << "The following is merged instructions" << std::endl;
				for (auto& m : mergedDict)
				{
					if (stackDict.find(m.first) == stackDict.end())
					{
						std::cerr << "(" << m.first << "->" << m.second << ") ";
						mergedDict_pair.emplace_back(m.first, m.second);
					}
				}
				std::cerr << std::endl << "------------------------------------" << std::endl;
			}

			container cont;
			topological_sort(dependencyG, std::back_inserter(cont));
			auto depGraphIndex = get(boost::vertex_name, dependencyG);
			std::set<int> hasVisited;
			std::vector<int> finalOrder;
			finalOrder.reserve(inputMaterials.size());

			for (auto ii = cont.rbegin(); ii != cont.rend(); ++ii)
			{
				const auto index = depGraphIndex[*ii];

				if (mapIntToUnionShared.find(index) != mapIntToUnionShared.end())
				{
					auto& elemSet = *(mapIntToUnionShared[index]);

					for (auto rit = elemSet.rbegin(); rit != elemSet.rend(); ++rit)
					{
						auto& rc = *rit;

						if (hasVisited.find(rc) == hasVisited.end() && mapIntToUnion.find(rc) != mapIntToUnion.end())
						{
							auto& elemSet_ = *(mapIntToUnion[rc]);
							for (auto rit_ = elemSet_.rbegin(); rit_ != elemSet_.rend(); ++rit_)
							{
								auto& rc_ = *rit_;

								//std::cerr << "(T4): " << rc_ << " ";
								finalOrder.push_back(rc_);
								hasVisited.insert(rc_);
							}
						}

						//std::cerr << "(T3): " << rc << " ";
						finalOrder.push_back(rc);
						hasVisited.insert(rc);
					}
				}

				if (hasVisited.find(index) == hasVisited.end() && mapIntToUnion.find(index) != mapIntToUnion.end())
				{
					auto& elemSet = *(mapIntToUnion[index]);

					for (auto rit = elemSet.rbegin(); rit != elemSet.rend(); ++rit)
					{
						auto& rc = *rit;

						//std::cerr << "(T2): " << rc << " ";
						finalOrder.push_back(rc);
						hasVisited.insert(rc);
					}
				}

				{
					//std::cerr << "(T1): " << index << " ";
					finalOrder.push_back(index);
				}
			}

			if (true)
			{
				if (finalOrder.size() != nProgs)
				{
					std::cerr << "finalOrder.size() != nProgs, exit" << std::endl;
					system("pause");
					exit(-1);
				}

				std::unordered_map<int, bool> tmpLinearity;
				for (int i = 0; i < finalOrder.size(); ++i)
				{
					auto& od = finalOrder[i];
					auto& mt = inputMaterials[od];

					if (PRICE_DICT.find(mt) == PRICE_DICT.end())
					{
						if (stockQuery.find(mt) == stockQuery.end())
						{
							std::cerr << "query error" << std::endl;
							std::cerr << mt << std::endl;
							system("pause");
						}
						int sr = stockQuery[mt];
						if (tmpLinearity.find(sr) == tmpLinearity.end())
						{
							std::cerr << "linearity problem, exit" << std::endl;
							system("pause");
							//exit(-1);
						}
					}
					tmpLinearity[od] = true;
				}
			}

			////////////////////////////////////////////////////////////////////
			//map stack and merge
			for (auto& o : stackDict_pair)
			{
				o.first = Math::Functs::VectorIndex(finalOrder, o.first);
				o.second = Math::Functs::VectorIndex(finalOrder, o.second);
			}
			auto components = Math::Functs::ConnectedComponents(finalOrder.size(),stackDict_pair);
			for (auto& o : mergedDict_pair)
			{
				o.first = Math::Functs::VectorIndex(finalOrder, o.first);
				o.second = Math::Functs::VectorIndex(finalOrder, o.second);
				o.first=components[Math::Functs::VectorIndex(components, o.first)].back();
				o.second = components[Math::Functs::VectorIndex(components, o.second)].front();
			}

			////////////////////////////////////////////////////////////////////

			if (GraphMap::get_instance().GetLogMerged())
			{
				GraphMap::get_instance().ImportOrders(finalOrder);
				GraphMap::get_instance().ImportMergedDicts(mergedDict_pair);
				GraphMap::get_instance().ImportStackDicts(stackDict_pair);

				// From order to output
				std::vector<std::tuple<int, int, int, int, int, int, double, double, double, std::string, std::string>> info;
				std::vector<std::string> strs;
				for (const auto& od : finalOrder)
				{
					strs.emplace_back("\n" + progs_[od]->_str);
					//eclassID enodeID equivalenceID lineID arrID cutLineID

					std::tuple<int, int, int, int, int, int, double, double, double, std::string, std::string>
						one_info(
							progs_[od]->eclassID,
							progs_[od]->enodeID,
							progs_[od]->equivalenceID,
							progs_[od]->lineID,
							progs_[od]->arr,
							progs_[od]->cutLineID,
							progs_[od]->ast.GetReference()->values[0],
							progs_[od]->ast.GetReference()->values[1],
							progs_[od]->ast.GetReference()->values[2],
							progs_[od]->ast.GetReferFace(),
							progs_[od]->ast.GetReferEdge());

					info.emplace_back(one_info);
					debug_nb++;
				}
				GraphMap::get_instance().ImportInfos(info);
				GraphMap::get_instance().ImportOutput(strs);
			}
		

			if(false) std::cerr << std::endl << "------------------------------------" << std::endl;
		}
		/****************************************************************************/

		double originalTime = 0;
		for (auto& ref : references)
		{
			switch (ref->get_type())
			{
			case TYPE_REFB:
				originalTime += TIME_BANDSAW;
				break;
			case TYPE_REFC:
				originalTime += TIME_CHOPSAW;
				break;
			case TYPE_REFJ:
				originalTime += TIME_JIGSAW;
				break;
			case TYPE_REFK:
				originalTime += TIME_TRACKSAW;
				break;
			case TYPE_REFD:
				originalTime += TIME_DRILL;
				break;
			default:
				break;
			}
		}

		double costPrecision = 0;
		for (auto i = 0; i < references.size(); ++i)
		{
			auto& r = references[i];
			if (r->get_type() == TYPE_REFC)
			{
				auto cr = dynamic_cast<Refc*>(r);
				costPrecision += PRECISION_CHOPSAW * (GetAnglePrecision(cr->values[0]) + GetAnglePrecision(cr->values[1]));
				costPrecision += PRECISION_CHOPSAW * (GetDimensionPrecision(cr->values[2]));
				costPrecision += PRECISION_CHOPSAW;
			}
			else if (r->get_type() == TYPE_REFB)
			{
				auto cr = dynamic_cast<Refb*>(r);

				double sumOfDimensions = 0;
				for (auto& v : cr->values)
				{
					sumOfDimensions += GetDimensionPrecision(v);
				}

				costPrecision += PRECISION_BANDSAW * sumOfDimensions;
				costPrecision += PRECISION_BANDSAW;
			}
			else if (r->get_type() == TYPE_REFJ)
			{
				auto cr = dynamic_cast<Refj*>(r);

				double sumOfDimensions = 0;
				for (auto& v : cr->values)
				{
					sumOfDimensions += GetDimensionPrecision(v);
				}

				costPrecision += PRECISION_JIGSAW * sumOfDimensions;
				costPrecision += PRECISION_JIGSAW;
			}
			else if (r->get_type() == TYPE_REFK)
			{
				auto cr = dynamic_cast<Refk*>(r);
				//costPrecision += PRECISION_TRACKSAW * (GetAnglePrecision(cr->values[0]) + GetAnglePrecision(cr->values[1]));
				costPrecision += PRECISION_TRACKSAW * (GetDimensionPrecision(cr->values[1]) + GetDimensionPrecision(cr->values[2]));
				costPrecision += PRECISION_TRACKSAW;
			}
			else if (r->get_type() == TYPE_REFD)
			{
				auto dr = dynamic_cast<Refd*>(r);
				costPrecision += PRECISION_DRILL * (GetDimensionPrecision(dr->values[1]) + GetDimensionPrecision(dr->values[2]));
				costPrecision += PRECISION_DRILL;
			}
		}

		//std::vector<double> evalResult = { materialCost, costPrecision, originalTime - result };
		std::vector<double> evalResult = { materialCost, costPrecision / static_cast<double>(nProgs), originalTime - result };
		// std::cerr << evalResult[0] << " " << evalResult[1] << " " << evalResult[2] << std::endl;

		evalResult[0] = (int)(evalResult[0] * 1000.0) / 1000.0;
		evalResult[1] = (int)(evalResult[1] * 1000.0) / 1000.0;
		evalResult[2] = (int)(evalResult[2] * 1000.0) / 1000.0;

		//evalResult[1] = 0.0;

		return evalResult;
	}

	// this function is an overloaded version for ``EncodedVec`` type
	std::vector<double> Evaluator::EvaluatePrograms(const EncodedVec& vec)
	{
		std::vector<ENode*> program;
		vec.GetAllPrograms(vec.root,program);

		return EvaluatePrograms(program);
	}
}