#include <optimizer.hpp>


#include "PreCompiled.h"
#include <QProgressDialog>
#include <sstream>
#include "TaskFeaturePick.h"
#include <QInputDialog>
#include <qmessagebox.h>



namespace OptCompiler
{
	namespace Heuristic
	{
		void UniformTree(std::shared_ptr<TreeNode> node, Atom type)
		{
			if (!node) return;

			if (node->GetN() != -1 && node->ec->ens[node->GetN()]->IsENodeProg())
			{
				// If its offsprings are programs
				std::vector<int> vecPossibleProgs;

				for (auto i = 0; i < node->ec->ens.size(); ++i)
				{
					auto& progEc = node->ec->ens[i];
					if (progEc->IsENodeProg() && progEc->ast.GetProcess()->get_type() == type)
					{
						vecPossibleProgs.push_back(i);
					}
				}

				if (!vecPossibleProgs.empty())
				{
					int rn = Config::get_instance().GetRandom(0, vecPossibleProgs.size() - 1);
					node->SetN(vecPossibleProgs[rn]);
				}
				else
				{
					node->SetN(0);
				}
			}

			for (auto& n : node->subNodes)
			{
				UniformTree(n, type);
			}
		}
	}

	int OptimizingCompiler::execute(const std::string& folder, const std::vector<std::vector<std::tuple<string, int, int>>>& encoded_vecs, const VectorPI2& encoded_tds)
	{
		auto OutputArrWLC = [&]()
		{
			//livefile
			std::ofstream file(Config::get_instance().LIVEFILE_FILE());

			auto& arr_wlcs = Config::get_instance().ecProgArrWLCs;
			for (auto& arr_wlc : arr_wlcs)
			{
				std::cerr << "Arr_WLC: " << arr_wlc.index << " wlc: " << arr_wlc.wlc << " arr: " << arr_wlc.arr << "  iters(" << arr_wlc.iters.size() << "): ";
				for (auto& iter : arr_wlc.iters)
					std::cerr << iter << " ";
				std::cerr << "\n";

				//file << arr_wlc.wlc << " " << arr_wlc.arr << " " << arr_wlc.iters.size() << " ";
				for (auto& iter : arr_wlc.iters)
				{
					//file << iter << " ";
					file << iter << " " << arr_wlc.arr << std::endl;
				}
				//file << "\n";
			}

			file.close();
		};

		GraphMap::get_instance().ClearResults();
		GraphMap::get_instance().SetLogMerged(false);//vis output results
		Config::get_instance().Clear();
		Config::get_instance().SetConfigs(folder);

		std::cerr << Config::get_instance().SUBMODULE_FILE() << std::endl;
		std::ifstream fin(Config::get_instance().SUBMODULE_FILE());
		if (!fin)
		{
			std::cerr << "can not open this file" << std::endl;
			system("pause");
		}

		const auto popSize = Config::get_instance().GetPopSize();
		const auto evoGens = Config::get_instance().GetEvoGens();
		const auto nIterations = Config::get_instance().GetIterations();
		const auto nTerminal = Config::get_instance().GetNTerminal();

		// 0 - Parse external files, and construct e-graph accordingly
		auto eg = new EGraph();
		std::unordered_map<int, EGraphParser::SubModule> sexpProgs;
		std::unordered_set<int> emptyProgs;
		EGraphParser::ParseXML(Config::get_instance().SUBMODULE_FILE(), sexpProgs, emptyProgs);
		ConstructENodes(eg, sexpProgs);
		auto root = EGraphParser::ParseEGraph(Config::get_instance().EGRAPHFILE_PATH(), eg);

		eg->ENodeNB();
		std::cerr << "Stock ENode NB: " << eg->stock_e_nodes << std::endl;
		std::cerr << "Union ENode NB: " << eg->union_e_nodes << std::endl;

		if (eg->stock_e_nodes == 0)
		{
			std::cerr << "There is no stock nodes in the current e-graph." << std::endl;
			std::cerr << "if (eg->stock_e_nodes == 0)" << std::endl;
			system("pause");
			return EXIT_FAILURE;
		}

		auto startTime = std::chrono::steady_clock::now();

		// 1 - Instantiate a pagmo problem constructing it from a UDP
		pagmo::problem prob{ pagmo::dtlz(1, 5, 3) };

		// 2 - Instantiate a pagmo algorithm
		pagmo::algorithm algo{ nsga2_egraph(evoGens) };
		//algo.set_verbosity(50);

		// 3 - Instantiate a population
		auto rootEN = new ENode();
		rootEN->ecs.push_back(root);

		// 4 - Start optimization
		double refPoint[3] = { -DBL_MAX, 6.25 * 1.2, -DBL_MAX };

		pagmo::archipelago archi{ Config::get_instance().GetNIslands(), algo, prob, popSize };

		std::ofstream initPareto(Config::get_instance().INIT_FILE());
		for (auto& island : archi)
		{
			pagmo::population pop{ prob, popSize, Config::get_instance().GetRandSeed() };

			//compute 
			pop.random_decision_tree(eg, rootEN, popSize, encoded_vecs, encoded_tds);;

			pop.ExportPareto(initPareto, pop);
			island.set_population(pop);

			auto m_pf = std::get<0>(fast_non_dominated_sorting(pop.get_f()))[0];

			auto& m_x = pop.get_x();
			auto& m_f = pop.get_f();

			for (population::size_type i = 0; i < m_x.size(); i++)
			{
				std::vector<ENode*> progs;
				m_x[i].GetAllPrograms(m_x[i].root, progs);

				if (std::find(m_pf.begin(), m_pf.end(), i) != m_pf.end())
				{
					for (auto& prog : progs)
						if (prog->ecs.empty())
							Config::get_instance().Push_Iter_ArrWLC(prog->arr, prog->wlc, 0);
				}

				auto curTime = 10.0 * progs.size();
				if (curTime > refPoint[2]) refPoint[2] = curTime;
			}

			for (auto& p : m_f)
			{
				auto curMat = 4.0 * p[0];
				if (curMat > refPoint[0])
				{
					refPoint[0] = curMat;
				}
			}
		}
		initPareto.close();

		OutputArrWLC();
		//return EXIT_SUCCESS;

		if (false)
		{

			auto iter = archi.begin();
			auto islandPop = iter->get_population();
			for (auto i = 0; i < islandPop.get_x().size(); ++i)
			{
				auto& p = islandPop.get_x()[i];
				Heuristic::UniformTree(p.root, TYPE_TRACKSAW);
				islandPop.set_xf(i, p, Evaluator::EvaluatePrograms(p));
			}
			iter->set_population(islandPop);

			++iter;

			islandPop = iter->get_population();
			for (auto i = 0; i < islandPop.get_x().size(); ++i)
			{
				auto& p = islandPop.get_x()[i];
				Heuristic::UniformTree(p.root, TYPE_CHOPSAW);
				islandPop.set_xf(i, p, Evaluator::EvaluatePrograms(p));
			}
			iter->set_population(islandPop);

			++iter;

			islandPop = iter->get_population();
			for (auto i = 0; i < islandPop.get_x().size(); ++i)
			{
				auto& p = islandPop.get_x()[i];
				Heuristic::UniformTree(p.root, TYPE_BANDSAW);
				islandPop.set_xf(i, p, Evaluator::EvaluatePrograms(p));
			}
			iter->set_population(islandPop);

			++iter;

			islandPop = iter->get_population();
			for (auto i = 0; i < islandPop.get_x().size(); ++i)
			{
				auto& p = islandPop.get_x()[i];
				Heuristic::UniformTree(p.root, TYPE_JIGSAW);
				islandPop.set_xf(i, p, Evaluator::EvaluatePrograms(p));
			}
			iter->set_population(islandPop);
		}

		std::cerr << "Refpoint = " << refPoint[0] << " " << refPoint[1] << " " << refPoint[2] << std::endl;

		auto& rp = CompilerConfig::Instance().GetRefPoint();

		refPoint[0] = rp[0];
		refPoint[1] = rp[1];
		refPoint[2] = rp[2];
		Config::get_instance().SetRefpoint(refPoint);

		//archi.evolve(nIterations);

		auto pop_ = archi[0].get_algorithm().evolve(archi[0].get_population());

		// 5 - Wait until finished
		//archi.wait_check();

		auto endTime = std::chrono::steady_clock::now();

		// 6 - Print the fitness of the best solution in each island
		if (true)
		{
			/*for (const auto& isl : archi) { print(isl.get_population());}*/
			print(pop_);
		}

		GraphMap::get_instance().ClearResults();
		GraphMap::get_instance().SetLogMerged(true);//vis output results

		for (int i = 0; i < pop_.get_x().size(); i++)
		{
			std::cerr << "#" << i << std::endl;
			auto score = Evaluator::EvaluatePrograms(pop_.get_x()[i]);

			std::vector<std::tuple<string, int, int>> tn;
			VectorPI1 td_edges;
			CollectTreeNodes(pop_.get_x()[i].root, tn, td_edges);
			GraphMap::get_instance().ImportTns(tn);
			GraphMap::get_instance().ImportTds(td_edges);
			GraphMap::get_instance().ImportScores(score);

			std::vector<std::tuple<string, int, int, Vector1i1>> tuples;
			EncodeTreeLeaf1(pop_.get_x()[i].root, tuples);
			GraphMap::get_instance().ImportTuples(tuples);

			tuples.clear();
			tn.clear();
			td_edges.clear();
			score.clear();
		}

		auto hyperVolume = Config::get_instance().PostprocessHypervolume();
		auto paretoFrontsNB = Config::get_instance().PostprocessParetoFrontsNB();

		GraphMap::get_instance().SaveToFile(Config::get_instance().RESULT_FILE(), hyperVolume[paretoFrontsNB.size() - 1]);

		auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();
		std::cerr << "Time elapsed in milliseconds = " << time_ms << " ms" << std::endl;

		std::ofstream timeStat(Config::get_instance().FOLDER_NAME() + std::to_string(time_ms));
		timeStat << std::to_string(time_ms);
		timeStat.close();


		auto uniqueParetoFrontsNB = Config::get_instance().PostprocessUniqueParetoFrontsNB();
		auto popDiffNb = Config::get_instance().PostprocessPopDiffNb();
		auto arrNB = Config::get_instance().PostprocessArrangesNB();

		std::ofstream hvfile(Config::get_instance().HYPERVOLUME_FILE());

		hvfile << "Ref: " << Config::get_instance().GetRefpoint()[0] << " " << Config::get_instance().GetRefpoint()[1] << " " << Config::get_instance().GetRefpoint()[2] << std::endl;

		for (int i = 0; i < paretoFrontsNB.size(); i++)
		{
			hvfile << "Iteration: " << i + 1 << " HV: " << hyperVolume[i] << " Unique ParetoFronts: " << uniqueParetoFrontsNB[i] << " ParetoFronts: " << paretoFrontsNB[i] << " PopDiffNb: " << popDiffNb[i] << " Arranges NB: " << arrNB[i] << std::endl;
		}
		hvfile.close();

		OutputArrWLC();

		return EXIT_SUCCESS;
	}

	FabInfo OptimizingCompiler::ReadParetoFronts(const std::string& path)
	{
		FabInfo fab_info;

		try
		{
			tinyxml2::XMLDocument doc;
			doc.LoadFile(path.c_str());

			if (doc.Error())
			{
				std::cerr << "Load xml error: " << path << std::endl;
				std::cerr << doc.ErrorStr() << std::endl;
				system("pause");
			}


			auto root = doc.FirstChildElement("root");

			fab_info.hv = atof(root->FirstChildElement("HyperVolume")->Attribute("value"));

			for (auto e = root->FirstChildElement("Individual"); e != nullptr; e = e->NextSiblingElement("Individual"))
			{
				FabInfo::Fab fab;
				fab.vis_path = path;

				fab.score = Math::Functs::SplitD(e->FirstChildElement("Score")->GetText(), ' ');
			
				if(fab.score[1]>0) fab.order = Math::Functs::SplitI(e->FirstChildElement("Order")->GetText(), ' ');

				fab.score[1] = fab.score[1] < 0.0 ? 0.0 : fab.score[1];//haisen
				auto tree = e->FirstChildElement("Tree");
				if (tree->GetLineNum() >= 1)
				{
					for (auto info = tree->FirstChildElement("Info"); info != nullptr; info = info->NextSiblingElement("Info"))
					{
						fab.EClassIDs.emplace_back(atoi(info->Attribute("EClassID")));
						fab.ArrIDs.emplace_back(atoi(info->Attribute("ArrID")));
						fab.ENodeIDs.emplace_back(atoi(info->Attribute("ENodeID")));
					}
				}


				auto lines = e->FirstChildElement("Lines");
				if (lines->GetLineNum() >= 1)
				{
					for (auto line = lines->FirstChildElement("Line"); line != nullptr; line = line->NextSiblingElement("Line"))
					{
						fab.referFE.emplace_back(
							line->Attribute("ReferFaceID"),
							line->Attribute("ReferEdgeID"),
							atof(line->Attribute("RefAngle1")),
							atof(line->Attribute("RefAngle2")),
							atof(line->Attribute("RefLength")));
					}
				}


				auto stack = e->FirstChildElement("OrderStack");
				if (stack->GetLineNum() >= 1)
				{
					for (auto pair = stack->FirstChildElement("Pair"); pair != nullptr; pair = pair->NextSiblingElement("Pair"))
					{
						//order_stacks
						auto v = Math::Functs::SplitI(pair->GetText(), ' ');
						fab.order_stack.emplace_back(v[0], v[1]);
					}
				}

				auto merge = e->FirstChildElement("OrderMerge");
				if (merge->GetLineNum() >= 1)
				{
					for (auto pair = merge->FirstChildElement("Pair"); pair != nullptr; pair = pair->NextSiblingElement("Pair"))
					{
						auto v = Math::Functs::SplitI(pair->GetText(), ' ');
						fab.order_merge.emplace_back(v[0], v[1]);
					}
				}


				auto encoded_vec = e->FirstChildElement("EncodedVec");
				if (encoded_vec->GetLineNum() >= 1)
				{
					for (auto tree_node = encoded_vec->FirstChildElement("TreeNode"); tree_node != nullptr; tree_node = tree_node->NextSiblingElement("TreeNode"))
					{
						std::tuple<string, int, int > t;
						std::get<0>(t) = tree_node->Attribute("ID");
						std::get<1>(t) = atoi(tree_node->Attribute("H"));
						std::get<2>(t) = atoi(tree_node->Attribute("N"));
						fab.encoded_vec.emplace_back(t);
					}
				}

				auto encoded_td = e->FirstChildElement("TDS");
				if (encoded_td->GetLineNum() >= 1)
				{
					for (auto pair = encoded_td->FirstChildElement("Pair"); pair != nullptr; pair = pair->NextSiblingElement("Pair"))
					{
						auto v = Math::Functs::SplitI(pair->GetText(), ' ');
						fab.encoded_td.emplace_back(v[0], v[1]);
					}
				}

				fab_info.fabs.emplace_back(fab);
			}
		}
		catch (...)
		{
			std::cerr << "auto ReadParetoFronts = [](const std::string& path)" << std::endl;
			system("pause");
		}

		return fab_info;
	}

	double OptimizingCompiler::ComputeHV(const std::vector<std::vector<double>>& hyper_volume_scores)
	{
		auto& rp = CompilerConfig::Instance().GetRefPoint();

		auto m_method = pagmo::hv3d().clone();
		auto hv_obj = pagmo::hypervolume(hyper_volume_scores, true);
		double hypvol = hv_obj.compute({ rp[0], rp[1], rp[2] }, *m_method);
		return hypvol;
	}

	std::vector<FabInfo::Fab> OptimizingCompiler::ParetoFronts_Score(const std::vector<FabInfo::Fab>& fabs, Vector1i1& fronts_index, int add_fronts)
	{
		if (fabs.size() <= 1) return fabs;

		Vector1d2 scores;
		for (auto& fab : fabs) scores.emplace_back(fab.score);

		auto ndf = std::get<0>(fast_non_dominated_sorting(scores)); // non dominated fronts [[0,3,2],[1,5,6],[4],...]
		for (auto front : ndf.front())
			fronts_index.emplace_back(front);

		std::vector<FabInfo::Fab> fabs_;
		for (auto& front_index : fronts_index)
			fabs_.emplace_back(fabs[front_index]);

		if (ndf.size()>1)
		{
			int nb = fabs_.size();
			for (int i = 1; i < ndf.size()&& fabs_.size() - nb<add_fronts; i++)
			{
				for (int j = 0; j < ndf[i].size() && fabs_.size() - nb < add_fronts; j++)
				{
					fabs_.emplace_back(fabs[ndf[i][j]]);
				}
			}
		}

		return fabs_;
		
	}

	int OptimizingCompiler::ParetoFronts_Score(const std::vector<FabInfo::Fab>& fabs, const FabInfo::Fab& fab)
	{
		Vector1d2 scores;
		for (auto& fab_ : fabs) scores.emplace_back(fab_.score);

		return ParetoFronts_Score(scores,fab.score);
	}

	int OptimizingCompiler::ParetoFronts_Score(const Vector1d2& fronts, const Vector1d1& score)
	{
		for (int j = 0; j < fronts.size(); j++)
		{
			auto front = fronts[j];

			if (front[0] == score[0] && front[1] == score[1] && front[2] == score[2])
			{
	
			}
			else
			{
				if (front[0] <= score[0] && front[1] <= score[1] && front[2] <= score[2])
				{
					return j;
				}
			}
		}
		return -1;
	}

	Vector1d2 OptimizingCompiler::ParetoFronts_Score(const Vector1d2& scores, Vector1i1& fronts_index)
	{
		std::vector<std::vector<double>> fronts;

		for (int i = 0; i < scores.size(); i++)
		{
			auto score = scores[i];

			bool b = true;
			std::vector<std::vector<double>> fronts_;
			std::vector<int> fronts_index_;

			for (int j = 0; j < fronts.size(); j++)
			{
				auto front = fronts[j];

				if (front[0] == score[0] && front[1] == score[1] && front[2] == score[2])
				{
					fronts_.emplace_back(front);
					fronts_index_.emplace_back(fronts_index[j]);
				}
				else
				{
					if (front[0] <= score[0] && front[1] <= score[1] && front[2] <= score[2])
					{
						b = false;
						fronts_.emplace_back(front);
						fronts_index_.emplace_back(fronts_index[j]);
					}
					else
					{
						if (front[0] >= score[0] && front[1] >= score[1] && front[2] >= score[2])
						{
							b = true;
						}
						else
						{
							fronts_.emplace_back(front);
							fronts_index_.emplace_back(fronts_index[j]);
						}
					}
				}
			}

			if (b)
			{
				fronts_.emplace_back(score);
				fronts_index_.emplace_back(i);
			}
			fronts = fronts_;
			fronts_index = fronts_index_;
		}

		return fronts;
	}

	const Vector1i2 FabInfo::GetArrIDs() const
	{
		Vector1i2 arrIDs;
		for (auto& fab : fabs)
			arrIDs.emplace_back(fab.ArrIDs);
		return arrIDs;
	}

	const Vector1d2 FabInfo::GetScores() const
	{
		Vector1d2 scores;
		for (auto& fab : fabs)
			scores.emplace_back(fab.score);
		return scores;
	}

	const VectorPI2 FabInfo::GetEncodedTds() const
	{
		VectorPI2 encoded_tds;
		for (auto& fab : fabs)
			encoded_tds.emplace_back(fab.encoded_td);
		return encoded_tds;
	}

	const VectorPI2 FabInfo::GetEgraphEncodes() const
	{
		VectorPI2 egraph_encodes;
		for (auto& fab : fabs)
			egraph_encodes.emplace_back(fab.joint_encode);
		return egraph_encodes;
	}

	const std::vector<std::vector<std::tuple<Base::string, int, int>>> FabInfo::GetEncodedVecs() const
	{
		std::vector<std::vector<std::tuple<string, int, int>>> encoded_vecs;
		for (auto& fab : fabs)
			encoded_vecs.emplace_back(fab.encoded_vec);
		return encoded_vecs;
	}

	const void FabInfo::OutputLog() const
	{
		for (int i = 0; i < fabs.size(); i++)
		{
			auto& fab = fabs[i];
			std::cerr << "Index: " << i << " score: "
				<< Math::Functs::IntString(fab.score, false, ",")
				<< " ArrIDs: " << Math::Functs::IntString(fab.ArrIDs, false, ",")
				<<" Joint_encode: "<< Math::Functs::IntString(fab.joint_encode, false, ",",";")
				<< " encoded_td: " << Math::Functs::IntString(fab.encoded_td, false, ",", ";") << std::endl;
		}
	}

	void FabInfo::Clear()
	{
		for (auto& fab : fabs)fab.Clear();
		std::vector<Fab>().swap(fabs);
	}

}
