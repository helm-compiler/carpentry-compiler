#pragma once
#include <iostream>
#include <nsgaii.hpp>
#include <evaluator.hpp>
#include <parsers.hpp>
#include <egraph.hpp>
#include <boost/filesystem.hpp>
#include <pagmo/problems/dtlz.hpp>
#include <pagmo/archipelago.hpp>
#include <random>
#include <config.hpp>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

namespace OptCompiler
{
	struct FabInfo
	{
	public:
		struct Fab
		{
			Vector1d1 score;//material/precision/time
			Vector1i1 EClassIDs;
			Vector1i1 ArrIDs;
			Vector1i1 ENodeIDs;

			string vis_path;

			std::vector<std::tuple<string, int, int>> encoded_vec;//???
			VectorPI1 encoded_td;//???

			Vector1i1 order;//cutting order
			VectorPI1 order_stack;//cutting order stacks
			VectorPI1 order_merge;//cutting order merges
			std::vector<std::tuple<std::string, std::string, double, double, double>> referFE;//cutting references 

			VectorPI1 joint_encode;

			Fab() {}
			Fab(const Vector1d1& score_, const Vector1i1& ArrIDs_):score(score_), ArrIDs(ArrIDs_){};

			void Clear()
			{
				Vector1d1().swap(score);
				Vector1i1().swap(EClassIDs);
				Vector1i1().swap(ArrIDs);
				Vector1i1().swap(ENodeIDs);

				std::vector<std::tuple<string, int, int>>().swap(encoded_vec);
				VectorPI1().swap(encoded_td);

				Vector1i1().swap(order);
				VectorPI1().swap(order_stack);
				VectorPI1().swap(order_merge);
				std::vector<std::tuple<std::string, std::string, double, double, double>>().swap(referFE);

				VectorPI1().swap(joint_encode);
			}
		};

		double hv;//hypervolume
		std::vector<Fab> fabs;

		int nsga_layer = -1;//ga nsga layer
		double nsga_crowding_distance = 0.0;//ga nsga crowding distance

		const Vector1i2 GetArrIDs() const;
		const Vector1d2 GetScores() const;
		const VectorPI2 GetEncodedTds() const;
		const VectorPI2 GetEgraphEncodes() const;
		const std::vector<std::vector<std::tuple<string, int, int>>> GetEncodedVecs() const;

		const void OutputLog() const;

		void Clear();
	};



	class OptimizingCompiler {
	public:
		OptimizingCompiler() = default;

		~OptimizingCompiler() = default;

		int execute(const std::string& folder, const std::vector<std::vector<std::tuple<string, int, int>>>& encoded_vecs = std::vector<std::vector<std::tuple<string, int, int>>>(), const VectorPI2& encoded_tds = VectorPI2());

		static FabInfo ReadParetoFronts(const std::string& path);

		static Vector1d2 ParetoFronts_Score(const Vector1d2& scores, Vector1i1& fronts_index = Vector1i1());

		static int ParetoFronts_Score(const Vector1d2& scores, const Vector1d1& score);
		static int ParetoFronts_Score(const std::vector<FabInfo::Fab>& fabs, const FabInfo::Fab& fab);

		static std::vector<FabInfo::Fab> ParetoFronts_Score(const std::vector<FabInfo::Fab>& fabs, Vector1i1& fronts_index = Vector1i1(), int add_fronts = 0);

		static double ComputeHV(const std::vector<std::vector<double>>& hyper_volume_scores);


	private:

		void ConstructENodes(EGraph* eg, std::unordered_map<int, EGraphParser::SubModule>& sexpProgs)
		{
			auto& cutLineId_cMaps = Config::get_instance().ecCutLineIdmaps;
			auto& enode_cMaps = Config::get_instance().ecENodeIdmaps;
			auto& arr_cMaps = Config::get_instance().ecArrIdmaps;
			auto& wlc_cMaps = Config::get_instance().ecWLCIdmaps;

			int eId = 10000;

			int eNodeId = 0;

			// Load all subprograms
			for (const auto& subMd : sexpProgs)
			{
				// loop: submodule
				auto eclassID = subMd.first;
				auto eclassID_str = Math::Functs::IntString(eclassID);
				auto ec = new EClass();//EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass EClass 
				ec->SetId(eclassID_str);

				auto& cutLineId_cMap = cutLineId_cMaps[eclassID_str];
				auto& enode_cMap = enode_cMaps[eclassID_str];
				auto& arr_cMap = arr_cMaps[eclassID_str];
				auto& wlc_cMap = wlc_cMaps[eclassID_str];

				for (auto i = 0; i < subMd.second.size(); ++i)
				{
					// loop: program
					ENode* enToProg = new ENode();//ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode 
					ec->AddEn(enToProg);

					auto enToProgId_str = eclassID_str + "_" + Math::Functs::IntString(i);

					enToProg->SetId(enToProgId_str);
					enToProg->eclassID = eclassID;
					enToProg->cutLinesIds = cutLineId_cMap[i];
					enToProg->enodeID = enode_cMap[i];
					enToProg->arr = arr_cMap[i];
					enToProg->wlc = wlc_cMap[i];

					auto& prog = subMd.second[i];
					// for each equivalence
					for(int j=0;j<prog.size();j++)
					{
						auto eq = prog[j];

						auto ecEquiv = new EClass();//EClass EClass EClass EClass EClass EClass EClass EClass EClass
						auto ecEquivId_str = enToProgId_str + "_" + Math::Functs::IntString(j);
						ecEquiv->SetId(ecEquivId_str);
						enToProg->ecs.push_back(ecEquiv);
						// loop: equivalence

						if (eq.size() != 1)
						{
							int dsad = 0;
						}

						for(int k=0;k<eq.size();k++)
						//for (auto& ln : eq)
						{
							auto ln = eq[k];
							// loop: line
							auto en = new ENode(ln.first, ln.second);//ENode ENode ENode ENode ENode ENode ENode
							auto enId_str = ecEquivId_str + "_" + Math::Functs::IntString(k);
							en->SetId(enId_str);
							en->eclassID = eclassID;
							en->cutLinesIds = cutLineId_cMap[i];
							en->enodeID = enode_cMap[i];
							en->arr = arr_cMap[i];
							en->wlc = wlc_cMap[i];
							en->cutLineID = cutLineId_cMap[i][j];
							en->equivalenceID = j;
							en->lineID = k;
							GraphMap::get_instance().cont[&(en->ast.g)] = ProgPointer(eclassID, i);
							ecEquiv->AddEn(en);
						}

						eg->AddEc(ecEquiv);
					}
				}
				eg->AddEc(ec);
			}
		}

		void ComputeReferencePoint(pagmo::population& pop)
		{
			double refPoint[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };

			for (auto& p : pop.get_f())
			{
				for (auto i = 0; i < 3; ++i)
				{
					if (p[i] > refPoint[i])
					{
						refPoint[i] = p[i];
					}
				}
			}
			std::cerr << "Ref Point = ";
			for (auto i = 0; i < 3; ++i)
				std::cerr << refPoint[i] << " ";
			std::cerr << std::endl;
		}

		void BuildTree(std::shared_ptr<TreeNode> node, std::deque<int>& ids)
		{
			if (!node) return;

			if (node->GetN() != -1)
			{
				node->SetN(ids.front());
				ids.erase(ids.begin());
			}

			for (auto& n : node->subNodes)
			{
				BuildTree(n, ids);
			}
		}

		void mutate(std::shared_ptr<TreeNode> node)
		{
			if (node->GetN() == -1)
			{
				for (auto& subNode : node->subNodes)
				{
					mutate(subNode);
				}
			}
			else
			{
				if (Config::get_instance().GetProb() < 0.02)
				{
					auto hb = static_cast<int>(node->ec->ens.size()) - 1;
					node->SetN(Config::get_instance().GetRandom(0, hb));
				}

				if (!node->IsProgram())
				{
					for (auto& subNode : node->subNodes)
					{
						mutate(subNode);
					}
				}
			}
		}

		void mutate(EncodedVec& child)
		{
			mutate(child.root);
		}

		std::string EncodeTreeLeafMain(const std::shared_ptr<TreeNode>& node)
		{
			if (!node) return "#";

			std::string str;

			if (node->val != "-1")
				str = node->val + "-" + std::to_string(node->GetN());

			if (!node->subNodes.empty())
			{
				for (auto& n : node->subNodes)
				{
					str += "," + EncodeTreeLeafMain(n);
				}
			}

			return str;
		}

		void ParseAndEvaluate(const std::string& file)
		{
			using namespace tinyxml2;
			std::string fsn = file;
			std::ofstream fileScore(fsn.append(".txt"));

			XMLDocument doc;
			doc.LoadFile(file.c_str());
			XMLElement* root = doc.FirstChildElement("root");
			for (XMLElement* e = root->FirstChildElement("Program"); e != nullptr; e = e->NextSiblingElement("Program"))
			{
				std::vector<std::string> oneProg;
				std::vector<ENode*> fullProg;
				for (auto ln = e->FirstChildElement("sexp"); ln != nullptr; ln = ln->NextSiblingElement("sexp"))
				{
					//std::cerr << ln->GetText() << std::endl;
					//oneProg.emplace_back(ln->GetText());
					std::string a1 = ln->GetText();
					std::string a2 = ln->GetText();
					auto en = new ENode(a1, a2);
					fullProg.push_back(en);
				}

				auto v = Evaluator::EvaluatePrograms(fullProg);
				std::cerr << v[0] << " " << v[1] << " " << v[2] << std::endl;
				fileScore << v[0] << " " << v[1] << " " << v[2] << std::endl;
				for (auto n : fullProg)
					delete n;
			}
			fileScore.close();
		}
	};
}