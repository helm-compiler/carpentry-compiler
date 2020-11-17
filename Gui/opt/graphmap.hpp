#pragma once
#include <iostream>
#include <string>
#include <parsers.hpp>
#include <map>
#include <config.hpp>
#include <fstream>
#include "tinyxml2.h"

namespace OptCompiler
{
	template <typename datum>
	using Vector2 = std::vector<std::vector<datum>>;

	template <typename datum>
	using Vector3 = std::vector<std::vector<std::vector<datum>>>;

	struct ProgPointer
	{
		ProgPointer() = default;

		ProgPointer(int m, int p) : moduleID(m), programID(p)
		{
		};
		int moduleID;
		int programID;
	};

	class GraphMap
	{
	public:
		static const std::string common_test_func_prefix;

		static GraphMap& get_instance()
		{
			static GraphMap instance;
			return instance;
		}

		GraphMap(const GraphMap& root) = delete;
		GraphMap& operator=(const GraphMap&) = delete;

		GraphMap()
		{
			logMerged = false;
			cont = std::unordered_map<SexpParser::Graph*, ProgPointer>();
		};

		void SetLogMerged(bool _i) { logMerged = _i; }

		bool GetLogMerged() const { return logMerged; }

		bool ImportOrders(const std::vector<int>& _o)
		{
			orders_.push_back(_o);
			return true;
		}

		//std::vector<std::unordered_map<int, int>> mergedDicts_;
		//std::vector<std::unordered_map<int, int>> StackDicts_;

		bool ImportMergedDicts(const std::vector<std::pair<int, int>>& mergedDict_)
		{
			mergedDicts_.emplace_back(mergedDict_);
			return true;
		}

		bool ImportStackDicts(const std::vector<std::pair<int, int>>& stackDict_)
		{
			stackDicts_.emplace_back(stackDict_);
			return true;
		}

		bool ImportScores(const std::vector<double>& _s)
		{
			scores_.push_back(_s);
			return true;
		}

		//std::vector<std::vector<std::tuple<int, int, int, Vector1i1>>> tuples_;

		void ImportTuples(const std::vector<std::tuple<string, int, int, Vector1i1>>& _ec)
		{
			tuples_.push_back(_ec);
		}

		void ImportTns(std::vector<std::tuple<string,int,int>>& tns_)
		{
			tns.emplace_back(tns_);
		}

		void ImportTds(VectorPI1& tds_)
		{
			tds.emplace_back(tds_);
		}

		void ImportInfos(const std::vector<std::tuple<int, int, int, int,int,int, double, double, double, std::string, std::string>>& _ec)
		{
			infos_.push_back(_ec);
		}

		void ImportTreeIds(const Vector1i1& _ec)
		{
			tree_ids_.push_back(_ec);
		}

		void ImportOutput(const std::vector<std::string>& _ec)
		{
			outputs_.push_back(_ec);
		}

		bool SaveToFile(const std::string& _f, const double &hv)
		{
			const auto declaration = R"(<?xml version="1.0" encoding="UTF-8"?>)";
			auto xmlDoc = std::make_shared<tinyxml2::XMLDocument>(_f.c_str());
			xmlDoc->Parse(declaration);

			const auto root = xmlDoc->NewElement("root");
			xmlDoc->InsertEndChild(root);

			if (orders_.size() != scores_.size())
			{
				std::cerr << "orders_.size() != scores_.size()" << std::endl;
				return false;
			}

			if (tuples_.size() != scores_.size())
			{
				std::cerr << "tuples_.size() != scores_.size()" << std::endl;
				return false;
			}

			//HyperVolume
			auto hyper_volume = xmlDoc->NewElement("HyperVolume");
			hyper_volume->SetAttribute("value", hv);
			root->InsertEndChild(hyper_volume);

			// Export scores and order
			for (auto i = 0; i < orders_.size(); ++i)
			{
				auto individual = xmlDoc->NewElement("Individual");
				individual->SetAttribute("ID", i);

				{
					auto score = xmlDoc->NewElement("Score");
					std::stringstream ss;
					for (int j = 0; j < scores_[i].size(); j++)
					{
						auto str = (j == scores_[i].size() - 1) ? "" : " ";
						ss << std::setprecision(10) << scores_[i][j] << str;
					}
					score->SetText(ss.str().c_str());
					individual->InsertEndChild(score);
				}

				{
					auto order = xmlDoc->NewElement("Order");
					std::stringstream ss;
					for (int j = 0; j < orders_[i].size(); j++)
					{
						auto str = (j == orders_[i].size() - 1) ? "" : " ";
						ss << std::setprecision(10) << orders_[i][j] << str;
					}
					order->SetText(ss.str().c_str());
					individual->InsertEndChild(order);
				}

				//std::vector<std::unordered_map<int, int>> mergedDicts_;
				//std::vector<std::unordered_map<int, int>> stackDicts_;

				{
					auto order = xmlDoc->NewElement("OrderStack");
					
					for (auto& o : stackDicts_[i])
					{
						std::stringstream ss;
						ss << o.first << " " << o.second;
						auto pair = xmlDoc->NewElement("Pair");
						pair->SetText(ss.str().c_str());
						order->InsertEndChild(pair);
					}
					individual->InsertEndChild(order);
				}

				{
					auto order = xmlDoc->NewElement("OrderMerge");
					for (auto& o : mergedDicts_[i])
					{
						std::stringstream ss;
						ss << o.first << " " << o.second;
						auto pair = xmlDoc->NewElement("Pair");
						pair->SetText(ss.str().c_str());
						order->InsertEndChild(pair);
					}
					individual->InsertEndChild(order);
				}

				{
					auto order = xmlDoc->NewElement("EncodedVec");
					for (auto& o : tns[i])
					{
						std::stringstream ss;
						ss << std::get<0>(o) << " " << std::get<1>(o) <<" "<< std::get<2>(o);
						auto tree_node = xmlDoc->NewElement("TreeNode");
						tree_node->SetAttribute("ID", std::get<0>(o).c_str());//e-class id
						tree_node->SetAttribute("H", std::get<1>(o));//height
						tree_node->SetAttribute("N", std::get<2>(o));//n
						//tree_node->SetText(ss.str().c_str());
						order->InsertEndChild(tree_node);

					
					}
					individual->InsertEndChild(order);
				}


				{
					auto order = xmlDoc->NewElement("TDS");
					for (auto& o : tds[i])
					{
						std::stringstream ss;
						ss << o.first << " " << o.second;
						auto pair = xmlDoc->NewElement("Pair");
						pair->SetText(ss.str().c_str());
						order->InsertEndChild(pair);
					}
					individual->InsertEndChild(order);
				}


				{
					auto tree = xmlDoc->NewElement("Tree");

					for (int j = 0; j < tuples_[i].size(); j++)
					{
						auto info = xmlDoc->NewElement("Info");
						info->SetAttribute("EClassID", std::string(std::get<0>(tuples_[i][j])).c_str());
						info->SetAttribute("ENodeID", std::string(std::to_string(std::get<1>(tuples_[i][j]))).c_str());
						info->SetAttribute("ArrID", std::string(std::to_string(std::get<2>(tuples_[i][j]))).c_str());
						info->SetAttribute("CutLinesID", std::string(Math::Functs::IntString(std::get<3>(tuples_[i][j]), false, " ")).c_str());
						tree->InsertEndChild(info);
					}
					individual->InsertEndChild(tree);
				}

				{
					auto Infos = xmlDoc->NewElement("Lines");

					//eclassID enodeID equivalenceID lineID arrID cutLineID referFace referEdge
					for (int j = 0; j < infos_[i].size(); j++)
					{
						auto info = xmlDoc->NewElement("Line");
						info->SetAttribute("EClassID", std::string(std::to_string(std::get<0>(infos_[i][j]))).c_str());
						info->SetAttribute("ENodeID", std::string(std::to_string(std::get<1>(infos_[i][j]))).c_str());
						info->SetAttribute("EquivalenceID", std::string(std::to_string(std::get<2>(infos_[i][j]))).c_str());
						info->SetAttribute("LineID", std::string(std::to_string(std::get<3>(infos_[i][j]))).c_str());
						info->SetAttribute("ArrID", std::string(std::to_string(std::get<4>(infos_[i][j]))).c_str());
						info->SetAttribute("CutLineID", std::string(std::to_string(std::get<5>(infos_[i][j]))).c_str());
						info->SetAttribute("RefAngle1", std::string(std::to_string(std::get<6>(infos_[i][j]))).c_str());
						info->SetAttribute("RefAngle2", std::string(std::to_string(std::get<7>(infos_[i][j]))).c_str());
						info->SetAttribute("RefLength", std::string(std::to_string(std::get<8>(infos_[i][j]))).c_str());
						info->SetAttribute("ReferFaceID", std::string(std::get<9>(infos_[i][j])).c_str());
						info->SetAttribute("ReferEdgeID", std::string(std::get<10>(infos_[i][j])).c_str());
						info->SetText(outputs_[i][j].c_str());
						Infos->InsertEndChild(info);
					}
					individual->InsertEndChild(Infos);
				}

				//{
				//	auto output = xmlDoc->NewElement("Output");
				//	output->SetText(outputs_[i].c_str());
				//	individual->InsertEndChild(output);
				//}
				root->InsertEndChild(individual);
			}


			xmlDoc->SaveFile(_f.c_str());
			return true;
		}

		void ClearResults()
		{
			orders_.clear();
			scores_.clear();
			cont.clear();
			tuples_.clear();
			tree_ids_.clear();
			outputs_.clear();
			mergedDicts_.clear();
			stackDicts_.clear();
			infos_.clear();
			tns.clear();
			tds.clear();
		}

		std::unordered_map<SexpParser::Graph*, ProgPointer> cont;

		bool logMerged;

		Vector2<int> orders_;
		Vector2<double> scores_;
		std::vector<std::vector<std::tuple<string, int, int, Vector1i1>>> tuples_;
		Vector1i2 tree_ids_;
		std::vector<std::vector<std::string>> outputs_;

		std::vector < std::vector<std::tuple<string, int, int>>>  tns;
		VectorPI2 tds;

		std::vector<std::vector<std::pair<int, int>>> mergedDicts_;
		std::vector<std::vector<std::pair<int, int>>> stackDicts_;

		//eclassID enodeID equivalenceID lineID arrID cutLineID
		//angle1,angle2, offset
		//referFace referEdge
		std::vector < std::vector < std::tuple<int, int, int, int, int, int,double,double,double, std::string, std::string>>> infos_;

	};
}