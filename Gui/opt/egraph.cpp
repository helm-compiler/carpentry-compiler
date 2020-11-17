#include "egraph.hpp"
#include <evaluator.hpp>
#include <sexp/parser.hpp>
#include <sexp/util.hpp>

namespace OptCompiler
{
	// deprecated
	std::string EncodeTree(const std::shared_ptr<TreeNode>& node)
	{
		if (!node) return "#";

		std::string str;
		if (node->val != "-1")
			str = node->val + "-" + std::to_string(node->GetN());

		for (auto& n : node->subNodes)
		{
			str += "," + EncodeTree(n);
		}

		return str;
	}

	std::string EncodeTreeLeafD(const std::shared_ptr<TreeNode>& node)
	{
		if (!node) return "#";

		std::string str;

		if (node->val != "-1")
			str = node->val + "-" + std::to_string(node->GetN());

		if (!node->subNodes.empty())
		{
			for (auto& n : node->subNodes)
			{
				str += "," + EncodeTreeLeafD(n);
			}
		}

		return str;
	}


	int CollectTreeNodes(const std::shared_ptr<TreeNode>& v, std::vector<std::tuple<string, int, int>>& tn, VectorPI1& edges)
	{
		int tn_id = tn.size();
		std::tuple<string, int, int> t;
		std::get<0>(t) = v->val;
		std::get<1>(t) = v->height;
		std::get<2>(t) = v->GetN();

		tn.push_back(t);

		for (auto& t : v->subNodes)
		{
			int t_id = CollectTreeNodes(t, tn, edges);
			edges.emplace_back(std::pair<int, int>(tn_id, t_id));
		}
		return tn_id;
	}

	void EncodeTreeLeaf1(const std::shared_ptr<TreeNode>& node,  std::vector<std::tuple<string, int, int, Vector1i1>> &tuples)
	{
		if (!node) return;

		if (node->val != "-1" && !node->IsProgram())
		{
			auto& cutLine_cMaps = Config::get_instance().ecCutLineIdmaps;
			auto& enode_cMaps = Config::get_instance().ecENodeIdmaps;
			auto& arr_cMaps = Config::get_instance().ecArrIdmaps;
			auto& wlc_cMaps = Config::get_instance().ecWLCIdmaps;
			if (cutLine_cMaps.find(node->val) != cutLine_cMaps.end() &&
				enode_cMaps.find(node->val) != enode_cMaps.end() &&
				arr_cMaps.find(node->val) != arr_cMaps.end() &&
				wlc_cMaps.find(node->val) != wlc_cMaps.end())
			{
				auto& cutline_cMap = cutLine_cMaps[node->val];
				auto& enode_cMap = enode_cMaps[node->val];
				auto& arr_cMap = arr_cMaps[node->val];
				auto& wlc_cMap = wlc_cMaps[node->val];
				if (cutline_cMap.find(node->GetN()) != cutline_cMap.end() &&
					enode_cMap.find(node->GetN()) != enode_cMap.end() &&
					arr_cMap.find(node->GetN()) != arr_cMap.end() &&
					wlc_cMap.find(node->GetN()) != wlc_cMap.end())
				{
					std::tuple<string, int, int, Vector1i1> eclass_enode_arrange_cutlines;
					std::get<0>(eclass_enode_arrange_cutlines) = node->val;
					std::get<1>(eclass_enode_arrange_cutlines) = enode_cMap[node->GetN()];
					std::get<2>(eclass_enode_arrange_cutlines) = arr_cMap[node->GetN()];
					std::get<3>(eclass_enode_arrange_cutlines) = cutline_cMap[node->GetN()];
					tuples.emplace_back(eclass_enode_arrange_cutlines);
				}
			}
			else
			{
				std::cerr << "error " << std::endl;
				exit(-1);
			}
		}

		if (!node->subNodes.empty())
		{
			for (auto& n : node->subNodes)
			{
				EncodeTreeLeaf1(n, tuples);
			}
		}
	}

	std::string EncodeTreeLeaf(const std::shared_ptr<TreeNode>& node)
	{
		if (!node) return "#";

		std::string str;

		if (node->val != "-1" && !node->IsProgram())
		{
			auto& cutLine_cMaps = Config::get_instance().ecCutLineIdmaps;
			auto& enode_cMaps = Config::get_instance().ecENodeIdmaps;
			auto& arr_cMaps = Config::get_instance().ecArrIdmaps;
			auto& wlc_cMaps = Config::get_instance().ecWLCIdmaps;
			if (cutLine_cMaps.find(node->val) != cutLine_cMaps.end() &&
				enode_cMaps.find(node->val) != enode_cMaps.end() && 
				arr_cMaps.find(node->val) != arr_cMaps.end() && 
				wlc_cMaps.find(node->val) != wlc_cMaps.end())
			{
				auto& cutline_cMap = cutLine_cMaps[node->val];
				auto& enode_cMap = enode_cMaps[node->val];
				auto& arr_cMap = arr_cMaps[node->val];
				auto& wlc_cMap = wlc_cMaps[node->val];
				if (cutline_cMap.find(node->GetN()) != cutline_cMap.end() &&
					enode_cMap.find(node->GetN()) != enode_cMap.end()&&
					arr_cMap.find(node->GetN()) != arr_cMap.end()&&
					wlc_cMap.find(node->GetN()) != wlc_cMap.end())
				{
					//str = std::to_string(node->val) + "-" +
					//	std::to_string(id_cMap[node->GetN()])+ "-"+
					//	std::to_string(arr_cMap[node->GetN()])+"-"+
					//	std::to_string(wlc_cMap[node->GetN()]) + ",";

					str = node->val + " " + //eclass id
						  std::to_string(arr_cMap[node->GetN()])+" "+ //arrange id
						  std::to_string(enode_cMap[node->GetN()]) + " ";//program id in each eclass

				}
			}
			else
			{
				std::cerr << "error " << std::endl;
				exit(-1);
			}
		}

		if (!node->subNodes.empty())
		{
			for (auto& n : node->subNodes)
			{
				str += EncodeTreeLeaf(n);
			}
		}

		return str;
	}

	std::ostream& operator<<(std::ostream& os, const EncodedVec& v)
	{
		std::vector<ENode*> progs;
		v.GetAllPrograms(v.root, progs);
		//std::cerr << "Size = " << progs.size() << std::endl; //debug

		/*
		for (auto& p : progs)
		{
			std::cerr << p->_str << std::endl;
		}
		*/
		// Evaluator::ExportProgramsToFile(progs, std::cerr);
		// std::cerr << std::endl;
		auto score= Evaluator::EvaluatePrograms(progs);
		//std::cerr << std::endl;//debug
		return os << EncodeTreeLeaf(v.root);
	}

	ENode::ENode() : ast()
	{
		ecs.reserve(0);
	}

	ENode::ENode(const std::string& str, bool isFile)
	{
		ecs.reserve(0);
		_str = str;
		if (isFile)
		{
			ast = SexpParser::ParseFile(str);
		}
		else
		{
			ast = SexpParser::Parse(str);
		}
	}

	ENode::ENode(const std::string& sexp, const std::string& ori)
	{
		ecs.reserve(0);
		_str = ori;
		if(!sexp.empty()) 
			ast = SexpParser::Parse(sexp);
	}

	void ENode::EnToEcs(const std::vector<EClass*>& es)
	{
		ecs = es;
	}

	bool ENode::IsENodeProg() const
	{
		return ecs.empty();
	}

	EClass::EClass() = default;

	EClass::EClass(const std::vector<ENode*>& es, string i)
	{
		ens = es;
		idstr = i; // for sub modules, the id is the number of the directory that contains them, for programs it can be -1.
	}

	void EClass::AddEn(ENode* en)
	{
		ens.push_back(en);
	}

	void EClass::SetId(string i)
	{
		idstr = i;
	}

	EGraph::EGraph() = default;

	EGraph::EGraph(const std::vector<EClass*> es)
	{
		ecs = es;
	}

	void EGraph::AddEc(EClass* ec)
	{
		ecs.push_back(ec);
		mapIdToIndex[ec->idstr] = ecs.size() - 1;
	}

	bool EGraph::CheckExistId(std::string id)
	{
		return mapIdToIndex.find(id) != mapIdToIndex.end();
	}

	EClass* EGraph::GetEcById(string id)
	{
		return ecs[mapIdToIndex[id]];
	}

	EClass* EGraph::GetEcByIndex(int index)
	{
		return ecs[index];
	}

	TreeNode::TreeNode() = default;

	TreeNode::TreeNode(std::string x, int h) : ec(nullptr), height(h), val(x), n(0)
	{
	}

	TreeNode::TreeNode(std::string x_, int h_, int n_) : ec(nullptr), height(h_), val(x_), n(n_)
	{

	}

	bool TreeNode::SetN(int _n)
	{
		if (_n == -1)
		{
			n = -1;
			return true;
		}

		if (_n >= ec->ens.size()) return false;

		n = _n;

		if (!subNodes.empty())
		{
			subNodes.clear();
		}

		if (!ec->ens[n]->IsENodeProg())
		{
			for (auto& sec : ec->ens[n]->ecs)
			{
				auto tn = std::make_shared<TreeNode>(sec->idstr, height);
				tn->ec = sec;
				tn->parentPtr = shared_from_this();
				const auto tmpN = Config::get_instance().GetRandom(0, int(sec->ens.size()) - 1);
				tn->SetN(tmpN);
				subNodes.push_back(tn);
			}
		}
		return true;
	}

	int TreeNode::GetN() const
	{
		return n;
	}

	void TreeNode::SetParent(std::shared_ptr<TreeNode> parPtr)
	{
		parentPtr = parPtr;
	}

	bool TreeNode::IsProgram() const
	{
		if (n == -1) return false;
		return ec->ens[n]->IsENodeProg();
	}

	std::shared_ptr<TreeNode> TreeNode::Clone() const
	{
		auto clone = std::make_shared<TreeNode>(val, height);
		clone->ec = ec;
		clone->n = n;
		clone->val = val;
		// avoid copy shared_ptr
		clone->subNodes.resize(subNodes.size());
		for (auto i = 0; subNodes.size() > i; ++i)
		{
			auto c_clone = subNodes[i]->Clone();
			c_clone->SetParent(clone);
			clone->subNodes[i] = c_clone;
		}
		return clone;
	}

	bool EncodedVec::RecursiveInitializeTree(std::shared_ptr<TreeNode>& root, ENode* r, int height) const
	{
		height++;
		for (auto& ec : r->ecs)
		{
			auto tn = std::make_shared<TreeNode>(ec->idstr, height);
			root->subNodes.push_back(tn);
			tn->ec = ec;
			tn->SetN(int(ec->ens.size()) - 1);
			//tn->SetN(0);
			tn->parentPtr = root;
		}
		return true;
	}

	EncodedVec::EncodedVec() : root(nullptr)
	{
	}

	EncodedVec::EncodedVec(ENode* prog)
	{
		// All vecs share the same root
		root = std::make_shared<TreeNode>("-1", 0);
		root->SetN(-1);
		const auto height = 1;
		RecursiveInitializeTree(root, prog, height);
	}

	EncodedVec::EncodedVec(EGraph* eg, ENode* prog, const std::vector<std::tuple<string, int, int>>& vecs, const VectorPI1& tds)
	{
		//eclass_id  
		//v->val,v->height,v->GetN();
		std::vector<std::shared_ptr<TreeNode>> nodes;
		for (auto& vec : vecs)
		{
			auto node = std::make_shared<TreeNode>(std::get<0>(vec), std::get<1>(vec), std::get<2>(vec));
			node->ec = eg->GetEcById(std::get<0>(vec));
			nodes.emplace_back(node);
		}
		for (auto& td : tds)
		{
			nodes[td.second]->SetParent(nodes[td.first]);
			nodes[td.first]->subNodes.emplace_back(nodes[td.second]);
		}

		root = nodes.front();
	}

	EncodedVec::EncodedVec(const EncodedVec& e)
	{
		if (e.root != nullptr)
			root = e.root->Clone();
	}

	std::string EncodedVec::EncodeStr(const EncodedVec& v)
	{
		std::vector<ENode*> progs;
		v.GetAllPrograms(v.root, progs);
		std::string node_ids="";
		for (auto& prog : progs)
			node_ids += prog->idstr +" ";
		return  node_ids;
		//return  EncodeTreeLeaf(v.root);
	}

	std::string EncodedVec::BuildSubTreeMapOld(std::shared_ptr<TreeNode> node,
		std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m)
	{
		if (!node) return "#";
		auto str = node->val;

		for (auto& n : node->subNodes)
		{
			str += "," + BuildSubTreeMapOld(n, m);
		}

		m[str] = node;

		return str;
	}

	auto EncodedVec::MatchCommonSubTreeOld(const std::shared_ptr<TreeNode>& node,
		std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m,
		std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>>&
		res) -> std::string
	{
		if (!node) return "#";
		auto str = node->val;

		for (auto& n : node->subNodes)
		{
			str += "," + MatchCommonSubTreeOld(n, m, res);
		}

		if (m.find(str) != m.end() && node->height != 0) res.emplace_back(node, m[str]);

		return str;
	}

	auto EncodedVec::BuildSubTreeMap(std::shared_ptr<TreeNode> node,std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m) -> int
	{
		if (node->subNodes.empty()) return 0;

		m[node->val] = node;

		for (auto& n : node->subNodes)
		{
			BuildSubTreeMap(n, m);
		}
		return 0;
	}

	auto EncodedVec::MatchCommonSubTree(const std::shared_ptr<TreeNode>& node,
		std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m,
		std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>>& res)
		-> bool
	{
		//if (node->subNodes.size() < 3) return true;
		if (node->subNodes.empty()) return true;

		if (m.find(node->val) != m.end() && node->height != 0 && node->val != "-1") res.emplace_back(node, m[node->val]);

		for (auto& n : node->subNodes)
		{
			MatchCommonSubTree(n, m, res);
		}
		return false;
	}

	bool EncodedVec::GetCommonSubTree(std::shared_ptr<TreeNode> a, std::shared_ptr<TreeNode> b,
		std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>>& res)
	{
		std::unordered_map<std::string, std::shared_ptr<TreeNode>> m;

		BuildSubTreeMap(a, m);

		MatchCommonSubTree(b, m, res);

		return !res.empty();
	}

	bool EncodedVec::GetCommonSubTree(std::shared_ptr<TreeNode> a, std::shared_ptr<TreeNode> b)
	{
		std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>> res;
		return GetCommonSubTree(a, b, res);
	}




	void EncodedVec::CheckError(const std::shared_ptr<TreeNode>& node)
	{
		if (node->GetN() != -1)
		{
			if (node->ec->idstr != node->val)
			{
				system("pause");
			}
		}

		for (auto& subNode : node->subNodes)
		{
			CheckError(subNode);
		}
	}

	void EncodedVec::CheckError2(const std::shared_ptr<TreeNode>& node, std::unordered_set<string>& already)
	{
		if (node->GetN() != -1)
		{
			if (already.find(node->val) != already.end())
			{
				system("pause");
			}

			else
			{
				already.insert(node->val);
			}
		}

		for (auto& subNode : node->subNodes)
		{
			CheckError2(subNode, already);
		}
	}

	bool EncodedVec::CrossOver(EncodedVec& v1, EncodedVec& v2, EncodedVec& s1, EncodedVec& s2)
	{
		//check_error(v1.root);
		//check_error(v2.root);
		s1 = EncodedVec(v1);
		s2 = EncodedVec(v2);
		std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>> res, rq;

		if (!GetCommonSubTree(v1.root, v2.root, rq)) return false;
		const auto coNode =
			Config::get_instance().GetRandom(0, int(rq.size()) - 1);

		auto& r = rq[coNode];
		if (r.first->parentPtr.lock() == nullptr || r.second->parentPtr.lock() == nullptr)
			return false;

		GetCommonSubTree(s1.root, s2.root, res);
		//if (GetCommonSubTree(s1.root, s2.root, res)) std::cerr << "seems like everything is fine" << std::endl;

		auto& nr = res[coNode];

		// step 1: find a ptr to its child
		unsigned a = -1, b = -1;

		for (auto i = 0; i < nr.first->parentPtr.lock()->subNodes.size(); ++i)
		{
			if (nr.first->parentPtr.lock()->subNodes[i] == nr.first)
			{
				a = i;
				// std::cerr << "find a" << std::endl;
			}
		}

		for (auto i = 0; i < nr.second->parentPtr.lock()->subNodes.size(); ++i)
		{
			if (nr.second->parentPtr.lock()->subNodes[i] == nr.second)
			{
				b = i;
				// std::cerr << "find b" << std::endl;
			}
		}

		// step 2: swap ptr to childs
		nr.first->parentPtr.lock()->subNodes[a].swap(nr.second->parentPtr.lock()->subNodes[b]);

		// step 3: swap ptr to parent
		nr.first->parentPtr.swap(nr.second->parentPtr);
		//std::unordered_set<int> c1, c2;
		//check_error2(s1.root, c1);
		//check_error2(s2.root, c2);

		return true;
	}

	void EncodedVec::GetAllPrograms(std::shared_ptr<TreeNode> node, std::vector<ENode*>& program)
	{
		if (!node) return;

		if (node->GetN() != -1&& node->ec->ens[node->GetN()]->IsENodeProg())
		{
			program.push_back(node->ec->ens[node->GetN()]); 
		}
		else
		{
			for (auto& n : node->subNodes)
			{
				GetAllPrograms(n, program);
			}
		}
	}

	std::size_t EncodedVec::size() const
	{
		std::vector<ENode*> program;
		GetAllPrograms(root, program);

		return program.size();
	}
}