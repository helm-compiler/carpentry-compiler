#pragma once
#include <parsers.hpp>
#include <random>
#include <graphmap.hpp>

namespace OptCompiler {

	class EClass;
	class TreeNode;

	void EncodeTreeLeaf1(const std::shared_ptr<TreeNode>& node, std::vector<std::tuple<string, int, int, Vector1i1>>& tuples);

	int CollectTreeNodes(const std::shared_ptr<TreeNode>& v, std::vector<std::tuple<string,int,int>>& tn, VectorPI1& edges);

	class ENode
	{
	public:
		ENode();

		ENode(const std::string& str, bool isFile = true);

		ENode(const std::string& sexp, const std::string& ori);

		~ENode()
		{
			ecs.clear();
			cutLinesIds.clear();
		}

		void EnToEcs(const std::vector<EClass*>& es);

		bool IsENodeProg() const;
		
		void SetId(string i) { idstr = i; };

		std::vector<EClass*> ecs;
		SexpParser::AST ast;
		std::string _str;
		int arr=-1;//arrange index in pce 
		int wlc=-1;//wlc index in wlcs
		int eclassID = -1;//eclass id 
		int enodeID=-1;//enode id in the current eclass; it's a local index number.
		int cutLineID = -1;//cutline id in the original cutlines of its arrange
		int equivalenceID = -1;// equivalence index, check the input xml format
		int lineID = -1; //line index in the current equivalence
		Vector1i1 cutLinesIds;	

		//int id{};
		std::string idstr{};
	};

	class EClass
	{
	public:

		EClass();

		EClass(const std::vector<ENode*>& es, string i);

		void AddEn(ENode* en);

		void SetId(string i);

		std::vector<ENode*> ens;
		//int id{};

		std::string idstr{};
	};


	class EGraph
	{
	public:
		EGraph();

		EGraph(std::vector<EClass*> es);

		void AddEc(EClass* ec);

		EClass* GetEcById(string id);

		bool CheckExistId(std::string id);

		EClass* GetEcByIndex(int index);

		std::vector<EClass*> ecs;
		std::unordered_map<string, int> mapIdToIndex;

		
		int stock_e_nodes=0;
		int union_e_nodes=0;


		void NB(EClass* ec)
		{
			for (auto& en : ec->ens)
			{
				if (en->ecs.empty())
				{
					stock_e_nodes++;
				}
				else
				{
					union_e_nodes++;
					for (auto& ec_ : en->ecs)
					{
						NB(ec_);
					}
				}
			}
		};

		void ENodeNB()
		{
			//for (int i = 0; i < ecs.size(); i++)
			//	NB(ecs[i]);
			NB(GetEcById("-1"));
		}
	};

	class TreeNode : public std::enable_shared_from_this<TreeNode>
	{
	public:
		TreeNode();

		TreeNode(std::string x, int h);
		TreeNode(std::string x_, int h_, int n_);

		bool SetN(int _n);

		int GetN() const;

		void SetParent(std::shared_ptr<TreeNode> parPtr);;

		bool IsProgram() const;

		std::shared_ptr<TreeNode> Clone() const;

	public:
		EClass* ec{};
		std::vector<std::shared_ptr<TreeNode>> subNodes;
		std::weak_ptr<TreeNode> parentPtr;
		int height{};
		std::string val{}; // id of e-class

	private:
		int n{}; // specific program n, n = -1 means it point to another ec
	};

	class EncodedVec
	{
		bool RecursiveInitializeTree(std::shared_ptr<TreeNode>& root, ENode* r, int height) const;

	public:
		EncodedVec();

		EncodedVec(ENode* prog);

		EncodedVec(EGraph* eg, ENode* prog, const std::vector<std::tuple<string,int,int>>& vecs, const VectorPI1& tds);

		EncodedVec(const EncodedVec& e);

		friend std::ostream& operator<<(std::ostream& os, const EncodedVec& v);

		static std::string EncodeStr(const EncodedVec& v);

		static std::string BuildSubTreeMapOld(std::shared_ptr<TreeNode> node,
			std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m);

		static std::string MatchCommonSubTreeOld(const std::shared_ptr<TreeNode>& node,
			std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m,
			std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>
			>& res);


		static int BuildSubTreeMap(std::shared_ptr<TreeNode> node, std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m);

		static bool MatchCommonSubTree(const std::shared_ptr<TreeNode>& node,
			std::unordered_map<std::string, std::shared_ptr<TreeNode>>& m,
			std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>>& res);

		static bool GetCommonSubTree(std::shared_ptr<TreeNode> a,
			std::shared_ptr<TreeNode> b,
			std::vector<std::pair<std::shared_ptr<TreeNode>, std::shared_ptr<TreeNode>>>& res);

		static bool GetCommonSubTree(std::shared_ptr<TreeNode> a, std::shared_ptr<TreeNode> b);

		static void CheckError(const std::shared_ptr<TreeNode>& node);

		static void CheckError2(const std::shared_ptr<TreeNode>& node, std::unordered_set<string>& already);

		static bool CrossOver(EncodedVec& v1, EncodedVec& v2, EncodedVec& s1, EncodedVec& s2);

		static void GetAllPrograms(std::shared_ptr<TreeNode> node, std::vector<ENode*>& program);

		std::size_t size() const;

		std::shared_ptr<TreeNode> root;
	};

}