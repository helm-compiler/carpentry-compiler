#include <parsers.hpp>
#include <sexp/parser.hpp>
#include <sexp/util.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <graphmap.hpp>

namespace OptCompiler
{
	namespace SexpParser
	{
		LLHelm* PatternMatching(const std::string& str)
		{
			const auto& cstr = str;
			if (cstr == "Seq") return new Seq();
			if (cstr == "Angle") return new Angle();
			if (cstr == "Lumber") return new Lumber();
			if (cstr == "Face") return new Face();
			if (cstr == "Edge") return new Edge();
			if (cstr == "Assign") return new Assign();
			if (cstr == "Refc") return new Refc();
			if (cstr == "Refb") return new Refb();
			if (cstr == "Refk") return new Refk();
			if (cstr == "Refj") return new Refj();
			if (cstr == "Refd") return new Refd();
			if (cstr == "Height") return new Height();
			if (cstr == "Stackable") return new Stackable();
			if (cstr == "Var") return new Var();
			if (cstr == "Tup") return new Tuple();
			if (cstr == "Assign") return new Assign();
			if (cstr == "Chopsaw") return new Chopsaw();
			if (cstr == "Length") return new Length();
			if (cstr == "Jigsaw") return new Jigsaw();
			if (cstr == "Tracksaw") return new Tracksaw();
			if (cstr == "XPrime") return new XPrime();
			if (cstr == "Path") return new Path();
			if (cstr == "Length") return new Length();
			if (cstr == "Bandsaw") return new Bandsaw();
			if (cstr == "Drill") return new Drill();

			return new LLHelm(cstr);
		}

		LLHelm* PatternMatching(const float num)
		{
			return new Float(num);
		}

		void GetLiterals(sexp::Value& val, std::vector<sexp::Value>& springs)
		{
			if (val.is_nil())
				return;
			if (val.is_cons())
			{
				springs.push_back(val.get_car());
				GetLiterals(val.get_cdr(), springs);
			}
		}

		void RecursiveParse(sexp::Value& value, AST& ast, GVertex vparent)
		{
			if (value.is_nil()) return;

			if (!value.is_cons())
			{
				auto v = boost::add_vertex(ast.g);
				boost::add_edge(vparent, v, ast.g);
				if (value.is_symbol())
				{
					ast.CreateVertex(v, value.as_string());
					ast.g[v] = PatternMatching(value.as_string());
					//std::cerr << value.as_string() << std::endl;
				}
				else if (value.is_real())
				{
					ast.CreateVertex(v, std::to_string(value.as_float()));
					ast.g[v] = PatternMatching(value.as_float());
					//std::cerr << value.as_float() << std::endl;
				}
				return;
			}

			auto v = boost::add_vertex(ast.g);
			boost::add_edge(vparent, v, ast.g);

			//std::cerr << value.get_car().as_string() << std::endl;
			ast.CreateVertex(v, value.get_car().as_string());
			ast.g[v] = PatternMatching(value.get_car().as_string());
			std::vector<sexp::Value> spt;
			GetLiterals(value.get_cdr(), spt);
			for (auto& l : spt)
			{
				RecursiveParse(l, ast, v);
			}
		}

		void AST::PreProcess()
		{
			GVertex vtxProcess;
			bfs_process_visitor vis(vtxProcess, outputParts_, inputMaterial_,referFace_,referEdge_);
			boost::depth_first_search(g, visitor(vis));
			process_ = g[vtxProcess];

			auto ref = std::pair<bool, Ref*>(false, nullptr);
			auto idmap = get(boost::vertex_index, g);
			std::vector<boost::default_color_type> vertex_color(num_vertices(g));
			auto vcmap = make_iterator_property_map(vertex_color.begin(), idmap);
			bfs_ref_visitor refVis(vcmap, ref);
			boost::breadth_first_search(g, vtxProcess, visitor(refVis));

			if (ref.first)
			{
				auto t = dynamic_cast<Ref*>(ref.second);
				ref_ = ref.second;
				ref_->FinalizeRef();
			}
		}
		int debug_debug_nb = 0;
		AST Parse(const std::string& str)
		{
			//std::cerr << str << std::endl;
			auto testExpr = sexp::Parser::from_string(str);
			AST ast;
			auto v = boost::add_vertex(ast.g);
			ast.g[v] = new LLHelm();

			RecursiveParse(testExpr, ast, v);

			ast.PreProcess();

			auto d= ast.GetReferFace();

			if (ast.GetReferFace().empty())
			{
				int d = 0;
			}
			// system("pause");
			/*
			std::ofstream file("dot.dot");
			boost::write_graphviz(file, ast.g);
			file.close();
			*/
			debug_debug_nb++;
			return ast;
		}

		AST ParseFile(const std::string& filename)
		{
			std::ifstream file(filename);
			auto testExpr = sexp::Parser::from_stream(file);
			file.close();

			AST ast;
			const auto v = boost::add_vertex(ast.g);
			ast.g[v] = new LLHelm();
			RecursiveParse(testExpr, ast, v);

			/*
			std::ofstream file("dot.dot");
			boost::write_graphviz(file, ast.g);
			file.close();
			*/
			return ast;
		}
	}

	void EGraphParser::ParseXML(const std::string& path, std::unordered_map<int, SubModule>& all_progs,
		std::unordered_set<int>& emptyProgs)
	{
		using namespace tinyxml2;
		all_progs.clear();

		XMLDocument doc;
		doc.LoadFile(path.c_str());
		XMLElement* root = doc.FirstChildElement("root");
		for (XMLElement* e = root->FirstChildElement("sub"); e != nullptr; e = e->NextSiblingElement("sub"))
		{
			std::cerr << e->Attribute("eclassID") << std::endl;
			int subId = atoi(e->Attribute("eclassID"));
			std::string subId_str = Math::Functs::IntString(subId);
			SubModule mod_progs;
			
			auto& cutLine_cMap = Config::get_instance().ecCutLineIdmaps[subId_str];
			auto& enode_cMap = Config::get_instance().ecENodeIdmaps[subId_str];
			auto& arr_cMap = Config::get_instance().ecArrIdmaps[subId_str];
			auto& wlc_cMap = Config::get_instance().ecWLCIdmaps[subId_str];

			//ecProgArrWLCs
			for (XMLElement* p = e->FirstChildElement("prog"); p != nullptr; p = p->NextSiblingElement("prog"))
			{
				Vector2<PairStrs> prog;
				Vector1i1 prog_cutLineIds;
				for (XMLElement* eq = p->FirstChildElement("equiv"); eq != nullptr; eq = eq->NextSiblingElement("equiv"))
				{
					std::vector<PairStrs> equiv_prog;
					for (auto ln = eq->FirstChildElement("line"); ln != nullptr; ln = ln->NextSiblingElement("line"))
					{
						const auto sexp = ln->FirstChildElement("sexp");
						const auto ori = ln->FirstChildElement("ori");

						if(!sexp->GetText())
							equiv_prog.emplace_back("", ori->GetText());
						else
							equiv_prog.emplace_back(sexp->GetText(), ori->GetText());
					}
					prog.push_back(equiv_prog);
					prog_cutLineIds.emplace_back(atoi(eq->Attribute("cutLineID")));
				}

				//atoi(p->Attribute("arr"));
				//atoi(p->Attribute("wlc"));
				cutLine_cMap[mod_progs.size()] = prog_cutLineIds;
				enode_cMap[mod_progs.size()] = atoi(p->Attribute("enodeID"));
				arr_cMap[mod_progs.size()] = atoi(p->Attribute("arrID"));
				wlc_cMap[mod_progs.size()] = atoi(p->Attribute("wlcID"));

				Config::get_instance().Push_ArrWLC(atoi(p->Attribute("arrID")), atoi(p->Attribute("wlcID")));

				mod_progs.push_back(prog);
			}

			all_progs[subId] = mod_progs;

			if (mod_progs.empty())
			{
				emptyProgs.insert(subId);
			}
		}
	}

	EClass* EGraphParser::ParseEGraph(const std::string& path, EGraph* eg)
	{
		using namespace tinyxml2;
		string rootID = "-1";
		XMLDocument doc;
		doc.LoadFile(path.c_str());
		XMLElement* root = doc.FirstChildElement("root");

		for (XMLElement* e = root->FirstChildElement("EClass"); e != nullptr; e = e->NextSiblingElement("EClass"))
		{
			std::cerr << "EClass ID = " << e->Attribute("ID") << std::endl;
			string ecID = e->Attribute("ID");
			//rootID = ecID;

			for (XMLElement* p = e->FirstChildElement("ENode"); p != nullptr; p = p->NextSiblingElement("ENode"))
			{
				ENode* en = new ENode();//ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode ENode 

				for (XMLElement* n = p->FirstChildElement("C"); n != nullptr; n = n->NextSiblingElement("C"))
				{
					std::string ecPtr = n->GetText();

					if (eg->CheckExistId(ecPtr))
					{
						en->ecs.push_back(eg->GetEcById(ecPtr));
					}
					else
					{
						//std::cerr << "if (eg->CheckExistId(ecPtr))" << std::endl;
						//system("pause");
					}
				}
				if (eg->CheckExistId(ecID)&& !en->ecs.empty())
					eg->GetEcById(ecID)->AddEn(en);
				else
				{
					//std::cerr << "if (eg->CheckExistId(ecID))" << std::endl;
					//system("pause");
				}
			}
		}

		std::cerr << "Root ID = " << rootID << std::endl;
		std::cerr << "Node in Root ID = " << eg->GetEcById(rootID)->ens.size() << std::endl;
		return eg->GetEcById(rootID);
	}
}