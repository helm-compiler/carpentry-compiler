#ifndef HL_PARSER_HEADER_
#define HL_PARSER_HEADER_

#include <string>
#include "value.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/edge_list.hpp>
#include "parser.hpp"
#include <iostream>

namespace HighLevelParser
{
	enum Atom
	{
		TYPE_HL_HELM,
		TYPE_Make_Stock,
		TYPE_Make_Sketch,
		TYPE_Query_Face_By_Closest_Point,
		TYPE_Query_Edge_By_Closest_Point,
		TYPE_Query_Arc_By_Closest_Center_And_Radius,
		TYPE_Geometry,
		TYPE_Constraint,
		TYPE_Start,
		TYPE_End,
		TYPE_Make_Cut,
		TYPE_Make_Hole,
		TYPE_PointOnObject,
		TYPE_Float,
		TYPE_Assign,
		TYPE_Line,
		TYPE_Angle,
		TYPE_Coincident,
		TYPE_Horizontal,
		TYPE_Vertical,
		TYPE_Distance,
		TYPE_DistanceX,
		TYPE_DistanceY,
		TYPE_Equal,
		TYPE_Parallel,
		TYPE_Circle,
		TYPE_Arc
	};

	class HL_HELM
	{
	public:
		HL_HELM() {};
		HL_HELM(const std::string& s) : str(s) {};
		virtual Atom GetType()
		{
			return TYPE_HL_HELM;
		}

		std::string str;
	};

	class Make_Stock : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Make_Stock; };
	};

	class Make_Sketch : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Make_Sketch; };
	};

	class Query_Face_By_Closest_Point : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Query_Face_By_Closest_Point; };
	};

	class Geometry : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Geometry; };
	};

	class Constraint : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Constraint; };
	};

	class Query_Edge_By_Closest_Point : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Query_Edge_By_Closest_Point; };
	};

	class Query_Arc_By_Closest_Center_And_Radius : public HL_HELM
	{
	public:
		Query_Arc_By_Closest_Center_And_Radius() {};
		Atom GetType() override { return TYPE_Query_Arc_By_Closest_Center_And_Radius; }
	};

	class Start : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Start; };
	};

	class End : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_End; };
	};

	class Make_Cut : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_Make_Cut; };

	};

	class PointOnObject : public HL_HELM
	{
	public:
		Atom GetType() override { return TYPE_PointOnObject; };
	};

	class Float : public HL_HELM
	{
	public:
		Float(float v) : value(v) {};
		Atom GetType() override { return TYPE_Float; }

		float value;
	};

	class Assign : public HL_HELM
	{
	public:
		Assign() {};
		Atom GetType() override { return TYPE_Assign; }
	};

	class Line : public HL_HELM
	{
	public:
		Line() { };
		Atom GetType() override { return TYPE_Line; }
	};

	class Angle : public HL_HELM
	{
	public:
		Angle() { };
		Atom GetType() override { return TYPE_Angle; }
	};

	class Coincident : public HL_HELM
	{
	public:
		Coincident() {};
		Atom GetType() override { return TYPE_Coincident; }
	};

	class Horizontal : public HL_HELM
	{
	public:
		Horizontal() {};
		Atom GetType() override { return TYPE_Horizontal; }
	};

	class Vertical : public HL_HELM
	{
	public:
		Vertical() {};
		Atom GetType() override { return TYPE_Vertical; }
	};

	class Distance : public HL_HELM
	{
	public:
		Distance() {};
		Atom GetType() override { return TYPE_Distance; }
	};

	class DistanceX : public HL_HELM
	{
	public:
		DistanceX() {};
		Atom GetType() override { return TYPE_DistanceX; }
	};

	class DistanceY : public HL_HELM
	{
	public:
		DistanceY() {};
		Atom GetType() override { return TYPE_DistanceY; }
	};

	class Equal : public HL_HELM
	{
	public:
		Equal() {};
		Atom GetType() override { return TYPE_Equal; }
	};

	class Parallel : public HL_HELM
	{
	public:
		Parallel() {};
		Atom GetType() override { return TYPE_Parallel; }
	};

	class Circle : public HL_HELM
	{
	public:
		Circle() {};
		Atom GetType() override { return TYPE_Circle; }
	};

	class Arc : public HL_HELM
	{
	public:
		Arc() {};
		Atom GetType() override { return TYPE_Arc; }
	};

	class Make_Hole : public HL_HELM
	{
	public:
		Make_Hole() {};
		Atom GetType() override { return TYPE_Make_Hole; }
	};

	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, HL_HELM*, boost::property<
		boost::vertex_index_t, std::size_t>,
		boost::no_property> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor GVertex;
	typedef boost::graph_traits<Graph>::edge_descriptor GEdge;
	typedef std::unordered_map<GVertex, std::string> VTS;


	// AST
	class AST
	{
	public:
		AST() = default;

		void CreateVertex(const GVertex& v, const std::string& str)
		{
			content[v] = str;
		}

		GVertex assign;
		GVertex operation;

		Graph g;
		VTS content;
	};

	class HLParser
	{
	public:
		HLParser();
		HLParser(const std::string& file);

		bool ParseFile(const std::string& file);

		HL_HELM* PatternMatching(const std::string& str)
		{
			const auto& cstr = str;
			if (cstr == "Assign") return new Assign();
			if (cstr == "Make_Sketch") return new Make_Sketch();
			if (cstr == "Make_Stock") return new Make_Stock();
			if (cstr == "Make_Cut") return new Make_Cut();
			if (cstr == "Make_Hole") return new Make_Hole();
			if (cstr == "Query_Face_By_Closest_Point") return new Query_Face_By_Closest_Point();
			if (cstr == "Query_Arc_By_Closest_Center_And_Radius") return new Query_Arc_By_Closest_Center_And_Radius();
			if (cstr == "Query_Edge_By_Closest_Point") return new Query_Edge_By_Closest_Point();
			if (cstr == "Geometry") return new Geometry();
			if (cstr == "Constraint") return new Constraint();
			if (cstr == "PointOnObject") return new PointOnObject();
			if (cstr == "Start") return new Start();
			if (cstr == "End") return new End();
			if (cstr == "Line") return new Line();
			if (cstr == "Circle") return new Circle();
			if (cstr == "Arc") return new Arc();
			if (cstr == "Angle") return new Angle();
			if (cstr == "Coincident") return new Coincident();
			if (cstr == "Horizontal") return new Horizontal();
			if (cstr == "Vertical") return new Vertical();
			if (cstr == "Distance") return new Distance();
			if (cstr == "DistanceX") return new DistanceX();
			if (cstr == "DistanceY") return new DistanceY();
			if (cstr == "Equal") return new Equal();
			if (cstr == "Parallel") return new Parallel();
			return new HL_HELM(cstr);
		}

		HL_HELM* PatternMatching(const float num)
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
					//std::cout << value.as_string() << std::endl;
				}
				else if (value.is_real())
				{
					ast.CreateVertex(v, std::to_string(value.as_float()));
					ast.g[v] = PatternMatching(value.as_float());
					//std::cout << value.as_float() << std::endl;
				}
				return;
			}

			auto v = boost::add_vertex(ast.g);
			boost::add_edge(vparent, v, ast.g);

			//std::cout << value.get_car().as_string() << std::endl;
			ast.CreateVertex(v, value.get_car().as_string());
			ast.g[v] = PatternMatching(value.get_car().as_string());
			std::vector<sexp::Value> spt;
			GetLiterals(value.get_cdr(), spt);
			for (auto& l : spt)
			{
				RecursiveParse(l, ast, v);
			}
		}

		void Parse(const std::string& str)
		{
			auto testExpr = sexp::Parser::from_string(str);

			assignAst.emplace_back();
			AST& ast = assignAst.back();

			auto v = boost::add_vertex(ast.g);
			ast.g[v] = new HL_HELM();

			RecursiveParse(testExpr, ast, v);

			int m = 0;
			for (auto& v : ast.g.m_vertices)
			{
				if (v.m_property->GetType() == TYPE_Assign)
				{
					for (auto& e : v.m_out_edges)
					{
						if (m == 0)
							ast.assign = e.get_target();
						else
							ast.operation = e.get_target();

						m++;
					}
					break;
				}
			}
		}

		std::vector<AST>& GetAst() { return assignAst; };

	private:
		std::vector<AST> assignAst;
	};
}
#endif