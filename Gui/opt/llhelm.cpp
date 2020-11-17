#include <llhelm.hpp>
#include <algorithm>
#include <set>

#define SCALE_FACTOR 1000000.

namespace OptCompiler
{

	Length::Length() = default;

	Length::Length(float f)
	{
		l = f;
	}

	int Length::cost_material(Length l)
	{
		return 0;
	}

	Atom Length::get_type()
	{
		return TYPE_LENGTH;
	}

	Angle::Angle() = default;

	Angle::Angle(float f)
	{
		a = f;
	}

	int Angle::cost_material(Angle a)
	{
		return 0;
	}

	Atom Angle::get_type()
	{
		return TYPE_ANGLE;
	}

	Lumber::Lumber() = default;

	Lumber::Lumber(std::string s)
	{
		uuid = s;
	}

	int Lumber::cost_material(Lumber l)
	{
		return 1;
	}

	Atom Lumber::get_type()
	{
		return TYPE_LUMBER;
	}

	XPrime::XPrime() = default;

	Atom XPrime::get_type()
	{
		return TYPE_XPRIME;
	}

	Face::Face() = default;

	Face::Face(std::string s)
	{
		uuid = s;
	}

	int Face::cost_material(Face l)
	{
		return 0;
	}

	Atom Face::get_type()
	{
		return TYPE_FACE;
	}

	Edge::Edge() = default;

	Edge::Edge(std::string s)
	{
		uuid = s;
	}

	int Edge::cost_material(Edge l)
	{
		return 0;
	}

	Atom Edge::get_type()
	{
		return TYPE_EDGE;
	}

	Refc::Refc() = default;

	int Refc::cost_material(Refc r)
	{
		return 0;
	}

	Atom Refc::get_type()
	{
		return TYPE_REFC;
	}

	void Refc::FinalizeRef()
	{
		intValues.reserve(values.size());
		for (auto& v : values)
		{
			intValues.push_back(SCALE_FACTOR * v);
		}
	}

	Refk::Refk() = default;

	int Refk::cost_material(Refk r)
	{
		return 0;
	}

	Atom Refk::get_type()
	{
		return TYPE_REFK;
	}

	void Refk::FinalizeRef()
	{
		intValues.reserve(values.size());
		for (auto& v : values)
		{
			intValues.push_back(SCALE_FACTOR * v);
		}
	}

	Refj::Refj() = default;

	int Refj::cost_material(Refj r)
	{
		return 0;
	}

	Atom Refj::get_type()
	{
		return TYPE_REFJ;
	}

	void Refj::FinalizeRef()
	{
		std::set<int> unqiueSet;
		intValues.reserve(values.size());
		for (auto& v : values)
		{
			const int scaledV = SCALE_FACTOR * v;
			if (unqiueSet.find(scaledV) == unqiueSet.end())
			{
				intValues.push_back(scaledV);
				unqiueSet.insert(scaledV);
			}
		}
	}

	Refb::Refb() = default;

	int Refb::cost_material(Refb r)
	{
		return 0;
	}

	Atom Refb::get_type()
	{
		return TYPE_REFB;
	}

	void Refb::FinalizeRef()
	{
		std::set<int> unqiueSet;
		intValues.reserve(values.size());
		for (auto& v : values)
		{
			const int scaledV = SCALE_FACTOR * v;
			if (unqiueSet.find(scaledV) == unqiueSet.end())
			{
				intValues.push_back(scaledV);
				unqiueSet.insert(scaledV);
			}
		}
	}

	Refd::Refd() = default;

	int Refd::cost_material(Refd r)
	{
		return 0;
	}

	Atom Refd::get_type()
	{
		return TYPE_REFD;
	}

	void Refd::FinalizeRef()
	{
		std::set<int> unqiueSet;
		intValues.reserve(values.size());
		for (auto& v : values)
		{
			const int scaledV = SCALE_FACTOR * v;
			if (unqiueSet.find(scaledV) == unqiueSet.end())
			{
				intValues.push_back(scaledV);
				unqiueSet.insert(scaledV);
			}
		}
	}

	Height::Height() = default;

	int Height::cost_material(Height h)
	{
		return 0;
	}

	Atom Height::get_type()
	{
		return TYPE_HEIGHT;
	}

	Pick::Pick() = default;

	Pick::Pick(int i)
	{
		id = i;
	}

	int Pick::getid()
	{
		return id;
	}

	int Pick::cost_material(Pick p)
	{
		return 10000000;
	}

	Atom Pick::get_type()
	{
		return TYPE_PICK;
	}

	Stackable::Stackable() = default;

	Stackable::Stackable(bool b)
	{
		stackable = b;
	}

	int Stackable::cost_material(Stackable b)
	{
		return 0;
	}

	Atom Stackable::get_type()
	{
		return TYPE_STACKABLE;
	}

	Bandsaw::Bandsaw() = default;

	Atom Bandsaw::get_type()
	{
		return TYPE_BANDSAW;
	}

	Tracksaw::Tracksaw() = default;

	Atom Tracksaw::get_type()
	{
		return TYPE_TRACKSAW;
	}

	Path::Path() = default;

	Atom Path::get_type()
	{
		return TYPE_PATH;
	}

	Chopsaw::Chopsaw() = default;

	Chopsaw::Chopsaw(Lumber l1, Face f1, Edge e1, Refc r, Stackable b, Height h1)
	{
		l = l1;
		f = f1;
		e = e1;
		rc = r;
		s = b;
		h = h1;
	}

	int Chopsaw::cost_material(Chopsaw cs)
	{
		return 0;
	}

	Atom Chopsaw::get_type()
	{
		return TYPE_CHOPSAW;
	}

	Var::Var() = default;

	Var::Var(std::string s)
	{
		name = s;
	}

	int Var::cost_material(Var v)
	{
		return 0;
	}

	Atom Var::get_type()
	{
		return TYPE_VAR;
	}

	Tuple::Tuple() = default;

	Tuple::Tuple(LLHelm l1, LLHelm r1)
	{
		l = l1;
		r = r1;
	}

	int Tuple::cost_material(Tuple t)
	{
		return 0;
	}

	Atom Tuple::get_type()
	{
		return TYPE_TUPLE;
	}

	Assign::Assign() = default;

	Assign::Assign(LLHelm l1, LLHelm r1)
	{
		l = l1;
		r = r1;
	}

	int Assign::cost_material(Assign a)
	{
		return 0;
	}

	Atom Assign::get_type()
	{
		return TYPE_ASSIGN;
	}

	Seq::Seq() = default;

	Seq::Seq(LLHelm l1, LLHelm r1)
	{
		l = l1;
		r = r1;
	}

	int Seq::cost_material(Seq s)
	{
		return 0;
	}

	Atom Seq::get_type()
	{
		return TYPE_SEQ;
	}

	Float::Float(float f) : num(f)
	{
	}

	Atom Float::get_type()
	{
		return TYPE_FLOAT;
	}

	Drill::Drill()
	{
	}

	Atom Drill::get_type()
	{
		return TYPE_DRILL;
	}

}