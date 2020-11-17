#pragma once
#include <string>
#include <vector>

namespace OptCompiler
{
	enum Atom
	{
		TYPE_LLHELM,
		TYPE_LENGTH,
		TYPE_ANGLE,
		TYPE_LUMBER,
		TYPE_FACE,
		TYPE_EDGE,
		TYPE_REFC,
		TYPE_REFB,
		TYPE_HEIGHT,
		TYPE_PICK,
		TYPE_STACKABLE,
		TYPE_CHOPSAW,
		TYPE_BANDSAW,
		TYPE_JIGSAW,
		TYPE_TRACKSAW,
		TYPE_DRILL,
		TYPE_XPRIME,
		TYPE_PATH,
		TYPE_TUPLE,
		TYPE_ASSIGN,
		TYPE_SEQ,
		TYPE_VAR,
		TYPE_FLOAT,
		TYPE_REFK,
		TYPE_REFJ,
		TYPE_REFD
	};

	class LLHelm
	{
	public:
		virtual ~LLHelm() = default;
		LLHelm() = default;

		LLHelm(const std::string& s) : content(s)
		{
		};

		virtual const std::string& ToString() { return content; };

		LLHelm eval(LLHelm t);

		int cost_material(LLHelm t);

		int cost_precsion(LLHelm t);

		int cost_time(LLHelm t);

		virtual Atom get_type()
		{
			return TYPE_LLHELM;
		}

	private:
		std::string content;
	};

	class Length : public LLHelm
	{
	public:
		float l{};

		Length();

		Length(float f);

		int cost_material(Length l);

		Atom get_type() override;
	};

	class Angle : public LLHelm
	{
	public:
		Angle();

		Angle(float f);

		int cost_material(Angle a);

		float a;

		Atom get_type() override;
	};

	class Lumber : public LLHelm
	{
	public:

		std::string uuid;

		Lumber();

		Lumber(std::string s);

		int cost_material(Lumber l);

		Atom get_type() override;
	};

	class XPrime : public LLHelm
	{
	public:
		XPrime();

		Atom get_type() override;
	};

	class Face : public LLHelm
	{
	public:

		std::string uuid;

		Face();

		Face(std::string s);

		int cost_material(Face l);

		Atom get_type() override;
	};

	class Edge : public LLHelm
	{
	public:

		std::string uuid;

		Edge();

		Edge(std::string s);

		int cost_material(Edge l);

		Atom get_type() override;
	};

	class Ref : public LLHelm
	{
	public:
		float height;
		std::vector<float> values;
		std::vector<int> intValues;
		bool hasHeight = false;
		bool hasValues = false;

		bool IsComplete() const { return hasHeight && hasValues; }
		virtual void FinalizeRef() {};
	};

	class Refc : public Ref
	{
	public:

		Refc();

		int cost_material(Refc r);

		Atom get_type() override;

		void FinalizeRef() override;
	};

	class Refk : public Ref
	{
	public:

		Refk();

		int cost_material(Refk r);

		Atom get_type() override;

		void FinalizeRef() override;
	};

	class Refj : public Ref
	{
	public:

		Refj();

		int cost_material(Refj r);

		Atom get_type() override;

		void FinalizeRef() override;
	};

	class Refb : public Ref
	{
	public:
		Refb();

		int cost_material(Refb r);

		std::vector<float> values;

		Atom get_type() override;

		void FinalizeRef() override;
	};

	class Refd : public Ref
	{
	public:
		Refd();

		int cost_material(Refd r);

		std::vector<float> values;

		Atom get_type() override;

		void FinalizeRef() override;
	};

	class Height : public LLHelm
	{
	public:
		Height();

		int cost_material(Height h);

		Atom get_type() override;
	};

	class Pick : public LLHelm
	{
	public:

		Pick();

		Pick(int i);

		int getid();

		int id;

		int cost_material(Pick p);

		Atom get_type() override;
	};

	class Stackable : public LLHelm
	{
	public:
		bool stackable;
		Stackable();

		Stackable(bool b);

		int cost_material(Stackable b);

		Atom get_type() override;
	};

	class Bandsaw : public LLHelm
	{
	public:
		Bandsaw();

		Atom get_type() override;
	};

	class Jigsaw : public LLHelm
	{
	public:
		Jigsaw() = default;

		Atom get_type() override
		{
			return TYPE_JIGSAW;
		}
	};

	class Tracksaw : public LLHelm
	{
	public:
		Tracksaw();

		Atom get_type() override;
	};

	class Drill : public LLHelm
	{
	public:
		Drill();

		Atom get_type() override;
	};

	class Path : public LLHelm
	{
	public:
		Path();

		Atom get_type() override;
	};

	class Chopsaw : public LLHelm
	{
	public:
		Lumber l;
		Face f;
		Edge e;
		Refc rc;
		Stackable s;
		Height h;

		Chopsaw();

		Chopsaw(Lumber l1, Face f1, Edge e1, Refc r, Stackable b, Height h1);

		int cost_material(Chopsaw cs);

		Atom get_type() override;
	};

	class Var : public LLHelm
	{
	public:
		std::string name;

		Var();

		Var(std::string s);

		int cost_material(Var v);

		Atom get_type() override;
	};

	class Tuple : public LLHelm
	{
	public:
		LLHelm l;
		LLHelm r;

		Tuple();

		Tuple(LLHelm l1, LLHelm r1);

		int cost_material(Tuple t);

		Atom get_type() override;
	};

	class Assign : public LLHelm
	{
	public:
		LLHelm l;
		LLHelm r;

		Assign();

		Assign(LLHelm l1, LLHelm r1);

		int cost_material(Assign a);

		Atom get_type() override;
	};

	class Seq : public LLHelm
	{
	public:
		LLHelm l;
		LLHelm r;

		Seq();

		Seq(LLHelm l1, LLHelm r1);

		int cost_material(Seq s);

		Atom get_type() override;
	};

	class Float : public LLHelm
	{
	public:
		Float(float f);;
		float num;

		Atom get_type() override;
	};
}