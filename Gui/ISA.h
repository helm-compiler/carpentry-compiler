#pragma once
#include <string>
#include <map>
#include <vector>
#include <iostream>

#define STRINGIFY( name ) # name
#define DEFINE_GET_SET_FUNCS double GetValue(const Property a) { return values[a]; }\
	double GetValue(const std::string& a) { if (xmap.find(a) != xmap.end()) return values[xmap[a]]; return 0.0; } \
	void SetValue(const Property a, const double v) { values[a] = v; }\
	void SetValue(const std::string& a, const double v) { if (xmap.find(a) != xmap.end()) values[xmap[a]] = v; }

class Process
{
public:
	virtual bool GetEnabled() const { return enabled; }
	virtual void SetEnabled(const bool e) { enabled = e; }
	virtual void SetValue(const std::string& a, const double v) = 0;
private:
	bool enabled = false;
};


struct ChopsawConfig : public Process
{
	enum Property
	{
		CHOPSAW_WIDTH_MIN = 0,
		CHOPSAW_WIDTH_MAX,
		CHOPSAW_HEIGHT,
		CHOPSAW_MAX_PERP_ANGLE,
		CHOPSAW_MIN_PERP_ANGLE,
		CHOPSAW_MAX_BEVE_ANGLE,
		CHOPSAW_MAX_HEIGHT,
		CHOPSAW_MAX_Y,
		CHOPSAW_MAX_X,
		CHOPSAW_MAX_X_PRIME,
		CHOPSAW_MIN_X_PRIME,
		CHOPSAW_MIN_Z
	};

	DEFINE_GET_SET_FUNCS

private:
	//std::vector<double> values = { 19.9, 40.1, 50.1, 60.01, -50.01, 45.01, 101.61, 152.41, 2438.41, 914.41, 148, 0.0 };
	std::vector<double> values = { 19.9, 40.1, 50.1, 90.01, -90.01, 89.01, 101.61, 152.41, 2438.41, 2500.0, 148, 0.0 };
	//Haisen
	std::map<std::string, Property> xmap = {
		{STRINGIFY(CHOPSAW_WIDTH_MIN), CHOPSAW_WIDTH_MIN},
		{STRINGIFY(CHOPSAW_WIDTH_MAX), CHOPSAW_WIDTH_MAX},
		{STRINGIFY(CHOPSAW_HEIGHT), CHOPSAW_HEIGHT},
		{STRINGIFY(CHOPSAW_MAX_PERP_ANGLE), CHOPSAW_MAX_PERP_ANGLE},
		{STRINGIFY(CHOPSAW_MIN_PERP_ANGLE), CHOPSAW_MIN_PERP_ANGLE},
		{STRINGIFY(CHOPSAW_MAX_BEVE_ANGLE), CHOPSAW_MAX_BEVE_ANGLE},
		{STRINGIFY(CHOPSAW_MAX_HEIGHT), CHOPSAW_MAX_HEIGHT},
		{STRINGIFY(CHOPSAW_MAX_Y), CHOPSAW_MAX_Y},
		{STRINGIFY(CHOPSAW_MAX_X), CHOPSAW_MAX_X},
		{STRINGIFY(CHOPSAW_MAX_X_PRIME), CHOPSAW_MAX_X_PRIME},
		{STRINGIFY(CHOPSAW_MIN_X_PRIME), CHOPSAW_MIN_X_PRIME},
		{STRINGIFY(CHOPSAW_MIN_Z), CHOPSAW_MIN_Z} };
};

struct BandsawConfig : public Process
{
	enum Property
	{
		BANDSAW_MAX_HEIGHT = 0,
		BANDSAW_MAX_X,
		BANDSAW_MAX_X_PRIME,
		BANDSAW_MAX_Y
	};

	DEFINE_GET_SET_FUNCS

private:
	std::vector<double> values = { 152.41, 660.41, 330.21, 609.61 };

	std::map<std::string, Property> xmap = {
		{STRINGIFY(BANDSAW_MAX_HEIGHT), Property::BANDSAW_MAX_HEIGHT},
		{STRINGIFY(BANDSAW_MAX_X), Property::BANDSAW_MAX_X },
		{STRINGIFY(BANDSAW_MAX_X_PRIME), Property::BANDSAW_MAX_X_PRIME},
		{STRINGIFY(BANDSAW_MAX_Y), Property::BANDSAW_MAX_Y }
	};
};

struct JigsawConfig : public Process
{
	enum Property
	{
		JIGSAW_MAX_HEIGHT = 0
	};

	DEFINE_GET_SET_FUNCS

private:
	std::vector<double> values = { 25.41 };

	std::map<std::string, Property> xmap = {
		{STRINGIFY(JIGSAW_MAX_HEIGHT), Property::JIGSAW_MAX_HEIGHT}
	};
};

struct TracksawConfig : public Process
{
	enum Property
	{
		TRACKSAW_MAX_HEIGHT = 0,
		TRACKSAW_MAX_X,
		TRACKSAW_MAX_X_PRIME,
		TRACKSAW_MAX_BEVE_ANGLE,
		TRACKSAW_MAX_Y
	};

	DEFINE_GET_SET_FUNCS

private:
	std::vector<double> values = { 25.41, 2438.41, 914.41, 45.01, 711.21 };

	std::map<std::string, Property> xmap = {
		{STRINGIFY(TRACKSAW_MAX_HEIGHT), Property::TRACKSAW_MAX_HEIGHT},
		{STRINGIFY(TRACKSAW_MAX_X), Property::TRACKSAW_MAX_X },
		{STRINGIFY(TRACKSAW_MAX_X_PRIME), Property::TRACKSAW_MAX_X_PRIME },
		{STRINGIFY(TRACKSAW_MAX_BEVE_ANGLE), Property::TRACKSAW_MAX_BEVE_ANGLE},
		{STRINGIFY(TRACKSAW_MAX_Y), Property::TRACKSAW_MAX_Y }
	};
};

struct DrillPressConfig : public Process
{
	enum Property
	{
	};

	DEFINE_GET_SET_FUNCS;

private:
	std::vector<double> values = {};

	std::map<std::string, Property> xmap;
};


template<typename T> std::shared_ptr<Process> _createProcess() { return std::make_shared<T>(); }


class ISA
{
public:
	ISA();
	~ISA() = default;
	bool AddProcess(const std::string& name, std::shared_ptr<Process> proc)
	{
		mapProcesses_[name] = proc;
		return true;
	}

	bool HasProcess(const std::string& name)
	{
		if (mapProcesses_.find(name) != mapProcesses_.end())
		{
			return mapProcesses_[name]->GetEnabled();
		}
		else
		{
			return false;
		}
	}

	std::shared_ptr<Process> GetProcess(const std::string& name)
	{
		if (HasProcess(name)) return mapProcesses_[name];
		else return nullptr;
	}

	
	typedef std::map<std::string, std::shared_ptr<Process>(*)()> map_type;

	static void InitISA()
	{
		if (!map) { map = new map_type; }
		map->insert(std::make_pair("Chopsaw", _createProcess<ChopsawConfig>));
		map->insert(std::make_pair("Bandsaw", _createProcess<BandsawConfig>));
		map->insert(std::make_pair("Jigsaw", _createProcess<JigsawConfig>));
		map->insert(std::make_pair("Tracksaw", _createProcess<TracksawConfig>));
		map->insert(std::make_pair("Drillpress", _createProcess<DrillPressConfig>));
	}

	static std::shared_ptr<Process> CreateProcess(std::string const& s) {
		map_type::iterator it = getMap()->find(s);
		if (it == getMap()->end())
			return 0;
		return it->second();
	}

protected:
	static map_type* getMap() {
		// never delete'ed. (exist until program termination)
		// because we can't guarantee correct destruction order 
		if (!map) { map = new map_type; }
		return map;
	}


private:
	std::map<std::string, std::shared_ptr<Process> > mapProcesses_;
	static map_type* map;
};