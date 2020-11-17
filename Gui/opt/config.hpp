#pragma once
#include <string>
#include <random>
#include <mutex>
#include <condition_variable>
#include <Mod/PartDesign/Gui/CompilerConfig.h>




class Config
{
public:
	static const std::string common_test_func_prefix;

	static Config& get_instance()
	{
		static Config instance;
		return instance;
	}

	Config(const Config& root) = delete;
	Config& operator=(const Config&) = delete;

	Config();

	int GetRandom(const int& _a, const int& _b)
	{
		std::uniform_int_distribution<> dis(_a, _b);
		return dis(gen);
	}

	double GetProb()
	{
		std::uniform_real_distribution<> drng(0., 1.); // to generate a number in [0, 1)
		return drng(gen);
	}

	std::string SUBMODULE_FILE() const { return submodule; }

	std::string FULLPROGS_PATH() const { return fullprogs; }

	std::string EGRAPHFILE_PATH() const { return egraphfile; }

	std::string ORIGINAL_FILE() const { return oriprogs; }

	std::string RESULT_FILE() const { return resultfile; }

	std::string HYPERVOLUME_FILE() const { return hypervolume; }

	//livefile
	std::string LIVEFILE_FILE() const { return livefile; }

	std::string INIT_FILE() const { return initfile; }

	std::string FOLDER_NAME() const { return foldername; }

	unsigned GetRandSeed() const { return rdSeed; }

	unsigned GetPopSize() const { return popSize; }

	unsigned GetNIslands() const { return nIslands; }

	unsigned GetEvoGens() const { return evoGens; }

	unsigned GetIterations() const { return nIterations; }

	int GetNTerminal() const { return nTerminal; }

	void AddHypervolume(const std::vector<float>& hv);
	void AddParetoFrontsNB(const std::vector<int>& nbs);
	void AddUniqueParetoFrontsNB(const std::vector<int>& nbs);
	void AddPopDifNb(const std::vector<int>& nbs);
	void AddArrangeNB(const std::vector<int>& nbs);

	void SetConfigs(const std::string taskName_="");

	std::vector<float> PostprocessHypervolume() const;

	std::vector<int> PostprocessParetoFrontsNB() const;
	std::vector<int> PostprocessUniqueParetoFrontsNB() const;
	std::vector<int> PostprocessPopDiffNb() const;

	std::vector<int> PostprocessArrangesNB() const;


	void SetRefpoint(double* p);

	const double* GetRefpoint() const;

	void UpdateHyperVolume(const std::vector<float>& updated);
	
	struct ArrWLC
	{
		int index = -1;
		int arr = -1;
		int wlc = -1;
		std::vector<int> iters;
		ArrWLC(int index_, int arr_, int wlc_):index(index_),arr(arr_), wlc(wlc_) {};
	};

	void Clear()
	{
		ecCutLineIdmaps.clear();
		ecENodeIdmaps.clear();
		ecArrIdmaps.clear();
		ecWLCIdmaps.clear();
		ecProgArrWLCs.clear();
		arrangesNB.clear();
		paretoFrontsNB.clear();
		hyperVolume.clear();

		//std::random_device rd; //Will be used to obtain a seed for the random number engine
		//std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()


		rdSeed = rd();
		popSize = 40;//ok
		nIslands = 1;//default(6) //ok
		evoGens = 30;
		nIterations = 1;//default(30)
		nTerminal = -1;//default(-1)
		lookUpInterval = 800;
		current = 0;
		//refPoint[3];
		currentHv = 0;
	}

	std::unordered_map<std::string, std::unordered_map<int, Vector1i1>> ecCutLineIdmaps;
	std::unordered_map<std::string, std::unordered_map<int, int>> ecENodeIdmaps;
	std::unordered_map<std::string, std::unordered_map<int, int>> ecArrIdmaps;
	std::unordered_map<std::string, std::unordered_map<int, int>> ecWLCIdmaps;

	std::vector<ArrWLC> ecProgArrWLCs;

	void Push_ArrWLC(int arr_, int wlc_)
	{
		for (auto& arr_wlc : ecProgArrWLCs)
			if (arr_wlc.arr == arr_ && arr_wlc.wlc == wlc_) return;
		ecProgArrWLCs.emplace_back(ArrWLC(ecProgArrWLCs.size(),arr_,wlc_));
	};



	bool Push_Iter_ArrWLC(int arr_, int wlc_, int iter_)
	{
		bool return_b = false;

		for (auto& arr_wlc : ecProgArrWLCs)
		{
			if (arr_wlc.arr == arr_ && arr_wlc.wlc == wlc_)
			{
				bool b = false;
				for (auto& iter : arr_wlc.iters)
				{
					if (iter == iter_)
					{
						b = true;
						break;
					}
				}

				if (!b)
				{
					arr_wlc.iters.emplace_back(iter_);
					return_b = true;
				}
				break;
			}
		}

		return return_b;
	}


private:
	std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
	std::string foldername = taskName + "/";
	std::string submodule = taskName + "/collapse/all_progs.xml";
	std::string fullprogs = taskName + "/subs_progs_only_progs/";
	std::string oriprogs = taskName + "/collapse/ori_progs.xml";
	std::string egraphfile = taskName + "/egraph.xml";
	std::string resultfile = taskName + "/result.xml";
	std::string hypervolume = taskName + "/collapse/hyervolume.txt";
	std::string livefile = taskName + "/collapse/livefile.txt";
	std::string initfile = taskName + "/collapse/init.txt";

	unsigned rdSeed = rd();
	unsigned popSize = 40;//ok
	unsigned nIslands = 1;//default(6) //ok
	unsigned evoGens = 30;
	unsigned nIterations = 1;//default(30)
	int nTerminal = -1;//default(-1)
	int lookUpInterval = 800;
	int current = 0;
	double refPoint[3];
	int currentHv = 0;

	
	Vector1i1 arrangesNB;
	Vector1i1 paretoFrontsNB;
	Vector1i1 uniqueParetoFrontsNB;
	Vector1i1 popDifNb;
	std::vector<float> hyperVolume;
	std::mutex mtx;
	std::condition_variable cv;

	std::string taskName = "bookcase-pd-0701";
};

inline Config::Config()  // NOLINT(cert-msc32-c, cert-msc51-cpp)
{
	taskName = CompilerConfig::Instance().GetEGraphOutputFolder();
	std::cerr << "Task name = " << taskName << std::endl;
	foldername = taskName + "/";
	submodule = taskName + "/collapse/all_progs.xml";
	fullprogs = taskName + "/subs_progs_only_progs/";
	oriprogs = taskName + "/collapse/ori_progs.xml";
	egraphfile = taskName + "/egraph.xml";
	resultfile = taskName + "/collapse/result.xml";
	hypervolume = taskName + "/collapse/hyervolume.txt";
	livefile = taskName + "/collapse/livefile.txt";
	initfile = taskName + "/collapse/init.txt";

	//hardware concurrency
	int processor_count = std::thread::hardware_concurrency();
	if (nIslands > processor_count) nIslands = processor_count;

	// random number engine
	rdSeed = rd();
	//rdSeed = rd();
	gen = std::mt19937(rdSeed);

	unsigned popSize = 40;//ok
	unsigned nIslands = 1;//default(6) //ok
	unsigned evoGens = 30;
	unsigned nIterations = 1;//default(30)

	popSize = CompilerConfig::Instance().GetOptPopSize();
	nIslands = CompilerConfig::Instance().GetOptNIslands();
	evoGens = CompilerConfig::Instance().GetOptEvoGens();
	nIterations = CompilerConfig::Instance().GetOptNIterations();
	nTerminal = CompilerConfig::Instance().GetOptNTerminal();
}

inline void Config::SetConfigs(const std::string taskName_)
{
	taskName = taskName_;
	if(taskName_.size()==0)
		taskName = CompilerConfig::Instance().GetEGraphOutputFolder();

	std::cerr << "Task name = " << taskName << std::endl;
	foldername = taskName + "/";
	submodule = taskName + "/collapse/all_progs.xml";
	fullprogs = taskName + "/subs_progs_only_progs/";
	oriprogs = taskName + "/collapse/ori_progs.xml";
	egraphfile = taskName + "/egraph.xml";
	resultfile = taskName + "/collapse/result.xml";
	hypervolume = taskName + "/collapse/hyervolume.txt";
	livefile = taskName + "/collapse/livefile.txt";
	initfile = taskName + "/collapse/init.txt";

	// random number engine
	rdSeed = 4;
	//rdSeed = rd();
	gen = std::mt19937(rdSeed);

	popSize = CompilerConfig::Instance().GetOptPopSize();
	nIslands = CompilerConfig::Instance().GetOptNIslands();
	evoGens = CompilerConfig::Instance().GetOptEvoGens();
	nIterations = CompilerConfig::Instance().GetOptNIterations();
	nTerminal = CompilerConfig::Instance().GetOptNTerminal();
}

inline void Config::UpdateHyperVolume(const std::vector<float>& updated)
{
	//std::unique_lock<std::mutex> lock(mtx);

	int ind = currentHv - GetEvoGens();

	for (auto i = 0; i < GetEvoGens()&&i< updated.size(); ++i)
	{
		hyperVolume[ind + i] = std::max(hyperVolume[ind + i], updated[i]);
	}
}
inline void Config::AddArrangeNB(const std::vector<int>& nbs)
{
	arrangesNB = nbs;
}

inline void Config::AddParetoFrontsNB(const std::vector<int>& nbs)
{
	//TODO: did not work  for parallel computation
	paretoFrontsNB = nbs;
}

inline void Config::AddUniqueParetoFrontsNB(const std::vector<int>& nbs)
{
	uniqueParetoFrontsNB = nbs;
}
inline void Config::AddPopDifNb(const std::vector<int>& nbs)
{
	popDifNb = nbs;
}

inline void Config::AddHypervolume(const std::vector<float>& hv)
{
	if (evoGens == 0) return;
	// wait untill all threads get here
	{
		std::unique_lock<std::mutex> lock(mtx);
		current++;
		lock.unlock();
	}

	// for each thread, add its
	{
		std::unique_lock<std::mutex> lck(mtx);
		if (current != nIslands) { cv.wait(lck); }
		
		if (current == nIslands)
		{
			auto curSize = hyperVolume.size();
			
			if (curSize > lookUpInterval)
			{
				if (std::abs(hyperVolume[curSize - 1]- hyperVolume[curSize - lookUpInterval - 1]) < 5)
				{
					evoGens = 0;
				}
			}

			if (!hyperVolume.empty())
				std::cerr << "Current lead hypervolume = " << hyperVolume.back() << std::endl;
			
			for (auto i = 0; i < evoGens; ++i)
			{
				hyperVolume.push_back(0);
			}
			currentHv += evoGens;
			std::cerr << hyperVolume.size() << std::endl;
		}

		UpdateHyperVolume(hv);
	
		current = 0;
		cv.notify_all();
	}
}

inline std::vector<float> Config::PostprocessHypervolume() const
{
	return hyperVolume;
}

inline std::vector<int> Config::PostprocessParetoFrontsNB() const
{
	return paretoFrontsNB;
}
inline std::vector<int> Config::PostprocessUniqueParetoFrontsNB() const
{
	return uniqueParetoFrontsNB;
}

inline std::vector<int> Config::PostprocessPopDiffNb() const
{
	return popDifNb;
}

inline std::vector<int> Config::PostprocessArrangesNB() const
{
	return arrangesNB;
}

inline void Config::SetRefpoint(double* p)
{
	for (int i = 0; i < 3; ++i)
	{
		refPoint[i] = p[i];
	}
}

inline const double* Config::GetRefpoint() const
{
	return refPoint;
}
