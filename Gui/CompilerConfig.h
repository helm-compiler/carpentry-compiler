#ifndef _CompilerConfig_H
#define _CompilerConfig_H

#pragma once
#include <string>
#include <boost/assign/list_of.hpp>
#include "PCE.h"
#include "ISA.h"


class CompilerConfig
{
	typedef std::vector<PartDesignGui::WLCollection> VecWLC;
	typedef std::vector<shared_ptr<PartDesignGui::PCEShape>> VecShape;

public:
	CompilerConfig();
	CompilerConfig(const CompilerConfig& root);
	CompilerConfig& operator=(const CompilerConfig&) = delete;
	static const std::string common_test_func_prefix;
	static CompilerConfig& Instance();

	//HMODULE hModule;
	HMODULE& GetHModule();

	void    LoadConfigFile(const std::string& file);
	void	SaveConfigFile(const std::string& file);
	bool    HasLoaded() const;
	double  GetKerf() const;
	double  GetKerf2() const;
	bool    IsRandomPruning() const;
	int     GetRandomSeed() const;
	int     GetMaximalPerm() const;
	void	WriteDefaultConfig(const std::string& filename);
	void	RegenerateTempFolder();
	double  GetPartMatchError() const;
	double  GetAngleMatchError() const;
	double  GetMaximalCuttingDistance() const;
	double	GetSafeDistance() const;
	int     GetTotalOrderNumber() const;
	int     GetPruningLimitation() const;
	VecWLC& GetVecWLC();
	VecShape& GetVecShape();
	const std::string& GetEGraphOutputFolder() const;
	const std::string& GetPackingOutputFolder() const;
	const std::string& GetTemporarilyFolder() const;
	std::shared_ptr<Process> GetProcess(const std::string& name);
	bool HasProcess(const std::string& name);

	int GetOptPopSize() const;
	int GetOptNIslands() const;
	int GetOptEvoGens() const;
	int GetOptNIterations() const;
	int GetOptNTerminal() const;


	const std::string& GetMainObj() const;

	const bool GetRunDesign() const;
	const bool GetRunPack() const;
	const bool GetRunOrder() const;
	const bool GetRunOpt() const;
	const bool GetSharedEgraph() const;

	const int GetInitGeneNb() const;
	const int GetExplorationNb() const;
	const int GetExplorationTerminalNb() const;
	const double GetExplorationM() const;
	const double GetExplorationCR() const;
	//double exploration_m_ = 0.4;
	//double exploration_cr_ = 0.95;

	const bool GetObjOutput() const;

	const Vector3d GetRefPoint() const;

	//t_branch_nb
	//t_bottom_up_select
	const int GetTBottomUpSelect() const;

	const int GetTBranchNb()const;
	const int GetTRandDelete()const;
	const int GetTShuffleSelection() const;
	const int GetTAddShuffle()const;
	const int GetTExactNSGANb() const;
	const int GetTAddFrontNb() const;
	const int GetTDAIFurther() const;

	const int GetPipelineSelection() const;
	const int GetTInitDVS() const;

	const bool GetTInitFab() const;
	const bool GetTInitDepth() const;
	const bool GetTInitSimilarEx() const;
	const bool GetTInitGa() const;
	//t_init_similar_ex

private:

	HMODULE hModule;

	bool   hasLoaded_ = false;
	double kerf_;
	double kerf2_;
	bool   randomPruning_;
	int    randomSeed_;
	int    maximalPerm_;
	double angleMatchError_{};
	double partMatchError_{};
	double maximalCuttingDistance_{};
	double safeDistance_{};
	int    totalOrderNumber_{};
	int    pruningLimitation_{};
	VecWLC wlcs_;
	VecShape shapes_;
	std::string egraph_o_Folder_;
	std::string packing_o_Folder_;
	std::string tempFolder;//TODO: remove it later maybe

	ISA isa_;

	//optimization
	int opt_popSize_;
	int opt_nIslands_;
	int opt_evoGens_;
	int opt_nIterations_;
	int opt_nTerminal_;

	//main_obj
	std::string main_obj_;

	bool run_design_=false;
	bool run_pack_=true;
	bool run_order_=true;
	bool run_opt_=true;

	//design exploration
	int init_gene_nb_;
	int exploration_nb_;
	int expoloration_terminal_nb_;
	double exploration_m_=0.4;
	double exploration_cr_ = 0.95;

	bool shared_egraph_ = true;

	bool obj_output = false;

	Vector3d ref_point_= Vector3d(-1.0,-1.0,-1.0);//hypervolume reference point

	//pipeline 3 4
	int t_bottom_up_select = 4;
	int t_branch_nb=2;
	int t_rand_delete = 2;
	int t_shuffle_selection = 2;
	int t_add_shuffle = 2;
	int t_exact_nsga_nb = 3;
	int t_add_front_nb = 2;
	int t_dai_further = 2;

	int t_init_dvs = 100;
	bool t_init_similar_ex = false;
	bool t_init_fab = true;
	bool t_init_depth = true;
	bool t_init_ga = false;

	//pipeline algorithm
	int pipeline_selection = 0;
};

#endif // PARTDESIGN_WORKBENCH_H