#include "CompilerConfig.h"
#include "tinyxml2.h"
#include <boost/filesystem.hpp>

CompilerConfig& CompilerConfig::Instance()
{
	static CompilerConfig instance;
	return instance;
}

CompilerConfig::CompilerConfig(const CompilerConfig& root) : kerf_(3.175), kerf2_(1.5875), randomPruning_(false),
randomSeed_(0),
maximalPerm_(10), opt_popSize_(40), opt_nIslands_(1), opt_evoGens_(30), opt_nIterations_(1),opt_nTerminal_(-1)
{
	RegenerateTempFolder();

	if (!hModule)
	{
		hModule = LoadLibrary(_T("carpentry_geom.dll"));
		if (!hModule) {
			DWORD dw = GetLastError(); // returns 0xc1 (193)
			std::cerr << "LoadLibrary failed with error code " << dw << "\n";
			system("pause");
		}
		else
			std::cerr << "LoadLibrary success\n";
	}
}

CompilerConfig::CompilerConfig() : kerf_(3.175), kerf2_(1.5875), randomPruning_(false), randomSeed_(0),
maximalPerm_(10), opt_popSize_(40), opt_nIslands_(1), opt_evoGens_(30), opt_nIterations_(1), opt_nTerminal_(-1)
{
	RegenerateTempFolder();

	if (!hModule)
	{
		hModule = LoadLibrary(_T("carpentry_geom.dll"));
		if (!hModule) {
			DWORD dw = GetLastError(); // returns 0xc1 (193)
			std::cerr << "LoadLibrary failed with error code " << dw << "\n";
			system("pause");
		}
		else
			std::cerr << "LoadLibrary success\n";
	}
}

HMODULE& CompilerConfig::GetHModule()
{
	return hModule;
}

void CompilerConfig::LoadConfigFile(const std::string& file)
{
	try
	{
		tinyxml2::XMLDocument doc;
		
		doc.LoadFile(file.c_str());

		if (doc.Error())
		{
			std::cerr << "Load xml error: "<<file << std::endl;
			std::cerr << doc.ErrorStr() << std::endl;
			system("pause");
		}

		auto root = doc.FirstChildElement("root");
		
		ISA::InitISA();
		for (auto e = root->FirstChildElement("Tool"); e != nullptr; e = e->NextSiblingElement("Tool"))
		{
			const auto toolType = e->Attribute("Type");
			auto process = ISA::CreateProcess(toolType);
			process->SetEnabled(true);
			for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
			{
				const auto paramName = p->Attribute("Name");
				process->SetValue(paramName, strtod(p->GetText(), nullptr));
			}

			isa_.AddProcess(toolType, process);

			
			if (strcmp("Bandsaw", toolType) == 0)
			{
				auto bandsaw = std::make_shared<BandsawConfig>();
				bandsaw->SetEnabled(true);
				for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
				{
					const auto paramName = p->Attribute("Name");
					bandsaw->SetValue(paramName, strtod(p->GetText(), nullptr));
				}
				isa_.AddProcess(toolType, bandsaw);
			}
			else if (strcmp("Chopsaw", toolType) == 0)
			{
				auto chopsaw = std::make_shared<ChopsawConfig>();
				chopsaw->SetEnabled(true);
				for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
				{
					const auto paramName = p->Attribute("Name");
					chopsaw->SetValue(paramName, strtod(p->GetText(), nullptr));
				}
				isa_.AddProcess(toolType, chopsaw);
			}
			else if (strcmp("Jigsaw", toolType) == 0)
			{
				auto jigsaw = std::make_shared<JigsawConfig>();
				jigsaw->SetEnabled(true);
				for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
				{
					const auto paramName = p->Attribute("Name");
					jigsaw->SetValue(paramName, strtod(p->GetText(), nullptr));
				}
				isa_.AddProcess(toolType, jigsaw);
			}
			else if (strcmp("Tracksaw", toolType) == 0)
			{
				auto tracksaw = std::make_shared<TracksawConfig>();
				tracksaw->SetEnabled(true);
				for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
				{
					const auto paramName = p->Attribute("Name");
					tracksaw->SetValue(paramName, strtod(p->GetText(), nullptr));
				}
				isa_.AddProcess(toolType, tracksaw);
			}
			else if (strcmp("Drillpress", toolType) == 0)
			{
				auto drillpress = std::make_shared<DrillPressConfig>();
				drillpress->SetEnabled(true);
				for (auto p = e->FirstChildElement("Para"); p != nullptr; p = p->NextSiblingElement("Para"))
				{
					const auto paramName = p->Attribute("Name");
					drillpress->SetValue(paramName, strtod(p->GetText(), nullptr));
				}
				isa_.AddProcess(toolType, drillpress);
			}
			
		}

		auto params = root->FirstChildElement("Params");

		if (params->FirstChildElement("RandomPruning"))
			randomPruning_ = atoi(params->FirstChildElement("RandomPruning")->GetText());
		if (params->FirstChildElement("RandomSeed"))
			randomSeed_ = atoi(params->FirstChildElement("RandomSeed")->GetText());
		if (params->FirstChildElement("MaximalPerm"))
		{
			maximalPerm_ = atoi(params->FirstChildElement("MaximalPerm")->GetText());
			if (maximalPerm_ < 0) maximalPerm_ = 1000;
		}
		if (params->FirstChildElement("Kerf"))
			kerf_ = atof(params->FirstChildElement("Kerf")->GetText());
		kerf2_ = kerf_ / 2.0;
		if (params->FirstChildElement("PartMatchError"))
			partMatchError_ = atof(params->FirstChildElement("PartMatchError")->GetText());
		if (params->FirstChildElement("AngleMatchError"))
			angleMatchError_ = atof(params->FirstChildElement("AngleMatchError")->GetText());
		if (params->FirstChildElement("MaximalCuttingDistance"))
			maximalCuttingDistance_ = atof(params->FirstChildElement("MaximalCuttingDistance")->GetText());
		if (params->FirstChildElement("SafeDistance"))
			safeDistance_ = atof(params->FirstChildElement("SafeDistance")->GetText());
		if (params->FirstChildElement("TotalOrderNumber"))
		{
			totalOrderNumber_ = atoi(params->FirstChildElement("TotalOrderNumber")->GetText());
			if (totalOrderNumber_ < 0)totalOrderNumber_ = 500;
		}
		if (params->FirstChildElement("PruningLimitation"))
		{
			pruningLimitation_ = atoi(params->FirstChildElement("PruningLimitation")->GetText());
			if (pruningLimitation_ < 0)pruningLimitation_ = 1000;
		}
		if (params->FirstChildElement("EGraphFolder"))
			egraph_o_Folder_ = params->FirstChildElement("EGraphFolder")->GetText();
		if (params->FirstChildElement("PackingFolder"))
			packing_o_Folder_ = params->FirstChildElement("PackingFolder")->GetText();
		if (params->FirstChildElement("MainObj"))
			main_obj_ = params->FirstChildElement("MainObj")->GetText();

		if (main_obj_.empty() || egraph_o_Folder_.empty())
		{
			main_obj_ = packing_o_Folder_ + "\\iges";
			egraph_o_Folder_ = packing_o_Folder_ + "\\e-graph";
			packing_o_Folder_ = packing_o_Folder_ + "\\packing";
		}

		//AngleMatchError
		//need to improve : detect valid 
		if (params->FirstChildElement("OptPopSize")) 
			opt_popSize_ = atoi(params->FirstChildElement("OptPopSize")->GetText());
		if (params->FirstChildElement("OptNIslands")) 
			opt_nIslands_ = atoi(params->FirstChildElement("OptNIslands")->GetText());
		if (params->FirstChildElement("OptEvoGens"))
		{
			opt_evoGens_ = atoi(params->FirstChildElement("OptEvoGens")->GetText());
			if (opt_evoGens_ < 0)opt_evoGens_ = 100;
		}
		if (params->FirstChildElement("OptNIterations")) 
			opt_nIterations_ = atoi(params->FirstChildElement("OptNIterations")->GetText());
		if (params->FirstChildElement("OptNTerminal"))
			opt_nTerminal_ = atoi(params->FirstChildElement("OptNTerminal")->GetText());

		if (params->FirstChildElement("RunDesign")) 
			run_design_ = atoi(params->FirstChildElement("RunDesign")->GetText());
		if (params->FirstChildElement("RunPack")) 
			run_pack_ = atoi(params->FirstChildElement("RunPack")->GetText());
		if (params->FirstChildElement("RunOrder"))
			run_order_ = atoi(params->FirstChildElement("RunOrder")->GetText());
		if (params->FirstChildElement("RunOpt")) 
			run_opt_ = atoi(params->FirstChildElement("RunOpt")->GetText());

		if (params->FirstChildElement("InitGeneNb"))
			init_gene_nb_ = atoi(params->FirstChildElement("InitGeneNb")->GetText());
		if (params->FirstChildElement("ExplorationNb"))
			exploration_nb_ = atoi(params->FirstChildElement("ExplorationNb")->GetText());
		if (params->FirstChildElement("ExplorationTerminalNb"))
			expoloration_terminal_nb_ = atoi(params->FirstChildElement("ExplorationTerminalNb")->GetText());

		if (params->FirstChildElement("ExplorationM"))
			exploration_m_ = atof(params->FirstChildElement("ExplorationM")->GetText());
		if (params->FirstChildElement("ExplorationCR"))
			exploration_cr_ = atof(params->FirstChildElement("ExplorationCR")->GetText());

		//expoloration_terminal_nb_

		if (params->FirstChildElement("SharedEGraph"))
			shared_egraph_ = atoi(params->FirstChildElement("SharedEGraph")->GetText());

		if (params->FirstChildElement("OBJOutput"))
			obj_output = atoi(params->FirstChildElement("OBJOutput")->GetText());

		//t_bottom_up_select
		if (params->FirstChildElement("TBottomUpSelect"))
			t_bottom_up_select = atoi(params->FirstChildElement("TBottomUpSelect")->GetText());

		if (params->FirstChildElement("TBranchNB"))
			t_branch_nb = atoi(params->FirstChildElement("TBranchNB")->GetText());
		if (params->FirstChildElement("TRandDelete"))
			t_rand_delete = atoi(params->FirstChildElement("TRandDelete")->GetText());
		if (params->FirstChildElement("TShuffleSelection"))
			t_shuffle_selection = atoi(params->FirstChildElement("TShuffleSelection")->GetText());
		if (params->FirstChildElement("TAddShuffle"))
			t_add_shuffle = atoi(params->FirstChildElement("TAddShuffle")->GetText());
		if (params->FirstChildElement("TExactNSGANB"))
			t_exact_nsga_nb = atoi(params->FirstChildElement("TExactNSGANB")->GetText());
		if (params->FirstChildElement("TAddFrontNB"))
			t_add_front_nb = atoi(params->FirstChildElement("TAddFrontNB")->GetText());

		if (params->FirstChildElement("TDAIFurther"))
			t_dai_further = atoi(params->FirstChildElement("TDAIFurther")->GetText());
		//t_dai_further

		if (params->FirstChildElement("TInitFab"))
			t_init_fab = atoi(params->FirstChildElement("TInitFab")->GetText());

		if (params->FirstChildElement("TInitDepth"))
			t_init_depth = atoi(params->FirstChildElement("TInitDepth")->GetText());

		if (params->FirstChildElement("TInitGa"))
			t_init_ga = atoi(params->FirstChildElement("TInitGa")->GetText());

		if (params->FirstChildElement("TInitSimilarEx"))
			t_init_similar_ex = atoi(params->FirstChildElement("TInitSimilarEx")->GetText());

		if (params->FirstChildElement("TInitDVS"))
			t_init_dvs = atoi(params->FirstChildElement("TInitDVS")->GetText());

		//pipeline_selection
		if (params->FirstChildElement("PipelineSelection"))
			pipeline_selection = atoi(params->FirstChildElement("PipelineSelection")->GetText());

		//ref_point_
		if (params->FirstChildElement("RefPoint"))
		{
			auto rp = Math::Functs::SplitD(params->FirstChildElement("RefPoint")->GetText(),' ');
			ref_point_ = Vector3d{ rp[0], rp[1], rp[2] };
		}

		for (auto e = root->FirstChildElement("Stock"); e != nullptr; e = e->NextSiblingElement("Stock"))
		{
			PartDesignGui::WLCollection wlc(wlcs_.size());
			
			wlcs_.emplace_back(wlc);
			for (auto p = e->FirstChildElement("WLC"); p != nullptr; p = p->NextSiblingElement("WLC"))
			{
				auto name = p->Attribute("Name");
				auto L = strtod(p->Attribute("L"), nullptr);
				auto W = strtod(p->Attribute("W"), nullptr);
				auto H = strtod(p->Attribute("H"), nullptr);
				auto UID = p->GetText();
				wlcs_.back().wls.emplace_back(PartDesignGui::WL(0, name, L, W, H, UID, kerf_));
			}
			wlcs_.back().SortWLS();
		}

		hasLoaded_ = true;
	}
	catch (...)
	{
		hasLoaded_ = false;
	}
}

void CompilerConfig::SaveConfigFile(const std::string& file)
{
}

bool CompilerConfig::HasLoaded() const
{
	return hasLoaded_;
}

double CompilerConfig::GetKerf() const
{
	return kerf_;
}

double CompilerConfig::GetKerf2() const
{
	return kerf2_;
}

bool CompilerConfig::IsRandomPruning() const
{
	return randomPruning_;
}

int CompilerConfig::GetRandomSeed() const
{
	return randomSeed_;
}

double CompilerConfig::GetPartMatchError() const
{
	return partMatchError_;
}

double CompilerConfig::GetAngleMatchError() const
{
	return angleMatchError_;
}

double CompilerConfig::GetMaximalCuttingDistance() const
{
	return maximalCuttingDistance_;
}

double CompilerConfig::GetSafeDistance() const
{
	return safeDistance_;
}

int CompilerConfig::GetTotalOrderNumber() const
{
	return totalOrderNumber_;
}

int CompilerConfig::GetPruningLimitation() const
{
	return pruningLimitation_;
}

CompilerConfig::VecWLC& CompilerConfig::GetVecWLC()
{
	return wlcs_;
}

CompilerConfig::VecShape& CompilerConfig::GetVecShape()
{
	return shapes_;
}

const std::string& CompilerConfig::GetEGraphOutputFolder() const
{
	return egraph_o_Folder_;
}

const std::string& CompilerConfig::GetPackingOutputFolder() const
{
	return packing_o_Folder_;
}

const std::string& CompilerConfig::GetTemporarilyFolder() const
{
	return tempFolder;
}

std::shared_ptr<Process> CompilerConfig::GetProcess(const std::string& name)
{
	return isa_.GetProcess(name);
}

bool CompilerConfig::HasProcess(const std::string& name)
{
	return isa_.HasProcess(name);
}

int CompilerConfig::GetMaximalPerm() const
{
	return maximalPerm_;
}

void CompilerConfig::WriteDefaultConfig(const std::string& filename)
{
	//auto fileContent = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><root><Tool Type = \"Bandsaw\"></Tool><Tool Type=\"Chopsaw\"><Para Name = \"CHOPSAW_MIN_X_PRIME\">-1e-4</Para></Tool><Tool Type = \"Jigsaw\"></Tool><Tool Type = \"Tracksaw\"></Tool><!-- <Tool Type = \"Drillpress\">	</Tool> -->	<Stock>	<WLC Name = \"wood_0.75x12x20\" L=\"304.8\" W=\"508.0\" H=\"19.05\">8538f8db-adb6-4b5d-8ec0-6d6bad325ac4</WLC>	<WLC Name = \"wood_0.75x24x20\" L=\"609.6\" W=\"508.0\" H=\"19.05\">786c393f-b7fa-4910-9fc7-1492021f1a77</WLC>	<!-- <WLC Name = \"wood_2x2\" L=\"762.0\" W=\"762.0\" H=\"19.05\">b024b206-900c-4e9a-a934-00293b4bd186</WLC> -->	</Stock>	<Stock>	<WLC Name = \"wood_0.5x12x20\" L=\"304.8\" W=\"609.6\" H=\"12.7\">af330a8f-db6b-46ee-9427-2af4ac227167</WLC>	<WLC Name = \"wood_0.5x24x20\" L=\"609.6\" W=\"609.6\" H=\"12.7\">44335d53-2916-4f71-81b3-2a5d76d71fff</WLC>	<!-- <WLC Name = \"wood_2x2\" L=\"1219.2\" W=\"609.6\" H=\"12.7\">b024b206-900c-4e9a-a934-00293b4bd186</WLC> -->	</Stock>	<Stock>	<WLC Name = \"wood_1x12x20\" L=\"304.8\" W=\"508.0\" H=\"38.1\">e49f6c42-16fe-4b9c-8198-0469f33182da</WLC>	<WLC Name = \"wood_1x24x20\" L=\"609.6\" W=\"508.0\" H=\"38.1\">aac77b14-f1cf-4950-b48c-c4d9e012fddb</WLC>	</Stock>	<Stock>	<WLC Name = \"lumber_24x2x2\" L=\"609.6\" W=\"38.1\" H=\"38.1\">15411ad4-b2b1-4dbe-b5f3-f6e5df128d31</WLC>	<WLC Name = \"lumber_48x2x2\" L=\"1219.2\" W=\"38.1\" H=\"38.1\">d01ecf45-d768-404c-8518-2dae6444cdff</WLC>	<WLC Name = \"lumber_96x2x2\" L=\"2438.4\" W=\"38.1\" H=\"38.1\">1df1a8cb-d30e-4c7f-93ab-7de601a610ad</WLC>	</Stock>	<Stock>	<WLC Name = \"lumber_24x2x4\" L=\"609.6\" W=\"88.9\" H=\"38.1\">b054fd41-d023-4cde-b5ca-cf6356ec6b14</WLC>	<WLC Name = \"lumber_48x2x4\" L=\"1219.2\" W=\"88.9\" H=\"38.1\">32b733c8-9bb2-4c2b-b3af-4844992d5807</WLC>	<WLC Name = \"lumber_96x2x4\" L=\"2438.4\" W=\"88.9\" H=\"38.1\">805da52d-c8e4-4ec6-8cfe-8a304f52ff24</WLC>	</Stock>	<Stock>	<WLC Name = \"lumber_24x2x4\" L=\"609.6\" W=\"88.9\" H=\"88.9\">764917a9-4b5a-4e52-a402-5cdfc7b0d055</WLC>	<WLC Name = \"lumber_48x2x4\" L=\"1219.2\" W=\"88.9\" H=\"88.9\">906eb3e2-1df2-4764-a150-e736a3a37483</WLC>	<WLC Name = \"lumber_96x2x4\" L=\"2438.4\" W=\"88.9\" H=\"88.9\">7eaa6d88-3d87-4afa-a5f5-2d0cae939480</WLC>	</Stock>	<Stock>	<WLC Name = \"lumber_24x2x8\" L=\"609.6\" W=\"184.15\" H=\"38.1\">652f5cd5-9e6e-48fa-b3f8-176f10fcdb56</WLC>	<WLC Name = \"lumber_48x2x8\" L=\"1219.2\" W=\"184.15\" H=\"38.1\">a032d738-d738-4880-8f50-36f4250db666</WLC>	<WLC Name = \"lumber_96x2x8\" L=\"2438.4\" W=\"184.15\" H=\"38.1\">b024b206-900c-4e9a-a934-00293b4bd186</WLC>	</Stock>	<Params>	<RandomPruning>0 </RandomPruning>	<RandomSeed>0 </RandomSeed>	<MaximalPerm>25 </MaximalPerm>	<Kerf>3.175 </Kerf>	<PartMatchError>0.01 </PartMatchError>	<MaximalCuttingDistance>200.0 </MaximalCuttingDistance>	<SafeDistance>2.54 </SafeDistance>	<TotalOrderNumber>50 </TotalOrderNumber>	<PruningLimitation>10</PruningLimitation><VisualizationFile>D:\\vis.txt</VisualizationFile><VisualizationOutputFile>D:\\vis\\vis</VisualizationOutputFile><EGraphFolder>D:\\e-graph</EGraphFolder><PackingFolder>D:\\packing</PackingFolder></Params></root>";
	std::ofstream file(filename);
	//file << fileContent;
	file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?><root>" << std::endl;
	file << "" << std::endl;
	file << "<Tool Type = \"Bandsaw\"></Tool>" << std::endl;
	file << "<Tool Type=\"Chopsaw\"><Para Name = \"CHOPSAW_MIN_X_PRIME\">-1e-4</Para></Tool>" << std::endl;
	file << "<Tool Type = \"Jigsaw\"></Tool>" << std::endl;
	file << "<Tool Type = \"Tracksaw\"></Tool>" << std::endl;
	file << "<!-- <Tool Type = \"Drillpress\">	</Tool> -->	" << std::endl;
	file << "" << std::endl;
	file << "<Stock>" << std::endl;
	file << "	<WLC Name = \"wood_0.75x12x20\" L=\"304.8\" W=\"508.0\" H=\"19.05\">8538f8db-adb6-4b5d-8ec0-6d6bad325ac4</WLC>	" << std::endl;
	file << "	<WLC Name = \"wood_0.75x24x20\" L=\"609.6\" W=\"508.0\" H=\"19.05\">786c393f-b7fa-4910-9fc7-1492021f1a77</WLC>	" << std::endl;
	file << "	<!-- <WLC Name = \"wood_2x2\" L=\"762.0\" W=\"762.0\" H=\"19.05\">b024b206-900c-4e9a-a934-00293b4bd186</WLC> -->" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "	" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"wood_0.5x12x20\" L=\"304.8\" W=\"609.6\" H=\"12.7\">af330a8f-db6b-46ee-9427-2af4ac227167</WLC>" << std::endl;
	file << "	<WLC Name = \"wood_0.5x24x20\" L=\"609.6\" W=\"609.6\" H=\"12.7\">44335d53-2916-4f71-81b3-2a5d76d71fff</WLC>	" << std::endl;
	file << "	<!-- <WLC Name = \"wood_2x2\" L=\"1219.2\" W=\"609.6\" H=\"12.7\">b024b206-900c-4e9a-a934-00293b4bd186</WLC> -->	" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "	" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"wood_1x12x20\" L=\"304.8\" W=\"508.0\" H=\"38.1\">e49f6c42-16fe-4b9c-8198-0469f33182da</WLC>" << std::endl;
	file << "	<WLC Name = \"wood_1x24x20\" L=\"609.6\" W=\"508.0\" H=\"38.1\">aac77b14-f1cf-4950-b48c-c4d9e012fddb</WLC>" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "	" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"lumber_24x2x2\" L=\"609.6\" W=\"38.1\" H=\"38.1\">15411ad4-b2b1-4dbe-b5f3-f6e5df128d31</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_48x2x2\" L=\"1219.2\" W=\"38.1\" H=\"38.1\">d01ecf45-d768-404c-8518-2dae6444cdff</WLC>" << std::endl;
	file << "	<WLC Name = \"lumber_96x2x2\" L=\"2438.4\" W=\"38.1\" H=\"38.1\">1df1a8cb-d30e-4c7f-93ab-7de601a610ad</WLC>" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"lumber_24x2x4\" L=\"609.6\" W=\"88.9\" H=\"38.1\">b054fd41-d023-4cde-b5ca-cf6356ec6b14</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_48x2x4\" L=\"1219.2\" W=\"88.9\" H=\"38.1\">32b733c8-9bb2-4c2b-b3af-4844992d5807</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_96x2x4\" L=\"2438.4\" W=\"88.9\" H=\"38.1\">805da52d-c8e4-4ec6-8cfe-8a304f52ff24</WLC>	" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"lumber_24x2x4\" L=\"609.6\" W=\"88.9\" H=\"88.9\">764917a9-4b5a-4e52-a402-5cdfc7b0d055</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_48x2x4\" L=\"1219.2\" W=\"88.9\" H=\"88.9\">906eb3e2-1df2-4764-a150-e736a3a37483</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_96x2x4\" L=\"2438.4\" W=\"88.9\" H=\"88.9\">7eaa6d88-3d87-4afa-a5f5-2d0cae939480</WLC>	" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "" << std::endl;
	file << "<Stock>	" << std::endl;
	file << "	<WLC Name = \"lumber_24x2x8\" L=\"609.6\" W=\"184.15\" H=\"38.1\">652f5cd5-9e6e-48fa-b3f8-176f10fcdb56</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_48x2x8\" L=\"1219.2\" W=\"184.15\" H=\"38.1\">a032d738-d738-4880-8f50-36f4250db666</WLC>	" << std::endl;
	file << "	<WLC Name = \"lumber_96x2x8\" L=\"2438.4\" W=\"184.15\" H=\"38.1\">b024b206-900c-4e9a-a934-00293b4bd186</WLC>" << std::endl;
	file << "</Stock>	" << std::endl;
	file << "" << std::endl;
	file << "<Params>	" << std::endl;
	file << "" << std::endl;
	file << "	<RandomPruning>0 </RandomPruning>	" << std::endl;
	file << "	<RandomSeed>0 </RandomSeed>	" << std::endl;
	file << "	<MaximalPerm>25 </MaximalPerm>" << std::endl;
	file << "	<Kerf>3.175 </Kerf>" << std::endl;
	file << "	<PartMatchError>0.01 </PartMatchError>" << std::endl;
	file << "	<AngleMatchError>0.08 </AngleMatchError>" << std::endl;
	file << "	<MaximalCuttingDistance>200.0 </MaximalCuttingDistance>" << std::endl;
	file << "	<SafeDistance>2.54 </SafeDistance>" << std::endl;
	file << "	<TotalOrderNumber>50 </TotalOrderNumber>" << std::endl;
	file << "	<PruningLimitation>10</PruningLimitation>" << std::endl;
	file << "	<OptPopSize>40</OptPopSize>" << std::endl;
	file << "	<OptNIslands>1</OptNIslands>" << std::endl;
	file << "	<OptEvoGens>30</OptEvoGens>" << std::endl;
	file << "	<OptNIterations>1</OptNIterations>" << std::endl;
	file << "	<VisualizationFile>D:\\vis.txt</VisualizationFile>" << std::endl;
	file << "	<VisualizationOutputFile>D:\\vis\\vis</VisualizationOutputFile>" << std::endl;
	file << "	<EGraphFolder>D:\\e-graph</EGraphFolder>" << std::endl;
	file << "	<PackingFolder>D:\\packing</PackingFolder>" << std::endl;
	file << "" << std::endl;
	file << "</Params>" << std::endl;
	file << "</root>" << std::endl;


	file.close();
}

void CompilerConfig::RegenerateTempFolder()
{
	boost::filesystem::path temp = boost::filesystem::unique_path();
	tempFolder = temp.string();
}

int CompilerConfig::GetOptPopSize() const
{
	return opt_popSize_;
}
int CompilerConfig::GetOptNIslands() const
{
	return opt_nIslands_;
}
int CompilerConfig::GetOptEvoGens() const
{
	return opt_evoGens_;
}
int CompilerConfig::GetOptNIterations() const
{
	return opt_nIterations_;
}

int CompilerConfig::GetOptNTerminal() const
{
	return opt_nTerminal_;
}

const std::string& CompilerConfig::GetMainObj() const
{
	return main_obj_;
}


const bool CompilerConfig::GetRunDesign() const
{
	return run_design_;
}

const bool CompilerConfig::GetObjOutput() const 
{
	return obj_output;
}

const Vector3d CompilerConfig::GetRefPoint() const
{
	return ref_point_;
}

const int CompilerConfig::GetTBottomUpSelect() const
{
	return t_bottom_up_select;
}

const bool CompilerConfig::GetRunPack() const
{
	return run_pack_;
}
const bool CompilerConfig::GetRunOrder() const
{
	return run_order_;
}
const bool CompilerConfig::GetRunOpt() const
{
	return run_opt_;
}

const bool CompilerConfig::GetSharedEgraph() const
{
	return shared_egraph_;
}

const int CompilerConfig::GetInitGeneNb() const
{
	return init_gene_nb_;
}

const int CompilerConfig::GetExplorationNb() const
{
	return exploration_nb_;
}

const int CompilerConfig::GetExplorationTerminalNb() const
{
	return expoloration_terminal_nb_;
}

const double CompilerConfig::GetExplorationM() const
{
	return exploration_m_;
}

const double CompilerConfig::GetExplorationCR() const
{
	return exploration_cr_;
}

const int CompilerConfig::GetTBranchNb()const
{
	return t_branch_nb;
}
const int CompilerConfig::GetTRandDelete()const
{
	return t_rand_delete;
}
const int CompilerConfig::GetTShuffleSelection() const
{
	return t_shuffle_selection;
}
const int CompilerConfig::GetTAddShuffle()const
{
	return t_add_shuffle;
}

const int CompilerConfig::GetTExactNSGANb() const
{
	return t_exact_nsga_nb;
}

const int CompilerConfig::GetTAddFrontNb() const
{
	return t_add_front_nb;
}

const int CompilerConfig::GetTDAIFurther() const
{
	return t_dai_further;
}

const int CompilerConfig::GetPipelineSelection() const
{
	return pipeline_selection;
}

const int CompilerConfig::GetTInitDVS() const
{
	return t_init_dvs;
}

const bool CompilerConfig::GetTInitFab() const
{
	return t_init_fab;
}

const bool CompilerConfig::GetTInitDepth() const
{
	return t_init_depth;
}

const bool CompilerConfig::GetTInitSimilarEx() const
{
	return t_init_similar_ex;
}

const bool CompilerConfig::GetTInitGa() const
{
	return t_init_ga;
}
