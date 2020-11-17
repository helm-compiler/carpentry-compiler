#ifndef PAGMO_ALGORITHMS_NSGA2_HPP
#define PAGMO_ALGORITHMS_NSGA2_HPP

#include <algorithm> // std::shuffle, std::transform
#include <iomanip>
#include <numeric> // std::iota, std::inner_product
#include <random>
#include <string>
#include <tuple>
#include <sstream>
#include <egraph.hpp>
#include <fstream>

#include <pagmo/algorithm.hpp> // needed for the cereal macro
#include <pagmo/exceptions.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/decompose.hpp>
#include <pagmo/rng.hpp>
#include <pagmo/utils/generic.hpp>         // uniform_real_from_range, some_bound_is_equal
#include <pagmo/utils/multi_objective.hpp> // crowding_distance, etc..
#include <pagmo/utils/hv_algos/hv_algorithm.hpp>
#include <pagmo/utils/hv_algos/hv_hv3d.hpp>
#include <pagmo/utils/hypervolume.hpp>
#include <evaluator.hpp>
#include <config.hpp>

namespace OptCompiler
{
	using namespace pagmo;

	struct RecordPop
	{
		std::vector<vector_double> m_f;
		std::vector<std::vector<OptCompiler::EncodedVec>> m_x;

		RecordPop(population pop)
		{
			UpdateRecordPop(pop);
		}

		void Clear()
		{
			std::vector<vector_double>().swap(m_f);
			std::vector<std::vector<OptCompiler::EncodedVec>>().swap(m_x);
		}

		int XNB() 
		{
			int nb = 0;
			for (auto& m_x_ : m_x)
			{
				nb += m_x_.size();
			}
			return nb;
		};

		void Push(const OptCompiler::EncodedVec& x, const vector_double& f)
		{
			for(int i=0;i<m_f.size();i++)
			{
				if (m_f[i][0] == f[0] && m_f[i][1] == f[1] && m_f[i][2] == f[2])
				{
					auto x_str = EncodedVec::EncodeStr(x);
					for (int j = 0; j < m_x[i].size(); j++)
						if (x_str == EncodedVec::EncodeStr(m_x[i][j])) return;
					m_x[i].emplace_back(x);
					return;
				}
			}
			m_f.emplace_back(f);
			m_x.emplace_back(std::vector<OptCompiler::EncodedVec>(1,x));
		}

		void UpdateRecordPop(const population& pop)
		{
			for (int i = 0; i < pop.get_x().size(); i++) Push(pop.get_x()[i], pop.get_f()[i]);

			if (m_f.size() >= 2)
			{
				auto ndf = std::get<0>(fast_non_dominated_sorting(m_f))[0];
				std::vector<vector_double> m_f_;
				std::vector<std::vector<OptCompiler::EncodedVec>> m_x_;
				for (auto& pf : ndf) m_f_.emplace_back(m_f[pf]);
				for (auto& pf : ndf) m_x_.emplace_back(m_x[pf]);
				m_f = m_f_;
				m_x = m_x_;
				std::vector<vector_double>().swap(m_f_);
				std::vector<std::vector<OptCompiler::EncodedVec>>().swap(m_x_);
			}
			//compute difference

		};

		int own_best_idx(population& old_pop, const population &new_pop, vector_double::size_type N)
		{
			UpdateRecordPop(new_pop);

			int xnb = XNB();

			vector_double cds(crowding_distance(m_f));

			population best_pop;

			if (xnb >= N)
			{
				//crowding distance
				vector_double cds(crowding_distance(m_f));

				// We now have front and crowding distance, we sort the front w.r.t. the crowding
				std::vector<vector_double::size_type> idxs(m_f.size());
				std::iota(idxs.begin(), idxs.end(), vector_double::size_type(0u));
				std::sort(idxs.begin(), idxs.end(), [&cds](vector_double::size_type idx1, vector_double::size_type idx2) {
					return detail::greater_than_f(cds[idx1], cds[idx2]);}); // Descending order1
		
				for (auto& idx : idxs)
				{
					for (auto& x : m_x[idx])
					{
						if (best_pop.size() >= N) goto XNBLabel;
						best_pop.push_back(x, m_f[idx]);
					}
					continue;
					XNBLabel: break;
				}
			}
			else
			{
				for (int i = 0; i < m_f.size(); i++)
					for (int j = 0; j < m_x[i].size(); j++)
						best_pop.push_back(m_x[i][j], m_f[i]);
				
				population temp_pop;
				for (int i = 0; i < old_pop.size(); i++)
				{
					auto temp_str = EncodedVec::EncodeStr(old_pop.get_x()[i]);
					for (int j = 0; j < best_pop.size(); j++)
						if (temp_str == EncodedVec::EncodeStr(best_pop.get_x()[j]))
							goto OldJumpLabel;
					temp_pop.push_back(old_pop.get_x()[i], old_pop.get_f()[i]);
					OldJumpLabel: continue;
				}
				for (int i = 0; i < new_pop.size(); i++)
				{
					auto temp_str = EncodedVec::EncodeStr(new_pop.get_x()[i]);
					for (int j = 0; j < best_pop.size(); j++) 
						if (temp_str == EncodedVec::EncodeStr(best_pop.get_x()[j]))
							goto NewJumpLabel;
					temp_pop.push_back(new_pop.get_x()[i], new_pop.get_f()[i]);
					NewJumpLabel: continue;
				}
	
				if (temp_pop.get_x().size() > N - xnb)
				{
					auto best_idx = select_best_N_mo(temp_pop.get_f(), N - xnb);
					for (int i = 0; i < best_idx.size(); i++)
						best_pop.push_back(temp_pop.get_x()[best_idx[i]], temp_pop.get_f()[best_idx[i]]);
				}
				else
				{
					for (int i = 0; i < temp_pop.size(); i++)
						best_pop.push_back(temp_pop.get_x()[i], temp_pop.get_f()[i]);

					while (true)
					{
						for (int i = 0; i < m_f.size(); i++)
						{
							for (int j = 0; j < m_x[i].size(); j++)
							{
								int d = best_pop.size();
								if (best_pop.size() >= N) goto JumpOut;
								best_pop.push_back(m_x[i][j], m_f[i]);
							}
						}
						continue;
						JumpOut: break;;
					}
				}
			}

			//compute different

			std::vector<std::string> old_pop_strs;
			for (population::size_type i = 0; i < N; ++i)
				old_pop_strs.emplace_back(old_pop.get_x()[i].EncodeStr(old_pop.get_x()[i]));

			int different_nb = 0;
			for (population::size_type i = 0; i < N; ++i)
			{
				auto str = best_pop.get_x()[i].EncodeStr(best_pop.get_x()[i]);
				if (std::find(old_pop_strs.begin(), old_pop_strs.end(), str) == old_pop_strs.end())
					different_nb++;
			}

			for (population::size_type i = 0; i < N; ++i)
			{
				old_pop.set_xf(i, best_pop.get_x()[i], best_pop.get_f()[i]);
			}

			return different_nb;
		};

	};

	/// Nondominated Sorting genetic algorithm II (NSGA-II)
	class nsga2_egraph
	{
	public:
		/// Single entry of the log (gen, fevals, ideal_point)
		typedef std::tuple<unsigned int, unsigned long long, vector_double> log_line_type;
		/// The log
		typedef std::vector<log_line_type> log_type;

		/// Constructor
		/**
		 * Constructs the NSGA II user defined algorithm.
		 *
		 * @param[in] gen Number of generations to evolve.
		 * @param[in] cr Crossover probability.
		 * @param[in] eta_c Distribution index for crossover.
		 * @param[in] m Mutation probability.
		 * @param[in] eta_m Distribution index for mutation.
		 * @param seed seed used by the internal random number generator (default is random)
		 * @throws std::invalid_argument if \p cr is not \f$ \in [0,1[\f$, \p m is not \f$ \in [0,1]\f$, \p eta_c is not in
		 * [1,100[ or \p eta_m is not in [1,100[.
		 */
		nsga2_egraph(unsigned gen = 1u, double cr = 0.95, double eta_c = 10., double m = 0.4, double eta_m = 50.,
		             unsigned seed = pagmo::random_device::next())
			: m_gen(gen), m_cr(cr), m_eta_c(eta_c), m_m(m), m_eta_m(eta_m), m_e(seed), m_seed(seed), m_verbosity(0u),
			  m_log()
		{
			seed = Config::get_instance().GetRandSeed();
			set_seed(seed);
			if (cr >= 1. || cr < 0.)
			{
				pagmo::detail::ex_thrower<std::invalid_argument>(__FILE__, __LINE__, __func__)(
					"The crossover probability must be in the [0,1] range, while a value of " + std::to_string(cr) +
					" was detected");
			}
			if (m < 0. || m > 1.)
			{
				pagmo::detail::ex_thrower<std::invalid_argument>(__FILE__, __LINE__, __func__)(
					"The mutation probability must be in the [0,1] range, while a value of " + std::to_string(cr) +
					" was detected");
			}
			if (eta_c < 1. || eta_c > 100.)
			{
				pagmo::detail::ex_thrower<std::invalid_argument>(__FILE__, __LINE__, __func__)(
					"The distribution index for crossover must be in [1, 100], while a value of " + std::
					to_string(eta_c) + " was detected");
			}
			if (eta_m < 1. || eta_m > 100.)
			{
				pagmo::detail::ex_thrower<std::invalid_argument>(__FILE__, __LINE__, __func__)(
					"The distribution index for mutation must be in [1, 100], while a value of " + std::to_string(eta_m)
					+ " was detected");
			}

			m_method = hv3d().clone();
		}

		/// Algorithm evolve method (juice implementation of the algorithm)
		/**
		 *
		 * Evolves the population for the requested number of generations.
		 *
		 * @param pop population to be evolved
		 * @return evolved population
		 * @throw std::invalid_argument if pop.get_problem() is stochastic, single objective or has non linear constraints.
		 * If \p int_dim is larger than the problem dimension. If the population size is smaller than 5 or not a multiple of
		 * 4.
		 */
		population evolve(population pop) const
		{
			// We store some useful variables
			const auto& prob = pop.get_problem();
			// This is a const reference, so using set_seed for example will not be
			// allowed
			// This getter does not return a const reference but a copy
			auto NP = pop.size();

			auto fevals0 = prob.get_fevals(); // discount for the fevals already made
			unsigned int count = 1u; // regulates the screen output

			std::vector<float> vechv;
			std::vector<int> pfnb;
			std::vector<int> popDifNb;
			std::vector<int> unipfnb;
			std::vector<int> arrnb;

			// PREAMBLE-------------------------------------------------------------------------------------------------
			// We start by checking that the problem is suitable for this
			// particular algorithm.
			if (NP < 5u || (NP % 4 != 0u))
			{
				pagmo::detail::ex_thrower<std::invalid_argument>(__FILE__, __LINE__, __func__)(
					"for NSGA-II at least 5 individuals in the population are needed and the "
					"population size must be a multiple of 4. Detected input population size is: " + std::
					to_string(NP));
			}
			// ---------------------------------------------------------------------------------------------------------

			// No throws, all valid: we clear the logs
			m_log.clear();


			int same_HV_tm=1;


			// Declarations
			std::vector<vector_double::size_type> best_idx(NP), shuffle1(NP), shuffle2(NP);
			vector_double::size_type parent1_idx, parent2_idx;
			EncodedVec child1, child2;

			std::iota(shuffle1.begin(), shuffle1.end(), 0u);
			std::iota(shuffle2.begin(), shuffle2.end(), 0u);

			RecordPop record_pop(pop);

			// Main NSGA-II loop

			int nT = Config::get_instance().GetNTerminal();
			auto nGen =  Config::get_instance().GetEvoGens();
			for (auto gen = 1u; gen <= nGen; gen++)
			{
				if (nT > 0 && same_HV_tm == nT)break;

				//std::cerr << "Opt: " << gen << " / " << nGen << std::endl;

				// 0 - Logs and prints (verbosity modes > 1: a line is added every m_verbosity generations)
				if (m_verbosity > 0u)
				{
					// Every m_verbosity generations print a log line
					if (gen % m_verbosity == 1u || m_verbosity == 1u)
					{
						// We compute the ideal point
						vector_double ideal_point = ideal(pop.get_f());
						// Every 50 lines print the column names
						if (count % 50u == 1u)
						{
							print("\n", std::setw(7), "Gen:", std::setw(15), "Fevals:");
							for (decltype(ideal_point.size()) i = 0u; i < ideal_point.size(); ++i)
							{
								if (i >= 5u)
								{
									print(std::setw(15), "... :");
									break;
								}
								print(std::setw(15), "ideal" + std::to_string(i + 1u) + ":");
							}
							print('\n');
						}
						print(std::setw(7), gen, std::setw(15), prob.get_fevals() - fevals0);
						for (decltype(ideal_point.size()) i = 0u; i < ideal_point.size(); ++i)
						{
							if (i >= 5u)
							{
								break;
							}
							print(std::setw(15), ideal_point[i]);
						}
						print('\n');
						++count;
						// Logs
						m_log.emplace_back(gen, prob.get_fevals() - fevals0, ideal_point);
					}
				}

				// At each generation we make a copy of the population into popnew
				population popnew;

				// We create some pseudo-random permutation of the population indexes
				std::shuffle(shuffle1.begin(), shuffle1.end(), m_e);
				std::shuffle(shuffle2.begin(), shuffle2.end(), m_e);

				// 1 - We compute crowding distance and non dominated rank for the current population
				auto fnds_res = fast_non_dominated_sorting(pop.get_f());
				auto ndf = std::get<0>(fnds_res); // non dominated fronts [[0,3,2],[1,5,6],[4],...]
				vector_double pop_cd(NP); // crowding distances of the whole population
				auto ndr = std::get<3>(fnds_res); // non domination rank [0,1,0,0,2,1,1, ... ]
				for (const auto& front_idxs : ndf)
				{
					if (front_idxs.size() == 1u)
					{
						// handles the case where the front has collapsed to one point
						pop_cd[front_idxs[0]] = std::numeric_limits<double>::infinity();
					}
					else
					{
						if (front_idxs.size() == 2u)
						{
							// handles the case where the front has collapsed to one point
							pop_cd[front_idxs[0]] = std::numeric_limits<double>::infinity();
							pop_cd[front_idxs[1]] = std::numeric_limits<double>::infinity();
						}
						else
						{
							std::vector<vector_double> front;
							for (auto idx : front_idxs)
							{
								front.push_back(pop.get_f()[idx]);
							}
							auto cd = crowding_distance(front);
							for (decltype(cd.size()) i = 0u; i < cd.size(); ++i)
							{
								pop_cd[front_idxs[i]] = cd[i];
							}
						}
					}
				}

				// 3 - We then loop through all individuals with increment 4 to select two pairs of parents that will
				// each create 2 new offspring
				for (decltype(NP) i = 0u; i < NP; i += 4)
				{
					// We create two offsprings using the shuffled list 1
					parent1_idx = tournament_selection(shuffle1[i], shuffle1[i + 1], ndr, pop_cd);
					parent2_idx = tournament_selection(shuffle1[i + 2], shuffle1[i + 3], ndr, pop_cd);
					crossover(child1, child2, parent1_idx, parent2_idx, pop);
					mutate(child1);
					mutate(child2);
					// we use prob to evaluate the fitness so
					// that its feval counter is correctly updated
					auto f1 = prob.fitness(child1);
					auto f2 = prob.fitness(child2);
					popnew.push_back(child1, f1);
					popnew.push_back(child2, f2);

					// We repeat with the shuffled list 2
					parent1_idx = tournament_selection(shuffle2[i], shuffle2[i + 1], ndr, pop_cd);
					parent2_idx = tournament_selection(shuffle2[i + 2], shuffle2[i + 3], ndr, pop_cd);
					crossover(child1, child2, parent1_idx, parent2_idx, pop);
					mutate(child1);
					mutate(child2);
					// we use prob to evaluate the fitness so
					// that its feval counter is correctly updated
					f1 = prob.fitness(child1);
					f2 = prob.fitness(child2);
					popnew.push_back(child1, f1);
					popnew.push_back(child2, f2);

				} // popnew now contains 2NP individuals

				int diff_nb = record_pop.own_best_idx(pop, popnew, NP);

				popDifNb.emplace_back(diff_nb);

				//hypervolume
				auto hv_obj = hypervolume(record_pop.m_f, true);
				const auto refPnt = Config::get_instance().GetRefpoint();
				float hypvol = hv_obj.compute({refPnt[0], refPnt[1], refPnt[2]}, *m_method); 

				if (!vechv.empty())
				{
					if (Math::Functs::IsAlmostZero(vechv.back() - hypvol))
						same_HV_tm++;
					else
						same_HV_tm = 1;
				}
	
				vechv.push_back(hypvol);


				//record the current iterations
				unipfnb.emplace_back(record_pop.m_f.size());
				pfnb.emplace_back(record_pop.XNB());
				
				arrnb.emplace_back(0);
				for (int i = 0; i < record_pop.m_x.size(); i++)
				{
					for (int j = 0; j < record_pop.m_x[i].size(); j++)
					{
						std::vector<ENode*> progs;
						record_pop.m_x[i][j].GetAllPrograms(record_pop.m_x[i][j].root, progs);
						for (auto& prog : progs)
							if (prog->ecs.empty() && Config::get_instance().Push_Iter_ArrWLC(prog->arr, prog->wlc, gen))
								arrnb.back()++;
					}
				}

				std::cerr << "iterating... " << gen << " / " << nGen << " hypvol: " << hypvol <<" unique_m_f: "<< record_pop.m_f.size()<<" unique_m_x: "<<record_pop.XNB()<<" pop_diff_nb: "<< diff_nb <<" arr_NB: "<< arrnb.back() <<" terminal: "<< same_HV_tm<<"/"<< nT << std::endl;

				////////////////////////////
			} // end of main NSGAII loop

			
			Config::get_instance().AddHypervolume(vechv);
			Config::get_instance().AddParetoFrontsNB(pfnb);
			Config::get_instance().AddUniqueParetoFrontsNB(unipfnb);
			Config::get_instance().AddPopDifNb(popDifNb);
			Config::get_instance().AddArrangeNB(arrnb);

			//popDifNb

			pop = population(pop.get_problem(), record_pop.XNB(), pop.get_seed());;

			for (int i = 0; i < record_pop.m_f.size(); i++)
				for (int j = 0; j < record_pop.m_x[i].size(); j++)
				{
					pop.push_back(record_pop.m_x[i][j], record_pop.m_f[i]);
					break;
				}

			record_pop.Clear();

			//std::cerr << "size of vechv = " << vechv.size() << std::endl;
			return pop;
		}

		/// Sets the seed
		/**
		 * @param seed the seed controlling the algorithm stochastic behavior
		 */
		void set_seed(unsigned int seed)
		{
			m_e.seed(seed);
			m_seed = seed;
		}
		
		/// Gets the seed
		/**
		 * @return the seed controlling the algorithm stochastic behavior
		 */
		unsigned int get_seed() const
		{
			return m_seed;
		}

		/// Sets the algorithm verbosity
		/**
		 * Sets the verbosity level of the screen output and of the
		 * log returned by get_log(). \p level can be:
		 * - 0: no verbosity
		 * - >0: will print and log one line each \p level generations.
		 *
		 * Example (verbosity 1):
		 * @code{.unparsed}
		 * Gen:        Fevals:        ideal1:        ideal2:        ideal3:
		 *   1              0      0.0257554       0.267768       0.974592
		 *   2             52      0.0257554       0.267768       0.908174
		 *   3            104      0.0257554       0.124483       0.822804
		 *   4            156      0.0130094       0.121889       0.650099
		 *   5            208     0.00182705      0.0987425       0.650099
		 *   6            260      0.0018169      0.0873995       0.509662
		 *   7            312     0.00154273      0.0873995       0.492973
		 *   8            364     0.00154273      0.0873995       0.471251
		 *   9            416    0.000379582      0.0873995       0.471251
		 *  10            468    0.000336743      0.0855247       0.432144
		 * @endcode
		 * Gen, is the generation number, Fevals the number of function evaluation used. The ideal point of the current
		 * population follows cropped to its 5th component.
		 *
		 * @param level verbosity level
		 */
		void set_verbosity(unsigned int level)
		{
			m_verbosity = level;
		};
		/// Gets the verbosity level
		/**
		 * @return the verbosity level
		 */
		unsigned int get_verbosity() const
		{
			return m_verbosity;
		}

		/// Algorithm name
		/**
		 * Returns the name of the algorithm.
		 *
		 * @return <tt> std::string </tt> containing the algorithm name
		 */
		std::string get_name() const
		{
			return "NSGA-II:";
		}

		/// Extra informations
		/**
		 * Returns extra information on the algorithm.
		 *
		 * @return an <tt> std::string </tt> containing extra informations on the algorithm
		 */
		std::string get_extra_info() const
		{
			std::ostringstream ss;
			stream(ss, "\tGenerations: ", m_gen);
			stream(ss, "\n\tCrossover probability: ", m_cr);
			stream(ss, "\n\tDistribution index for crossover: ", m_eta_c);
			stream(ss, "\n\tMutation probability: ", m_m);
			stream(ss, "\n\tDistribution index for mutation: ", m_eta_m);
			stream(ss, "\n\tSeed: ", m_seed);
			stream(ss, "\n\tVerbosity: ", m_verbosity);
			return ss.str();
		}

		/// Get log
		/**
		 * A log containing relevant quantities monitoring the last call to evolve. Each element of the returned
		 * <tt>std::vector</tt> is a nsga2::log_line_type containing: Gen, Fevals, ideal_point
		 * as described in nsga2::set_verbosity
		 * @return an <tt>std::vector</tt> of nsga2::log_line_type containing the logged values Gen, Fevals,
		 * ideal_point
		 */
		const log_type& get_log() const
		{
			return m_log;
		}

		/// Object serialization
		/**
		 * This method will save/load \p this into the archive \p ar.
		 *
		 * @param ar target archive.
		 *
		 * @throws unspecified any exception thrown by the serialization of the UDP and of primitive types.
		 */
		template <typename Archive>
		void serialize(Archive& ar)
		{
			ar(m_gen, m_cr, m_eta_c, m_m, m_eta_m, m_e, m_seed, m_verbosity, m_log);
		}

	private:
		vector_double::size_type tournament_selection(vector_double::size_type idx1, vector_double::size_type idx2,
		                                              const std::vector<vector_double::size_type>& non_domination_rank,
		                                              std::vector<double>& crowding_d) const
		{
			if (non_domination_rank[idx1] < non_domination_rank[idx2]) return idx1;
			if (non_domination_rank[idx1] > non_domination_rank[idx2]) return idx2;
			if (crowding_d[idx1] > crowding_d[idx2]) return idx1;
			if (crowding_d[idx1] < crowding_d[idx2]) return idx2;
			std::uniform_real_distribution<> drng(0., 1.); // to generate a number in [0, 1)
			return ((drng(m_e) > 0.5) ? idx1 : idx2);
		}

		void crossover(EncodedVec& child1, EncodedVec& child2, vector_double::size_type parent1_idx,
		               vector_double::size_type parent2_idx, const population& pop) const
		{
			// Parents decision vectors
			EncodedVec parent1 = pop.get_x()[parent1_idx];
			EncodedVec parent2 = pop.get_x()[parent2_idx];
			child1 = EncodedVec(parent1);
			child2 = EncodedVec(parent2);

			std::uniform_real_distribution<> drng(0., 1.);
			if (drng(m_e) <= m_cr)
			{
				EncodedVec::CrossOver(parent1, parent2, child1, child2);
			}
		}

		void mutate(std::shared_ptr<TreeNode> node) const
		{
			if (node->GetN() == -1)
			{
				for (auto& subNode : node->subNodes)
				{
					mutate(subNode);
				}
			}
			else
			{
				std::uniform_real_distribution<> drng(0., 1.); // to generate a number in [0, 1)
				if (drng(m_e) <= m_m)
				{
					auto hb = static_cast<int>(node->ec->ens.size()) - 1;
					node->SetN(uniform_integral_from_range_cm(0, hb, m_e));
				}

				if (!node->IsProgram())
				{
					for (auto& subNode : node->subNodes)
					{
						mutate(subNode);
					}
				}
			}
		}

		void mutate(EncodedVec& child) const
		{
			mutate(child.root);
		}

		void check_error(std::shared_ptr<TreeNode> node) const
		{
			if (node->GetN() == -1)
			{
				for (auto& subNode : node->subNodes)
				{
					check_error(subNode);
				}
			}
			else
			{
				if (node->ec->idstr != node->val)
				{
					system("pause");
				}
			}
		}

		unsigned int m_gen;
		double m_cr;
		double m_eta_c;
		double m_m;
		double m_eta_m;
		mutable detail::random_engine_type m_e;
		unsigned int m_seed;
		unsigned int m_verbosity;
		mutable log_type m_log;
		std::shared_ptr<hv_algorithm> m_method;
	};
} // namespace pagmo


#endif
