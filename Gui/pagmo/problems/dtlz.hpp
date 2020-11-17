/* Copyright 2017-2018 PaGMO development team

This file is part of the PaGMO library.

The PaGMO library is free software; you can redistribute it and/or modify
it under the terms of either:

  * the GNU Lesser General Public License as published by the Free
	Software Foundation; either version 3 of the License, or (at your
	option) any later version.

or

  * the GNU General Public License as published by the Free Software
	Foundation; either version 3 of the License, or (at your option) any
	later version.

or both in parallel, as here.

The PaGMO library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the PaGMO library.  If not,
see https://www.gnu.org/licenses/. */

#ifndef PAGMO_PROBLEM_DTLZ_HPP
#define PAGMO_PROBLEM_DTLZ_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <pagmo/detail/constants.hpp>
#include <pagmo/exceptions.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>

#include "egraph.hpp"

// MINGW-specific warnings.
#if defined(__GNUC__) && defined(__MINGW32__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=pure"
#endif



namespace pagmo
{
	class dtlz
	{
	public:
		/** Constructor
		 *
		 * Will construct a problem of the DTLZ test-suite.
		 *
		 * @param prob_id problem id
		 * @param dim the problem dimension (size of the decision vector)
		 * @param fdim number of objectives
		 * @param alpha controls density of solutions (used only by DTLZ4)
		 *
		 * @throw std::invalid_argument if the prob_id is not in [1 .. 7], if fdim is less than 2 or if fdim or dim_param
		 * are larger than an implementation defiend value
		 *
		 */
		dtlz(unsigned int prob_id = 1u, vector_double::size_type dim = 5u, vector_double::size_type fdim = 3u,
			unsigned int alpha = 100u)
			: m_prob_id(prob_id), m_alpha(alpha), m_dim(dim), m_fdim(fdim)
		{

			if (prob_id == 0u || prob_id > 7u) {
				pagmo_throw(std::invalid_argument, "DTLZ test suite contains seven (prob_id = [1 ... 7]) problems, prob_id="
					+ std::to_string(prob_id) + " was detected");
			}
			if (fdim < 2u) {
				pagmo_throw(std::invalid_argument, "DTLZ test problem have a minimum of 2 objectives: fdim="
					+ std::to_string(fdim) + " was detected");
			}
			// We conservatively limit these dimensions to avoid checking overflows later
			if (fdim > std::numeric_limits<decltype(fdim)>::max() / 3u) {
				pagmo_throw(std::invalid_argument, "The number of objectives is too large");
			}
			if (dim > std::numeric_limits<decltype(dim)>::max() / 3u) {
				pagmo_throw(std::invalid_argument, "The problem dimension is too large");
			}
			if (dim <= fdim) {
				pagmo_throw(std::invalid_argument, "The problem dimension has to be larger than the number of objectives.");
			}
		}
		/// Fitness computation
		/**
		 * Computes the fitness for this UDP
		 *
		 * @param x the decision vector.
		 *
		 * @return the fitness of \p x.
		 */
		vector_double fitness(const OptCompiler::EncodedVec& x) const
		{
			return OptCompiler::Evaluator::EvaluatePrograms(x);
// 			vector_double retval;
// 			retval.push_back(0.3);
// 			retval.push_back(0.7);
// 			retval.push_back(0.9);
// 			return retval;
		}
		/// Number of objectives
		/**
		 *
		 * It returns the number of objectives.
		 *
		 * @return the number of objectives
		 */
		vector_double::size_type get_nobj() const
		{
			return m_fdim;
		}
		/// Box-bounds
		/**
		 *
		 * It returns the box-bounds for this UDP, [0,1] for each component
		 *
		 * @return the lower and upper bounds for each of the decision vector components
		 */
		std::pair<vector_double, vector_double> get_bounds() const
		{
			return { vector_double(m_dim, 0.), vector_double(m_dim, 1.) };
		}

		/// Problem name
		/**
		 *
		 *
		 * @return a string containing the problem name
		 */
		std::string get_name() const
		{
			return "DTLZ" + std::to_string(m_prob_id);
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
			ar(m_prob_id, m_dim, m_fdim, m_alpha);
		}

	private:
		// Problem dimensions
		unsigned int m_prob_id;
		// used only for DTLZ4
		unsigned int m_alpha;
		// dimension parameter
		vector_double::size_type m_dim;
		// number of objectives
		vector_double::size_type m_fdim;
	};
} // namespace pagmo

PAGMO_REGISTER_PROBLEM(pagmo::dtlz)

#if defined(__GNUC__) && defined(__MINGW32__)
#pragma GCC diagnostic pop
#endif

#endif