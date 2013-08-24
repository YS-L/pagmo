/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/con2uncon.h"
#include "../types.h"
#include "base.h"
#include "cstrs_core.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an immune system constraints handling algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method. Its number of
 * generations should be set to 1.
 * @param[in] original_algo_immune pagmo::algorithm to use as 'original' optimization method
 * for the immune system
 * @param[in] gen number of generations.
 * @param[in] select_method the method used for selecting the antibodies.
 * @param[in] inject_method the method used for reinjecting the antibodies.
 * @paran[in] phi the feasible fraction selection to compute the mean value
 * @paran[in] gamma number of antigens selected / number of total antigens
 * @paran[in] sigma number of antibodies / number of antigens
 * @param[in] ftol stopping criteria on the f tolerance
 * @param[in] xtol stopping criteria on the x tolerance
 * @throws value_error if gen is negative
 */
cstrs_core::cstrs_core(const base &original_algo, int gen,
					   int repair_iter,
					   double repair_tolerance,
					   double repair_step_size,
					   double ftol, double xtol):
	base(),m_gen(gen),m_repair_iter(repair_iter),m_repair_tolerance(repair_tolerance),m_repair_step_size(repair_step_size),m_ftol(ftol),m_xtol(xtol)
{
	m_original_algo = original_algo.clone();

	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor.
cstrs_core::cstrs_core(const cstrs_core &algo):
	base(algo),m_original_algo(algo.m_original_algo->clone()),m_gen(algo.m_gen),
	m_repair_iter(algo.m_repair_iter),m_repair_tolerance(algo.m_repair_tolerance),m_repair_step_size(algo.m_repair_step_size),m_ftol(algo.m_ftol), m_xtol(algo.m_xtol)
{}

/// Clone method.
base_ptr cstrs_core::clone() const
{
	return base_ptr(new cstrs_core(*this));
}

/// Evolve implementation.
/**
 * Run the CORE algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_core::evolve(population &pop) const
{	
	// store useful variables
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for CORE
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and CORE is not suitable to solve it");
	}
	if(prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,"The problem is multiobjective and CORE is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if(pop_size == 0) {
		return;
	}

	// generates the unconstrained problem
	problem::con2uncon prob_unconstrained(prob);

	// associates the population to this problem
	population pop_uncon(prob_unconstrained);

	// Main CORE loop
	for(int k=0; k<m_gen; k++) {
		if(k%10) {
			for(population::size_type i=0; i<pop_size; i++) {
				pop.repair(i, m_repair_iter, m_repair_tolerance, m_repair_step_size);
			}
		}

		// the population is repaired, it can be now used in the new unconstrained population
		pop_uncon.clear();
		// the initial popluation contains the initial random population
		for(population::size_type i=0; i<pop_size; i++) {
			pop_uncon.push_back(pop.get_individual(i).cur_x);
		}

		m_original_algo->evolve(pop_uncon);

		// push back the population in the main problem
		pop.clear();
		for(population::size_type i=0; i<pop_size; i++) {
			pop.push_back(pop_uncon.get_individual(i).cur_x);
		}

		// Check the exit conditions (every 40 generations, just as DE)
		if(k % 40 == 0) {
			decision_vector tmp(prob_dimension);

			double dx = 0;
			for(decision_vector::size_type i=0; i<prob_dimension; i++) {
				tmp[i] = pop.get_individual(pop.get_worst_idx()).best_x[i] - pop.get_individual(pop.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}

			if(dx < m_xtol ) {
				if (m_screen_output) {
					std::cout << "Exit condition -- xtol < " << m_xtol << std::endl;
				}
				break;
			}

			double mah = std::fabs(pop.get_individual(pop.get_worst_idx()).best_f[0] - pop.get_individual(pop.get_best_idx()).best_f[0]);

			if(mah < m_ftol) {
				if(m_screen_output) {
					std::cout << "Exit condition -- ftol < " << m_ftol << std::endl;
				}
				break;
			}

			// outputs current values
			if(m_screen_output) {
				std::cout << "Generation " << k << " ***" << std::endl;
				std::cout << "    Best global fitness: " << pop.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}
		}
	}

//	// store the final population in the main population
//	pop.clear();
//	for(population::size_type i=0; i<pop_size; i++) {
//		pop.push_back(pop_mixed.get_individual(i).cur_x);
//	}
}

/// Algorithm name
std::string cstrs_core::get_name() const
{
	return m_original_algo->get_name() + "[CORE]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_core::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_core::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_core::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " ";
	s << "\n\tConstraints handled with CORE algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_core);
