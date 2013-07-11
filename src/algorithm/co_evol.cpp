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
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/death_penalty.h"
#include "../types.h"
#include "base.h"
#include "co_evol.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an self adaptive algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @throws value_error if stop is negative
 */
co_evol::co_evol(const base &original_algo, int gen, double pop_ratio):base(),m_gen(gen),m_pop_ratio(pop_ratio)
{
	m_original_algo = original_algo.clone();
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if((pop_ratio >= 1.) || (pop_ratio <= 0.)) {
		pagmo_throw(value_error,"the population ratio pop_2/pop_1 must be between 0 and 1");
	}
}

/// Copy constructor.
co_evol::co_evol(const co_evol &algo):base(algo),m_original_algo(algo.m_original_algo->clone()),m_gen(algo.m_gen),m_pop_ratio(algo.m_pop_ratio)
{}

/// Clone method.
base_ptr co_evol::clone() const
{
	return base_ptr(new co_evol(*this));
}

/// Evolve implementation.
/**
 * Run the Self-Adaptive algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void co_evol::evolve(population &pop) const
{	
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for Self-Adaptive
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and Self-Adaptive is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	// split the population into two populations
	// the population P1 associated to the modified problem with penalized fitness
	// the reason P2 is not a population type is that only one iteration of prob2
	// is necessary and it is not easy to locally update the fitness

	// creates the problem 1
	// TODO: CHANGE WITH CO EVOL
	//problem::co_evol prob_1(prob);
	problem::death_penalty prob_1(prob);

	//population pop_2;

	// sub population size
	population::size_type sub_pop_1_size = population::size_type(pop_size*m_pop_ratio);
	population::size_type sub_pop_2_size = population::size_type(pop_size/m_pop_ratio);

	// initialize the population 1
	population sub_pop_1(prob_1,0);
//	population sub_pop_1(sub_pop_1_size);
	sub_pop_1.clear();
	for(population::size_type i=0; i<sub_pop_1_size; i++) {
		sub_pop_1.push_back(pop.get_individual(i).cur_x);
	}

	// initialize the population 2, composed of individuals only
//	std::vector<population::individual_type> pop_2(sub_pop_2_size);
//	for(population::size_type i=0; i<sub_pop_2_size; i++) {
//		pop_2.push_back(population::individual_type());
//		pop_2.back().cur_x.resize(2);
//		pop_2.back().cur_v.resize(2);
//		pop_2.back().cur_c.resize(0);
//		pop_2.back().cur_f.resize(1);
//		pop_2.back().best_x.resize(2);
//		pop_2.back().best_c.resize(0);
//		pop_2.back().best_f.resize(1);
//	}

	std::vector<decision_vector> sub_pop_2_x(sub_pop_2_size);
	std::vector<fitness_vector> sub_pop_2_f(sub_pop_2_size);
	for(population::size_type i=0; i<sub_pop_2_size; i++) {
		sub_pop_2_x[i] = decision_vector(2);
		sub_pop_2_f[i] = fitness_vector(1);
	}

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {
		std::cout << "current generation: " << k << std::endl;

		// for each individuals of pop 2

		for(int j=0; j<sub_pop_2_size; j++) {
			// set decision vector encoding w1 and w2 in prob 1
			//prob_1.set_co_decision(sub_pop_2_x[i]);

			// the problem 1 depends on the population 2
			// sub_pop_1.problem().set_population(pop_2);

			// creates an instance of the population 1
			population pop_1_instance(sub_pop_1);


			// evolve the population 1 instance
			m_original_algo->evolve(pop_1_instance);

			// store indexes of feasible individuals
			std::vector<population::size_type> feasible_idx(0);

			for(population::size_type i=0; i<sub_pop_1_size; i++) {
				const population::individual_type &current_individual = pop_1_instance.get_individual(i);

				if(prob.feasibility_x(current_individual.cur_x)) {
					feasible_idx.push_back(i);
				}
			}

			// computes the average fitness
			int number_of_feasible=feasible_idx.size();
			double average_fitness=0.;
			for(int i=0; i<feasible_idx.size();i++) {
				const population::individual_type &current_individual = pop_1_instance.get_individual(feasible_idx.at(i));

				average_fitness += current_individual.cur_f.at(0);
			}
			average_fitness = (average_fitness / number_of_feasible) + number_of_feasible;

			// assuming minimization (convert to objective function)
			average_fitness = - average_fitness;

			// store the average fitness to the individual j
			sub_pop_2_f[j][0] = average_fitness;
		}

		// need to evolve the population 2 now

//		// update the population pop
//		pop.clear();
//		for(pagmo::population::size_type i=0; i<pop_new.size(); i++) {
//			pop.push_back(pop_new.get_individual(i).cur_x);
//		}


//		std::cout << pop.champion() << std::endl;
	}
}

/// Algorithm name
std::string co_evol::get_name() const
{
	return m_original_algo->get_name() + "[Co-Evolution]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr co_evol::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void co_evol::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string co_evol::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm: " << m_original_algo->get_name() << ' ';
	s << "\n\tConstraints handled with co-evolution algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::co_evol);
