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
#include "../problem/cstrs_co_evolution.h"
#include "../types.h"
#include "base.h"
#include "cstrs_co_evolution.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs a co-evolution adaptive algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @param[in] original_algo_2 pagmo::algorithm to use as 'original' optimization method for population 2
 * @param[in] pop_2_size population size for the penalty encoding population.
 * @param[in] gen number of generations.
 * @param[in] problem::cstrs_co_evolution::method_type the method used for the population 2.
 * Three posssibililties are available: SIMPLE, SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS.
 * The simple one is the original version of the Coello/He implementation. The SPLIT_NEQ_EQ,
 * splits the equalities and inequalities constraints in two different sets for the
 * penalty weigths, containing respectively inequalities and equalities weigths. The
 * SPLIT_CONSTRAINTS splits the constraints in M set of weigths with M the number of
 * constraints.
 * @param[in] pen_lower_bound the lower boundary used for penalty.
 * @param[in] pen_upper_bound the upper boundary used for penalty.
 * @throws value_error if stop is negative
 */
cstrs_co_evolution::cstrs_co_evolution(const base &original_algo, const base &original_algo_2, int pop_2_size,
									   int gen,method_type method, double pen_lower_bound,
									   double pen_upper_bound):
	base(),m_gen(gen),m_pop_2_size(pop_2_size),m_method(method),
	m_pen_lower_bound(pen_lower_bound),m_pen_upper_bound(pen_upper_bound)
{
	m_original_algo = original_algo.clone();
	m_original_algo_2 = original_algo_2.clone();
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if(pop_2_size <= 0) {
		pagmo_throw(value_error,"the population of the population 2 be greater than 0");
	}
}

/// Copy constructor.
cstrs_co_evolution::cstrs_co_evolution(const cstrs_co_evolution &algo):
	base(algo),m_original_algo(algo.m_original_algo->clone()),
	m_original_algo_2(algo.m_original_algo_2->clone()),m_gen(algo.m_gen),
	m_pop_2_size(algo.m_pop_2_size),m_method(algo.m_method),
	m_pen_lower_bound(algo.m_pen_lower_bound),m_pen_upper_bound(algo.m_pen_upper_bound)
{}

/// Clone method.
base_ptr cstrs_co_evolution::clone() const
{
	return base_ptr(new cstrs_co_evolution(*this));
}

/// Evolve implementation.
/**
 * Run the co-evolution algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_co_evolution::evolve(population &pop) const
{	
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for co-evolution
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and co-evolution is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if(pop_size == 0) {
		return;
	}

	// split the population into two populations
	// the population P1 associated to the modified problem with penalized fitness
	// the reason P2 is not a population type is that only one iteration of prob2
	// is necessary and it is not easy to locally update the fitness

	// creates the problem 1 and 2
	problem::cstrs_co_evolution prob_1(prob, m_method);
	problem::cstrs_co_evolution_2 prob_2(prob, prob_1.get_expected_penalty_coeff_size());

	// sub population size
	population::size_type sub_pop_1_size = pop_size;
	population::size_type sub_pop_2_size = m_pop_2_size;

	// set up the structure for the population 1 and 2
	std::vector< std::vector<decision_vector> > sub_pop_1_x_vector(sub_pop_2_size);
	std::vector< std::vector<fitness_vector> > sub_pop_1_f_vector(sub_pop_2_size);

	for(population::size_type i=0; i<sub_pop_2_size; i++) {
		sub_pop_1_x_vector[i] = std::vector<decision_vector>(sub_pop_1_size);
		sub_pop_1_f_vector[i] = std::vector<fitness_vector>(sub_pop_1_size);
	}

	// initialize sub_pop_1 vectors
	for(population::size_type j=0; j<sub_pop_2_size; j++) {
		std::vector<decision_vector> &sub_pop_1_x = sub_pop_1_x_vector.at(j);
		std::vector<fitness_vector> &sub_pop_1_f = sub_pop_1_f_vector.at(j);

		for(population::size_type i=0; i<sub_pop_1_size; i++) {
			sub_pop_1_x[i] = pop.get_individual(i).cur_x;
			sub_pop_1_f[i] = pop.get_individual(i).cur_f;
		}
	}

	std::vector<decision_vector> sub_pop_2_x(sub_pop_2_size);
	std::vector<fitness_vector> sub_pop_2_f(sub_pop_2_size);

	// population 2 problem related variables
	// decision vector is depends on the method we use for 1
	int pop_2_x_dimension = prob_1.get_expected_penalty_coeff_size();

	decision_vector lb(pop_2_x_dimension);
	decision_vector ub(pop_2_x_dimension);

	std::fill(lb.begin(),lb.end(),m_pen_lower_bound);
	std::fill(ub.begin(),ub.end(),m_pen_upper_bound);

	prob_2.set_bounds(lb,ub);

	// population 2 initialization
	// might be done with the problem_2 definition, but didn't want to create a
	// sub_pop_2 here...
	for(population::size_type j=0; j<sub_pop_2_size; j++) {
		sub_pop_2_x[j] = decision_vector(pop_2_x_dimension,0.);

		// choose random coefficients between lower bound and upper bound
		for(population::size_type i=0; i<pop_2_x_dimension;i++) {
			sub_pop_2_x[j][i] = boost::uniform_real<double>(lb[i],ub[i])(m_drng);
		}
		sub_pop_2_f[j] = fitness_vector(1);
	}

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {
		// for each individuals of pop 2, evolve the current population,
		// and store the position of the feasible idx
		for(population::size_type j=0; j<sub_pop_2_size; j++) {

			std::vector<decision_vector> &sub_pop_1_x = sub_pop_1_x_vector.at(j);
			std::vector<fitness_vector> &sub_pop_1_f = sub_pop_1_f_vector.at(j);

			// modify the problem by settin decision vector encoding penalty
			// coefficients w1 and w2 in prob 1

			prob_1.set_penalty_coeff(sub_pop_2_x.at(j));

			// creating the subpopulation pop 1 instance based on the
			// updated prob 1
			population sub_pop_1_instance(prob_1,0);

			// initialize the instance with pop_1 chromosomes
			for(population::size_type i=0; i<sub_pop_1_size; i++) {
				sub_pop_1_instance.push_back(sub_pop_1_x[i]);
			}

			// evolve the population 1 instance
			m_original_algo->evolve(sub_pop_1_instance);

			// store the new chromosomes
			for(population::size_type i=0; i<sub_pop_1_size; i++) {
				sub_pop_1_x[i] = sub_pop_1_instance.get_individual(i).cur_x;
				sub_pop_1_f[i] = sub_pop_1_instance.get_individual(i).cur_f;
			}
		}

		// set up penalization variables needs for the population 2
		// the constraints has not been evaluated yet.
		prob_2.update_penalty_coeff(sub_pop_2_x, sub_pop_1_x_vector, sub_pop_1_f_vector);
		population sub_pop_2(prob_2,0);

		// initialize the instance with pop_1 chromosomes
		for(population::size_type i=0; i<sub_pop_2_size; i++) {
			sub_pop_2.push_back(sub_pop_2_x[i]);
		}

		m_original_algo_2->evolve(sub_pop_2);

		// store the new chromosomes
		for(population::size_type i=0; i<sub_pop_2_size; i++) {
			sub_pop_2_x[i] = sub_pop_2.get_individual(i).cur_x;
			sub_pop_2_f[i] = sub_pop_2.get_individual(i).cur_f;
		}
	}

	// store the best fitness population in the final pop
	population::size_type best_idx = 0;
	for(population::size_type j=1; j<sub_pop_2_size; j++) {
		if(sub_pop_2_f[j][0] < sub_pop_2_f[best_idx][0]) {
			best_idx = j;
		}
	}

	// store the final population in the main population
	pop.clear();
	for(population::size_type i=0; i<sub_pop_1_size; i++) {
		pop.push_back(sub_pop_1_x_vector[best_idx][i]);
	}

	std::cout << pop.champion() << std::endl;
}

/// Algorithm name
std::string cstrs_co_evolution::get_name() const
{
	return m_original_algo->get_name() + "[Co-Evolution]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_co_evolution::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_co_evolution::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_co_evolution::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " - " << m_original_algo_2->get_name() << " ";
	s << "\n\tConstraints handled with co-evolution algorithm";
	return s.str();
}
}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_co_evolution);
