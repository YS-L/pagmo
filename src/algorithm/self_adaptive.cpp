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
#include "../types.h"
#include "base.h"
#include "self_adaptive.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an self adaptive algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @throws value_error if stop is negative or perturb is not in [0,1]
 */
self_adaptive::self_adaptive(const base &original_algo, int gen):base(),m_gen(gen)
{
	m_original_algo = original_algo.clone();
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor.
self_adaptive::self_adaptive(const self_adaptive &algo):base(algo),m_original_algo(algo.m_original_algo->clone()),m_gen(algo.m_gen)
{}

/// Clone method.
base_ptr self_adaptive::clone() const
{
	return base_ptr(new self_adaptive(*this));
}

/// Evolve implementation.
/**
 * Run the Self-Adaptive algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void self_adaptive::evolve(population &pop) const
{
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type dimension = prob.get_dimension(), prob_i_dimension = prob.get_i_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type dimension_continuous = dimension - prob_i_dimension;

	// get the constraints dimension
	constraint_vector c(prob.get_c_dimension(),0);
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints = prob.get_c_dimension() -
			prob.get_ic_dimension();

	//We perform some checks to determine wether the problem/population are suitable for Self-Adaptive
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and Self-Adaptive is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	const std::vector<double> &c_tol = prob.get_c_tol();

	std::vector<pagmo::population::size_type> feasible_idx;
	std::vector<pagmo::population::size_type> infeasible_idx;

	// Main Self-Adaptive loop
	for (int k=0; k < m_gen; k++) {
		// check if the population contains at least one individual that is non feasible

		feasible_idx.clear();
		infeasible_idx.clear();

		// evaluates the scaling factor
		for(pagmo::population::size_type i=0; i < pop_size; i++) {
			const population::individual_type &current_individual = pop.get_individual(i);
			if(prob.feasibility_x(current_individual.cur_x)) {
				feasible_idx.push_back(i);
			} else {
				infeasible_idx.push_back(i);
			}
		}

		// test a faire sur non feasible !
		if(feasible_idx.size() == 0) {
			// use real fitness as fitness evolution without touching the initial population
			// TODO! Need to convert the problem to unconstrained problem and evove the population
			// against this newly defined problem

			//m_original_algo->evolve(pop);

		} else {
			// scaling factor is the max constraint violation value in the current population
			constraint_vector c_scaling(prob.get_c_dimension(),0);

			// evaluates the scaling factor
			for(pagmo::population::size_type i=0; i<pop_size; i++) {

				// updates the current constraint vector
				const population::individual_type &current_individual = pop.get_individual(i);
				prob.compute_constraints(c,current_individual.cur_x);

				// sets the right definition of the constraints (can be in base problem? currently used
				// by con2mo as well)
				for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
					c[j] = std::abs(c.at(j)) - c_tol.at(j);
				}
				for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
					c[j] = std::max(0.,c.at(j));
				}

				// computes scaling
				for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
					c_scaling[j] = std::max(c_scaling[j],c[j]);
				}
			}

			// evaluate solutions infeasibility
			std::vector<double> solution_infeasibility(pop_size,0.);

			for(pagmo::population::size_type i=0; i<pop_size; i++) {

				// updates the current constraint vector
				const population::individual_type &current_individual = pop.get_individual(i);
				prob.compute_constraints(c,current_individual.cur_x);

				// sets the right definition of the constraints (can be in base problem? currently used
				// by con2mo as well)
				for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
					c[j] = std::abs(c.at(j)) - c_tol.at(j);
				}
				for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
					c[j] = std::max(0.,c.at(j));
				}

				// computes soution infeasibility
				for(problem::base::c_size_type j=0; j < prob_c_dimension; j++) {
					solution_infeasibility[i] += c[j]/c_scaling[j];
				}
				solution_infeasibility[i] /= prob_c_dimension;
			}

			// identification of bounding solutions

			population::size_type hat_down_idx = 0;
			population::size_type hat_up_idx = 0;

			// first case
			if(feasible_idx.size() > 0) {
				// x_hat_down = feasible individual with lowest objective value in p
				for(pagmo::population::size_type i=0; i < feasible_idx.size(); i++) {
					const population::individual_type &current_individual = pop.get_individual(feasible_idx.at(i));

					if(prob.compare_fitness(current_individual.cur_f, pop.get_individual(hat_down_idx).cur_f)) {
						hat_down_idx = i;
					}
				}

				// hat down is now availlable
				decision_vector x_hat_down = pop.get_individual(hat_down_idx).cur_x;
				fitness_vector f_hat_down =  pop.get_individual(hat_down_idx).cur_f;

				// x_hat_up
				bool pop_contains_infeasible_f_better_x_hat_down = false;
				for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
					const population::individual_type &current_individual = pop.get_individual(infeasible_idx.at(i));

					if(prob.compare_fitness(current_individual.cur_f, f_hat_down)) {
						pop_contains_infeasible_f_better_x_hat_down = true;
						// stores this individual as the first possible individual
						hat_up_idx = i;
						break;
					}
				}

				if(pop_contains_infeasible_f_better_x_hat_down) {
					// gets the individual with maximum infeasibility and objfun lower than f_hat_down
					for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
						const population::individual_type &current_individual =
								pop.get_individual(infeasible_idx.at(i));

						if( prob.compare_fitness(current_individual.cur_f, f_hat_down) &&
							(solution_infeasibility.at(infeasible_idx.at(i)) >= solution_infeasibility.at(hat_up_idx)) ) {

							if(solution_infeasibility.at(infeasible_idx.at(i)) == solution_infeasibility.at(hat_up_idx)) {
								if(prob.compare_fitness(current_individual.cur_f, pop.get_individual(hat_up_idx).cur_f)) {
									hat_up_idx = i;
								}
							} else {
								hat_up_idx = i;
							}
						}
					}

					// apply penalty 1


				} else { // the worst is the one that has the maximum infeasibility
					for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
						const population::individual_type &current_individual =pop.get_individual(infeasible_idx.at(i));

						if( (solution_infeasibility.at(infeasible_idx.at(i)) >= solution_infeasibility.at(hat_up_idx)) ) {
							if(solution_infeasibility.at(infeasible_idx.at(i)) == solution_infeasibility.at(hat_up_idx)) {
								if(!prob.compare_fitness(current_individual.cur_f, pop.get_individual(hat_up_idx).cur_f)) {
									hat_up_idx = i;
								}
							} else {
								hat_up_idx = i;
							}
						}
					}


				// do not apply penalty 1

				}
				// hat up is now availlable

			} else { // case where there is no feasible solution in the population
				// best is the individual with the lowest infeasibility
				for(pagmo::population::size_type i=0; i<pop_size; i++) {
					if(solution_infeasibility.at(i) < solution_infeasibility.at(hat_down_idx)) {
						hat_down_idx = i;
					}
				}
				// worst individual
				for(pagmo::population::size_type i=0; i<pop_size; i++) {
					if(solution_infeasibility.at(i) > solution_infeasibility.at(hat_up_idx)) {
						hat_up_idx = i;
					}
				}

				// apply penalty 1 to the population
				// I think I killed the possibility of using maximization here...
				for(pagmo::population::size_type i=0; i<pop_size; i++) {
					// updates the fitness of the population
					// f_dot[i] =
				}
			}
		}
	}
}

/// Algorithm name
std::string self_adaptive::get_name() const
{
	return m_original_algo->get_name() + "[Self-Adaptive]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr self_adaptive::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void self_adaptive::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string self_adaptive::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm: " << m_original_algo->get_name() << ' ';
	s << "\n\tConstraints handled with Self-Adaptive algorithm";
	return s.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::self_adaptive);
