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

#include <cmath>
#include <algorithm>

#include "../exceptions.h"
#include "../types.h"
#include "../population.h"
#include "base.h"
#include "self_adaptive.h"

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a self-adaptive
 * as constraints handling technique.
 *
 */
self_adaptive::self_adaptive(const base &problem, const population &pop):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_pop(pop),
	m_solution_infeasibility(pop.size(),0.)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
	set_population(pop);
}

/// Copy Constructor. Performs a deep copy
self_adaptive::self_adaptive(const self_adaptive &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_pop(prob.m_pop),
	m_solution_infeasibility(m_pop.size(),0.)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
	set_population(m_pop);
}

/// Clone method.
base_ptr self_adaptive::clone() const
{
	return base_ptr(new self_adaptive(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
void self_adaptive::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	constraint_vector c(m_original_problem->get_c_dimension(),0);
	m_original_problem->compute_constraints(c,x);

	if(m_original_problem->feasibility_c(c)) {
		m_original_problem->objfun(f, x);
	} else {
		// finds the position of the asked decision_vector

		pagmo::population::size_type idx_x = -1;

		for(pagmo::population::size_type i=0; i<m_pop.size(); i++) {

			const population::individual_type &current_individual = m_pop.get_individual(i);
			if(current_individual.cur_x == x) {
				idx_x = i;
				break;
			}
		}
		if(idx_x == -1) {
			// if the original problem could not be found, it means that we are evaluating a new population
			// in that case, we use the original problem
			m_original_problem->objfun(f, x);
			return;
		}

//		if(m_apply_penalty_1) {
//			double inf_tilde =
//					(m_solution_infeasibility.at(idx_x) - m_solution_infeasibility.at(m_hat_down_idx)) /
//					(m_solution_infeasibility.at(m_hat_up_idx) - m_solution_infeasibility.at(m_hat_down_idx));

		//	m_f_dot[idx_x][0] = m_pop.get_individual(idx_x).cur_f[0] + inf_tilde *
		//			(m_original_problem->objfun(m_pop.get_individual(m_hat_down_idx).cur_x)[0] -
		//			m_original_problem->objfun(m_pop.get_individual(m_hat_up_idx).cur_x)[0]) ;

//		}
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string self_adaptive::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with self-adaptive method ";
	oss << std::endl;
	return oss.str();
}

std::string self_adaptive::get_name() const
{
	return m_original_problem->get_name() + " [self_adaptive]";
}

/// Sets the population and updates penalty coefficients.
/**
 * Will update penalty coefficients based on the given popuation.
 */
void self_adaptive::set_population(const population &pop)
{
	m_pop = pop;

	// Let's store some useful variables.
	const problem::base::size_type dimension = m_original_problem->get_dimension(), prob_i_dimension = m_original_problem->get_i_dimension();
	const decision_vector &lb = m_original_problem->get_lb(), &ub = m_original_problem->get_ub();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type dimension_continuous = dimension - prob_i_dimension;

	// get the constraints dimension
	constraint_vector c(m_original_problem->get_c_dimension(),0);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	std::vector<pagmo::population::size_type> feasible_idx;
	std::vector<pagmo::population::size_type> infeasible_idx;

	// Main Self-Adaptive loop
	//for (int k=0; k < m_gen; k++) {
	// check if the population contains at least one individual that is non feasible

	feasible_idx.clear();
	infeasible_idx.clear();

	// evaluates the scaling factor
	for(pagmo::population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);
		if(m_original_problem->feasibility_x(current_individual.cur_x)) {
			feasible_idx.push_back(i);
		} else {
			infeasible_idx.push_back(i);
		}
	}

	// evaluate solutions infeasibility
	//std::vector<double> solution_infeasibility(pop_size,0.);
	compute_solution_infeasibility(m_solution_infeasibility, pop);

	// identification of bounding solutions
	m_hat_down_idx = 0;
	m_hat_up_idx = 0;

	// first case
	if(feasible_idx.size() > 0) {
		// x_hat_down = feasible individual with lowest objective value in p
		for(pagmo::population::size_type i=0; i < feasible_idx.size(); i++) {
			const population::individual_type &current_individual = pop.get_individual(feasible_idx.at(i));

			if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_down_idx).cur_f)) {
				m_hat_down_idx = feasible_idx.at(i);
			}
		}

		// hat down is now availlable
		decision_vector x_hat_down = pop.get_individual(m_hat_down_idx).cur_x;
		fitness_vector f_hat_down =  pop.get_individual(m_hat_down_idx).cur_f;

		// x_hat_up
		bool pop_contains_infeasible_f_better_x_hat_down = false;
		for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
			const population::individual_type &current_individual = pop.get_individual(infeasible_idx.at(i));

			if(m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down)) {
				pop_contains_infeasible_f_better_x_hat_down = true;
				// stores this individual as the first possible individual
				m_hat_up_idx = infeasible_idx.at(i);
				break;
			}
		}

		if(pop_contains_infeasible_f_better_x_hat_down) {
			// gets the individual with maximum infeasibility and objfun lower than f_hat_down
			for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
				const population::individual_type &current_individual =
						pop.get_individual(infeasible_idx.at(i));

				if( m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down) &&
					(m_solution_infeasibility.at(infeasible_idx.at(i)) >= m_solution_infeasibility.at(m_hat_up_idx)) ) {

					if(m_solution_infeasibility.at(infeasible_idx.at(i)) == m_solution_infeasibility.at(m_hat_up_idx)) {
						if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_up_idx).cur_f)) {
							m_hat_up_idx = infeasible_idx.at(i);
						}
					} else {
						m_hat_up_idx = infeasible_idx.at(i);
					}
				}
			}

			// apply penalty 1
			m_apply_penalty_1 = true;

		} else { // the worst is the one that has the maximum infeasibility
			for(pagmo::population::size_type i=0; i < infeasible_idx.size(); i++) {
				const population::individual_type &current_individual =pop.get_individual(infeasible_idx.at(i));

				if( (m_solution_infeasibility.at(infeasible_idx.at(i)) >= m_solution_infeasibility.at(m_hat_up_idx)) ) {
					if(m_solution_infeasibility.at(infeasible_idx.at(i)) == m_solution_infeasibility.at(m_hat_up_idx)) {
						if(!m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_up_idx).cur_f)) {
							m_hat_up_idx = infeasible_idx.at(i);
						}
					} else {
						m_hat_up_idx = infeasible_idx.at(i);
					}
				}

				// do not apply penalty 1
				m_apply_penalty_1 = false;
			}
		}
		// hat up is now availlable

	} else { // case where there is no feasible solution in the population
		// best is the individual with the lowest infeasibility
		for(pagmo::population::size_type i=0; i < pop_size; i++) {
			if(m_solution_infeasibility.at(i) < m_solution_infeasibility.at(m_hat_down_idx)) {
				m_hat_down_idx = i;
			}
		}
		// worst individual
		for(pagmo::population::size_type i=0; i<pop_size; i++) {
			if(m_solution_infeasibility.at(i) > m_solution_infeasibility.at(m_hat_up_idx)) {
				m_hat_up_idx = i;
			}
		}

		// apply penalty 1 to the population
		m_apply_penalty_1 = true;
	}

	// stores the hat round idx, i.e. the solution with highest objective
	// function value in the population
	m_hat_round_idx = 0;
	for(pagmo::population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_round_idx).cur_f)) {
			m_hat_round_idx = i;
		}
	}

	fitness_vector f_hat_round = pop.get_individual(m_hat_round_idx).cur_f;
	fitness_vector f_hat_down = pop.get_individual(m_hat_down_idx).cur_f;
	fitness_vector f_hat_up = pop.get_individual(m_hat_up_idx).cur_f;

	// evaluates scaling factor
	if(f_hat_up[0] <= f_hat_down[0]) {
		m_scaling_factor = (f_hat_round[0] - f_hat_down[0]) / f_hat_down[0];
	} else if(f_hat_up[0] == f_hat_round[0]) {
		m_scaling_factor = 0.;
	} else if(f_hat_up[0] > f_hat_down[0]) {
		m_scaling_factor = (f_hat_round[0] - f_hat_up[0]) / f_hat_up[0];
	}

	// evaluates the f_dot for all the current population
	m_f_dot.resize(pop_size);

	for(pagmo::population::size_type i=0; i<infeasible_idx.size(); i++) {

		if(m_apply_penalty_1) {
			double inf_tilde =
					(m_solution_infeasibility.at(infeasible_idx.at(i)) - m_solution_infeasibility.at(m_hat_down_idx)) /
					(m_solution_infeasibility.at(m_hat_up_idx) - m_solution_infeasibility.at(m_hat_down_idx));

//			m_f_dot[infeasible_idx.at(i)][0] = m_pop.get_individual(infeasible_idx.at(i)).cur_f[0] + inf_tilde *
//					(m_original_problem->objfun(m_pop.get_individual(m_hat_down_idx).cur_x)[0] -
//					m_original_problem->objfun(m_pop.get_individual(m_hat_up_idx).cur_x)[0]) ;

		}
	}
}

void self_adaptive::compute_solution_infeasibility(std::vector<double> &solution_infeasibility, const population &pop)
{
	// Let's store some useful variables.
	const problem::base::size_type dimension = m_original_problem->get_dimension(), prob_i_dimension = m_original_problem->get_i_dimension();
	const decision_vector &lb = m_original_problem->get_lb(), &ub = m_original_problem->get_ub();
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	constraint_vector c(m_original_problem->get_c_dimension(),0);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	constraint_vector c_scaling(m_original_problem->get_c_dimension(),0);

	// evaluates the scaling factor
	for(pagmo::population::size_type i=0; i<pop_size; i++) {

		// updates the current constraint vector
		const population::individual_type &current_individual = pop.get_individual(i);
		m_original_problem->compute_constraints(c,current_individual.cur_x);

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
	solution_infeasibility.resize(pop_size,0.);

	for(pagmo::population::size_type i=0; i<pop_size; i++) {

		// updates the current constraint vector
		const population::individual_type &current_individual = pop.get_individual(i);
		m_original_problem->compute_constraints(c,current_individual.cur_x);

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
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::self_adaptive);

