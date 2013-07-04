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
	m_pop(pop)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with self adaptive meta problem.");
	}

	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
	update_fitness(m_pop);
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
	m_pop(prob.m_pop)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
	update_fitness(m_pop);
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
	// finds the position of the asked decision_vector
	population::size_type idx_x = -1;

	for(population::size_type i=0; i<m_pop.size(); i++) {
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

		//std::cout << "could not find fitness: " << f << std::endl;
	} else {
		f = m_fitness[idx_x];
		//std::cout << "found fitness: " << f << std::endl;
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

/// Sets the population.
/**
 * Will update fitness as well based on the given popuation.
 */
void self_adaptive::set_population(const population &pop)
{
	m_pop.clear();
	for(population::size_type i=0; i<pop.size(); i++) {
		// Evaluate according to the new fitness;
		m_pop.push_back(pop.get_individual(i).cur_x);
	}

	update_fitness(m_pop);
}

void self_adaptive::update_fitness(const population &pop)
{
	// Let's store some useful variables.
	const population::size_type pop_size = pop.size();

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	std::vector<population::size_type> feasible_idx(0);
	std::vector<population::size_type> infeasible_idx(0);

	// store indexes of feasible and non feasible individuals
	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		if(m_original_problem->feasibility_x(current_individual.cur_x)) {
			feasible_idx.push_back(i);
		} else {
			infeasible_idx.push_back(i);
		}
	}

	// evaluates the fitness for the whole current population
	m_fitness.resize(pop_size);
	for(population::size_type i=0; i<pop_size; i++) {
		m_fitness[i] = pop.get_individual(i).cur_f;

	}

	// if the population is only feasible, then nothing is done
	if(infeasible_idx.size() == 0) {
		return;
	}

	bool apply_penalty_1 = false;
	double scaling_factor = 0.;

	std::vector<double> solution_infeasibility(pop_size);
	std::fill(solution_infeasibility.begin(),solution_infeasibility.end(),0.);

	population::size_type m_hat_down_idx = -1;
	population::size_type m_hat_up_idx = -1;
	population::size_type m_hat_round_idx = -1;

	// evaluate solutions infeasibility
	compute_solution_infeasibility(solution_infeasibility, pop);

	// first case, the population contains at least one feasible solution
	if(feasible_idx.size() > 0) {
		// initialize hat_down_idx
		m_hat_down_idx = feasible_idx.at(0);

		// x_hat_down = feasible individual with lowest objective value in p
		for(population::size_type i=0; i<feasible_idx.size(); i++) {
			const population::size_type current_idx = feasible_idx.at(i);
			const population::individual_type &current_individual = pop.get_individual(current_idx);

			if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_down_idx).cur_f)) {
				m_hat_down_idx = current_idx;
			}
		}

		// hat down is now availlable
		fitness_vector f_hat_down = pop.get_individual(m_hat_down_idx).cur_f;

		// x_hat_up value depends if the population contains infeasible individual with objective
		// function better than f_hat_down
		bool pop_contains_infeasible_f_better_x_hat_down = false;
		for(population::size_type i=0; i<infeasible_idx.size(); i++) {
			const population::size_type current_idx = infeasible_idx.at(i);
			const population::individual_type &current_individual = pop.get_individual(current_idx);

			if(m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down)) {
				pop_contains_infeasible_f_better_x_hat_down = true;

				// initialize m_hat_up_idx
				m_hat_up_idx = current_idx;

				break;
			}
		}

		if(pop_contains_infeasible_f_better_x_hat_down) {
			// m_hat_up_idx is already initizalized

			// gets the individual with maximum infeasibility and objfun lower than f_hat_down
			for(population::size_type i=0; i<infeasible_idx.size(); i++) {
				const population::size_type current_idx = infeasible_idx.at(i);
				const population::individual_type &current_individual = pop.get_individual(current_idx);

				if(m_original_problem->compare_fitness(current_individual.cur_f, f_hat_down) &&
					(solution_infeasibility.at(current_idx) >= solution_infeasibility.at(m_hat_up_idx)) ) {

					if(solution_infeasibility.at(current_idx) == solution_infeasibility.at(m_hat_up_idx)) {
						if(m_original_problem->compare_fitness(current_individual.cur_f, pop.get_individual(m_hat_up_idx).cur_f)) {
							m_hat_up_idx = current_idx;
						}
					} else {
						m_hat_up_idx = current_idx;
					}
				}
			}

			// apply penalty 1
			apply_penalty_1 = true;

		} else {
			// all the infeasible soutions have an objective function value greater than f_hat_down
			// the worst is the one that has the maximum infeasibility
			// initialize m_hat_up_idx
			m_hat_up_idx = infeasible_idx.at(0);

			for(population::size_type i=0; i<infeasible_idx.size(); i++) {
				const population::size_type current_idx = infeasible_idx.at(i);
				const population::individual_type &current_individual = pop.get_individual(current_idx);

				if(solution_infeasibility.at(current_idx) >= solution_infeasibility.at(m_hat_up_idx)) {
					if(solution_infeasibility.at(current_idx) == solution_infeasibility.at(m_hat_up_idx)) {
						if(m_original_problem->compare_fitness(pop.get_individual(m_hat_up_idx).cur_f, current_individual.cur_f)) {
							m_hat_up_idx = current_idx;
						}
					} else {
						m_hat_up_idx = current_idx;
					}
				}
			}

			// do not apply penalty 1
			apply_penalty_1 = false;
		}

	} else { // case where there is no feasible solution in the population
		// best is the individual with the lowest infeasibility
		m_hat_down_idx = 0;
		m_hat_up_idx = 0;

		// case of equality? what do we do?
		for(population::size_type i=0; i<pop_size; i++) {
			if(solution_infeasibility.at(i) <= solution_infeasibility.at(m_hat_down_idx)) {
				m_hat_down_idx = i;
			}
		}
		// worst individual
		for(population::size_type i=0; i<pop_size; i++) {
			if(solution_infeasibility.at(i) >= solution_infeasibility.at(m_hat_up_idx)) {
				m_hat_up_idx = i;
			}
		}

		// apply penalty 1 to the population
		apply_penalty_1 = true;
	}

	// stores the hat round idx, i.e. the solution with highest objective
	// function value in the population
	m_hat_round_idx = 0;
	for(population::size_type i=0; i<pop_size; i++) {
		const population::individual_type &current_individual = pop.get_individual(i);

		if(m_original_problem->compare_fitness(pop.get_individual(m_hat_round_idx).cur_f, current_individual.cur_f)) {
			m_hat_round_idx = i;
		}
	}

	// get the objective function values of the three individuals
	fitness_vector f_hat_round = pop.get_individual(m_hat_round_idx).cur_f;//m_original_problem->objfun(pop.get_individual(m_hat_round_idx).cur_x);//
	fitness_vector f_hat_down =  pop.get_individual(m_hat_down_idx).cur_f;//m_original_problem->objfun(pop.get_individual(m_hat_down_idx).cur_x);//
	fitness_vector f_hat_up = pop.get_individual(m_hat_up_idx).cur_f;//m_original_problem->objfun(pop.get_individual(m_hat_up_idx).cur_x);//

	// apply penalty
	if(apply_penalty_1) {
		double inf_tilde = 0.;
		for(population::size_type i=0; i<infeasible_idx.size(); i++) {
			const population::size_type current_idx = infeasible_idx.at(i);
			inf_tilde = (solution_infeasibility.at(current_idx) - solution_infeasibility.at(m_hat_down_idx)) /
					(solution_infeasibility.at(m_hat_up_idx) - solution_infeasibility.at(m_hat_down_idx));

			m_fitness[current_idx][0] += inf_tilde * (f_hat_down[0] - f_hat_up[0]);
		}
	}

	scaling_factor = 0.;
	// evaluates scaling factor
	if(m_original_problem->compare_fitness(f_hat_up, f_hat_down)) {
	//if(f_hat_up[0] <= f_hat_down[0]) {
		scaling_factor = (f_hat_round[0] - f_hat_down[0]) / f_hat_down[0];
	}
	if(f_hat_up[0] == f_hat_round[0]) {
		scaling_factor = 0.;
	}
	if(m_original_problem->compare_fitness(f_hat_down, f_hat_up)) {
	//if(f_hat_up[0] > f_hat_down[0]) {
		scaling_factor = (f_hat_round[0] - f_hat_up[0]) / f_hat_up[0];
	}

	// apply penalty 2
	for(population::size_type i=0; i<infeasible_idx.size(); i++) {
		const population::size_type current_idx = infeasible_idx.at(i);

		m_fitness[current_idx][0] = m_fitness[current_idx][0] +
				scaling_factor * std::fabs(m_fitness[current_idx][0]) *
				( (std::exp(2. * solution_infeasibility.at(current_idx)) - 1.) / (std::exp(2.) - 1.) );
	}

	// this should be handled by the original algorithm itself

//	// sets the final fitness to the whole population (even the feasible?)
//	// scale all fitness values from 0 (worst) to absolute value of the best fitness
//	fitness_vector worstfit = m_fitness[0];
//	for(population::size_type i=0; i<pop_size; i++) {
//		if(m_original_problem->compare_fitness(worstfit, m_fitness[i])) {
//			worstfit = m_fitness[i];
//		}
//	}
//	for(population::size_type i=0; i<pop_size; i++) {
//		m_fitness[i][0] = worstfit[0] - m_fitness[i][0];
//	}
}

void self_adaptive::compute_solution_infeasibility(std::vector<double> &solution_infeasibility, const population &pop)
{
	// Let's store some useful variables.
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	constraint_vector c(m_original_problem->get_c_dimension(), 0.);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	constraint_vector c_scaling(m_original_problem->get_c_dimension(),0.);

	// evaluates the scaling factor
	for(population::size_type i=0; i<pop_size; i++) {
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
	solution_infeasibility.resize(pop_size);
	std::fill(solution_infeasibility.begin(),solution_infeasibility.end(),0.);

	for(population::size_type i=0; i<pop_size; i++) {
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
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			// test needed otherwise the c_scaling can be 0, and division by 0 occurs
			if(c_scaling[j] > 0.) {
				solution_infeasibility[i] += c[j]/c_scaling[j];
			}
		}
		solution_infeasibility[i] /= prob_c_dimension;
	}
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::self_adaptive);

