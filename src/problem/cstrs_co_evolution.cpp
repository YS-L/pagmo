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
#include "cstrs_co_evolution.h"

namespace pagmo { namespace problem {
/**
 * Constructor of co-evolution problem using initial constrained problem
 *
 * Note: This problem is not inteded to be used by itself. Instead use the
 * co-evolution algorithm if you want to solve constrained problems.
 *
 * @param[in] problem base::problem to be modified to use a co-evolution
 * as constraints handling technique.
 * @param[in] method method_type to used for the co-evolution constraints
 * handling technique. Three posssibililties are available: SIMPLE,
 * SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS. The simple one is the original
 * version of the Coello/He implementation. The SPLIT_NEQ_EQ, splits the
 * equalities and inequalities constraints in two different sets for the
 * penalty weigths, containing respectively inequalities and equalities
 * weigths. The SPLIT_CONSTRAINTS splits the constraints in M set of weigths
 * with M the number of constraints.
 *
 */
cstrs_co_evolution::cstrs_co_evolution(const base &problem, const method_type method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_method(method)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co-evolution meta problem.");
	}

	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());

	m_penalty_coeff.resize(2);
	std::fill(m_penalty_coeff.begin(),m_penalty_coeff.end(),0.);
}

/// Copy Constructor. Performs a deep copy
cstrs_co_evolution::cstrs_co_evolution(const cstrs_co_evolution &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_penalty_coeff(prob.m_penalty_coeff),
	m_method(prob.m_method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Clone method.
base_ptr cstrs_co_evolution::clone() const
{
	return base_ptr(new cstrs_co_evolution(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the penalized fitness if the decision vector penalize with the penalty
 *  coefficient given.
 */
void cstrs_co_evolution::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_original_problem->objfun(f, x);

	std::vector<double> sum_viol;
	std::vector<int> num_viol;

	compute_penalty(sum_viol,num_viol,x);

	// assuming minimization
	switch(m_method)
	{
	case SIMPLE:
	{
		f[0] += sum_viol.at(0) * m_penalty_coeff.at(0) + double(num_viol.at(0)) * m_penalty_coeff.at(1);
		break;
	}
	case SPLIT_NEQ_EQ:
	{
		f[0] += sum_viol.at(0) * m_penalty_coeff.at(0) + double(num_viol.at(0)) * m_penalty_coeff.at(1);
		f[0] += sum_viol.at(1) * m_penalty_coeff.at(2) + double(num_viol.at(1)) * m_penalty_coeff.at(3);
		break;
	}
	case SPLIT_CONSTRAINTS:
	{
		int c_dimension = m_original_problem->get_c_dimension();
		for(int i=0; i<c_dimension; i++) {
			f[0] += sum_viol.at(0+i) * m_penalty_coeff.at(0+i*2) + double(num_viol.at(0+i)) * m_penalty_coeff.at(1+i*2);
		}
		break;
	}
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool cstrs_co_evolution::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string cstrs_co_evolution::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with co-evolution method ";
	oss << std::endl;
	return oss.str();
}

std::string cstrs_co_evolution::get_name() const
{
	return m_original_problem->get_name() + " [cstrs_co_evolution]";
}

/// Updates the fitness information based on the population.
/**
 *  By calling this method, penalties coefficients and
 *  fitnesses for the whole given population are computed.
 */
void cstrs_co_evolution::set_penalty_coeff(const std::vector<double> &penalty_coeff)
{
	switch(m_method)
	{
	case SIMPLE:
	{
		if(penalty_coeff.size() != 2) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 2.");
		}
		break;
	}
	case SPLIT_NEQ_EQ:
	{
		if(penalty_coeff.size() != 4) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 4.");
		}
		break;
	}
	case SPLIT_CONSTRAINTS:
		if(penalty_coeff.size() != 2 * m_original_problem->get_c_dimension()) {
			pagmo_throw(value_error,"The size of the penalty coefficient vector is not 2*number of constraints");
		}
		break;
	}
	m_penalty_coeff = penalty_coeff;
}

/// Returns the size of the penalty coefficient the problem expects
/// depending on the method used.
int cstrs_co_evolution::get_expected_penalty_coeff_size() {
	switch(m_method)
	{
	case SIMPLE:
	{
		return 2;
		break;
	}
	case SPLIT_NEQ_EQ:
	{
		return 4;
		break;
	}
	case SPLIT_CONSTRAINTS:
		return 2*m_original_problem->get_c_dimension();
		break;
	}
}

/// Computes the penalty depending on the provided penalty coefficient.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in,out] std::vector<double> sum_viol sum of violation vector.
 * @param[in,out] std::vector<double> num_viol number of violation vector.
 * @param[in] decision_vector x.
 */
void cstrs_co_evolution::compute_penalty(std::vector<double> &sum_viol, std::vector<int> &num_viol, const decision_vector &x) const
{
	// get the constraints dimension
	constraint_vector c(m_original_problem->get_c_dimension(), 0.);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	// updates the current constraint vector
	m_original_problem->compute_constraints(c,x);

	// sets the right definition of the constraints
	for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
		c[j] = std::abs(c.at(j)) - c_tol.at(j);
	}
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		c[j] = std::max(0.,c.at(j));
	}

	// updates the vectors depending on the method
	switch(m_method)
	{
	case SIMPLE:
	{
		sum_viol.resize(1);
		num_viol.resize(1);
		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			sum_viol[0] += c.at(j);
		}

		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[0] += 1;
			}
		}
		break;
	}
	case SPLIT_NEQ_EQ:
	{
		sum_viol.resize(2);
		num_viol.resize(2);
		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			sum_viol[0] += c.at(j);
		}
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			sum_viol[1] += c.at(j);
		}

		for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[0] += 1;
			}
		}
		for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[1] += 1;
			}
		}
		break;
	}
	case SPLIT_CONSTRAINTS:
	{
		sum_viol.resize(prob_c_dimension);
		num_viol.resize(prob_c_dimension);
		std::fill(sum_viol.begin(),sum_viol.end(),0.);
		std::fill(num_viol.begin(),num_viol.end(),0.);

		// update sum_num_viol
		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			sum_viol[j] += c.at(j);
		}

		for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
			if(!m_original_problem->test_constraint(c, j)) {
				num_viol[j] += 1;
			}
		}
		break;
	}
	}
}


/**
 * Constructor of co-evolution problem using initial constrained problem
 *
 * Note: This problem is not inteded to be used by itself. Instead use the
 * co-evolution algorithm if you want to solve constrained problems.
 *
 * @param[in] problem base::problem to be modified to use a co-evolution
 * as constraints handling technique.
 *
 */
cstrs_co_evolution_2::cstrs_co_evolution_2(const base &problem, int dimension):
	base(dimension,
		 problem.get_i_dimension(),
		 1,
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_sub_pop_2_x_vector(std::vector<decision_vector>(0)),
					   m_sub_pop_1_f_vector(std::vector< std::vector<double> >(0)),
					   m_feasible_count_vector(std::vector<int>(0)),
					   m_feasible_fitness_sum_vector(std::vector<double>(0)),
					   m_max_feasible_fitness(0.),
					   m_total_sum_viol(0.),
					   m_total_num_viol(0)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with co evolution meta problem.");
	}

	set_bounds(0.,10000.);
}

/// Copy Constructor. Performs a deep copy
cstrs_co_evolution_2::cstrs_co_evolution_2(const cstrs_co_evolution_2 &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_sub_pop_2_x_vector(prob.m_sub_pop_2_x_vector),
	m_sub_pop_1_f_vector(prob.m_sub_pop_1_f_vector),
	m_feasible_count_vector(prob.m_feasible_count_vector),
	m_feasible_fitness_sum_vector(prob.m_feasible_fitness_sum_vector),
	m_max_feasible_fitness(prob.m_max_feasible_fitness),
	m_total_sum_viol(prob.m_total_sum_viol),
	m_total_num_viol(prob.m_total_num_viol)
{
	set_bounds(prob.get_lb(),prob.get_ub());
}

/// Clone method.
base_ptr cstrs_co_evolution_2::clone() const
{
	return base_ptr(new cstrs_co_evolution_2(*this));
}

/// Implementation of the objective function.
/// (Wraps over the original implementation)
/**
 *  Returns the penalized fitness if the decision vector penalize with the penalty
 *  coefficient given.
 */
void cstrs_co_evolution_2::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	// search where the x vector is located
	int position=-1;
	for(int i=0; i<m_sub_pop_2_x_vector.size(); i++) {
		if(m_sub_pop_2_x_vector.at(i) == x) {
			position = i;
			break;
		}
	}

	if(position != -1) {
		// computes the fitness for population 2
		// case the solution containts at least one feasible individual
		if(m_feasible_count_vector.at(position) > 0) {
			f[0] = m_feasible_fitness_sum_vector[position] / double(m_feasible_count_vector[position]) - double(m_feasible_count_vector[position]);
		} else {
			f[0] = m_max_feasible_fitness +
					m_total_sum_viol[position]/double(m_total_num_viol[position]) -
					double(m_total_num_viol[position]);
		}
	}
	else {
		// what to do?
		f[0] = m_max_feasible_fitness;
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool cstrs_co_evolution_2::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string cstrs_co_evolution_2::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with co-evolution method ";
	oss << std::endl;
	return oss.str();
}

std::string cstrs_co_evolution_2::get_name() const
{
	return m_original_problem->get_name() + " [cstrs_co_evolution_2]";
}

/// Updates the fitness information based on the population.
/**
 *  By calling this method, penalties coefficients and
 *  fitnesses for the whole given population are computed.
 */
void cstrs_co_evolution_2::update_penalty_coeff(const std::vector<decision_vector> &pop_2_x,
												const std::vector< std::vector<decision_vector> > &sub_pop_1_x_vector,
												const std::vector< std::vector<fitness_vector> > &sub_pop_1_f_vector)
{
	m_sub_pop_2_x_vector = pop_2_x;

	population::size_type sub_pop_1_size = sub_pop_1_x_vector.at(0).size();
	population::size_type sub_pop_2_size = pop_2_x.size();

	m_feasible_count_vector.resize(sub_pop_2_size);
	m_feasible_fitness_sum_vector.resize(sub_pop_2_size);

	m_max_feasible_fitness = 0.;
	// computes the average fitness for the full population
	for(population::size_type j=0; j<sub_pop_2_size; j++) {
		const std::vector<decision_vector> &sub_pop_1_x = sub_pop_1_x_vector.at(j);
		const std::vector<fitness_vector> &sub_pop_1_f = sub_pop_1_f_vector.at(j);

		// computes the number of feasible solutions and their sum for the current population
		int feasible_count = 0;
		double feasible_fitness_sum = 0.;

		for(population::size_type i=0; i<sub_pop_1_size; i++) {
			const decision_vector &current_x = sub_pop_1_x[i];
			const fitness_vector &current_f = sub_pop_1_f[i];

			if(m_original_problem->feasibility_x(current_x)) {
				feasible_count++;
				feasible_fitness_sum += current_f[0];

				// computes max feasible value
				if(i==0) {
					m_max_feasible_fitness = current_f[0];
				} else {
					if(m_max_feasible_fitness < current_f[0]) {
						m_max_feasible_fitness = current_f[0];
					}
				}
			}
		}
		m_feasible_count_vector[j] = feasible_count;
		m_feasible_fitness_sum_vector[j] = feasible_fitness_sum;
	}

	m_total_sum_viol.resize(sub_pop_2_size);
	m_total_num_viol.resize(sub_pop_2_size);

	// computes the total sum_viol and num_viol
	for(population::size_type j=0; j<sub_pop_2_size; j++) {
		m_total_sum_viol[j] = 0.;
		m_total_num_viol[j] = 0;

		double sum_viol_temp;
		int num_viol_temp;
		for(population::size_type i=0; i<sub_pop_1_size; i++) {
			const std::vector<decision_vector> &sub_pop_1_x = sub_pop_1_x_vector.at(j);

			// the compute penalty needs to be recomputed (we can't use
			// the same for both pop1 and pop2 evolution)
			// we really should consider a way to store violation information
			// in the population as well :(
			compute_penalty(sum_viol_temp,num_viol_temp,sub_pop_1_x[i]);
			m_total_sum_viol[j] += sum_viol_temp;
			m_total_num_viol[j] += num_viol_temp;
		}
	}
}

/// Computes the penalty depending on the provided penalty coefficient.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in,out] std::vector<double> sum_viol sum of violation vector.
 * @param[in,out] std::vector<double> num_viol number of violation vector.
 * @param[in] decision_vector x.
 */
void cstrs_co_evolution_2::compute_penalty(double &sum_viol, int &num_viol, const decision_vector &x) const
{
	// get the constraints dimension
	constraint_vector c(m_original_problem->get_c_dimension(), 0.);
	problem::base::c_size_type prob_c_dimension = m_original_problem->get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	// updates the current constraint vector
	m_original_problem->compute_constraints(c,x);

	// sets the right definition of the constraints
	for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
		c[j] = std::abs(c.at(j)) - c_tol.at(j);
	}
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		c[j] = std::max(0.,c.at(j));
	}

	// update sum_num_viol
	sum_viol = 0.;
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		sum_viol += c.at(j);
	}

	num_viol = 0;
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		if(!m_original_problem->test_constraint(c, j)) {
			num_viol += 1;
		}
	}
}


}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_co_evolution);
BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_co_evolution_2);
