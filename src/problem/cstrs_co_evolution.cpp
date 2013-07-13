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
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a self-adaptive
 * as constraints handling technique.
 *
 */
cstrs_co_evolution::cstrs_co_evolution(const base &problem):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone())
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	// check that the dimension of the problem is 1
	if (m_original_problem->get_f_dimension() != 1) {
		pagmo_throw(value_error,"The original fitness dimension of the problem must be one, multi objective problems can't be handled with self adaptive meta problem.");
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
	m_penalty_coeff(prob.m_penalty_coeff)
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
 *  Returns the penalized fitness if the decision vector is found in the
 *  given population or the non penalized objective function of the underlying
 *  problem otherwise.
 */
void cstrs_co_evolution::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	m_original_problem->objfun(f, x);

	double coeff = 0.;
	int viol = 0;

	compute_penalty(coeff,viol,x);

	// assuming minimization
	f[0] += coeff * m_penalty_coeff.at(0) + double(viol) * m_penalty_coeff.at(1);
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
	m_penalty_coeff = penalty_coeff;
}

/// Computes the solution infeasibility measure.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in,out] std::vector<double solution infeasibility vector to update.
 * @param[in] population pop.
 */
void cstrs_co_evolution::compute_penalty(double &coeff, int &viol, const decision_vector &x) const
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

	// update coeff
	coeff = 0.;
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		coeff += c.at(j);
	}

	viol = 0;
	for(problem::base::c_size_type j=0; j<prob_c_dimension; j++) {
		if(!m_original_problem->test_constraint(c, j)) {
			viol += 1;
		}
	}
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cstrs_co_evolution);

