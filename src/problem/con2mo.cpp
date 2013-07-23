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
#include "base.h"
#include "con2mo.h"

static int __mo_dimension__(const pagmo::problem::base &original_problem,
							const pagmo::problem::con2mo::method_type method)
{
	if( method > 2 || method < 0) {
		pagmo_throw(value_error, "The constrained to multi-objective method must be set to OBJ_CSTRS for Coello type constrained to multi-objective, to OBJ_CSTRSVIO for COMOGO type constrained to nobj+1 objectives problem or to OBJ_EQVIO_INEQVIO for COMOGO type constrained to nobj+2 objectives problem.");
	}

	if(original_problem.get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	switch(method)
	{
	case pagmo::problem::con2mo::OBJ_CSTRS:
	{
		return original_problem.get_f_dimension() + original_problem.get_c_dimension();
		break;
	}
	case pagmo::problem::con2mo::OBJ_CSTRSVIO:
	{
		return original_problem.get_f_dimension() + 1;
		break;
	}
	case pagmo::problem::con2mo::OBJ_EQVIO_INEQVIO:
	{
		return original_problem.get_f_dimension() + 2;
		break;
	}
	default:
	{
		return 0.;
		break;
	}
	}
}

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a constrained to
 * multi-objective handling technique.
 * @param[in] method method_type to be modified to use a single constrained
 * to multi-objective approach defined with OBJ_CSTRS, OBJ_CSTRSVIO or OBJ_EQVIO_INEQVIO
 *
 */
con2mo::con2mo(const base &problem, const method_type method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 __mo_dimension__(problem, method),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_method(method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Copy Constructor. Performs a deep copy
con2mo::con2mo(const con2mo &prob):
	base((int)prob.get_dimension(),
		 prob.get_i_dimension(),
		 prob.get_f_dimension(),
		 prob.get_c_dimension(),
		 prob.get_ic_dimension(),
		 prob.get_c_tol()),
	m_original_problem(prob.m_original_problem->clone()),
	m_method(prob.m_method)
{
	set_bounds(m_original_problem->get_lb(),m_original_problem->get_ub());
}

/// Clone method.
base_ptr con2mo::clone() const
{
	return base_ptr(new con2mo(*this));
}

/// Implementation of the objective functions.
/// (Wraps over the original implementation)
void con2mo::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	constraint_vector c(m_original_problem->get_c_dimension(),0.);
	m_original_problem->compute_constraints(c,x);

	// modify equality constraints to behave as inequality constraints:
	c_size_type number_of_eq_constraints = m_original_problem->get_c_dimension() -
			m_original_problem->get_ic_dimension();

	const std::vector<double> &c_tol = m_original_problem->get_c_tol();

	for(c_size_type i=0; i<number_of_eq_constraints; i++) {
		c[i] = std::abs(c.at(i)) - c_tol.at(i);
	}

	decision_vector original_f(m_original_problem->get_f_dimension(),0.);
	m_original_problem->objfun(original_f,x);

	fitness_vector::size_type original_nbr_obj = original_f.size();
	constraint_vector::size_type number_of_constraints = c.size();

	// clean the fitness vector
	for(f_size_type i=0; i<f.size(); i++) {
		f[i] = 0.;
	}

	// in all cases, the first objectives holds the initial objectives
	for(f_size_type i=0; i<original_nbr_obj; i++) {
		f[i] = original_f.at(i);
	}

	switch(m_method)
	{
	case OBJ_CSTRS:
	{
		constraint_vector::size_type number_of_violated_constraints = 0;

		// computes the number of satisfied constraints
		for(c_size_type i=0; i<number_of_constraints; i++){
			if(!m_original_problem->test_constraint(c,i))
				number_of_violated_constraints += 1;
		}

		for(c_size_type i=0; i<number_of_constraints; i++) {
			if(c.at(i) > 0.) {
				f[original_nbr_obj+i] = c.at(i);
			} else if(number_of_violated_constraints != 0) {
				f[original_nbr_obj+i] = number_of_violated_constraints;
			} else {
				f[original_nbr_obj+i] = 0.;
				for(f_size_type j=0; j<original_nbr_obj; j++) {
					f[original_nbr_obj+i] += original_f.at(j);
				}
			}
		}
		break;
	}
	case OBJ_CSTRSVIO:
	{
		for(c_size_type i=0; i<number_of_constraints; i++) {
			if(c.at(i) > 0.) {
				f[original_nbr_obj] += c.at(i);
			}
		}
		break;
	}
	case OBJ_EQVIO_INEQVIO:
	{
		// treating equality constraints
		for(c_size_type i=0; i<number_of_eq_constraints; i++) {
			if(c.at(i) > 0.) {
				f[original_nbr_obj] += c.at(i);
			}
		}

		for(c_size_type i=number_of_eq_constraints; i<number_of_constraints; i++) {
			if(c.at(i) > 0.) {
				f[original_nbr_obj+1] += c.at(i);
			}
		}
		break;
	}
	default:
		pagmo_throw(value_error, "Error: There are only 3 methods for the constrained to multi-objective!");
		break;
	}
}

/// Implementation of fitness vectors comparison.
/**
 * @brief compare_fitness_impl calls the compare_fitness method of the original problem.
 * @return true if v_f1 is dominating v_f2, false otherwise.
 */
bool con2mo::compare_fitness_impl(const fitness_vector &v_f1, const fitness_vector &v_f2) const
{
	return m_original_problem->compare_fitness(v_f1,v_f2);
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string con2mo::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with constrained to multi-objective, method ";
	switch(m_method){
	case OBJ_CSTRS: {
		oss << "OBJ_CSTRS ";
		break;
	}
	case OBJ_CSTRSVIO: {
		oss << "OBJ_CSTRSVIO ";
		break;
	}
	case OBJ_EQVIO_INEQVIO: {
		oss << "OBJ_EQVIO_INEQVIO ";
		break;
	}
	}
	oss << std::endl;
	return oss.str();
}

std::string con2mo::get_name() const
{
	std::string method;

	switch(m_method){
	case OBJ_CSTRS: {
		method = "OBJ_CSTRS ";
		break;
	}
	case OBJ_CSTRSVIO: {
		method = "OBJ_CSTRSVIO ";
		break;
	}
	case OBJ_EQVIO_INEQVIO: {
		method = "OBJ_EQVIO_INEQVIO ";
		break;
	}
	}
	return m_original_problem->get_name() + " [con2mo, method_" + method + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::con2mo);

