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

namespace pagmo { namespace problem {

/**
 * Constructor using initial constrained problem
 *
 * @param[in] problem base::problem to be modified to use a constrained to
 * multi-objective handling technique.
 * @param[in] method method_type to be modified to use a simple constrained
 * to multi-objective if defined with SIMPLE and a Coello constrained to
 * multi-objective with COELLO.
 *
 */
con2mo::con2mo(const base &problem, const method_type method):
	base((int)problem.get_dimension(),
		 problem.get_i_dimension(),
		 problem.get_f_dimension() + problem.get_c_dimension(),
		 0,
		 0,
		 0.),
	m_original_problem(problem.clone()),
	m_method(method)
{
	if(m_original_problem->get_c_dimension() <= 0){
		pagmo_throw(value_error,"The original problem has no constraints.");
	}

	if( method > 1 || method < 0) {
		pagmo_throw(value_error, "the constrained to multi-objective method must be set to 0 for simple constrained to multi-objective or to 1 for Coello constrained to multi-objective problem.");
	}

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

	decision_vector original_f(m_original_problem->get_f_dimension(),0.);
	m_original_problem->objfun(original_f,x);

	fitness_vector::size_type original_nbr_obj = original_f.size();
	constraint_vector::size_type number_of_constraints = c.size();

	// in all cases, the first objectives holds the initial objectives
	for(f_size_type i=0; i<original_nbr_obj; i++) {
		f[i] = original_f.at(i);
	}

	switch(m_method)
	{
	case SIMPLE:
	{
		for(c_size_type i=0; i<number_of_constraints; i++) {
			f[original_nbr_obj+i] = std::max(0., c.at(i));
		}
		break;
	}
	case COELLO:
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
	default:
		pagmo_throw(value_error, "Error: There are only 2 methods for the constraints to multi-objective!");
		break;
	}
}

/// Extra human readable info for the problem.
/**
 * Will return a formatted string containing the type of constraint handling
 */
std::string con2mo::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	oss << "\n\tConstraints handled with constraints to multi-objective, method ";
	switch(m_method){
	case SIMPLE: {
		oss << "SIMPLE ";
		break;
	}
	case COELLO: {
		oss << "COELLO ";
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
	case SIMPLE: {
		method = "SIMPLE ";
		break;
	}
	case COELLO: {
		method = "COELLO ";
		break;
	}
	}
	return m_original_problem->get_name() + " [con2mo, method_" + method + "]";
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::con2mo);

