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

#ifndef PAGMO_PROBLEM_CSTRS_CO_EVOLUTION_H
#define PAGMO_PROBLEM_CSTRS_CO_EVOLUTION_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Constrainted co evolution meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in self adaptive constraints handling.
 *
 * The key idea of this constraint handling technique is to .
 *
 * @see R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE cstrs_co_evolution : public base
{
public:
	/// Type of co-evolution.
	/**
	* Definition of three types of co-evolution: SIMPLE, SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS.
	* The SIMPLE, is co-evolution defined by COELLO. The SPLIT_NEQ_EQ, splits equalities and
	* inequalities constraints (4 penalty coefficients). The SPLIT_CONSTRAINTS split the
	* number of coefficients upon the number of penlaty coefficients (2 * c_dimension).
	*/
	// co-evolution simple, split_neq_eq, split_constraints
	enum method_type {SIMPLE = 0, SPLIT_NEQ_EQ = 1, SPLIT_CONSTRAINTS = 2};

	//constructors
	cstrs_co_evolution(const base & = cec2006(4), const method_type = SIMPLE);

	//copy constructor
	cstrs_co_evolution(const cstrs_co_evolution &);
	base_ptr clone() const;
	std::string get_name() const;

	void set_penalty_coeff(const std::vector<double> &);
	int get_expected_penalty_coeff_size();

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;

private:
	void compute_penalty(std::vector<double> &, std::vector<int> &, const decision_vector &) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
		ar & m_penalty_coeff;
		ar & const_cast<method_type &>(m_method);
	}
	base_ptr m_original_problem;
	std::vector<double> m_penalty_coeff;

	const method_type m_method;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::cstrs_co_evolution);

#endif // PAGMO_PROBLEM_cstrs_co_evolution_H
