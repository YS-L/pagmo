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

#ifndef PAGMO_PROBLEM_ANTIBODIES_PROBLEM_H
#define PAGMO_PROBLEM_ANTIBODIES_PROBLEM_H

#include <string>
#include <boost/functional/hash.hpp>
#include <boost/serialization/map.hpp>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"
#include "../algorithm/immune_system.h"

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

class __PAGMO_VISIBLE antibodies_problem : public base
{
public:
	// hamming distance, euclidean distance
	enum method_type {HAMMING = 0, EUCLIDEAN = 1};

public:
	//constructors
	antibodies_problem(const base & = cec2006(4), const method_type = HAMMING);

	//copy constructor
	antibodies_problem(const antibodies_problem &);
	base_ptr clone() const;
	std::string get_name() const;

	void set_antigens(const std::vector<decision_vector> &);

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;
	bool compare_fitness_impl(const fitness_vector &, const fitness_vector &) const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
		ar & m_pop_antigens;
		ar & const_cast<method_type &>(m_method);
		ar & m_bit_encoding;
		ar & m_max_encoding_integer;
	}

	base_ptr m_original_problem;
	std::vector<decision_vector> m_pop_antigens;

	const method_type m_method;

	// encoding size
	int m_max_encoding_integer;
	int m_bit_encoding;

	// function to compute the distance
	double compute_distance(const decision_vector &x) const;

	// function for the hamming distance
	std::vector<int> double_to_binary(const double &number, const double &lb, const double &ub) const;

};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::antibodies_problem);

#endif // PAGMO_PROBLEM_antibodies_problem_H
