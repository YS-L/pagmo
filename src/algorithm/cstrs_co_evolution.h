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

#ifndef PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H
#define PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "base.h"
#include "jde.h"
#include "sga.h"
#include "../problem/cstrs_co_evolution.h"

namespace pagmo { namespace algorithm {

/// Self-Adaptive Fitness constraints handling meta-algorithm
/**
 *
 * Seld-Adaptive Fitness constraints handling is a meta-algorithm that allow
 * to solve constrained optimization problems. The key idea of this constraint
 * handling technique is to represent the constraint violation by a single
 * infeasibility measure, and to adapt dynamically the penalization of infeasible solutions.
 *
 * This meta-algorithm is based on the problem self-adaptive.
 *
 * @see Farmani, R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
		
class __PAGMO_VISIBLE cstrs_co_evolution: public base
{
public:
	cstrs_co_evolution(const base & = jde(), const base & = sga(1), int pop_2_size = 30, int gen = 1,
					   problem::cstrs_co_evolution::method_type method = problem::cstrs_co_evolution::SIMPLE,
					   double pen_lower_bound = 0., double pen_upper_bound = 100000.);
	cstrs_co_evolution(const cstrs_co_evolution &);
	base_ptr clone() const;

public:
	void evolve(population &) const;
	std::string get_name() const;
	base_ptr get_algorithm() const;
	void set_algorithm(const base &);

protected:
	std::string human_readable_extra() const;

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_algo;
		ar & m_original_algo_2;
		ar & const_cast<int &>(m_gen);
		ar & m_method;
		ar & m_pop_2_size;
		ar & m_pen_lower_bound;
		ar & m_pen_upper_bound;
	}
	base_ptr m_original_algo;
	base_ptr m_original_algo_2;
	//Number of generations
	const int m_gen;
	// population 2 size
	int m_pop_2_size;
	// problem 2 variables
	problem::cstrs_co_evolution::method_type m_method;
	double m_pen_lower_bound;
	double m_pen_upper_bound;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::cstrs_co_evolution);

#endif // PAGMO_ALGORITHM_CSTRS_CO_EVOLUTION_H
