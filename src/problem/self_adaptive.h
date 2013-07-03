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

#ifndef PAGMO_PROBLEM_SELF_ADAPTIVE_H
#define PAGMO_PROBLEM_SELF_ADAPTIVE_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Constrainted self adaptive meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in death penalty constraints handling.
 *
 * The key idea of this constraint handling technique is to represent the
 * constraint violation by a single infeasibility measure, and to adapt
 * dynamically the penalization of infeasible solutions.
 *
 * @see R., & Wright, J. A. (2003). Self-adaptive fitness formulation for constrained optimization.
 * Evolutionary Computation, IEEE Transactions on, 7(5), 445-455 for the paper introducing the method.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE self_adaptive : public base
{
public:
	//constructors
	self_adaptive(const base & = cec2006(4), const population &pop = population(cec2006(4)));

	//copy constructor
	self_adaptive(const self_adaptive &);
	base_ptr clone() const;
	std::string get_name() const;

	//update
	void set_population(const population &pop);

protected:
	std::string human_readable_extra() const;
	void objfun_impl(fitness_vector &, const decision_vector &) const;

private:
	void compute_solution_infeasibility(std::vector<double> &solution_infeasibility, const population &pop);
	void update_fitness(const population &);

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
		ar & m_original_problem;
		ar & m_pop;
	}
	base_ptr m_original_problem;
	population m_pop;

	bool m_apply_penalty_1;
	double m_scaling_factor;
	std::vector<double> m_solution_infeasibility;

	std::vector<fitness_vector> m_fitness;

	population::size_type m_hat_down_idx;
	population::size_type m_hat_up_idx;
	population::size_type m_hat_round_idx;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::self_adaptive);

#endif // PAGMO_PROBLEM_self_adaptive_H
