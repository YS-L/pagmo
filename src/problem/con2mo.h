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

#ifndef PAGMO_PROBLEM_CON2MO_H
#define PAGMO_PROBLEM_CON2MO_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "cec2006.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// Constrainted to Multi-Objective meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in multi-objective problem.
 *
 * Two implementations of the constrained to multi-objective are available. For a problem with m constraints,
 * m+1 objective functions are defined with the objective function as first objective and clipped positive
 * constraints values for others. The second one is the constrained to nulti-objective defined by
 * Coello Coello, where the objectives due to constraints includes number of violated constraints and
 * objective function as well.
 *
 * @see Coello Coello, C. A. (2002). Theoretical and numerical constraint-handling techniques used with evolutionary algorithms: a survey of the state of the art. Computer methods in applied mechanics and engineering, 191(11), 1245-1287.
 *
 * @see Coello, C. A. C. (2000). Treating constraints as objectives for single-objective evolutionary optimization. Engineering Optimization+ A35, 32(3), 275-308.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE con2mo : public base
{
	public:
		/// Type of constraints to multi-objective.
		/**
		* Definition of two types of constrained to multi-objective.
		* Simple distributes the objective and constraints to m+1 objective functions
		* defined with the objective function as first objective and clipped positive
		* constraints values for others. COELLO one is the same as Simple but
		* objectives due to constraints includes number of violated constraints and
		* objective function as well
		*/
		//death penalty type simple or kuri
		enum method_type {SIMPLE = 0, COELLO = 1};

		//constructors
		con2mo(const base & = cec2006(4), const method_type = SIMPLE);

		//copy constructor
		con2mo(const con2mo &);
		base_ptr clone() const;
		std::string get_name() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(fitness_vector &, const decision_vector &) const;

	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_original_problem;
			ar & const_cast<method_type &>(m_method);
		}
		base_ptr m_original_problem;

		const method_type m_method;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::con2mo);

#endif // PAGMO_PROBLEM_CON2MO_H
