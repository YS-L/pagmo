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

#ifndef PAGMO_ALGORITHM_VEGA_H
#define PAGMO_ALGORITHM_VEGA_H

#include <string>

#include "../config.h"
#include "../population.h"
#include "../serialization.h"
#include "sga.h"
#include "base.h"

namespace pagmo { namespace algorithm {

/// VEGA based multi-objective meta-algorithm
/**
 * Implements a meta-algorithm class that wraps some other algorithms,
 * resulting in multi-objective algorithm based on VEGA.
 *
 * @see Schaffer, J. D. (1985, July). Multiple objective optimization with vector evaluated genetic algorithms. In Proceedings of the 1st international Conference on Genetic Algorithms (pp. 93-100). L. Erlbaum Associates Inc.
 *
 * @see Deb Kalyanmoy (2001, June). Multi-Objective Optimization Using Evolutionary Algorithms.
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */

class __PAGMO_VISIBLE vega: public base
{
	public:
		// constructors
		vega(const base & = sga());

		// copy constructors
		vega(const vega &);
		base_ptr clone() const;
		std::string get_name() const;

		void evolve(population &) const;

	protected:
		std::string human_readable_extra() const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_mo_algo;
		}
		base_ptr m_mo_algo;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::algorithm::vega);

#endif // PAGMO_ALGORITHM_VEGA_H
