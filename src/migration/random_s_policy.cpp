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

#include <algorithm>
#include <vector>

#include "../population.h"
#include "base.h"
#include "base_s_policy.h"
#include "random_s_policy.h"
#include "../exceptions.h"

namespace pagmo { namespace migration {

/// Constructor from migration rate and type.
/**
 * @param[in] rate migration rate.
 * @param[in] type migration rate type.
 *
 * @see base_s_policy::base_s_policy.
 */
random_s_policy::random_s_policy(const double &rate, rate_type type):base_s_policy(rate,type) {}

base_s_policy_ptr random_s_policy::clone() const
{
	return base_s_policy_ptr(new random_s_policy(*this));
}

std::vector<population::individual_type> random_s_policy::select(population &pop) const
{
	pagmo_assert(get_n_individuals(pop) <= pop.size());
	// Gets the number of individuals to select
	const population::size_type migration_rate = get_n_individuals(pop);
	
	// Create a temporary array of individuals.
	std::vector<population::individual_type> result;
	
	// Create an array of indices
	std::vector<population::size_type> candidates_idx(boost::numeric_cast<std::vector<population::size_type>::size_type>(pop.size()));
	
	// Fill it with indices
	iota(candidates_idx.begin(),candidates_idx.end(),population::size_type(0));

	// shuffle
	random_shuffle(candidates_idx.begin(),candidates_idx.end());
	
	// Selects the n first individuals
	for (population::size_type i = 0; i< migration_rate; ++i) {
		result.push_back(pop.get_individual(candidates_idx[i]));
	}
	
	return result;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::migration::random_s_policy);
