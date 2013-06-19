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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../types.h"
#include "base.h"
#include "vega.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs a vega algorithm
 *
 * @param[in] algo pagmo::algorithm to use as multi-objective algorithm
 */
vega::vega(const base &algo):base(),
	m_mo_algo(algo.clone())
{

}

/// Copy constructor.
vega::vega(const vega &algo):base(algo),m_mo_algo(algo.clone())
{

}

/// Clone method.
base_ptr vega::clone() const
{
	return base_ptr(new vega(*this));
}

/// Evolve implementation.
/**
 * Run the vega algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void vega::evolve(population &pop) const
{
	// variables
	const problem::base &prob = pop.problem();

}

/// Algorithm name
std::string vega::get_name() const
{
	return m_mo_algo->get_name() + " [VEGA]";
}

/// Extra human readable algorithm info.
/**
				   * @return a formatted string displaying the parameters of the algorithm.
				   */
std::string vega::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_mo_algo->human_readable_extra() << std::endl;
	oss << "\n\tVEGA multi-objective algorithm";
	oss << std::endl;
	return oss.str();
}

}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::vega);
