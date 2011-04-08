/*****************************************************************************
 *   Copyright (C) 2004-2009 The PaGMO development team,                     *
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

#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "kur.h"

namespace pagmo { namespace problem {

/**
 * Will construct Kursawe's study problem.
 *
 * @see problem::base constructors.
 */
kur::kur():base(3,0,2)
{
	// Set bounds.
	set_lb(-5.0);
	set_ub(5.0);
}

/// Clone method.
base_ptr kur::clone() const
{
	return base_ptr(new kur(*this));
}

/// Implementation of the objective function.
void kur::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	pagmo_assert(f.size() == 2);
	pagmo_assert(x.size() == 3);

	f[0] = 0;
	f[1] = 0;

	f[0] += -10 * exp( -0.2 * sqrt(x[0]*x[0] - x[1]*x[1]));
	f[1] += pow(x[0],0.8) + 5 * sin(pow(x[0],3));
	for(problem::base::size_type i = 1; i < 3; ++i) {
		f[0] += -10 * exp( -0.2 * sqrt(x[i-1]*x[i-1] - x[i]*x[i]));
		f[1] += pow(x[i],0.8) + 5 * sin(pow(x[i],3));
	}
}

std::string kur::get_name() const
{
	return "Kursawe's study";
}
}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::kur);
