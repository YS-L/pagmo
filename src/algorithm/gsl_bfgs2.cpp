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

#include <gsl/gsl_multimin.h>
#include <string>

#include "../population.h"
#include "gsl_bfgs2.h"
#include "gsl_gradient.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Will invoke internally the constructor from algorithm::gsl_gradient with the specified parameters.
 *
 * @see gsl_gradient::gsl_gradient().
 */
gsl_bfgs2::gsl_bfgs2(int max_iter, const double &grad_tol, const double &numdiff_step_size, const double &step_size, const double &tol):
	gsl_gradient(max_iter,grad_tol,numdiff_step_size,step_size,tol) {}

/// Clone method.
base_ptr gsl_bfgs2::clone() const
{
	return base_ptr(new gsl_bfgs2(*this));
}

const gsl_multimin_fdfminimizer_type *gsl_bfgs2::get_gsl_minimiser_ptr() const
{
	return gsl_multimin_fdfminimizer_vector_bfgs2;
}

}}

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::gsl_bfgs2);
