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

#ifndef PAGMO_PROBLEM_SNOPT_TOYPROB_H
#define PAGMO_PROBLEM_SNOPT_TOYPROB_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo { namespace problem {

/// Test problem from SNOPT.
/**
 * Simple constrained minimisation test problem. Minimise
 * \f[
 * 	f(x,y) = y,
 * \f]
 * subject to
 * <center>
 * \f$
 * 	\begin{array}{rcl}
 * 		x^2 + 4y^2 &\leq & 4, \\
 * 		\left(x - 2\right)^2 + y^2 &\leq & 5.
 * 	\end{array}
 * \f$
 * </center>
 * The solution of this problem is in \f$ \left( x,y \right) = \left( 0,-1 \right) \f$. Search bounds are set to \f$ x \in \left[ 0,10 \right] \f$ and
 * \f$ y \in \left[ -10,10 \right] \f$ upon problem construction.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE snopt_toyprob: public base
{
	public:
		snopt_toyprob();
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		void set_sparsity(int &, std::vector<int>&, std::vector<int>&) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

}}

BOOST_CLASS_EXPORT_KEY(pagmo::problem::snopt_toyprob);

#endif
