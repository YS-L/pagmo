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

#ifndef PAGMO_PROBLEM_ZDT3_H
#define PAGMO_PROBLEM_ZDT3_H

#include <string>

#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// ZDT3 problem
/**
 *
 * This is a box-constrained continuous 30-dimension multi-objecive problem.
 * \f[
 *	g\left(x\right) = 1 + 9 \left(\sum_{i=2}^{30} x_i \right) / \left( n-1 \right)
 * \f]
 * \f[
 * 	F_1 \left(x\right) = x_1
 * \f]
 * \f[
 *      F_2 \left(x\right) = g(x) \left[ 1 - \sqrt{x_1 / g(x)} - x_1/g(x) \sin(10 \pi x_1) \right]  x \in \left[ 0,1 \right].
 *
 * \f]
 *
 * @see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.18.4257&rep=rep1&type=pdf
 * @author Andrea Mambrini (andrea.mambrini@gmail.com)
 */

class __PAGMO_VISIBLE zdt3 : public base
{
	public:
		zdt3();
		base_ptr clone() const;
		std::string get_name() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_pi;
		}
		double m_pi;
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::zdt3);

#endif // PAGMO_PROBLEM_ZDT3_H
