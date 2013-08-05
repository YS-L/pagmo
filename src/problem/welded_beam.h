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

#ifndef PAGMO_PROBLEM_WELDED_BEAM_H
#define PAGMO_PROBLEM_WELDED_BEAM_H

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "base.h"

namespace pagmo{ namespace problem {

/// The welded beam design problem: Constrained Real-Parameter Optimization
/**
 *
 * This class instanciates the welded beam design problem used to test
 * constrained problems. It is constrained, continuous, single objective problems
 *
 * @see http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page4679.htm
 *
 * @author Jeremie Labroquere (jeremie.labroquere@gmail.com)
 */
class __PAGMO_VISIBLE welded_beam : public base
{
public:
	welded_beam();
	base_ptr clone() const;
	std::string get_name() const;

protected:
	void objfun_impl(fitness_vector &, const decision_vector &) const;
	void compute_constraints_impl(constraint_vector &, const decision_vector &) const;

private:
	void initialize_best(void);

	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<base>(*this);
	}
};

}} //namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::welded_beam);

#endif
