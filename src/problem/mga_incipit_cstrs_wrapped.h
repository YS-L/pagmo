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

#ifndef PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_WRAPPED_H
#define PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_WRAPPED_H

#include <string>

#include "../config.h"
#include "../serialization.h"
#include "../types.h"
#include "../types.h"
#include "base.h"
#include "../keplerian_toolbox/planet_js.h"
#include "../keplerian_toolbox/epoch.h"


namespace pagmo{ namespace problem {

/// The beginning of the GTOC6 Jupiter Capture Trajectory

/**
 *
 * A PyGMO global optimization problem (box-bounded, continuous) representing a capture
 * in the Jupiter system. Two constraints are considered: 1. Closest distance, 2. Time of flight
 *
 * Decision vector:
 * [t0,u,v,T0] + [sin(beta1), cos(beta1), rp1/rP1, eta1,T1] + .... 
 * 
 * @author Dario Izzo (dario.izzo@esa.int)
 */
class __PAGMO_VISIBLE mga_incipit_cstrs_wrapped: public base
{
	public:
		mga_incipit_cstrs_wrapped(const std::vector<kep_toolbox::planet_ptr> = construct_default_sequence(),
			 const kep_toolbox::epoch t0_l = kep_toolbox::epoch(7305.0),
			 const kep_toolbox::epoch t0_u = kep_toolbox::epoch(11323.0),
			 const std::vector<std::vector<double> > tof = construct_default_tofs(),
			 double Tmax = 365.25,
			 double Dmin = 0.2
			 );
		mga_incipit_cstrs_wrapped(const mga_incipit_cstrs_wrapped&);
		base_ptr clone() const;
		
		std::string get_name() const;
		std::string pretty(const std::vector<double> &x) const;
		void set_tof(const std::vector<std::vector<double> >&);
		const std::vector<std::vector<double> >& get_tof() const;
		std::vector<kep_toolbox::planet_ptr> get_sequence() const;
	protected:
		void objfun_impl(fitness_vector &, const decision_vector &) const;
		void compute_constraints_impl(constraint_vector &, const decision_vector &) const;
		std::string human_readable_extra() const;
		
	private:
		static const std::vector<kep_toolbox::planet_ptr> construct_default_sequence() {
			std::vector<kep_toolbox::planet_ptr> retval;
			retval.push_back(kep_toolbox::planet_js("io").clone());
			retval.push_back(kep_toolbox::planet_js("io").clone());
			retval.push_back(kep_toolbox::planet_js("europa").clone());
			return retval;
		}
		static const std::vector<std::vector<double> > construct_default_tofs() {
			std::vector<std::vector<double> > retval;
			std::vector<double> dumb(2);
			dumb[0] = 100;dumb[1] = 200;
			retval.push_back(dumb);
			dumb[0] = 3;dumb[1] = 200;
			retval.push_back(dumb);
			dumb[0] = 4;dumb[1] = 100;
			retval.push_back(dumb);
			return retval;
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & m_seq;
			ar & m_tof;
			ar & const_cast<double &>(m_tmax);
			ar & const_cast<double &>(m_dmin);
		}
		std::vector<kep_toolbox::planet_ptr> m_seq;
		std::vector<std::vector<double> > m_tof;
		const double m_tmax;
		const double m_dmin;
};

}} // namespaces

BOOST_CLASS_EXPORT_KEY(pagmo::problem::mga_incipit_cstrs_wrapped)
#endif // PAGMO_PROBLEM_MGA_INCIPIT_CSTRS_WRAPPED_H
