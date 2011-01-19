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
 *   the Free Software Foundation; either version 3 of the License, or       *
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

#ifndef PAGMO_ALGORITHM_PYTHON_BASE_H
#define PAGMO_ALGORITHM_PYTHON_BASE_H

#include <boost/python/class.hpp>
#include <string>

#include "../../src/algorithm/base.h"
#include "../../src/config.h"
#include "../../src/exceptions.h"
#include "../../src/population.h"
#include "../../src/python_locks.h"
#include "../../src/serialization.h"
#include "../utils.h"

namespace pagmo { namespace algorithm {

// Wrapper for implementing algorithms in Python.
class __PAGMO_VISIBLE python_base: public base, public boost::python::wrapper<base>
{
	public:
		python_base():base() {}
		python_base(const base &p):base(p) {}
		base_ptr clone() const
		{
			gil_state_lock lock;
			base_ptr retval = this->get_override("__copy__")();
			if (!retval) {
				pagmo_throw(std::runtime_error,"algorithms's __copy__() method returns a NULL pointer, please check the implementation");
			}
			return retval;
		}
		std::string human_readable_extra() const
		{
			gil_state_lock lock;
			if (boost::python::override f = this->get_override("human_readable_extra")) {
				return f();
			}
			return base::human_readable_extra();
		}
		std::string default_human_readable_extra() const
		{
			return this->base::human_readable_extra();
		}
		std::string get_name() const
		{
			gil_state_lock lock;
			if (boost::python::override f = this->get_override("get_name")) {
				return f();
			}
			return base::get_name();
		}
		std::string default_get_name() const
		{
			return this->base::get_name();
		}
		bool is_thread_safe() const
		{
			// All calls re-implementable from Python are protected by a mutex and hence thread-safe.
			// TODO: correct this to make it match with the comment above :)
			return false;
		}
		void evolve(population &p) const
		{
			p = py_evolve(p);
		}
		// Changed implementations from Python.
		population py_evolve(const population &p) const
		{
			gil_state_lock lock;
			if (boost::python::override f = this->get_override("evolve")) {
				const population retval = f(p);
				return retval;
			}
			pagmo_throw(not_implemented_error,"the __evolve__() method has not been implemented");
		}
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
			ar & boost::serialization::base_object<boost::python::wrapper<base> >(*this);
		}
};

} }

BOOST_CLASS_EXPORT(pagmo::algorithm::python_base);

#endif
