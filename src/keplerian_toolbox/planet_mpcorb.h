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

#ifndef PLANET_MPCORB_H
#define PLANET_MPCORB_H

#include <string>

#include "../serialization.h"
#include "planet.h"


namespace kep_toolbox{

/// Minor Planet (keplerian)
/**
 * This class derives from the planet class and allow to instantiate planets of
 * from the MPCORB database using their names or row id. The file MPCORB.DAT is searched
 * in the current directory.
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class planet_mpcorb : public planet
{
public:
	/**
	 * Construct a minor planet from a line of the MPCORB.DAT file. Default value is the MPCORB.DAT line
	 * for the dwarf planet Ceres.
	 * \param[in] name a string containing one line of MPCORB.DAT
	 */
	planet_mpcorb(const std::string & = "00001    3.34  0.12 K107N 113.41048   72.58976   80.39321   10.58682  0.0791382  0.21432817   2.7653485  0 MPO110568  6063  94 1802-2006 0.61 M-v 30h MPCW       0000      (1) Ceres              20061025");
	planet_ptr clone() const;
	static epoch packed_date2epoch(std::string);

private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & boost::serialization::base_object<planet>(*this);
	}
	static int packed_date2number(char c);
};


} /// End of namespace kep_toolbox

BOOST_CLASS_EXPORT(kep_toolbox::planet_mpcorb);

#endif // PLANET_MPCORB_H
