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

#ifndef KEPLERIAN_TOOLBOX_LAMBERT_PROBLEM_H
#define KEPLERIAN_TOOLBOX_LAMBERT_PROBLEM_H

#include<vector>

#include "astro_constants.h"


namespace kep_toolbox {

/// Lambert Problem (3 dimensional)
/**
 * This class represent a Lambert's problem. When instantiated it assumes a prograde orbit (unless otherwise stated)
 * and evaluates all the solutions to the problem calling accordingly the lambert_3D routine with all possible rev numbers
 * After the object is instantiated the solutions can be retreived using the appropriate getters. Note that the
 * number of solutions will be N_max*2 + 1, where N_max is the maximum number of revolutions as evaluated by the
 * lambert_find_N routine.
 *
 * NOTE: The class has been tested extensively via monte carlo runs checked with numerical propagation. After 500000
 * lambert solvers a maximum error of 10^-8 has been found defined as the magnitude of the difference between r2 and
 * the same vector as evaluated propagating numerically with the v1 found. Average number of iterations is 7. Speed is
 * 3.2 sec to call lambert_3d 500374 times. Tests have been made on a 2.16Ghz Intel Core Duo
 *
 * @see kep_toolbox::lambert_2d
 * @see kep_toolbox::lambert_3d
 * @see kep_toolbox::lambert_find_N
 *
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class lambert_problem
{
public:
	friend std::ostream &operator<<(std::ostream &, const lambert_problem &);
	lambert_problem(const array3D &r1, const array3D &r2, const double &tof, const double& mu = 1., const int &cw = 0);
	const std::vector<array3D>& get_v1() const;
	const std::vector<array3D>& get_v2() const;
	const std::vector<double>& get_a() const;
	const std::vector<double>& get_p() const;
	const bool& is_reliable() const;
	int get_Nmax() const;
private:
	const array3D m_r1, m_r2;
	const double m_tof;
	const double m_mu;
	int m_lw;
	std::vector<array3D> m_v1;
	std::vector<array3D> m_v2;
	std::vector<int> m_iters;
	std::vector<double> m_a;
	std::vector<double> m_p;
	double m_s,m_c;
	int m_Nmax;
	bool m_has_converged;

};
std::ostream &operator<<(std::ostream &, const lambert_problem &);
} //namespaces

#endif // KEPLERIAN_TOOLBOX_LAMBERT_PROBLEM_H
