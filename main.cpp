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

#include <iostream>
#include <iomanip>
#include "src/pagmo.h"

using namespace pagmo;

// Example in C++ of the use of PaGMO 1.1.4

int main()
{
	pagmo::problem::cec2006 prob_constrained(7);

	//pagmo::algorithm::monte_carlo algo(1); //only one generation for the algo!
	//pagmo::algorithm::sga algo(1); //only one generation for the algo!
	pagmo::algorithm::sga_gray algo(1,0.9,0.04, 1000,
							   algorithm::sga_gray::mutation::GAUSSIAN, 0.1,
							   algorithm::sga_gray::selection::ROULETTE,
							   algorithm::sga_gray::crossover::EXPONENTIAL); //only one generation for the algo!
	//pagmo::algorithm::cmaes algo(1); //only one generation for the algo!
	//pagmo::algorithm::de algo(1); //only one generation for the algo!
	//pagmo::algorithm::pso algo(1); //only one generation for the algo!
	pagmo::algorithm::self_adaptive algo_constrained(algo, 5000);

	std::cout << algo_constrained;

	for (size_t i=0; i<20; ++i) {
		pagmo::population pop(prob_constrained,70);
		algo_constrained.evolve(pop);
		std::cout << pop.champion();
	}

	std::cout << algo_constrained << std::endl;
	std::cout << prob_constrained << std::endl;

//	pagmo::island isl = island(algo_constrained, prob_constrained, 70);

//	for (size_t i=0; i<20; ++i){
//		isl.evolve(1);
//		std::cout << isl.get_population().champion();
//	}

	return 0;
}
