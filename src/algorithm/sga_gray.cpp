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

#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/round.hpp>
#include <string>
#include <vector>
#include <algorithm>

#include "../exceptions.h"
#include "../population.h"
#include "../types.h"
#include "base.h"
#include "sga_gray.h"
#include "../problem/base_stochastic.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Allows to specify in detail all the parameters of the algorithm.
 *
 * @param[in] gen Number of generations to evolve.
 * @param[in] cr Crossover probability (of each allele if binomial crossover)
 * @param[in] m Mutation probability (of each allele)
 * @param[in] elitism The best individual is reinserted in the population each elitism generations
 * @param[in] mut Mutation type. One of sga_gray::mutation::GAUSSIAN, sga_gray::mutation::RANDOM
 * @param[in] width Mutation width. When gaussian mutation is selected is the width of the mutation
 * @param[in] sel Selection type. One of sga_gray::selection::BEST20, sga_gray::selection::ROULETTE
 * @param[in] cro Crossover type. One of sga_gray::crossover::BINOMIAL, sga_gray::crossover::EXPONENTIAL
 * @throws value_error if gen is negative, crossover probability is not \f$ \in [0,1]\f$, mutation probability is not \f$ \in [0,1]\f$,
 * elitism is <=0
 *
 */
sga_gray::sga_gray(int gen, const double &cr, const double &m, int elitism, mutation::type mut, double width, selection::type sel, crossover::type cro)
	:base(),m_gen(gen),m_cr(cr),m_m(m),m_elitism(elitism),m_mut(mut,width),m_sel(sel),m_cro(cro)
{
	if (gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if (cr > 1 || cr < 0) {
		pagmo_throw(value_error,"crossover probability must be in the [0,1] range");
	}
	if (m < 0 || m > 1) {
		pagmo_throw(value_error,"mutation probability must be in the [0,1] range");
	}
	if (elitism < 1) {
		pagmo_throw(value_error,"elitisim must be greater than zero");
	}
	if (width <0 || width >1) {
		pagmo_throw(value_error,"mutation width must be in the [0,1] range");
	}

	m_bit_encoding = 4;
	m_max_encoding_integer = int(std::pow(2., m_bit_encoding));

}

/// Clone method.
base_ptr sga_gray::clone() const
{
	return base_ptr(new sga_gray(*this));
}

/// Evolve implementation.
/**
 * Run the simple genetic algorithm for the number of generations specified in the constructors.
 * At each improvment velocity is also updated.
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */

void sga_gray::evolve(population &pop) const
{
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << double_to_binary(0,0,15) << std::endl;
//	std::cout << double_to_binary(1,0,15) << std::endl;
//	std::cout << double_to_binary(2,0,15) << std::endl;
//	std::cout << double_to_binary(3,0,15) << std::endl;
//	std::cout << double_to_binary(4,0,15) << std::endl;
//	std::cout << double_to_binary(5,0,15) << std::endl;
//	std::cout << double_to_binary(6,0,15) << std::endl;
//	std::cout << double_to_binary(7,0,15) << std::endl;
//	std::cout << double_to_binary(8,0,15) << std::endl;
//	std::cout << double_to_binary(9,0,15) << std::endl;
//	std::cout << double_to_binary(10,0,15) << std::endl;
//	std::cout << double_to_binary(11,0,15) << std::endl;
//	std::cout << double_to_binary(12,0,15) << std::endl;
//	std::cout << double_to_binary(13,0,15) << std::endl;
//	std::cout << double_to_binary(14,0,15) << std::endl;
//	std::cout << double_to_binary(15,0,15) << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << binary_to_double(double_to_binary(0,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(1,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(2,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(3,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(4,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(5,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(6,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(7,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(8,0,15),0,15) << std::endl;
//	std::cout << binary_to_double(double_to_binary(9,0,15),0,15) << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << binary_to_double(double_to_binary(0,0,4),0,4) << std::endl;
//	std::cout << binary_to_double(double_to_binary(1,0,4),0,4) << std::endl;
//	std::cout << binary_to_double(double_to_binary(2,0,4),0,4) << std::endl;
//	std::cout << binary_to_double(double_to_binary(3,0,4),0,4) << std::endl;
//	std::cout << binary_to_double(double_to_binary(4,0,4),0,4) << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << double_to_binary(-1,-1,1) << std::endl;
//	std::cout << double_to_binary(-0.5,-1,1) << std::endl;
//	std::cout << double_to_binary(0,-1,1) << std::endl;
//	std::cout << double_to_binary(0.5,-1,1) << std::endl;
//	std::cout << double_to_binary(1,-1,1) << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << "------------  GRAY  ----------------" << std::endl;
//	std::cout << "------------------------------------" << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(0,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(1,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(2,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(3,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(4,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(5,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(6,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(7,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(8,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(9,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(10,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(11,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(12,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(13,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(14,0,15))) << std::endl;
//	std::cout << gray_to_binary(binary_to_gray(double_to_binary(15,0,15))) << std::endl;
//	std::cout << "------------------------------------" << std::endl;

	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const problem::base::size_type D = prob.get_dimension(), Di = prob.get_i_dimension(), prob_c_dimension = prob.get_c_dimension(), prob_f_dimension = prob.get_f_dimension();
	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
	const population::size_type NP = pop.size();
	const problem::base::size_type Dc = D - Di;

	//We perform some checks to determine wether the problem/population are suitable for sga_gray
	if ( prob_c_dimension != 0 ) {
		pagmo_throw(value_error,"The problem is not box constrained and sga_gray is not suitable to solve it");
	}

	if ( prob_f_dimension != 1 ) {
		pagmo_throw(value_error,"The problem is not single objective and sga_gray is not suitable to solve it");
	}

	if (NP < 5) {
		pagmo_throw(value_error,"for sga_gray at least 5 individuals in the population are needed");
	}

	// Get out if there is nothing to do.
	if (m_gen == 0) {
		return;
	}
	// Some vectors used during evolution are allocated here.
	decision_vector dummy(D,0);			//used for initialisation purposes
	std::vector<decision_vector > X(NP,dummy), Xnew(NP,dummy);

	std::vector<fitness_vector > fit(NP);		//fitness

	fitness_vector bestfit;
	decision_vector bestX(D,0);

	std::vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);
	std::vector <int> selection(NP);

	int tempID;
	std::vector<int> fitnessID(NP);

	// Initialise the phenotypes and their fitness to that of the initial deme
	for (pagmo::population::size_type i = 0; i<NP; i++ ) {
		X[i]	=	pop.get_individual(i).cur_x;
		fit[i]	=	pop.get_individual(i).cur_f;
	}

	// Find the best member and store in bestX and bestfit
	double bestidx = pop.get_best_idx();
	bestX = pop.get_individual(bestidx).cur_x;
	bestfit = pop.get_individual(bestidx).cur_f;


	// Main sga_gray loop
	for (int j = 0; j<m_gen; j++) {

		switch (m_sel) {
		case selection::BEST20: { //selects the best 20% and puts multiple copies in Xnew
			//Sort the individuals according to their fitness
			for (pagmo::population::size_type i=0; i<NP; i++) fitnessID[i]=i;
			for (pagmo::population::size_type i=0; i < (NP-1); ++i) {
				for (pagmo::population::size_type j=i+1; j<NP; ++j) {
					if ( prob.compare_fitness(fit[j],fit[i]) ) {
						//swap fitness values
						fit[i].swap(fit[j]);
						//swap id's
						tempID = fitnessID[i];
						fitnessID[i] = fitnessID[j];
						fitnessID[j] = tempID;
					}
				}
			}
			int best20 = NP/5;
			for (pagmo::population::size_type i=0; i<NP; ++i) {
				selection[i] = fitnessID[i % best20];
			}
			break;
		}

		case selection::ROULETTE: {
			//We scale all fitness values from 0 (worst) to absolute value of the best fitness
			fitness_vector worstfit=fit[0];
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if (prob.compare_fitness(worstfit,fit[i])) worstfit=fit[i];
			}

			for (pagmo::population::size_type i = 0; i < NP; i++) {
				selectionfitness[i] = fabs(worstfit[0] - fit[i][0]);
			}

			// We build and normalise the cumulative sum
			cumsumTemp[0] = selectionfitness[0];
			for (pagmo::population::size_type i = 1; i< NP; i++) {
				cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
			}
			for (pagmo::population::size_type i = 0; i < NP; i++) {
				cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
			}

			//we throw a dice and pick up the corresponding index
			double r2;
			for (pagmo::population::size_type i = 0; i < NP; i++) {
				r2 = m_drng();
				for (pagmo::population::size_type j = 0; j < NP; j++) {
					if (cumsum[j] > r2) {
						selection[i]=j;
						break;
					}
				}
			}
			break;
		}
		}

		//Xnew stores the new selected generation of genotypes
		for (pagmo::population::size_type i = 0; i < NP; i++) {
			Xnew[i]=X[selection[i]];
		}


//		std::vector<int> dummy(Xnew.at(0).size() * m_bit_encoding, 0);
//		std::vector< std::vector<int> > encoded_decision_vector(Xnew.size(),dummy);

//		for(std::vector::size_type i=0; i<Xnew.size(); i++) {
//			decision_vector &current_x = Xnew.at(i);
//			encoded_decision_vector[i] = encode(current_x, lb, ub);
//		}

//		// encode Xnew
//		encode(Xnew);

		//2 - Crossover
		{
			int r1,L;
			decision_vector  member1,member2;

			for (pagmo::population::size_type i=0; i< NP; i++) {
				//for each chromosome selected i.e. in Xnew
				member1 = Xnew[i];
				//we select a mating patner different from the self (i.e. no masturbation)
				do {
					r1 = boost::uniform_int<int>(0,NP - 1)(m_urng);
				} while ( r1 == boost::numeric_cast<int>(i) );
				member2 = Xnew[r1];
				//and we operate crossover
				switch (m_cro) {
				//0 - binomial crossover
				case crossover::BINOMIAL: {
					size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
					for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
						if ((m_drng() < m_cr) || L + 1 == D) { /* change at least one parameter */
							member1[n] = member2[n];
						}
						n = (n+1)%D;
					}
					break; }
					//1 - exponential crossover
				case crossover::EXPONENTIAL: {
					size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
					L = 0;
					do {
						member1[n] = member2[n];
						n = (n+1) % D;
						L++;
					}  while ( (m_drng() < m_cr) && (L < boost::numeric_cast<int>(D)) );
					break; }
				}
				Xnew[i] = member1;

			} }

		//3 - Mutation
		switch (m_mut.m_type) {
		case mutation::GAUSSIAN: {
			boost::normal_distribution<double> dist;
			boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > delta(m_drng,dist);
			for (pagmo::problem::base::size_type k = 0; k < Dc;k++) { //for each continuous variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = (delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
					}
				}
			}
			for (pagmo::problem::base::size_type k = Dc; k < D;k++) { //for each integer variable
				double std = (ub[k]-lb[k]) * m_mut.m_width;
				for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
					if (m_drng() < m_m) {
						double mean = Xnew[i][k];
						double tmp = boost::math::iround(delta() * std + mean);
						if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) Xnew[i][k] = tmp;
					}
				}
			}
			break;
		}
		case mutation::RANDOM: {
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				for (pagmo::problem::base::size_type j = 0; j < Dc;j++) { //for each continuous variable
					if (m_drng() < m_m) {
						Xnew[i][j] = boost::uniform_real<double>(lb[j],ub[j])(m_drng);
					}
				}
				for (pagmo::problem::base::size_type j = Dc; j < D;j++) {//for each integer variable
					if (m_drng() < m_m) {
						Xnew[i][j] = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
					}
				}
			}
			break;
		}
		}

		// decode
		decode(Xnew);

		// If the problem is a stochastic optimization chage the seed and re-evaluate taking care to update also best and local bests
		try
		{
			//4 - Evaluate the new population (stochastic problem)
			dynamic_cast<const pagmo::problem::base_stochastic &>(prob).set_seed(m_urng());
			pop.clear(); // Removes memory based on different seeds (champion and best_x, best_f, best_c)
			
			// We re-evaluate the best individual (for elitism)
			prob.objfun(bestfit,bestX);
			// Re-evaluate wrt new seed the particle position and memory
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				// We evaluate here the new individual fitness
				prob.objfun(fit[i],Xnew[i]);
				// We update the velocity (in case coupling with PSO via archipelago)
				//dummy = Xnew[i];
				//std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				///We now set the cleared pop. cur_x is the best_x, re-evaluated with new seed.
				pop.push_back(Xnew[i]);
				//pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		catch (const std::bad_cast& e)
		{
			//4 - Evaluate the new population (deterministic problem)
			for (pagmo::population::size_type i = 0; i < NP;i++) {
				prob.objfun(fit[i],Xnew[i]);
				dummy = Xnew[i];
				std::transform(dummy.begin(), dummy.end(), pop.get_individual(i).cur_x.begin(), dummy.begin(),std::minus<double>());
				//updates x and v (cache avoids to recompute the objective function and constraints)
				pop.set_x(i,Xnew[i]);
				pop.set_v(i,dummy);
				if (prob.compare_fitness(fit[i], bestfit)) {
					bestfit = fit[i];
					bestX = Xnew[i];
				}
			}
		}
		
		//5 - Reinsert best individual every m_elitism generations
		if (j % m_elitism == 0) {
			int worst=0;
			for (pagmo::population::size_type i = 1; i < NP;i++) {
				if ( prob.compare_fitness(fit[worst],fit[i]) ) worst=i;
			}
			Xnew[worst] = bestX;
			fit[worst] = bestfit;
			dummy = Xnew[worst];
			std::transform(dummy.begin(), dummy.end(), pop.get_individual(worst).cur_x.begin(), dummy.begin(),std::minus<double>());
			//updates x and v (cache avoids to recompute the objective function)
			pop.set_x(worst,Xnew[worst]);
			pop.set_v(worst,dummy);
		}
		X = Xnew;
	} // end of main sga_gray loop
}

/// Algorithm name
std::string sga_gray::get_name() const
{
	return "A Simple Genetic Algorithm";
}

/// Extra human readable algorithm info.
/**
 * Will return a formatted string displaying the parameters of the algorithm.
 */
std::string sga_gray::human_readable_extra() const
{
	std::ostringstream s;
	s << "gen:" << m_gen << ' ';
	s << "CR:" << m_cr << ' ';
	s << "M:" << m_m << ' ';
	s << "elitism:" << m_elitism << ' ';
	s << "mutation:";
	switch (m_mut.m_type) {
	case mutation::RANDOM: {
		s << "RANDOM ";
		break;
	}
	case mutation::GAUSSIAN: {
		s << "GAUSSIAN (" << m_mut.m_width << ") ";
		break;
	}
	}
	s << "selection:";
	switch (m_sel) {
	case selection::BEST20: {
		s << "BEST20 ";
		break;
	}
	case selection::ROULETTE: {
		s << "ROULETTE ";
		break;
	}
	}
	s << "crossover:";
	switch (m_cro) {
	case crossover::EXPONENTIAL: {
		s << "EXPONENTIAL ";
		break;
	}
	case crossover::BINOMIAL: {
		s << "BINOMIAL ";
		break;
	}
	}

	return s.str();
}

std::vector<int> sga_gray::selection(const std::vector<fitness_vector> &pop_f, const problem::base &prob) const
{
	int method = 1;

	population::size_type NP = pop_f.size();

	std::vector<int> selection(NP);

	switch (method) {
	case 0: {
		//selects the best 20% and puts multiple copies in Xnew
		// in practice, for performance reasons, the selection, fitnessID should not be created each time
		// thus we might prefer to use a class of operators
		int tempID;
		std::vector<int> fitnessID(NP);

		for (population::size_type i=0; i<NP; i++) {
			fitnessID[i]=i;
		}
		for (population::size_type i=0; i < (NP-1); ++i) {
			for (population::size_type j=i+1; j<NP; ++j) {
				//if ( prob.compare_fitness(pop_f[fitnessID[j]],pop_f[fitnessID[i]]) ) {
				if ( pop_f[fitnessID[j]]<pop_f[fitnessID[i]]) {
					//swap fitness values
					// fit[i].swap(fit[j]);
					//swap id's
					tempID = fitnessID[i];
					fitnessID[i] = fitnessID[j];
					fitnessID[j] = tempID;
				}
			}
		}

		// fitnessID now contains the position of individuals ranked from best to worst
		int best20 = NP/5;
		for (pagmo::population::size_type i=0; i<NP; ++i) {
			selection[i] = fitnessID[i % best20]; // multiple copies
		}
		break;
	}
	case 1: { // ROULETTE
		std::vector<double> selectionfitness(NP), cumsum(NP), cumsumTemp(NP);

		//We scale all fitness values from 0 (worst) to absolute value of the best fitness
		fitness_vector worstfit=pop_f[0];
		for (population::size_type i = 1; i < NP;i++) {
			//if (prob.compare_fitness(worstfit,pop_f[i])) worstfit=pop_f[i];
			if (worstfit<pop_f[i]) worstfit=pop_f[i];
		}

		for (population::size_type i = 0; i < NP; i++) {
			selectionfitness[i] = fabs(worstfit[0] - pop_f[i][0]);
		}

		// We build and normalise the cumulative sum
		cumsumTemp[0] = selectionfitness[0];
		for (population::size_type i = 1; i< NP; i++) {
			cumsumTemp[i] = cumsumTemp[i - 1] + selectionfitness[i];
		}
		for (population::size_type i = 0; i < NP; i++) {
			cumsum[i] = cumsumTemp[i]/cumsumTemp[NP-1];
		}

		//we throw a dice and pick up the corresponding index
		double r2;
		for (population::size_type i = 0; i < NP; i++) {
			r2 = m_drng();
			for (population::size_type j = 0; j < NP; j++) {
				if (cumsum[j] > r2) {
					selection[i]=j;
					break;
				}
			}
		}
		break;
		}
	}
	return selection;
}

void sga_gray::crossover(std::vector<decision_vector> &pop_x) const
{
	int cross_method=0;

	// crossover probability
	const double crossover_rate = 0.8;

	population::size_type NP = pop_x.size();
	fitness_vector::size_type D = pop_x.at(0).size();

	int r1,L;
	decision_vector  member1,member2;

	for (population::size_type i=0; i<NP; i++) {
		//for each chromosome selected i.e. in pop_x
		member1 = pop_x[i];

		//we select a mating patner different from the self (i.e. no masturbation)
		do {
			r1 = boost::uniform_int<int>(0,NP - 1)(m_urng);
		} while ( r1 == boost::numeric_cast<int>(i) );
		member2 = pop_x[r1];

		// and we operate crossover
		switch (cross_method) {
		//0 - binomial crossover
		case 0: {
			size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
			for (size_t L = 0; L < D; ++L) { /* perform D binomial trials */
				if ((m_drng() < crossover_rate) || L + 1 == D) { /* change at least one parameter */
					member1[n] = member2[n];
				}
				n = (n+1)%D;
			}
			break;
		}
		// 1 - exponential crossover
		case 1: {
			size_t n = boost::uniform_int<int>(0,D-1)(m_urng);
			L = 0;
			do {
				member1[n] = member2[n];
				n = (n+1) % D;
				L++;
			}  while ( (m_drng() < crossover_rate) && (L < boost::numeric_cast<int>(D)) );
			break;
		}
		}
		pop_x[i] = member1;
	}
}

void sga_gray::mutate(std::vector<decision_vector> &pop_x, const problem::base &prob) const
{
	int method=0;

	double width=0.1;
	double mutation_rate = 0.02;

	const problem::base::size_type D = pop_x.at(0).size();

	const problem::base::size_type Di = 0.;
	const double lb[] = {1.,1000};
	const double ub[] = {1.,1000};

	const problem::base::size_type Dc = D - Di;

	const population::size_type NP = pop_x.size();

//	const problem::base::size_type Di = prob.get_i_dimension();
//	const decision_vector &lb = prob.get_lb(), &ub = prob.get_ub();
//	const population::size_type NP = pop_x.size();
//	const problem::base::size_type Dc = D - Di;

	switch (method) {
	case 0: { // GAUSSIAN
		boost::normal_distribution<double> dist;
		boost::variate_generator<boost::lagged_fibonacci607 &, boost::normal_distribution<double> > delta(m_drng,dist);
		for (problem::base::size_type k = 0; k < Dc;k++) { //for each continuous variable
			double std = (ub[k]-lb[k]) * width;
			for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
				if (m_drng() < mutation_rate) {
					double mean = pop_x[i][k];
					double tmp = (delta() * std + mean);
					if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) pop_x[i][k] = tmp;
				}
			}
		}
		for (problem::base::size_type k = Dc; k < D;k++) { //for each integer variable
			double std = (ub[k]-lb[k]) * width;
			for (pagmo::population::size_type i = 0; i < NP;i++) { //for each individual
				if (m_drng() < mutation_rate) {
					double mean = pop_x[i][k];
					double tmp = boost::math::iround(delta() * std + mean);
					if ( (tmp < ub[k]) &&  (tmp > lb[k]) ) pop_x[i][k] = tmp;
				}
			}
		}
		break;
		}
	case 1: { // RANDOM
		for (population::size_type i = 0; i < NP;i++) {
			for (pagmo::problem::base::size_type j = 0; j < Dc;j++) { //for each continuous variable
				if (m_drng() < mutation_rate) {
					pop_x[i][j] = boost::uniform_real<double>(lb[j],ub[j])(m_drng);
				}
			}
			for (problem::base::size_type j = Dc; j < D;j++) {//for each integer variable
				if (m_drng() < mutation_rate) {
					pop_x[i][j] = boost::uniform_int<int>(lb[j],ub[j])(m_urng);
				}
			}
		}
		break;
		}
	}
}



std::vector<int> sga_gray::double_to_binary(const double &number, const double &lb, const double &ub) const
{
	// initialize the vector of size m_bit_encoding
	std::vector<int> binary(m_bit_encoding, 0);

	// convert the current number into its binary representation considering the domain
	// available
	int temp_number = (number - lb) * (m_max_encoding_integer - 1) / (ub - lb) ;

	// store the binary number
	int position = 0;
	while (temp_number!=0)
	{
		if( (temp_number % 2) == 0 ) {
			binary[position] = 0;
		} else {
			binary[position] = 1;
		}
		temp_number = (int)std::floor(temp_number/2);
		position++;
	}
	// reverse the order as this algorithm provides the reverse binary reprentation
	std::reverse(binary.begin(),binary.end());

	return binary;
}

double sga_gray::binary_to_double(const std::vector<int> &binary, const double &lb, const double &ub) const
{
	// find the representation of the binary number in the integer domain
	int temp_number = 0;
	for(int i=0; i<m_bit_encoding; i++) {
		temp_number += binary.at(m_bit_encoding - 1 - i) * std::pow(2.,i);
	}

	// rescaling back into the domain double domain
	double number = temp_number * (ub - lb) / (m_max_encoding_integer - 1) + lb;

	return number;
}


std::vector<int> sga_gray::gray_to_binary(const std::vector<int> &gray) const
{
	// uses the XOR table
	int length = gray.size();

	std::vector<int> binary = gray;

	// the MSB is the same
	for(int i=1; i<length; i++) {
		binary[i] = (binary[i-1]+gray[i]) % 2;
	}

	return binary;
}

std::vector<int> sga_gray::binary_to_gray(const std::vector<int> &binary) const
{
	int length = binary.size();
	std::vector<int> gray = binary;

	// the MSB is the same
	// uses the XOR table
	for(int i = 1; i < length; i++) {
		gray[i] = binary[i]^binary[i-1];
	}

	return gray;
}

std::vector<int> sga_gray::encode_decision(const decision_vector &x, const decision_vector &lb, const decision_vector &ub) const
{
	std::vector<int> encoded_x(x.size() * m_bit_encoding, 0);

	for(decision_vector::size_type i=0; i<x.size(); i++) {
		std::vector<int> encoded_gene = binary_to_gray(double_to_binary(x.at(i),lb.at(i),ub.at(i)));
		// copy the gene at the right location
		for(decision_vector::size_type j=0; j<m_bit_encoding; j++) {
			encoded_x[i*m_bit_encoding + j] = encoded_gene[j];
		}
	}

	return encoded_x;
}
void sga_gray::decode(std::vector<decision_vector> &x_vector) const
{

}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::sga_gray);
