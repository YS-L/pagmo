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
#include <string>
#include <vector>

#include "../exceptions.h"
#include "../population.h"
#include "../problem/base.h"
#include "../problem/cstrs_co_evolution.h"
#include "../types.h"
#include "base.h"
#include "cstrs_co_evolution.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an self adaptive algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @throws value_error if stop is negative
 */
cstrs_co_evolution::cstrs_co_evolution(const base &original_algo, int gen, int pop_2_size):base(),m_gen(gen),m_pop_2_size(pop_2_size)
{
	m_original_algo = original_algo.clone();
	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
	if( pop_2_size <= 0 ) {
		pagmo_throw(value_error,"the population of the population 2 be greater than 0");
	}
}

/// Copy constructor.
cstrs_co_evolution::cstrs_co_evolution(const cstrs_co_evolution &algo):base(algo),m_original_algo(algo.m_original_algo->clone()),m_gen(algo.m_gen),m_pop_2_size(algo.m_pop_2_size)
{}

/// Clone method.
base_ptr cstrs_co_evolution::clone() const
{
	return base_ptr(new cstrs_co_evolution(*this));
}

/// Evolve implementation.
/**
 * Run the Self-Adaptive algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void cstrs_co_evolution::evolve(population &pop) const
{	
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for Self-Adaptive
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and Self-Adaptive is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if (pop_size == 0) {
		return;
	}

	// split the population into two populations
	// the population P1 associated to the modified problem with penalized fitness
	// the reason P2 is not a population type is that only one iteration of prob2
	// is necessary and it is not easy to locally update the fitness

	// creates the problem 1
	// TODO: CHANGE WITH CO EVOL
	problem::cstrs_co_evolution prob_1(prob);

	//population pop_2;

	// sub population size
	population::size_type sub_pop_1_size = pop_size;
	population::size_type sub_pop_2_size = m_pop_2_size;

	// initialize the population 1
	population sub_pop_1(prob_1,0);
//	population sub_pop_1(sub_pop_1_size);
	sub_pop_1.clear();
	for(population::size_type i=0; i<sub_pop_1_size; i++) {
		sub_pop_1.push_back(pop.get_individual(i).cur_x);
	}

	// initialize the population 2, composed of individuals only
	std::vector<decision_vector> sub_pop_2_x(sub_pop_2_size);
	std::vector<fitness_vector> sub_pop_2_f(sub_pop_2_size);

	// population 2 decision vector is 2*number of constraints
	int pop_2_x_size = 2;
	// int pop_2_x_size = 2 * prob_c_dimension;

	decision_vector dummy(pop_2_x_size,10.); //used for initialisation purposes
	for(population::size_type i=0; i<sub_pop_2_size; i++) {
		sub_pop_2_x[i] = decision_vector(pop_2_x_size,0.);

		// choose random coefficients between 1 and 1000
		for(int j = 0; j<pop_2_x_size;j++)
		{
			sub_pop_2_x[i][j] = boost::uniform_real<double>(1.,1000.)(m_drng);
		}
		sub_pop_2_f[i] = fitness_vector(1);
	}

	std::vector<decision_vector> sub_pop_2_x_new(sub_pop_2_size,dummy);

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {
		std::cout << "current generation: " << k << std::endl;

		// for each individuals of pop 2

		for(population::size_type j=0; j<sub_pop_2_size; j++) {
			// set decision vector encoding penalty coefficients
			// w1 and w2 in prob 1
			prob_1.set_penalty_coeff(sub_pop_2_x.at(j));

			// reevaluates the population 1 with these new coefficients
			prob_1.reset_caches();
			for(population::size_type i=0; i<sub_pop_2_size; i++) {
				sub_pop_1.set_x(i,sub_pop_1.get_individual(i).cur_x);
			}

			// use an instance of the population
			population pop_1_instance(sub_pop_1);

			// evolve the population 1 instance
			// DO WE NEED TO CLONE THE ORIGINAL
			// ALGO TO AVOID KEEPING THE ALGO HISTORY AS WELL?
			m_original_algo->evolve(pop_1_instance);

			// store indexes of feasible individuals
			std::vector<population::size_type> feasible_idx;

			feasible_idx.clear();
			for(population::size_type i=0; i<sub_pop_1_size; i++) {
				const population::individual_type &current_individual = pop_1_instance.get_individual(i);

				if(prob.feasibility_x(current_individual.cur_x)) {
					feasible_idx.push_back(i);
				}
			}

			// computes the average fitness
			int feasible_size = feasible_idx.size();
			double average_fitness=0.;
			for(int i=0; i<feasible_size;i++) {
				const population::individual_type &current_individual = pop_1_instance.get_individual(feasible_idx.at(i));

				average_fitness += current_individual.cur_f.at(0);
			}
			average_fitness /= feasible_size;

			// average scaling (NOT DETAILED IN COELLO PAPER)
//			for(population::size_type i=1; i<sub_pop_1_size; i++) {
//				const population::individual_type &current_individual = pop_1_instance.get_individual(i);

//				if(prob.compare_fitness(worst_fitness , current_individual.cur_f)) {
//					feasible_idx.push_back(i);
//				}
//			}

			double average_scaling = 1.;
			average_fitness = average_fitness * average_scaling + feasible_size;

			// assuming minimization (convert to objective function)
			average_fitness = - average_fitness;

			// store the average fitness to the individual j
			sub_pop_2_f[j][0] = average_fitness;

			// how do we select the population pop ?
			std::cout << pop_1_instance.champion() << std::endl;
		}

		// need to evolve the population 2 now
		// selection process, returns a vector
		// of position containing the selected members
		std::vector<int> sub_pop_2_s_idx = selection(sub_pop_2_x, sub_pop_2_f, prob);

		// sub_pop_2_x_new stores the new selected generation of chromosomes
		for (pagmo::population::size_type i = 0; i < sub_pop_2_size; i++) {
			sub_pop_2_x_new[i] = sub_pop_2_x[sub_pop_2_s_idx[i]];
		}

		// crossover the new population
		crossover(sub_pop_2_x_new);

		// mutate the new population
		mutate(sub_pop_2_x_new,prob);

		// elitism?

		// we don't need to reevaluate the population 2 as only the decision vector is used
		sub_pop_2_x = sub_pop_2_x_new;

		for(population::size_type i=0; i < sub_pop_2_size; i++)
			std::cout << "sub_pop_2_x"<< i << " " << sub_pop_2_x[i];

//		// updates the population 1
//		for(int j=0; j<sub_pop_1_size; j++) {
//			sub_pop_1[]
//		}
	}

	// we evaluate the population pop

	std::cout << pop.champion() << std::endl;
}

/// Algorithm name
std::string cstrs_co_evolution::get_name() const
{
	return m_original_algo->get_name() + "[Co-Evolution]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr cstrs_co_evolution::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void cstrs_co_evolution::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string cstrs_co_evolution::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithm: " << m_original_algo->get_name() << ' ';
	s << "\n\tConstraints handled with co-evolution algorithm";
	return s.str();
}


std::vector<int> cstrs_co_evolution::selection(const std::vector<decision_vector> &pop_x, const std::vector<fitness_vector> &pop_f, const problem::base &prob) const
{
	int method = 1;

	population::size_type NP = pop_x.size();

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

void cstrs_co_evolution::crossover(std::vector<decision_vector> &pop_x) const
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

void cstrs_co_evolution::mutate(std::vector<decision_vector> &pop_x, const problem::base &prob) const
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

// Coello is calling fixed point representation what others
// call


//// --------------------------------------------
//function y = cma_genophenotransform(gp, x)
//// input argument: tlist (an object), x vector to be transformed
//// elements in gp :
////         typical_x = Nx1 vector of middle of domains, default == 0
////         scaling = Nx1 vector, default == 1
//  if gp.skip
//    y = x;
//  else
//    y = zeros(gp.typical_x); // only to set the size
//    y(gp.scaling > 0) = x;
//    y = gp.typical_x + y .* gp.scaling;
//  end
//endfunction

//// --------------------------------------------
//function x = cma_phenogenotransform(gp, x)
//// input argument: tlist (or mlist) according to Gilles Pujol
//  if ~gp.skip
//    idx = gp.scaling > 0;
//    x = (x(idx) - gp.typical_x(idx)) ./ gp.scaling(idx);
//  end
//endfunction


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::cstrs_co_evolution);
