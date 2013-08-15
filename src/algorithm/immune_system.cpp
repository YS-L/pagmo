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
#include "../problem/antibodies_problem.h"
#include "../problem/unconstrain.h"
#include "../types.h"
#include "base.h"
#include "immune_system.h"

namespace pagmo { namespace algorithm {

/// Constructor.
/**
 * Constructs an immune system constraints handling algorithm
 *
 * @param[in] original_algo pagmo::algorithm to use as 'original' optimization method
 * @param[in] original_algo_penalties pagmo::algorithm to use as 'original' optimization method
 * for population representing the penaltiesweights
 * @param[in] pop_penalties_size population size for the penalty encoding population.
 * @param[in] gen number of generations.
 * @param[in] problem::immune_system::method_type the method used for the population 2.
 * Three posssibililties are available: SIMPLE, SPLIT_NEQ_EQ and SPLIT_CONSTRAINTS.
 * The simple one is the original version of the Coello/He implementation. The SPLIT_NEQ_EQ,
 * splits the equalities and inequalities constraints in two different sets for the
 * penalty weigths, containing respectively inequalities and equalities weigths. The
 * SPLIT_CONSTRAINTS splits the constraints in M set of weigths with M the number of
 * constraints.
 * @param[in] pen_lower_bound the lower boundary used for penalty.
 * @param[in] pen_upper_bound the upper boundary used for penalty.
 * @throws value_error if stop is negative
 */
immune_system::immune_system(const base &original_algo, const base &original_algo_immune,
							 int gen, method_type method,double ftol, double xtol):
	base(),m_gen(gen),m_method(method),m_ftol(ftol),m_xtol(xtol)
{
	m_original_algo = original_algo.clone();
	m_original_algo_immune = original_algo_immune.clone();

	if(gen < 0) {
		pagmo_throw(value_error,"number of generations must be nonnegative");
	}
}

/// Copy constructor.
immune_system::immune_system(const immune_system &algo):
	base(algo),m_original_algo(algo.m_original_algo->clone()),
	m_original_algo_immune(algo.m_original_algo_immune->clone()),m_gen(algo.m_gen),
	m_method(algo.m_method),m_ftol(algo.m_ftol),m_xtol(algo.m_xtol)
{}

/// Clone method.
base_ptr immune_system::clone() const
{
	return base_ptr(new immune_system(*this));
}

/// Evolve implementation.
/**
 * Run the co-evolution algorithm
 *
 * @param[in,out] pop input/output pagmo::population to be evolved.
 */
void immune_system::evolve(population &pop) const
{	
	// Let's store some useful variables.
	const problem::base &prob = pop.problem();
	const population::size_type pop_size = pop.size();
	const problem::base::size_type prob_dimension = prob.get_dimension();

	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = prob.get_c_dimension();

	//We perform some checks to determine wether the problem/population are suitable for co-evolution
	if(prob_c_dimension < 1) {
		pagmo_throw(value_error,"The problem is not constrained and co-evolution is not suitable to solve it");
	}
	if(prob.get_f_dimension() != 1) {
		pagmo_throw(value_error,"The problem is multiobjective and co-evolution is not suitable to solve it");
	}

	// Get out if there is nothing to do.
	if(pop_size == 0) {
		return;
	}

	// generates the unconstrained problem
	problem::unconstrain prob_unconstrained(prob);

	// associates the population to this problem
	population pop_mixed(prob_unconstrained);
	std::vector<decision_vector> pop_mixed_c(pop_size);

	// initializaton of antigens vector and antibodies population

	// TODO: CHANGE
	//population::size_type pop_antigens_size = 30; // typically half of the feasible solution
	//population::size_type pop_antibodies_size = 10; // typically 1/3 of the antigens

	double phi = 0.5; // feasible fraction selection to compute the mean value
	double gamma = 0.5; // number of antigens selected / number of feasibles
	double sigma = 1./3.; // number of antibodies / number of antigens

	// antigens population
	std::vector<decision_vector> pop_antigens;

	// vector containing the best antigens position
	std::vector<population::size_type> pop_antigens_pool;
	std::vector<population::size_type> pop_antibodies_pool;

	// Main Co-Evolution loop
	for(int k=0; k<m_gen; k++) {

		std::cout << "iter: " << k << std::endl;

		pop_mixed.clear();
		// the initial popluation contains the initial random population
		for(population::size_type i=0; i<pop_size; i++) {
			pop_mixed.push_back(pop.get_individual(i).cur_x);
		}

		pop_antigens.clear();

		// clearing the pools
		pop_antigens_pool.clear();
		pop_antibodies_pool.clear();

		// first of all we compute the constraints
		for(population::size_type i=0; i<pop_size; i++) {
			const population::individual_type &current_individual = pop_mixed.get_individual(i);
			pop_mixed_c[i] = prob.compute_constraints(current_individual.cur_x);
		}

		// we find if there are feasible individuals
		bool has_feasible = false;

		for(population::size_type i=0; i<pop_size; i++) {

			if(prob.feasibility_c(pop_mixed_c[i])) {
				has_feasible = true;
				break;
			}
		}

		// if feasible solutions exist in the population, the
		// antigens are selected by being close to the average fitness
		// of a sub population

		if(has_feasible) {
			// the pop_antigens_pool is based on the feasible population
			for(population::size_type i=0; i<pop_size; i++) {
				const population::individual_type &current_individual = pop.get_individual(i);

				if(prob.feasibility_c(current_individual.cur_c)) {
					pop_antigens_pool.push_back(i);
				} else {
					pop_antibodies_pool.push_back(i);
				}
			}

			population::size_type pop_antigens_pool_size = pop_antigens_pool.size();

			// we sort the pop_antigens_pool according to the fitness
			for(population::size_type i=0; i<pop_antigens_pool_size-1; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				const population::individual_type &current_individual = pop.get_individual(current_idx);

				for(population::size_type j=i+1; j<pop_antigens_pool_size; j++) {
					const population::size_type &current_second_idx = pop_antigens_pool.at(j);
					const population::individual_type &current_second_individual = pop.get_individual(current_second_idx);

					if(prob.compare_fitness(current_second_individual.cur_f, current_individual.cur_f)) {
						std::swap(pop_antigens_pool[i],pop_antigens_pool[j]);
					}
				}
			}

			// a subset of the best antigens in the pool is selected to compute the fitness average
			int pop_antigens_subset_size = std::max((int)(phi * pop_antigens_pool_size),1);

			// we compute the mean fitness value from the subset of the best individuals in the population
			double mean_fitness = 0.;
			for(population::size_type i=0; i<pop_antigens_subset_size; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				mean_fitness += pop.get_individual(current_idx).cur_f.at(0);
			}
			mean_fitness /= pop_antigens_subset_size;

			// the population around this mean fitness is selected to fill the antigens population
			// finds the position of the individual closest to the mean fitness

			// from the antigens pool we select a fraction from ones that are closest to the mean value
			// according to Hjale and Lee. The idea here is to get as much information about the
			// feasible domain as possible.

			int mean_position_idx=0;
			for(population::size_type i=0; i<pop_antigens_pool_size; i++) {
				const population::size_type &current_idx = pop_antigens_pool.at(i);
				if(mean_fitness <= pop.get_individual(current_idx).cur_f.at(0)) {
					mean_position_idx = i;
					break;
				}
			}

			// generates the antigens population
			int pop_antigens_size = std::max((int)(gamma * pop_antigens_pool_size), 1);

			int begin_antigen_idx = mean_position_idx - pop_antigens_size/2;
			int end_antigen_idx = mean_position_idx + pop_antigens_size/2;

			// move the selection range depending on the boundaries
			if(mean_position_idx - pop_antigens_size/2 < 0) {
				begin_antigen_idx = 0;
				end_antigen_idx = std::min((int)(end_antigen_idx + (pop_antigens_size/2 - mean_position_idx)),(int)pop_antigens_pool_size);
			}
			if(mean_position_idx + pop_antigens_size/2 >= pop_antigens_size) {
				begin_antigen_idx = std::max((int)(begin_antigen_idx - (mean_position_idx + pop_antigens_size/2 - pop_antigens_size)),0);
				end_antigen_idx = pop_antigens_size;
			}

			// we select individuals around the mean
			for(population::size_type i=begin_antigen_idx; i<end_antigen_idx; i++) {
				population::size_type current_individual_idx = pop_antigens_pool.at(i);
				pop_antigens.push_back(pop.get_individual(current_individual_idx).cur_x);
			}

			if(pop_antigens.size() == 0) {
				pop_antigens.push_back(pop.get_individual(pop_antigens_pool.at(mean_position_idx)).cur_x);
			}

		} else {
			// no feasible founds, the antigen population
			// is selected accordingly to selected method

			switch(m_method) {
			case(BEST_ANTIBODY): {
				// Coello method the antigen population contains only one antigen which is the
				// individual with lowest constraints violation
				population::size_type best_idx = 0;
				for(population::size_type i=1; i<pop_size; i++) {

					const population::individual_type &best_individual = pop.get_individual(best_idx);
					const population::individual_type &current_individual = pop.get_individual(i);

					if(prob.compare_constraints(current_individual.cur_c, best_individual.cur_c)) {
						best_idx = i;
					}
				}
				pop_antigens_pool.push_back(best_idx);
				pop_antigens.push_back(pop.get_individual(best_idx).cur_x);

				// antibodies
				for(population::size_type i=0; i<pop_size; i++) {
					if(i != best_idx){
						pop_antibodies_pool.push_back(i);
					}
				}
				break;
			}
			case(INFEASIBILITY): {
				// Personal method where the infeasibility is used
				// to select the antigen population
				for(population::size_type i=0; i<pop_size; i++) {
					pop_antigens_pool.push_back(i);
				}

				population::size_type pop_antigens_pool_size = pop_antigens_pool.size();

				// we sort the pop_antigens_pool according to the fitness
				for(population::size_type i=0; i<pop_antigens_pool_size-1; i++) {
					const population::size_type &current_idx = pop_antigens_pool.at(i);
					const population::individual_type &current_individual = pop.get_individual(current_idx);

					for(population::size_type j=i+1; j<pop_antigens_pool_size; j++) {
						const population::size_type &current_second_idx = pop_antigens_pool.at(j);
						const population::individual_type &current_second_individual = pop.get_individual(current_second_idx);

						if(prob.compare_constraints(current_second_individual.cur_c, current_individual.cur_c)) {
							std::swap(pop_antigens_pool[i],pop_antigens_pool[j]);
						}
					}
				}

				// fill the antigens with the best ones
				int pop_antigens_size = std::max((int)(gamma * pop_antigens_pool_size), 1);

				for(population::size_type i=0; i<pop_antigens_size; i++) {
					population::size_type current_individual_idx = pop_antigens_pool.at(i);
					pop_antigens.push_back(pop.get_individual(current_individual_idx).cur_x);
				}

				// antibodies
				// not the best implementation, trying to find a better one...
				for(population::size_type i=0; i<pop_size; i++) {
					if(std::find(pop_antigens_pool.begin(), pop_antigens_pool.begin()+pop_antigens_size, i) == pop_antigens_pool.begin()+pop_antigens_size) {
						pop_antibodies_pool.push_back(i);
					}
				}
				break;
			}
			}
		}

		population::size_type initial_pop_antibodies_pool_size = pop_antibodies_pool.size();

		if(initial_pop_antibodies_pool_size != 0) {

			// initial pool size
			// random shuffle of the antibodies as the antibodies population
			// are randomly selected without replacement from the antibodies pool
			if(initial_pop_antibodies_pool_size > 1) {
				for (population::size_type i=0; i<initial_pop_antibodies_pool_size; i++)
				{
					int j = boost::uniform_int<int>(0, initial_pop_antibodies_pool_size - 1)(m_urng);
					std::swap(pop_antibodies_pool[i], pop_antibodies_pool[j]);
				}
			}

			// select the antibodies population with the requested size:
			// antibodies population size = 1/3 antigens population size

			population::size_type pop_antigens_size = pop_antigens.size();

			int min_individual_for_algo = 8;

			population::size_type pop_antibodies_size = std::max( (int)(sigma * pop_antigens_size), min_individual_for_algo);
			pop_antibodies_size = std::min( pop_antibodies_size, initial_pop_antibodies_pool_size);

			//population::size_type pop_antibodies_size = std::max((int)(0.5 * pop_antigens_size), 6);
			//population::size_type pop_antibodies_size = std::max((int)(pop_antibodies_pool.size()), 6);

			// the problem can be updated with antigenes, need to be done here to avoid a cast
			problem::antibodies_problem prob_antibodies(prob, problem::antibodies_problem::EUCLIDEAN);
			prob_antibodies.set_antigens(pop_antigens);

			// immune system initialization
			population pop_antibodies(prob_antibodies);
			pop_antibodies.clear();

			for(population::size_type i=0; i<pop_antibodies_size; i++) {
				pop_antibodies.push_back(pop.get_individual(pop_antibodies_pool.at(i)).cur_x);
			}

			// ensure that the antibodies population has at least 6 individuals for de, sga, 8 for jde...
			int extra_antibodies_size = min_individual_for_algo - pop_antibodies_size;

			// this test is needed as for some reason i<extra_antibodies_size can provoque an
			// infinite loop if extra_antibodies_size < 0! Due to population::size_type?
			if(extra_antibodies_size > 0) {
				// add extra needed by randomly selecting the individuals in the pool
				for(population::size_type i=0; i<extra_antibodies_size; i++) {
					if(initial_pop_antibodies_pool_size > 1) {
						int j = boost::uniform_int<int>(0, initial_pop_antibodies_pool_size - 1)(m_urng);

						pop_antibodies.push_back(pop.get_individual(pop_antibodies_pool.at(j)).cur_x);
					} else {
						pop_antibodies.push_back(pop.get_individual(pop_antibodies_pool.at(0)).cur_x);
					}
				}
			}

			pop_antigens_size = pop_antigens.size();
			pop_antibodies_size = pop_antibodies.size();

			// run the immune system
			m_original_algo_immune->evolve(pop_antibodies);

			// if using this version instead, this is not working (need to find why):

			// sets the mixed population with all current best designs
			// and the copies of constraint conditioned designs
			for(population::size_type i=0; i<initial_pop_antibodies_pool_size; i++) {
				const population::size_type &current_idx = pop_antibodies_pool.at(i);
				// multiple copies of the best antibody
				// should be an option (with copies of the 25% best)
				pop_mixed.set_x(current_idx, pop_antibodies.champion().x);
			}
		}

		// for individuals, evolve the mixed population
		// which is an unconstrained problem
		// only one iteration should be done...
		m_original_algo->evolve(pop_mixed);

		// we update the population pop with the new population

		// store the final population in the main population, just for checking the right convergence
		pop.clear();
		for(population::size_type i=0; i<pop_size; i++) {
			pop.push_back(pop_mixed.get_individual(i).cur_x);
		}

		std::cout << pop.champion() << std::endl;

		// Check the exit conditions (every 40 generations, just as DE)
		if(k % 40 == 0) {
			decision_vector tmp(prob_dimension);

			double dx = 0;
			for(decision_vector::size_type i=0; i<prob_dimension; i++) {
				tmp[i] = pop_mixed.get_individual(pop_mixed.get_worst_idx()).best_x[i] - pop_mixed.get_individual(pop_mixed.get_best_idx()).best_x[i];
				dx += std::fabs(tmp[i]);
			}

			if(dx < m_xtol ) {
				if (m_screen_output) {
					std::cout << "Exit condition -- xtol < " << m_xtol << std::endl;
				}
				return;
			}

			double mah = std::fabs(pop_mixed.get_individual(pop_mixed.get_worst_idx()).best_f[0] - pop_mixed.get_individual(pop_mixed.get_best_idx()).best_f[0]);

			if(mah < m_ftol) {
				if(m_screen_output) {
					std::cout << "Exit condition -- ftol < " << m_ftol << std::endl;
				}
				return;
			}

			// outputs current values
			if(m_screen_output) {
				std::cout << "Generation " << k << " ***" << std::endl;
				std::cout << "    Best global fitness: " << pop.champion().f << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
				std::cout << "    xtol: " << dx << ", ftol: " << mah << std::endl;
			}
		}
	}
}

/// Algorithm name
std::string immune_system::get_name() const
{
	return m_original_algo->get_name() + "[Immune]";
}

/// Get a copy of the internal local algorithm.
/**
 * @return algorithm::base_ptr to a copy of the internal local algorithm.
 */
base_ptr immune_system::get_algorithm() const
{
	return m_original_algo->clone();
}

/// Set algorithm.
/**
 * A copy of the input algorithm will be set as the internal local algorithm.
 *
 * @param[in] algo algorithm to be set as local algorithm.
 */
void immune_system::set_algorithm(const base &algo)
{
	m_original_algo = algo.clone();
}

/// Extra human readable algorithm info.
/**
 * @return a formatted string displaying the parameters of the algorithm.
 */
std::string immune_system::human_readable_extra() const
{
	std::ostringstream s;
	s << "algorithms: " << m_original_algo->get_name() << " - " << m_original_algo_immune->get_name() << " ";
	s << "\n\tConstraints handled with immune system algorithm";
	return s.str();
}

/// Computes the solution infeasibility measure for the given constraint,
/// need the constraints scaling to be updated before calling this method.
/**
 * Updates the solution infeasibility vector with the population given.
 * @param[in] constraint_vector c.
 * @param[out] solution infeasibility.
 */
double immune_system::compute_solution_infeasibility(const constraint_vector &c, const constraint_vector &c_scaling, const problem::base& original_prob) const
{
	// get the constraints dimension
	problem::base::c_size_type prob_c_dimension = original_prob.get_c_dimension();
	problem::base::c_size_type number_of_eq_constraints =
			original_prob.get_c_dimension() -
			original_prob.get_ic_dimension();

	double solution_infeasibility = 0.;

	const std::vector<double> &c_tol = original_prob.get_c_tol();

	// computes solution infeasibility with the right definition of the constraints (can be in base problem? currently used
	// by con2mo as well)
	for(problem::base::c_size_type j=0; j<number_of_eq_constraints; j++) {
		// test needed otherwise the c_scaling can be 0, and division by 0 occurs
		if(c_scaling[j] > 0.) {
			solution_infeasibility += std::max(0.,(std::abs(c.at(j)) - c_tol.at(j))) / c_scaling[j];
		}
	}
	for(problem::base::c_size_type j=number_of_eq_constraints; j<prob_c_dimension; j++) {
		if(c_scaling[j] > 0.) {
			solution_infeasibility += std::max(0.,c.at(j)) / c_scaling[j];
		}
	}

	solution_infeasibility /= prob_c_dimension;

	return solution_infeasibility;
}


}} //namespaces

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::algorithm::immune_system);
