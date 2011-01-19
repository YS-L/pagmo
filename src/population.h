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

#ifndef PAGMO_POPULATION_H
#define PAGMO_POPULATION_H

#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "config.h"
#include "problem/base.h"
#include "rng.h"
#include "serialization.h"
#include "types.h"

namespace pagmo
{

// Forward declarations.
class base_island;
struct population_access;

/// Population class.
/**
 * This class contains an instance of an optimisation problem and a group of candidate solutions represented by the class individual_type. On creation,
 * the population is associated to a problem and initialised with random decision vectors. An instance of champion_type automatically keeps track of the best solution
 * ever appeared in the population. Methods are offered to get and manipulate the single individuals.
 *
 * Additionally, the population class holds for each individual I a "domination list", constituted by the list of individuals (identified by their positional index
 * in the population) which are dominated by I. Individual I1 is dominated by individual I2 if problem::base::compare_fc on the fitness and constraints vectors
 * of I1 and I2 respectively returns true. The best individual in the population is the one dominating the highest number of individuals.
 *
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
class __PAGMO_VISIBLE population
{
		friend class base_island;
		friend struct population_access;
	public:
		/// Individuals stored in the population.
		/**
		 * Individuals store the current decision and velocity vectors, the current constraint vector and the current fitness vector. They also
		 * keep memory of the best decision, constraint and fitness vectors "experienced" so far by the individual.
		 */
		struct individual_type
		{
				/// Current decision vector.
				decision_vector		cur_x;
				/// Current velocity vector.
				decision_vector		cur_v;
				/// Current constraint vector.
				constraint_vector	cur_c;
				/// Current fitness vector.
				fitness_vector		cur_f;
				/// Best decision vector so far.
				decision_vector		best_x;
				/// Best constraint vector so far.
				constraint_vector	best_c;
				/// Best fitness vector so far.
				fitness_vector		best_f;
				/// Human-readable representation.
				/**
				* @return formatted string containing the values of the data members.
				*/
				std::string human_readable() const
				{
					std::ostringstream oss;
					oss << "\tDecision vector:\t\t" << cur_x << '\n';
					oss << "\tVelocity vector:\t\t" << cur_v << '\n';
					oss << "\tConstraint vector:\t\t" << cur_c << '\n';
					oss << "\tFitness vector:\t\t\t" << cur_f << '\n';
					oss << "\tBest decision vector:\t\t" << best_x << '\n';
					oss << "\tBest constraint vector:\t\t" << best_c << '\n';
					oss << "\tBest fitness vector:\t\t" << best_f << '\n';
					return oss.str();
				}
			private:
				friend class boost::serialization::access;
				template <class Archive>
				void save(Archive &ar, const unsigned int version) const
				{
					custom_vector_double_save(ar,cur_x,version);
					custom_vector_double_save(ar,cur_v,version);
					custom_vector_double_save(ar,cur_c,version);
					custom_vector_double_save(ar,cur_f,version);
					custom_vector_double_save(ar,best_x,version);
					custom_vector_double_save(ar,best_c,version);
					custom_vector_double_save(ar,best_f,version);
				}
				template <class Archive>
				void load(Archive &ar, const unsigned int version)
				{
					custom_vector_double_load(ar,cur_x,version);
					custom_vector_double_load(ar,cur_v,version);
					custom_vector_double_load(ar,cur_c,version);
					custom_vector_double_load(ar,cur_f,version);
					custom_vector_double_load(ar,best_x,version);
					custom_vector_double_load(ar,best_c,version);
					custom_vector_double_load(ar,best_f,version);
				}
				template <class Archive>
				void serialize(Archive &ar, const unsigned int version)
				{
					boost::serialization::split_member(ar,*this,version);
				}
		};
		/// Population champion.
		/**
		 * A champion is the best individual that ever lived in the population. It is defined by a decision vector, a constraint vector and a fitness vector.
		 */
		struct champion_type
		{
				/// Decision vector.
				decision_vector		x;
				/// Constraint vector.
				constraint_vector	c;
				/// Fitness vector.
				fitness_vector		f;
				/// Human-readable representation.
				/**
				* @return formatted string containing the values of the data members.
				*/
				std::string human_readable() const
				{
					std::ostringstream oss;
					oss << "\tDecision vector:\t" << x << '\n';
					oss << "\tConstraints vector:\t" << c << '\n';
					oss << "\tFitness vector:\t\t" << f << '\n';
					return oss.str();
				}
			private:
				friend class boost::serialization::access;
				template <class Archive>
				void save(Archive &ar, const unsigned int version) const
				{
					custom_vector_double_save(ar,x,version);
					custom_vector_double_save(ar,c,version);
					custom_vector_double_save(ar,f,version);
				}
				template <class Archive>
				void load(Archive &ar, const unsigned int version)
				{
					custom_vector_double_load(ar,x,version);
					custom_vector_double_load(ar,c,version);
					custom_vector_double_load(ar,f,version);
				}
				template <class Archive>
				void serialize(Archive &ar, const unsigned int version)
				{
					boost::serialization::split_member(ar,*this,version);
				}
		};
		/// Underlying container type.
		typedef std::vector<individual_type> container_type;
		/// Population size type.
		typedef container_type::size_type size_type;
		/// Const iterator.
		typedef container_type::const_iterator const_iterator;
		explicit population(const problem::base &, int = 0);
		population(const population &);
		population &operator=(const population &);
		const individual_type &get_individual(const size_type &) const;
		const std::vector<size_type> &get_domination_list(const size_type &) const;
		const problem::base &problem() const;
		const champion_type &champion() const;
		std::string human_readable_terse() const;
		std::string human_readable() const;
		size_type get_best_idx() const;
		size_type get_worst_idx() const;
		void set_x(const size_type &, const decision_vector &);
		void set_v(const size_type &, const decision_vector &);
		void push_back(const decision_vector &);
		size_type size() const;
		const_iterator begin() const;
		const_iterator end() const;
		size_type n_dominated(const individual_type &) const;
		void reinit(const size_type &);
		void reinit();
	private:
		void init_velocity(const size_type &);
		void update_champion(const size_type &);
		void update_dom_list(const size_type &);
		struct domination_comp {
			domination_comp(const population &pop):m_pop(pop) {}
			bool operator()(const individual_type &i1, const individual_type &i2) const
			{
				pagmo_assert(&i1 >= &m_pop.m_container.front() && &i1 <= &m_pop.m_container.back());
				pagmo_assert(&i2 >= &m_pop.m_container.front() && &i2 <= &m_pop.m_container.back());
				const size_type idx1 = &i1 - &m_pop.m_container.front(), idx2 = &i2 - &m_pop.m_container.front();
				return m_pop.m_dom_list[idx1].size() > m_pop.m_dom_list[idx2].size();
			}
			const population &m_pop;
		};
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & m_prob;
			ar & m_container;
			ar & m_dom_list;
			ar & m_champion;
			ar & m_drng;
			ar & m_urng;
		}
		// Data members.
		// Problem.
		problem::base_ptr			m_prob;
		// Container of individuals.
		container_type				m_container;
		// List of dominated individuals.
		std::vector<std::vector<size_type> >	m_dom_list;
		// Population champion.
		champion_type				m_champion;
		// Double precision random number generator.
		mutable	rng_double			m_drng;
		// uint32 random number generator.
		mutable	rng_uint32			m_urng;
};

__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const population &);
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const population::individual_type &);
__PAGMO_VISIBLE_FUNC std::ostream &operator<<(std::ostream &, const population::champion_type &);

struct population_access
{
	static problem::base_ptr &get_problem_ptr(population &);
};

}

namespace boost { namespace serialization {

template <class Archive>
inline void save_construct_data(Archive &ar, const pagmo::population *pop, const unsigned int)
{
	// Save data required to construct instance.
	pagmo::problem::base_ptr prob = pop->problem().clone();
	ar << prob;
}

template <class Archive>
inline void load_construct_data(Archive &ar, pagmo::population *pop, const unsigned int)
{
	// Retrieve data from archive required to construct new instance.
	pagmo::problem::base_ptr prob;
	ar >> prob;
	// Invoke inplace constructor to initialize instance of the population.
	::new(pop)pagmo::population(*prob);
}

}} //namespaces

#endif
