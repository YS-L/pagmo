/*****************************************************************************
 *   Copyright (C) 2008, 2009 Advanced Concepts Team (European Space Agency) *
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

// 06/02/2009: Initial version by Francesco Biscani.

#ifndef PAGMO_ONE_WAY_RING_TOPOLOGY_H
#define PAGMO_ONE_WAY_RING_TOPOLOGY_H

#include "../../config.h"
#include "base_topology.h"
#include "graph_topology.h"
#include "island.h"

class __PAGMO_VISIBLE one_way_ring_topology: public base_topology, public graph_topology {
	public:
		one_way_ring_topology(const double &);
		one_way_ring_topology(const one_way_ring_topology &);
		virtual one_way_ring_topology *clone() const {return new one_way_ring_topology(*this);}
		virtual void push_back(const island &);
		virtual void reset();
		virtual void pre_evolution(island &);
		virtual void post_evolution(island &);
	private:
		one_way_ring_topology &operator=(const one_way_ring_topology &);
		size_t	m_first;
		size_t	m_last;
};

#endif