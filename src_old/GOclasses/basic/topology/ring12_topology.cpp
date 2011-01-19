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

// 13/01/2009: Initial version by Marek Ruciński.

#include "ring12_topology.h"

ring12_topology::ring12_topology():graph_topology(), a(0), b(0), c(0), d(0) {}

ring12_topology::ring12_topology(const ring12_topology &r)
		:graph_topology(r), a(r.a), b(r.b), c(r.c), d(r.d) {}

ring12_topology &ring12_topology::operator=(const ring12_topology &)
{
	pagmo_assert(false);
	return *this;
}

void ring12_topology::push_back(const size_t& id)
{
	// Store frequently-used variables.
	const size_t t_size = get_number_of_nodes();

	//Add the new node to the graph
	add_node(id);

	switch (t_size) {
	case 0:
		// If topology is empty, update the id of the first element.
		a = id;
		break;

	case 1:
		// Add connections to the only existing element.
		add_edge(id, a);
		add_edge(a, id);
		b = id;
		break;

	case 2:
		// Add new connections
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		c = id;
		break;

	case 3:
		// Add new connections
		add_edge(a, id);
		add_edge(id, a);
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);
		d = id;
		break;

	case 4:
		// Add new connections and update c,d

		// First, put the node in between a and d. Former +1 link a-d becomes now a +2 link
		add_edge(d, id);
		add_edge(id, d);
		add_edge(a, id);
		add_edge(id, a);

		// Also add connections id-b and id-c. Nothing needs to be dropped.
		add_edge(b, id);
		add_edge(id, b);
		add_edge(c, id);
		add_edge(id, c);

		// Update nodes
		c = d;
		d = id;

		// Now we have a pentagram :)
		break;

	default:
		// Insert the new node between d and a.

		// First, put the node in between a and d. Former +1 link a-d becomes now a +2 link
		add_edge(d, id);
		add_edge(id, d);
		add_edge(a, id);
		add_edge(id, a);

		// Now drop connections b-d and c-a, as they are no longer +2 but +3, replacing them by edges b-id and c-id
		remove_edge(a, c);
		remove_edge(c, a);
		add_edge(id, c);
		add_edge(c, id);

		remove_edge(b, d);
		remove_edge(d, b);
		add_edge(id, b);
		add_edge(b, id);

		// Finally, update the nodes
		c = d;
		d = id;
		break;
	}
}
