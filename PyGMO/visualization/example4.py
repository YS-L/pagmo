#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# Copyright (C) 2004-2009 The PaGMO development team,
# Advanced Concepts Team (ACT), European Space Agency (ESA)
# http://apps.sourceforge.net/mediawiki/pagmo
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers
# http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits
# act@esa.int
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


from PyGMO import visualization


"""
Example of an Earth gravity assist with deep space maneuvers.
"""


# Run some tests
if __name__ == "__main__":

   # Create the engine
   traj = visualization.Trajectory3D( "EarthEarthJupiter.csv",
         640, 480,
         24.*3600., 1000., 1000., 1000. ) # Unit conversions: days->s, km->m
   traj.addPlanets( "EarthEarthJupiter_flybyinfo.txt" )
   traj.setUnits( "km", "d", "km/s" )

   # Create some stuff
   traj.repeat( True )
   traj.axes( True )

   # Start the engine
   traj.start()






