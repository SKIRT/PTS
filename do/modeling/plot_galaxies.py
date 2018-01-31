#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_galaxies Plot the positions of the galaxies in the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.services.s4g import get_galaxy_names, has_galaxy
from pts.core.basics.map import Map
from pts.core.basics.log import log
from pts.core.basics.table import SmartTable

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

#
definition.add_optional("ngalaxies", "positive_integer", "max number of galaxies")

# Get configuration
config = parse_arguments("plot_galaxies", definition)

# -----------------------------------------------------------------

# Create the database instance
database = DustPediaDatabase()
#username, password = get_account()
#database.login(username, password)

# -----------------------------------------------------------------

sample = DustPediaSample()
galaxy_names = sample.get_names()

# -----------------------------------------------------------------



# -----------------------------------------------------------------

# Create table
table = SmartTable.from_dictionaries(*galaxies)
print(table)

# -----------------------------------------------------------------

from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random

fig = pyplot.figure()
ax = Axes3D(fig)

sequence_containing_x_vals = list(range(0, 100))
sequence_containing_y_vals = list(range(0, 100))
sequence_containing_z_vals = list(range(0, 100))

random.shuffle(sequence_containing_x_vals)
random.shuffle(sequence_containing_y_vals)
random.shuffle(sequence_containing_z_vals)

ax.scatter(sequence_containing_x_vals, sequence_containing_y_vals, sequence_containing_z_vals)
pyplot.show()

# -----------------------------------------------------------------
