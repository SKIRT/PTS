#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.check_attenuation Check the attenuation corrected of the prepared images.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.preparation.preparer import load_statistics
from pts.core.filter.filter import parse_filter
from pts.magic.services.extinction import GalacticExtinction
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Create the configuration
config = parse_arguments("check_attenuation", definition)

# -----------------------------------------------------------------

# Modeling path
modeling_path = fs.cwd()

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

properties = get_galaxy_properties(modeling_path)

# -----------------------------------------------------------------

#extinction = GalacticExtinction(environment.galaxy_name)
extinction = GalacticExtinction(properties.center)

# -----------------------------------------------------------------

print("")

# Loop over the names
for prep_name in environment.preparation_names:

    # Load the statistics
    statistics = load_statistics(modeling_path, prep_name)

    # Determine filter
    fltr = parse_filter(prep_name)

    # Get the extinction
    ext = extinction.extinction_for_filter(fltr)

    if ext == 0.0:
        if statistics.attenuation != 0.0: log.warning(prep_name + ": extinction is zero but preparation extinction value was " + str(statistics.attenuation))
        continue

    # Ratio
    ratio = statistics.attenuation / ext
    rel = abs((statistics.attenuation - ext)/ext)

    #print(prep_name, statistics.attenuation, ext, ratio, rel * 100)

    print(prep_name)
    print("")
    print(" - preparation: " + str(statistics.attenuation))
    print(" - real: " + str(ext))
    print(" - ratio: " + str(ratio))
    print(" - rel difference: " + str(rel * 100) + "%")
    print("")

# -----------------------------------------------------------------
