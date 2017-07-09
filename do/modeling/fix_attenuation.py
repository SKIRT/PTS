#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.fix_attenuation Fix the attenuation correction of the prepared images.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.preparation.preparer import load_statistics
from pts.core.filter.filter import parse_filter
from pts.magic.services.attenuation import GalacticAttenuation
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties
from pts.core.tools.stringify import tostr
from pts.core.tools import numbers

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Create the configuration
config = parse_arguments("fix_attenuation", definition)

# -----------------------------------------------------------------

# Modeling path
modeling_path = fs.cwd()

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# -----------------------------------------------------------------

properties = get_galaxy_properties(modeling_path)

# -----------------------------------------------------------------

attenuation = GalacticAttenuation(properties.center)

# -----------------------------------------------------------------

fix = dict()

# -----------------------------------------------------------------

# Loop over the names
for prep_name in environment.preparation_names:

    # Info
    log.info("Checking the '" + prep_name + "' image ...")

    # Load the preparation statistics
    statistics = load_statistics(modeling_path, prep_name)

    # Determine filter
    fltr = parse_filter(prep_name)

    # Get the extinction
    att = attenuation.extinction_for_filter(fltr)

    # Attenuation zero but still corrected for attenuation that is nonzero
    if att == 0.0 and statistics.attenuation != 0.0:

        # Give warning and add
        log.warning(prep_name + ": attenuation is zero but preparation attenuation value was " + str(statistics.attenuation))
        fix[prep_name] = (att, statistics.attenuation)

    # Attenuation nonzero but not corrected
    if att != 0.0 and statistics.attenuation == 0.0:

        # Give warning and add
        log.warning(prep_name + ": attenuation is nonzero but not corrected")
        fix[prep_name] = (att, statistics.attenuation)

    # Attenuations not equal
    elif not numbers.is_close(att, statistics.attenuation):

        # Give warning and add
        log.warning(prep_name + ": attenuation of " + tostr(att) + " does not correspond to value of " + tostr(statistics.attenuation) + " used for correction")
        fix[prep_name] = (att, statistics.attenuation)

# -----------------------------------------------------------------

print(fix)

# -----------------------------------------------------------------