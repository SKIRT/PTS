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
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.core.environment import GalaxyModelingEnvironment
from pts.modeling.preparation.preparer import load_statistics, get_statistics_path
from pts.core.filter.filter import parse_filter
from pts.magic.services.attenuation import GalacticAttenuation
from pts.core.tools.logging import log
from pts.modeling.component.galaxy import get_galaxy_properties
from pts.core.tools.stringify import tostr
from pts.core.tools import numbers
from pts.modeling.core.environment import verify_modeling_cwd
from pts.modeling.preparation.preparer import get_step_image_paths_with_cached
from pts.core.launch.pts import execute_pts_remote, execute_pts_local
from pts.core.tools import filesystem as fs
from pts.core.tools.serialization import write_dict

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()
config = parse_arguments("fix_attenuation", definition)

# -----------------------------------------------------------------

# Modeling path
modeling_path = verify_modeling_cwd()

# Load the modeling environment
environment = GalaxyModelingEnvironment(modeling_path)

# Load the caching remote
remote = environment.cache_host_id

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

# Remember the correction factors
correction_factors = dict()

# -----------------------------------------------------------------

# Loop over the prep names that need to be fixed
for prep_name in fix:

    # Get the image paths after and including the extinction correction step
    paths = get_step_image_paths_with_cached(modeling_path, prep_name, environment.cache_host_id, after_step="extinction", inclusive=True)

    # Determine the correction factor
    actual = float(fix[prep_name][0])
    mistaken = float(fix[prep_name][1])
    factor = 10**(actual - mistaken)

    # Debugging
    log.debug("The correction factor for the " + prep_name + " image is " + str(factor))
    correction_factors[prep_name] = factor

    #print(paths)
    #continue

    # Loop over the images
    for step in paths:

        # Inform the user
        log.info("Correcting the " + prep_name + " image after the " + step + " step ...")

        # Get the path
        path = paths[step]

        # Correct all the images by multiplying with the factor
        # Local
        if fs.is_file(path):
            log.debug("Fixing locally ...")
            execute_pts_local("multiply", path, factor, backup=True, debug=True)

        # Remote
        else:
            log.debug("Fixing remotely ...")
            execute_pts_remote(remote, "multiply", path, factor, backup=True, debug=True)

    # CHANGE ATTENUATION VALUE IN PREPARATION STATISTICS

    # Load the preparation statistics
    statistics_path = get_statistics_path(modeling_path, prep_name)
    statistics = load_statistics(modeling_path, prep_name)

    # BACKUP THE STATISTICS
    backup_statistics_path = fs.appended_filepath(statistics_path, "_backup")
    statistics.saveto(backup_statistics_path)

    # Set new attenuation value
    statistics.attenuation = actual

    # Save
    statistics.saveto(statistics_path)

# -----------------------------------------------------------------

# Write the correction factors
factors_path = fs.join(environment.prep_path, "attenuation_correction_factors.dat")
write_dict(correction_factors, factors_path)

# -----------------------------------------------------------------
