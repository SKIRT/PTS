#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.estimate Estimate the resource requirements for a certain ski file

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.advanced.resources import ResourceEstimator
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter, initialize_log

# -----------------------------------------------------------------

# Create the command-line parser
definition = ConfigurationDefinition()
definition.add_required("skifile", "file_path", "name of the ski file")
definition.add_required("nthreads", "positive_integer", "number of parallel threads")
definition.add_required("nprocesses", "positive_integer", "number of parallel processes")

# Parse
setter = ArgumentConfigurationSetter("estimate")
config = setter.run(definition)

# -----------------------------------------------------------------

# Initialize the logger
log = initialize_log(config)
log.info("Starting estimate ...")

# -----------------------------------------------------------------

# Determine the full path to the parameter file
ski_path = config.skifile

# Create and run a ResourceEstimator oject
estimator = ResourceEstimator()
estimator.run(ski_path, config.nprocesses, config.nthreads)

# Inform the user about the resource requirements
log.info("This simulation requires " + estimator.memory + " GB of virtual memory")
#log.info("This simulation requires a walltime of " + estimator.walltime)

# -----------------------------------------------------------------
