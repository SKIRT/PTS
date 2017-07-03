#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.dustpedia.core.sdss import sdss_bands
from pts.core.tools.parallelization import ncores

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")
definition.add_optional("bands", "string_list", "the bands (u/g/r/i/z)", choices=sdss_bands, default=sdss_bands)

# Output
definition.add_optional("output", "string", "the name of the output directory", default="out")

# Advanced options
definition.add_optional("max_nobservations_u", "positive_integer", "maximum number of u band observations")
definition.add_optional("max_nobservations_g", "positive_integer", "maximum number of g band observations")
definition.add_optional("max_nobservations_r", "positive_integer", "maximum number of r band observations")
definition.add_optional("max_nobservations_i", "positive_integer", "maximum number of i band observations")
definition.add_optional("max_nobservations_z", "positive_integer", "maximum number of z band observations")

# Parallelization
definition.add_optional("nprocesses", "positive_integer", "number of processes for parallel execution", max(8, ncores()))

# Advanced
definition.add_optional("fields_directories", "string_string_dictionary", "dictionary of fields directories where the keys are the bands")

# -----------------------------------------------------------------
