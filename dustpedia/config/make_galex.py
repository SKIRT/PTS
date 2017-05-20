#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.dustpedia.core.galex import galex_bands
from pts.core.tools.parallelization import ncores

# -----------------------------------------------------------------

# Configuration
definition = ConfigurationDefinition()

# Galaxy name
definition.add_required("galaxy_name", "string", "the name of the galaxy")
definition.add_optional("bands", "string_list", "the bands (NUV/FUV)", choices=galex_bands, default=galex_bands)

# Output
definition.add_optional("output", "string", "the name of the output directory", default="out")

# Advanced
definition.add_optional("max_nobservations_fuv", "positive_integer", "limit the number of FUV observations")
definition.add_optional("max_nobservations_nuv", "positive_integer", "limit the number of NUV observations")

# Parallelization
definition.add_optional("nprocesses", "positive_integer", "number of processes for parallel execution", max(8, ncores()))

# Advanced
definition.add_optional("download_directories", "string_string_dictionary", "dictionary of download directories where the keys are the bands")

definition.add_optional("manual_selection", "string_string_list_dictionary", "dictionary of selected observations to use (no filtering will be done then)")

# -----------------------------------------------------------------
