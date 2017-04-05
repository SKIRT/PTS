#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids
from pts.modeling.tests.base import possible_free_parameters, default_free_parameters
from pts.modeling.fitting.configuration import genetic_definition

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

definition.add_flag("plot", "do plotting", False)

definition.add_flag("rotate", "rotation", True)

definition.add_optional("nrandom_sources", "positive_integer", "number of random sources", 100)

definition.add_flag("add_catalogued_sources", "add catalogued sources", False)

definition.add_optional("point_source_catalogs", "string_list", "point source catalogs", ["II/246"])

definition.add_optional("nfilters_stars", "positive_integer", "number of filters to use where stars are visible", 5)
definition.add_optional("nfilters_extra", "positive_integer", "number of filters to use where stars are not visible", 2)

definition.add_optional("psf_model", "string", "model to use for the PSF", "airydisk", choices=["gaussian", "airydisk"])

definition.add_optional("noise_stddev", "real", "stddev of noise", 2.)

# -----------------------------------------------------------------
