#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

definition.add_flag("plot", "do plotting", False)

definition.add_flag("rotate", "rotation", True)

definition.add_optional("nrandom_sources", "positive_integer", "number of random sources", 100)

definition.add_flag("vary_fwhm", "vary the FWHM slightly for each star", True)

definition.add_flag("add_catalogued_sources", "add catalogued sources", False)

definition.add_optional("point_source_catalogs", "string_list", "point source catalogs", ["II/246"])

definition.add_optional("nfilters_stars", "positive_integer", "number of filters to use where stars are visible", 4)
definition.add_optional("nfilters_extra", "positive_integer", "number of filters to use where stars are not visible", 2)

definition.add_optional("psf_model", "string", "model to use for the PSF", "airydisk", choices=["gaussian", "airydisk"])

definition.add_optional("noise_stddev", "real", "stddev of noise", 5.)

definition.add_flag("only_local", "only evaluate PSFs locally", True)

definition.add_optional("galaxy_sersic_index", "real", "sersic index", 3.)
definition.add_optional("galaxy_relative_asymptotic_radius", "positive_real", "radius of the galaxy mask relative to the effective radius", 4.)
definition.add_flag("limit_galaxy", "limit the rendering of the galaxy to a certain ellipse", False)

# Parallel and remote
definition.add_optional("nprocesses", "positive_integer", "number of parallel processes to use", 1)
definition.add_optional("remote", "string", "remote host ID", choices=find_host_ids())

# -----------------------------------------------------------------
