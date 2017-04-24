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

# Number of frames
definition.add_optional("nframes", "positive_integer", "number of frames", 2)

# Size of frames
definition.add_optional("npixels", "positive_integer", "number of pixels of the frames", 500)

# Number of sources
definition.add_optional("nrandom_sources", "positive_integer", "number of point sources", 100)

# Flags
definition.add_flag("vary_fwhm", "vary the FWHM", False)

# PSF model
definition.add_fixed("psf_model", "psf model", "gaussian")

definition.add_optional("noise_stddev", "real", "stddev of noise", 5.)

# -----------------------------------------------------------------
