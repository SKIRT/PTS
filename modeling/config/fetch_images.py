#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from .fetch import definition
from ...core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Add required settings
definition.add_required("remote", "string", "remote host to use for creating the GALEX and SDSS data", choices=find_host_ids(schedulers=False))

# Add flags
definition.add_flag("errors", "also download the error frames from the DustPedia archive")

# H alpha
definition.add_required("halpha_url", "url", "url of the h-alpha image")
definition.add_optional("halpha_flux", "quantity", "flux of the h-alpha image to which the image should be rescaled")
# The total H-alpha flux (reference: FAR-ULTRAVIOLET AND Ha IMAGING OF NEARBY SPIRAL GALAXIES: THE OB STELLAR,
# POPULATION IN THE DIFFUSE IONIZED GAS (Hoopes et. al 2001)

# -----------------------------------------------------------------
