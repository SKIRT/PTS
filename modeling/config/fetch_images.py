#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.fetch import definition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Add required settings
definition.add_positional_optional("remote", "string", "remote host to use for creating the GALEX and SDSS data", choices=find_host_ids(schedulers=False))
definition.add_flag("attached", "run remote in attached mode", False)

# Add flags
definition.add_flag("errors", "also download the error frames from the DustPedia archive")

# H alpha
definition.add_optional("halpha_url", "url", "url of the h-alpha image") # optional for when fetch_images is relaunched, and the Halpha image is already present
definition.add_optional("halpha_flux", "quantity", "flux of the h-alpha image to which the image should be rescaled")
# The total H-alpha flux (reference: FAR-ULTRAVIOLET AND Ha IMAGING OF NEARBY SPIRAL GALAXIES: THE OB STELLAR,
# POPULATION IN THE DIFFUSE IONIZED GAS (Hoopes et. al 2001)

# Flags
definition.add_flag("make_poisson", "perform the poisson error mosaicing", True)

# Number of processes for Poisson error map making
definition.add_optional("nprocesses", "positive_integer", "number of parallel processes", 1)

# Advanced
definition.add_optional("max_nobservations_mosaic", "positive_integer", "maximum number of observations to use for the mosaic and poisson frames (GALEX and SDSS)")

# Other images
definition.add_optional("other_urls", "url_list", "urls of other images to fetch")

# -----------------------------------------------------------------
