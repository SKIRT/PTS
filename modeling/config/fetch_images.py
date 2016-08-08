#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from .fetch import definition
from ...core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Add required settings
definition.add_required("remote", "string", "the remote host to use for creating the GALEX and SDSS data", choices=find_host_ids())

# Add flags
definition.add_flag("errors", "also download the error frames from the DustPedia archive")

# H alpha
definition.add_required("halpha_url", "url", "the url of the h-alpha image", choices=["https://ned.ipac.caltech.edu/img/2001ApJ...559..878H/MESSIER_081:I:Ha:hwb2001.fits.gz"]) # this link is for M81
definition.add_optional("halpha_flux", "quantity", "the flux of the h-alpha image to which the image should be rescaled", default="7.8e40 erg/s", convert_default=True) # this flux is for M81
# The total H-alpha flux (reference: FAR-ULTRAVIOLET AND Ha IMAGING OF NEARBY SPIRAL GALAXIES: THE OB STELLAR,
# POPULATION IN THE DIFFUSE IONIZED GAS (Hoopes et. al 2001)

# -----------------------------------------------------------------
