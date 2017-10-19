#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.modeling.config.maps import definition

# -----------------------------------------------------------------

# Remake?
definition.add_flag("remake", "remake already existing maps", False)

# REplot?
definition.add_flag("replot", "remake already existing plots", False)

# -----------------------------------------------------------------

# Smoothing
definition.add_flag("smooth", "smooth the maps by convolving them with a gaussian kernel", True)
definition.add_optional("smoothing_factor", "positive_real", "factor that determines the smoothing kernel FWHM based on the original FWHMs of the maps", 2.)

# -----------------------------------------------------------------

# Plot
definition.add_flag("plot", "plotting", False)

# -----------------------------------------------------------------
