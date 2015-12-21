#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.eagle Processing and visualizing EAGLE output data with SKIRT.
#
# The eagle package provides Python modules that help processing and visualizing EAGLE output data with SKIRT.
#
# EAGLE is a cosmological hydrodynamical simulation conducted by the VIRGO consortium. The simulation
# includes detailed recipes for star formation, gas cooling, stellar evolution, and feedback from supernovae and AGN.
# Using more than 10 billion particles, it numerically resolves thousands of galaxies in a representative cosmological
# volume that also contains groups and clusters.
#
# The EAGLE simulation output is stored in a (large) set of data files in the HDF5 format, documented at the
# <a href="http://www.hdfgroup.org/HDF5/">HFD5 home page</a>. The output is organized in \em snapshots, where
# each snapshot represents the state of the universe at a particular time (or equivalently, redshift).
#
# The classes and functions in this package allow:
#  - extracting relevant information from the EAGLE output in a format that can be used in SKIRT
#  - scheduling large numbers of SKIRT simulations and managing/retrieving the results
#  - visualizing the SKIRT results in meaningful ways
#
# SKIRT is a Monte Carlo radiative transfer code developed by the Astronomical Observatory at the Ghent University.
# This package also uses the functionality offered by PTS, the Python toolkit for SKIRT.
#
