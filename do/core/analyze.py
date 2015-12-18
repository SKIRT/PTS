#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.analyze Perform analysis for some specific SPH data file(s).
#
# This script performs analysis for some specific SPH data file(s).
# The in/out filenames and other parameters are hardcoded in the script.
# Thus the script mainly serves as an example of how to use the pts.sphanalyze module.
#

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.prep.sphanalyze import *

# -----------------------------------------------------------------

inpath = "/Users/pcamps/EAGLE/Snapshot50/dat/"

# get a list of galaxies in order of increasing nr of particles
galaxies = np.loadtxt(inpath + "contents.dat", dtype=np.int, usecols=(0,1))[::-1]

# make a circle plot for each one
for group,subgroup in galaxies:
    name = "galaxy_{0}_{1}_".format(group, subgroup)
    analyze(inpath+name+"stars.dat", plotfile=name+"stars.png", scale=0.001, circleplot=True)
    analyze(inpath+name+"gas.dat", plotfile=name+"gas.png", scale=0.001, circleplot=True)

# -----------------------------------------------------------------
