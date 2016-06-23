#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.evolve This is the main module of the pyevolve,
#  every other module is above this namespace
#
# This package ...
#

# -----------------------------------------------------------------

__all__ = ["Consts", "Crossovers", "DBAdapters", "FunctionSlot",
           "G1DBinaryString", "G1DList", "G2DBinaryString",
           "G2DList", "GAllele", "GenomeBase", "GPopulation",
           "GSimpleGA", "GTree", "Initializators",
           "Migration", "Mutators", "Network", "Scaling", "Selectors",
           "Statistics", "Util"]

__version__ = '0.6'
__author__ = 'Christian S. Perone'

import Consts
import sys

if sys.version_info[:2] < Consts.CDefPythonRequire:
   raise Exception("Python 2.5+ required, the version %s was found on your system !" % (sys.version_info[:2],))

del sys

# -----------------------------------------------------------------
