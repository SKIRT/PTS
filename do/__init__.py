#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
#  Package initialization file
# -----------------------------------------------------------------

## \package pts.do Python scripts exposing PTS functionality to the command line.
#
# The scripts residing in this directory are intended to be invoked from the Terminal command line, so strictly
# speaking they are not part of a package since they cannot be imported as a module.
# Most scripts support one or more command line arguments and call on functions in the Python Toolkit for SKIRT (PTS),
# effectively exposing PTS functionality to the command line.
#
# To facilitate things even further, the \c __main__ module in this package allows executing any of the scripts
# residing in this directory from the command line, without having to specify its location.
#
