#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.analyse.statistics Contains functions for analysing the statistics file

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import tables

# -----------------------------------------------------------------

# COLUMNS:
# line = [self.getIdentify(), generation]
# line.extend(stats.asTuple())

# Statistics:
# "rawMax": "Maximum raw score",
# "rawMin": "Minimum raw score",
# "rawAve": "Average of raw scores",
# "rawDev": "Standard deviation of raw scores",
# "rawVar": "Raw scores variance",
# "fitMax": "Maximum fitness",
# "fitMin": "Minimum fitness",
# "fitAve": "Fitness average",

# -----------------------------------------------------------------

column_names = ["rawMax", "rawMin", "rawAve", "rawDev", "rawVar", "fitMax", "fitMin", "fitAve"]

# -----------------------------------------------------------------

def load_statistics(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    # Open the file
    statistics = tables.from_file(path, format="ascii")

    # Set column names
    for index in range(len(column_names)): statistics.rename_column("col" + str(index+1), column_names[index])

    # Return the statistics table
    return statistics

# -----------------------------------------------------------------
