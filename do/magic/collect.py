#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.collect Collect ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.train import Collector

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("mode", type=str, help="the collection mode (star or saturation)")
parser.add_argument("--debug", type=str, help="debug mode")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create the collector
collector = Collector.from_arguments(arguments)

# Run the collector
collector.run()

# -----------------------------------------------------------------
