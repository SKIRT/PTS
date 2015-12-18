#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.install

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.prep.update import SkirtUpdater

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("--remote", type=str, help="update SKIRT on a remote system")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output for the update procedure")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a SkirtUpdater instance and run it
updater = SkirtUpdater.from_arguments(arguments)
updater.run()

# -----------------------------------------------------------------
