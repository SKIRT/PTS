#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.core.install

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.prep.installation import SkirtInstaller

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("--remote", type=str, help="install SKIRT remotely")
parser.add_argument("--private", action="store_true", help="use the private SKIRT repository")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output for the installation procedure")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Create a SkirtInstaller instance and run it
installer = SkirtInstaller.from_arguments(arguments)
installer.run()

# -----------------------------------------------------------------
