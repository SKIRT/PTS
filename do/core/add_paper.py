#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.add_paper Add a paper to the local collection of papers.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools.papers import Papers

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Arguments
parser.add_argument("label", type=str, help="label")
parser.add_argument("url", type=str, help="url")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

papers = Papers()

papers.add_entry(arguments.label, arguments.url)

papers.save()

# -----------------------------------------------------------------
