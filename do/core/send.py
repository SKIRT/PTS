#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.send Send a file or directory to one of the configured remotes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

parser.add_argument("filename", type=str, help="the name (or path) of the file or directory to send")
parser.add_argument("remote", type=str, help="the remote host to send to")
parser.add_argument("--path", type=str, help="the path to the remote directory to send to")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------


# -----------------------------------------------------------------
