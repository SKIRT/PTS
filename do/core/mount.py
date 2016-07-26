#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.mount Mount a remote configured in PTS into a local directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()


parser.add_argument("remote", type=str, help="the remote host to send to")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------


# -----------------------------------------------------------------
