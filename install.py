#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.install Install PTS after obtaining the source code.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from core.prep.installation import PTSInstaller

# -----------------------------------------------------------------

# Create the PTS installer
installer = PTSInstaller()

# Run the installer
installer.run()

# -----------------------------------------------------------------
