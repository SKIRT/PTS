#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

path = fs.cwd()

# -----------------------------------------------------------------

# Create the DustPedia sample
sample = DustPediaSample()
galaxy_name = sample.get_name("M81")

# Create the database
database = DustPediaDatabase()

# Login
username, password = get_account()
database.login(username, password)

# Get the galaxy info
info = database.get_galaxy_info(galaxy_name)

print(info)

# -----------------------------------------------------------------
