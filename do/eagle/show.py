#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.show Show selected records in the EAGLE SKIRT-runs database.
#
# This script prints the contents of selected records in the EAGLE SKIRT-runs database.
# The script expects a single command-line argument specifying the selection criteria in the form of
# an SQL expression that can be used in a WHERE clause. If no argument is present, all records are shown.

# -----------------------------------------------------------------

# Import standard modules
import sys

# Import the relevant PTS classes and modules
from pts.eagle import database

# -----------------------------------------------------------------

# get the command-line argument (default is a value that converts to "true" in SQL, meaning "all")
selection = sys.argv[1] if len(sys.argv) > 1 else "1"

# show the database records
db = database.Database()
db.show(db.select(selection))
db.close()

# -----------------------------------------------------------------
