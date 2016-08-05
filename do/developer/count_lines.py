#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.count_lines Count the number of lines of code in the PTS project. 

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from inspect import getmembers, isfunction, getdoc

# Import the relevant PTS classes and modules
from pts.core.tools import parsing
from pts.core.tools import formatting as fmt
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

nlines = 0

nlines_per_file = dict()

for path in fs.files_in_path(introspection.pts_package_dir, extension="py", recursive=True):

    # Open the module file
    with open(path, 'r') as pyfile:

        nlines_per_file[path] = 0

        for line in pyfile:

            line = line.rstrip("\n")

            if line == "": continue
            if line.startswith("#"): continue

            nlines += 1

            nlines_per_file[path] += 1

# State the number of lines
print("PTS contains " + str(nlines) + " lines of python code")

# State the number of lines for each module
for path in nlines_per_file:
    print(" - " + path + ": " + str(nlines_per_file[path]) + " lines")

# -----------------------------------------------------------------
