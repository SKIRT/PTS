#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_datacubes Show datacubes in current working directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.misc.datacubes import load_datacubes_in_cwd
from pts.magic.tools import plotting
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

definition = ConfigurationDefinition()
definition.add_flag("share_normalization", "share normalization between the frames", True)
definition.add_flag("show_axes", "show axes", True)
config = parse_arguments("show_datacubes", definition)

# -----------------------------------------------------------------

# Load
datacubes = load_datacubes_in_cwd()

# -----------------------------------------------------------------

# Loop over the datacubes
for instr_name in datacubes:
    for contribution in datacubes[instr_name]:

        # Set title
        title = instr_name + ": " + contribution

        # Plot
        datacube = datacubes[instr_name][contribution]
        plotting.plot_datacube(datacube, title=title, share_normalization=config.share_normalization, show_axes=config.show_axes)

# -----------------------------------------------------------------
