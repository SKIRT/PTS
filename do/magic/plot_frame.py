#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.plot_frame Plot a frame.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.core.frame import Frame
from pts.magic.tools import plotting

# -----------------------------------------------------------------

definition = ConfigurationDefinition()

definition.add_required("frame", "file_path", "path of the frame")

# Get the configuration
config = parse_arguments("plot_frame", definition)

# -----------------------------------------------------------------

# Open the frame
frame = Frame.from_file(config.frame)

frame.replace_nans_by_zeroes()
frame.replace_infs_by_zeroes()

# -----------------------------------------------------------------

plotting.plot_frame(frame, interval=[0,1], scale="linear")

# -----------------------------------------------------------------
