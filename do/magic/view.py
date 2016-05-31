#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.view Open a FITS image with the magic viewer.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

from astropy.io import fits

# Import the relevant PTS classes and modules
from pts.magic.view import MagicViewer
from pts.core.tools import time
from pts.magic.core.frame import Frame
from pts.magic.core.image import Image

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# The filename
parser.add_argument("filename", type=str, help="the name (or path) of the FITS file")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# control_panels.py : <-------

#from __future__ import absolute_import
#from ztv.source_panel import SourcePanel
#from ztv.plot_panel import PlotPanel
#from ztv.phot_panel import PhotPanel
#from ztv.stats_panel import StatsPanel
#from ztv.color_panel import ColorPanel
#from ztv_examples.fits_faker_panel.fits_faker_panel import FitsFakerPanel

#control_panels_to_load = [("Source", SourcePanel),
#                          ("Color", ColorPanel),
#                          ("Plot", PlotPanel),
#                          ("Stats", StatsPanel),
#                          ("Phot", PhotPanel),
#                          ("Faker", FitsFakerPanel)
#                          ]

# ----->

#viewer = MagicViewer(control_panels_module_path='control_panels')


# Load the image first
image = Image.from_file(arguments.filename)

# Create the viewer
viewer = MagicViewer(title="Magic viewer")

# Load the image into the viewer
viewer.load_image(image)

# Set options
viewer.control_panel('Color')
viewer.cmap("hot")
viewer.scaling('Log')

# Wait for the GUI to be closed
viewer.wait()

# -----------------------------------------------------------------
