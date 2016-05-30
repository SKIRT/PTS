#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.magic.view

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.magic.view import MagicViewer

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

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

viewer = MagicViewer()

# -----------------------------------------------------------------
