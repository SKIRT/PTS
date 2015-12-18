#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plottimeline Plot a timeline for the different simulation phases of a SKIRT simulation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.extract.timeline import TimeLineExtractor
from pts.core.plot.timeline import TimeLinePlotter

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted timeline information
timeline_file_path = os.path.join(os.getcwd(), "timeline.dat")
if os.path.isfile(timeline_file_path):

    # Create a TimeLineExtractor instance from the saved timeline data
    extractor = TimeLineExtractor.open_table(timeline_file_path)

# If extracted timeline information is not present, first perform the extraction
else:

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Create a new TimeLineExtractor instance
    extractor = TimeLineExtractor()

    # Run the extractor
    extractor.run(simulation)

# Determine the path to the plotting directory
plot_path = os.path.join(os.getcwd())

# Create a MemoryPlotter instance
plotter = TimeLinePlotter()

# Run the memory plotter
plotter.run(extractor.table, plot_path)

# -----------------------------------------------------------------
