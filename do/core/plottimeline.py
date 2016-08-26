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

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.extract.timeline import TimeLineExtractor, TimeLineTable
from pts.core.plot.timeline import TimeLinePlotter
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("table", "save the extracted timeline table")

# Get configuration
setter = ArgumentConfigurationSetter("plottimeline")
config = setter.run(definition)

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted timeline information
timeline_table_path = fs.join(fs.cwd(), "timeline.dat")
if fs.is_file(timeline_table_path): table = TimeLineTable.from_file(timeline_table_path)

# If extracted timeline information is not present, first perform the extraction
else:

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Create a new TimeLineExtractor instance
    extractor = TimeLineExtractor()

    # Run the extractor and get the timeline table
    table = extractor.run(simulation)

# -----------------------------------------------------------------

if config.table and not fs.is_file(timeline_table_path): table.saveto(timeline_table_path)

# -----------------------------------------------------------------

# Determine the path to the plotting directory
plot_path = fs.join(fs.cwd())

# Create a TimeLinePlotter instance
plotter = TimeLinePlotter()

# Run the timeline plotter
plotter.run(table, plot_path)

# -----------------------------------------------------------------
