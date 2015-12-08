#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.plotmemory Make plots of the memory usage for a SKIRT simulation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from pts.core.simulation.simulation import createsimulations
from pts.core.extract.memory import MemoryExtractor
from pts.core.plot.memory import MemoryPlotter

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted memory information
memory_file_path = os.path.join(os.getcwd(), "memory.dat")
if os.path.isfile(memory_file_path):

    # Create a MemoryExtractor instance from the saved memory data
    extractor = MemoryExtractor.open_table(memory_file_path)

# If extracted memory information is not present, first perform the extraction
else:

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Create a new MemoryExtractor instance
    extractor = MemoryExtractor()

    # Run the extractor
    extractor.run(simulation)

# Determine the path to the plot file
plot_path = os.path.join(os.getcwd(), "memory.pdf")

# Create a MemoryPlotter instance
plotter = MemoryPlotter()

# Run the memory plotter
plotter.run(extractor.table, plot_path)

# -----------------------------------------------------------------
