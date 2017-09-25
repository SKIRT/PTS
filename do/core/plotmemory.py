#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plotmemory Make plots of the memory usage for a SKIRT simulation.
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.extract.memory import extract_memory_cwd, MemoryUsageTable
from pts.core.plot.memory import MemoryPlotter
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("table", "save the extracted memory table")

# Get the configuration
config = parse_arguments("plotmemory", definition)

# -----------------------------------------------------------------

# Look for a file in the current working directory that contains extracted memory information
memory_table_path = fs.join(fs.cwd(), "memory.dat")
if fs.is_file(memory_table_path): table = MemoryUsageTable.from_file(memory_table_path)

# If extracted memory information is not present, first perform the extraction
else: table = extract_memory_cwd()

# -----------------------------------------------------------------

if config.table and not fs.is_file(memory_table_path): table.saveto(memory_table_path)

# -----------------------------------------------------------------

# Determine the path to the plotting directory
plot_path = fs.join(fs.cwd())

# Create a MemoryPlotter instance
plotter = MemoryPlotter()

# Run the memory plotter
plotter.run(table, plot_path)

# -----------------------------------------------------------------
