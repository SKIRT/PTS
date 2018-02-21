#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.core.launch.options import LoggingOptions, AnalysisOptions
from pts.core.simulation.output import output_type_choices

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Input and output
definition.add_optional("simulation_input", "directory_path", "input directory for the simulation(s)", letter="i")
definition.add_optional("simulation_output", "directory_path", "output directory for the simulation(s)", fs.cwd(), letter="o", convert_default=True)

# Various flags
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")
definition.add_flag("emulate", "emulate the simulation while limiting computation")

# Other
definition.add_flag("keep", "keep remote input and output")
definition.add_flag("keep_input", "keep remote input specifically")
definition.add_optional("retrieve_types", "string_list", "types of output files that have to be retrieved/retained (None means everything)", choices=output_type_choices)

# Special things
definition.add_flag("dry", "dry run (don't actually launch the simulations)", False)
definition.add_flag("attached", "launch the simulations in attached mode (only works if remotes without scheduling system are used)")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Analysis options
definition.import_section_from_composite_class("analysis", "simulation analysis options", AnalysisOptions)

# Show stuff
definition.add_flag("show", "show", True)

# Show the output of finished simulations
definition.add_flag("show_finished", "show the output of finished simulations", False)

# -----------------------------------------------------------------

# Retrieve and analyse finished simulations?
definition.add_flag("retrieve", "retrieve finished simulations", True)
definition.add_flag("analyse", "analyse retrieved simulations", True)

# -----------------------------------------------------------------
