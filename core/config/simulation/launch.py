#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Input and output
definition.add_optional("input", "directory_path", "input directory for the simulation(s)", letter="i")
definition.add_optional("output", "directory_path", "output directory for the simulation(s)", letter="o")

# Various flags
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")
definition.add_flag("brief", "enable brief console logging", letter="b")
definition.add_flag("verbose", "enable verbose logging", letter="v")
definition.add_flag("memory", "enable memory logging", letter="m")
definition.add_flag("allocation", "enable memory (de)allocation logging", letter="a")
definition.add_flag("emulate", "emulate the simulation while limiting computation", letter="e")

# Other
definition.add_flag("keep", "keep remote input and output")
retrieve_type_choices = ["isrf", "abs", "temp", "sed", "image", "image-total", "image-direct", "image-transparent", "image-scattered", "image-dust", "image-dustscattered", "celltemp", "log", "wavelengths", "grid", "grho", "trho", "convergence"]
definition.add_optional("retrieve_types", "string_list", "types of output files that have to be retrieved (None means everything should be downloaded)", choices=retrieve_type_choices)

definition.add_flag("dry", "dry run (don't actually launch the simulations)", False)
definition.add_flag("attached", "launch the simulations in attached mode (only works if remotes without scheduling system are used)")

# -----------------------------------------------------------------
