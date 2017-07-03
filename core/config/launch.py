#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.tools import filesystem as fs
from pts.core.launch.options import LoggingOptions, AnalysisOptions

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Input and output
definition.add_optional("input", "directory_path", "input directory for the simulation(s)", letter="i")
definition.add_optional("output", "directory_path", "output directory for the simulation(s)", fs.cwd(), letter="o", convert_default=True)

# Various flags
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")
definition.add_flag("emulate", "emulate the simulation while limiting computation")

# Other
definition.add_flag("keep", "keep remote input and output")
retrieve_type_choices = dict()
retrieve_type_choices["isrf"] = "interstellar radiation field strength"
retrieve_type_choices["abs"] = "absorption luminosities"
retrieve_type_choices["temp"] = "temperature"
retrieve_type_choices["sed"] = "SEDs"
retrieve_type_choices["image"] = "all datacubes"
retrieve_type_choices["image-total"] = "datacubes of total emission"
retrieve_type_choices["image-direct"] = "datacubes of direct emission"
retrieve_type_choices["image-transparent"] = "datacubes of transparent emission"
retrieve_type_choices["image-scattered"] = "datacubes of scattered emission"
retrieve_type_choices["image-dust"] = "datacubes of dust emission"
retrieve_type_choices["image-dustscattered"] = "datacubes of scattered dust emission"
retrieve_type_choices["celltemp"] = "temperature per dust cell"
retrieve_type_choices["log"] = "log files"
retrieve_type_choices["wavelengths"] = "wavelength files"
retrieve_type_choices["grid"] = "grid files"
retrieve_type_choices["grho"] = "grid dust density"
retrieve_type_choices["trho"] = "theoretical dust density"
retrieve_type_choices["convergence"] = "convergence file"
definition.add_optional("retrieve_types", "string_list", "types of output files that have to be retrieved/retained (None means everything)", choices=retrieve_type_choices)

# Special things
definition.add_flag("dry", "dry run (don't actually launch the simulations)", False)
definition.add_flag("attached", "launch the simulations in attached mode (only works if remotes without scheduling system are used)")

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Analysis options
definition.import_section_from_composite_class("analysis", "simulation analysis options", AnalysisOptions)

# -----------------------------------------------------------------
