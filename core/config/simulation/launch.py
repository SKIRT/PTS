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
definition.add_flag("emulate", "emulate the simulation while limiting computation")

# Other
definition.add_flag("keep", "keep remote input and output")
retrieve_type_choices = ["isrf", "abs", "temp", "sed", "image", "image-total", "image-direct", "image-transparent", "image-scattered", "image-dust", "image-dustscattered", "celltemp", "log", "wavelengths", "grid", "grho", "trho", "convergence"]
retrieve_type_descriptions = dict()
retrieve_type_descriptions["isrf"] = "interstellar radiation field strength"
retrieve_type_descriptions["abs"] = "absorption luminosities"
retrieve_type_descriptions["temp"] = "temperature"
retrieve_type_descriptions["sed"] = "SEDs"
retrieve_type_descriptions["image"] = "all datacubes"
retrieve_type_descriptions["image-total"] = "datacubes of total emission"
retrieve_type_descriptions["image-direct"] = "datacubes of direct emission"
retrieve_type_descriptions["image-transparent"] = "datacubes of transparent emission"
retrieve_type_descriptions["image-scattered"] = "datacubes of scattered emission"
retrieve_type_descriptions["image-dust"] = "datacubes of dust emission"
retrieve_type_descriptions["image-dustscattered"] = "datacubes of scattered dust emission"
retrieve_type_descriptions["celltemp"] = "temperature per dust cell"
retrieve_type_descriptions["log"] = "log files"
retrieve_type_descriptions["wavelengths"] = "wavelength files"
retrieve_type_descriptions["grid"] = "grid files"
retrieve_type_descriptions["grho"] = "grid dust density"
retrieve_type_descriptions["trho"] = "theoretical dust density"
retrieve_type_descriptions["convergence"] = "convergence file"
definition.add_optional("retrieve_types", "string_list", "types of output files that have to be retrieved/retained (None means everything)", choices=retrieve_type_choices, choice_descriptions=retrieve_type_descriptions)

definition.add_flag("dry", "dry run (don't actually launch the simulations)", False)
definition.add_flag("attached", "launch the simulations in attached mode (only works if remotes without scheduling system are used)")

# -----------------------------------------------------------------
