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

# Add required arguments
definition.add_required("filename", "file_path", "the name/path of the ski file")

# Add positional arguments
definition.add_positional_optional("remote", "string", "the remote host on which to run the simulation (if none is specified, the simulation is run locally")
definition.add_optional("input", "directory_path", "the simulation input directory", letter="i")
definition.add_optional("output", "directory_path", "the simulation output directory", letter="o")
definition.add_optional("cluster", "string", "the name of the cluster", letter="c")
definition.add_optional("parallel", "integer_tuple", "the parallelization scheme (processes, threads)", letter="p")
definition.add_optional("walltime", "duration", "an estimate for the walltime of the simulation for the specified parallelization scheme")

# Various flags
definition.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")
definition.add_flag("brief", "enable brief console logging", letter="b")
definition.add_flag("verbose", "enable verbose logging", letter="v")
definition.add_flag("memory", "enable memory logging", letter="m")
definition.add_flag("allocation", "enable memory (de)allocation logging", letter="a")
definition.add_flag("emulate", "emulate the simulation while limiting computation", letter="e")

# Extraction settings
definition.add_section("extraction", "settings for extraction data from the simulation's log files")
definition.sections["extraction"].add_flag("progress", "extract information about the progress in the different simulation phases")
definition.sections["extraction"].add_flag("timeline", "extract timeline information for the different simulation phases on the different processes")
definition.sections["extraction"].add_flag("memory", "extract information about the memory usage during the simulation")

definition.add_section("plotting", "settings for plotting simulation data")
definition.sections["plotting"].add_flag("seds", "make plots of the simulated SEDs")
definition.sections["plotting"].add_flag("grids", "make plots of the dust grid")
definition.sections["plotting"].add_flag("progress", "make plots of the progress of the simulation phases as a function of time")
definition.sections["plotting"].add_flag("timeline", "plot the timeline for the different processes")
definition.sections["plotting"].add_flag("memory", "plot the memory consumption as a function of time")

definition.add_section("misc", "settings for creating data of various types from the simulation output")
definition.sections["misc"].add_flag("wave", "make a wavelength movie through the simulated datacube(s)")
definition.sections["misc"].add_flag("rgb", "make RGB images from the simulated datacube(s)")
definition.sections["misc"].add_flag("fluxes", "calculate observed fluxes from the SKIRT output SEDs")
definition.sections["misc"].add_flag("images", "make observed images form the simulated datacube(s)")
definition.sections["misc"].add_optional("refsed", "file_path", "the path to a reference SED file against which the simulated SKIRT SEDs should be plotted")
definition.sections["misc"].add_optional("filters", "string_list", "the names of the filters for which to recreate the observations")
definition.sections["misc"].add_optional("instruments", "string_list", "the names of the instruments for which to recreate the observations")
definition.sections["misc"].add_optional("wcs", "file_path", "the path to the FITS file for which the WCS should be set as the WCS of the recreated observed images")
definition.sections["misc"].add_optional("unit", "string", "the unit to which the recreated observed images should be converted")

# OTHER SETTINGS
definition.add_flag("keep", "keep remote input and output")
retrieve_type_choices = ["isrf", "abs", "temp", "sed", "image", "image-total", "image-direct", "image-transparent", "image-scattered", "image-dust", "image-dustscattered", "celltemp", "log", "wavelengths", "grid", "grho", "trho", "convergence"]
definition.add_optional("retrieve", "string_list", "the types of output files that have to be retrieved", choices=retrieve_type_choices)

# -----------------------------------------------------------------
