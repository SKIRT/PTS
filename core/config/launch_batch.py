#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.basics.host import find_host_ids

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition()

definition.add_optional("remotes", "string_list", "remote host IDs to use", choices=find_host_ids(), default=find_host_ids())
definition.add_optional("extra_remote", "string", "remote host ID to use for the extra simulations", choices=find_host_ids())

definition.add_optional("input", "directory_path", "input directory for the simulations")
definition.add_optional("output", "directory_path", "output directory for the simulations")

definition.add_flag("keep", "keep remote input and output", False)
definition.add_flag("shared_input", "whether the different simulations share their input folder", False)

definition.add_optional("retrieve_types", "string_list", "The types of output that should be retrieved for the simulation (None means everything should be downloaded)")

definition.add_optional("cores_per_process", "integer", "the number of cores to use per process (if the parallelization is not specified by the user for a certain remote host)")

definition.add_flag("data_parallel", "enable data parallelization mode (if the parallelization is not specified by the user of the batch launcher for a certain remote host)")

definition.add_flag("group_simulations", "group multiple simulations in one job", False)

definition.add_optional("group_walltime", "real", "preferred walltime per job of grouped simulations")

definition.add_optional("timing_table_path", "file_path", "path to the timing table")
definition.add_optional("memory_table_path", "file_path", "path to the memory table")

definition.add_flag("dry", "dry run (don't actually launch the simulations)", False)

definition.add_flag("attached", "launch the simulations in attached mode (only works if remotes without scheduling system are used)")

definition.add_section("logging", "logging options")

definition.sections["logging"].add_flag("brief", "brief console logging", False)
definition.sections["logging"].add_flag("verbose", "verbose logging mode", False)
definition.sections["logging"].add_flag("memory", "memory logging", False)
definition.sections["logging"].add_flag("allocation", "memory (de)allocation logging", False)
definition.sections["logging"].add_optional("allocation_limit", "real", "memory (de)allocation logging lower limit (GB)", 1e-5)

definition.add_section("analysis", "simulation analysis options")

definition.sections["analysis"].add_flag("relative", "Treat the specified paths for the simulation output directories and directories for analysis (plotting, extraction, misc) as relative to the ski file (the simulation directory)", False)

definition.sections["analysis"].add_section("extraction", "extraction options")
definition.sections["analysis"].sections["extraction"].add_optional("path", "directory_path", "path of the directory where the extracted data should be placed")
definition.sections["analysis"].sections["extraction"].add_flag("progress", "extract progress", False)
definition.sections["analysis"].sections["extraction"].add_flag("timeline", "extract timeline", False)
definition.sections["analysis"].sections["extraction"].add_flag("memory", "extract memory", False)

definition.sections["analysis"].add_section("plotting", "plotting options")
definition.sections["analysis"].sections["plotting"].add_optional("path", "directory_path", "path of the directory where the plots should be placed")
definition.sections["analysis"].sections["plotting"].add_flag("seds", "plot seds")
definition.sections["analysis"].sections["plotting"].add_flag("grids", "plot grids")
definition.sections["analysis"].sections["plotting"].add_flag("progress", "plot progress")
definition.sections["analysis"].sections["plotting"].add_flag("timeline", "plot timeline")
definition.sections["analysis"].sections["plotting"].add_flag("memory", "plot memory")
definition.sections["analysis"].sections["plotting"].add_optional("reference_sed", "file_path", "reference SED file path")

definition.sections["analysis"].add_section("misc", "miscellaneous options")
definition.sections["analysis"].sections["misc"].add_flag("rgb", "make rgb image(s)")
definition.sections["analysis"].sections["misc"].add_flag("wave", "make wave movie")
definition.sections["analysis"].sections["misc"].add_flag("fluxes", "calculate observed fluxes")
definition.sections["analysis"].sections["misc"].add_flag("images", "create observed images")
definition.sections["analysis"].sections["misc"].add_optional("observation_filters", "string_list", "observation filters")
definition.sections["analysis"].sections["misc"].add_optional("observation_instruments", "string_list", "simulation instruments for which to create observations")
definition.sections["analysis"].sections["misc"].add_optional("make_images_remote", "string", "remote host ID to create the observed images on", choices=find_host_ids())
definition.sections["analysis"].sections["misc"].add_optional("images_wcs", "file_path", "FITS file with WCS info for the observations")
definition.sections["analysis"].sections["misc"].add_optional("images_unit", "string", "unit for the images")
definition.sections["analysis"].sections["misc"].add_optional("images_kernels", "dictionary", "kernel paths for convolving the observed images")

definition.sections["analysis"].add_optional("timing_table_path", "file_path", "timing table path for recording the simulation timing results")
definition.sections["analysis"].add_optional("memory_table_path", "file_path", "memory table path for recording the simulation memory usage")

# -----------------------------------------------------------------
