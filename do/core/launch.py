#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.launch Launch a SKIRT simulation locally or remotely with the best performance, based on the
#  current load of the system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse

# Import the relevant PTS classes and modules
from pts.core.launch.launcher import SkirtLauncher
from pts.core.launch.remotelauncher import SkirtRemoteLauncher
from pts.core.tools import logging, time, parsing
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import Configuration

# -----------------------------------------------------------------

# TODO: work on this further

# Create the configuration
config = Configuration()

# Add required arguments
config.add_required("filename", "absolute_path", "the name/path of the ski file")

# Add positional arguments
config.add_positional_optional("remote", str, "the remote host on which to run the simulation (if none is specified, the simulation is run locally")
config.add_optional("input", "absolute_path", "the simulation input directory", letter="i")
config.add_optional("output", "absolute_path", "the simulation output directory", letter="o")
config.add_optional("cluster", str, "the name of the cluster", letter="c")
config.add_optional("parallel", "int_tuple", "the parallelization scheme (processes, threads)", letter="p")
config.add_optional("walltime", "duration", "an estimate for the walltime of the simulation for the specified parallelization scheme")

# Flags
config.add_flag("relative", "treats the given input and output paths as being relative to the ski/fski file")
config.add_flag("brief", "enable brief console logging", letter="b")
config.add_flag("verbose", "enable verbose logging", letter="v")
config.add_flag("memory", "enable memory logging", letter="m")
config.add_flag("allocation", "enable memory (de)allocation logging", letter="a")
config.add_flag("emulate", "emulate the simulation while limiting computation", letter="e")

config.add_section("extraction")

config.add_section("plotting")

config.add_section("misc")

# Read the configuration settings from the provided command-line arguments
config.read()

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="the name of the ski/fski file")
parser.add_argument("-i", "--input", type=str, help="the simulation input path")
parser.add_argument("-o", "--output", type=str, help="the simulation output path")
parser.add_argument("-r", "--remote", type=str, help="run the simulation remotely")
parser.add_argument("-c", "--cluster", type=str, help="add the name of the cluster if different from the default")
parser.add_argument("-p", "--parallel", type=parsing.int_tuple, help="specify the parallelization scheme (processes, threads per process)")
parser.add_argument("-t", "--walltime", type=parsing.duration, help="specify an estimate for the walltime of the simulation for the specified parallelization scheme")
parser.add_argument("--relative", action="store_true", help="treats the given input and output paths as being relative to the ski/fski file")
parser.add_argument("--brief", action="store_true", help="enable brief console logging")
parser.add_argument("-v", "--verbose", action="store_true", help="enable verbose logging mode")
parser.add_argument("-m", "--memory", action="store_true", help="enable memory logging mode")
parser.add_argument("-a", "--allocation", action="store_true", help="enable memory (de)allocation logging mode")
parser.add_argument("-e", "--emulate", action="store_true", help="emulate the simulation while limiting computation")
parser.add_argument("--extractprogress", action="store_true", help="extract the progress from the log files")
parser.add_argument("--extracttimeline", action="store_true", help="extract the timeline from the log files")
parser.add_argument("--extractmemory", action="store_true", help="extract the memory usage from the log files")
parser.add_argument("--plotseds", action="store_true", help="make plots of the output SEDs")
parser.add_argument("--plotgrids", action="store_true", help="make plots of the dust grid")
parser.add_argument("--plotprogress", action="store_true", help="make plots of the progress of the different processes as a function of time")
parser.add_argument("--plottimeline", action="store_true", help="make a plot of the timeline for the different processes")
parser.add_argument("--plotmemory", action="store_true", help="make a plot of the memory consumption as a function of time")
parser.add_argument("--refsed", type=str, help="specify the path to a reference SED file against which the simulated SKIRT SEDs should be plotted")
parser.add_argument("--makergb", action="store_true", help="add this option to make RGB images from the SKIRT output")
parser.add_argument("--makewave", action="store_true", help="add this option to make a wave movie from the SKIRT output")
parser.add_argument("--fluxes", action="store_true", help="add this option to calculate observed fluxes from the SKIRT output SEDs")
parser.add_argument("--images", action="store_true", help="add this option to make observed images from the SKIRT output datacubes")
parser.add_argument("--filters", parsing.string_list, help="the names of the filters for which to recreate the observations (seperated by commas)")
parser.add_argument("--instruments", parsing.string_list, help="the names of the instruments for which to recreate the observations (seperated by commas)")
parser.add_argument("--wcs", type=str, help="the path to the FITS file for which the WCS should be set as the WCS of the recreated observed images")
parser.add_argument("--unit", type=str, help="the unit to which the recreated observed images should be converted")
parser.add_argument("--debug", action="store_true", help="add this option to enable debug output")
parser.add_argument('--report', action='store_true', help='write a report file')
parser.add_argument("--keep", action="store_true", help="add this option to keep the remote input and output")
parser.add_argument("--retrieve", type=parsing.string_list, help="specify the types of output files that have to be retrieved")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Determine the full path to the parameter file
arguments.filepath = fs.absolute(arguments.file)

# Determine the full path to the input and output directories
if arguments.input is not None: arguments.input_path = fs.absolute(arguments.input)
if arguments.output is not None: arguments.output_path = fs.absolute(arguments.output)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), time.unique_name("launch") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting launch ...")

# -----------------------------------------------------------------

# Either create a SkirtRemoteLauncher or a SkirtLauncher
if arguments.remote: launcher = SkirtRemoteLauncher.from_arguments(arguments)
else: launcher = SkirtLauncher.from_arguments(arguments)

# Run the launcher (the simulation is performed locally or remotely depending on which launcher is used)
launcher.run()

# -----------------------------------------------------------------
