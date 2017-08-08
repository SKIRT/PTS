#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.scaling_plots Make the scaling plots from Verstocken et al., 2017, 
#  provided that you have the scaling test suite directory in your current directory

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
from subprocess import call, STDOUT

# Import the relevant PTS classes and modules
from pts.core.basics.log import log
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools import filesystem as fs
from pts.core.tools import time
from pts.core.tools import introspection
#from pts.core.tools.strings import split_except_within_single_quotes

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required
definition.add_required("suite_name", "string", "name of the scaling test suite directory (must be inside current directory)")

# Add flag
definition.add_flag("small", "make plots smaller so that they fit in one column of a publication and have a suiting resolution")

# Get configuration
config = parse_arguments("scaling_plots", definition)

# -----------------------------------------------------------------

# Set figsize
if config.small:
    figsize = "8,6"
    figsize_timelines = "8,8"
else:
    figsize = "12,9"
    figsize_timelines = "12,12"

# -----------------------------------------------------------------

# Locate the scaling test suite directory
suite_path = fs.join(fs.cwd(), config.suite_name)
if not fs.is_directory(suite_path): raise ValueError("The directory '" + suite_path + "' does not exist")

# -----------------------------------------------------------------

# Make directory for output
output_path = fs.create_directory_in(fs.cwd(), time.unique_name("scaling_plots"))

# Make subdirectories
single_node_path = fs.create_directory_in(output_path, "Single-node comparison")
multi_node_path = fs.create_directory_in(output_path, "Load balancing and multi-node scaling")
communication_path = fs.create_directory_in(output_path, "Communication")
hybridization_path = fs.create_directory_in(output_path, "Hybridization")
photon_packages_path = fs.create_directory_in(output_path, "Increased number of photon packages")
memory_path = fs.create_directory_in(output_path, "Memory scaling")

# -----------------------------------------------------------------

pts_path = introspection.pts_executable_path

# -----------------------------------------------------------------

FNULL = open(os.devnull, 'w')

# -----------------------------------------------------------------

# Navigate to the pure scaling directory
path = fs.join(suite_path, "1.pureScaling")
fs.change_cwd(path)

# Inform the user
log.info("Plotting the single-node scaling ...")

# Make plot of single node scaling
command = pts_path + " plot_scaling speedup total --debug --output '" + single_node_path + "' --plot/not_add_border --plot/not_add_legend_border --not_fit --plot/figsize " + figsize
#command = split_string_except_within_single_quotes(command)
#print(command) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------

# Navigate to the multi node scaling directory
path = fs.join(suite_path, "4.multinodeScaling")
fs.change_cwd(path)

# Inform the user
log.info("Plotting the multi-node scaling and timelines ...")

# Make plots of multi node scaling
command = pts_path + " plot_scaling efficiency,speedup,runtime total --debug --output '" + multi_node_path + "' --plot/not_add_border --plot/not_add_legend_border --not_fit --plot/figsize " + figsize
#call(split_string_except_within_single_quotes(command)) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# Make timeline plots
command = pts_path + " plot_timelines --debug --output '" + multi_node_path + "' --figsize " + figsize_timelines
#call(split_string_except_within_single_quotes(command))
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------

# Navigate to the communication directory
path = fs.join(suite_path, "2.communication", "vsStandard")
fs.change_cwd(path)

# Inform the user
log.info("Plotting the communication times ...")

# Make plots of communication time
command = pts_path + " plot_scaling runtime communication --debug --output '" + communication_path  + "' --hetero --split_communication " \
          "--communication_subphases 'stellar absorption communication,emission spectra communication' --not_fit " \
          "--plot/not_add_border --plot/not_add_legend_border --plot/legend_below --plot/figsize " + figsize
#call(split_string_except_within_single_quotes(command)) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------

# Navigate to the single node directory
path = fs.join(suite_path, "3.nodeScaling")
fs.change_cwd(path)

# Inform the user
log.info("Plotting hybridization on a single node ...")

# Make plots of the hybridization scaling
command = pts_path + " plot_scaling runtime total --debug --hybridisation --output '" + hybridization_path + "' --normalize_runtimes --plot/not_ylog --plot/not_add_border --plot/not_add_legend_border --plot/figsize " + figsize + " --input run1"
#call(split_string_except_within_single_quotes(command)) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------

# Navigate to the suite directory
fs.change_cwd(suite_path)

# Inform the user
log.info("Plotting scaling with increased number of photon packages ...")

# Make plots
command = pts_path + " plot_scaling --all_properties --input '4.multinodeScaling,5.multinodeScalingx10' --output '" + photon_packages_path + "' --debug --hetero --extrapolation/timing/npackages --fitting/pure_scaling_behaviour/communication 0,1,2 --not_fit --plot/not_add_border --plot/not_add_legend_border --plot/figsize 8,6 --plot/legend_below"
#call(split_string_except_within_single_quotes(command)) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------

# Navigate to the multi node scaling directory
path = fs.join(suite_path, "4.multinodeScaling")
fs.change_cwd(path)

# Inform the user
log.info("Plotting memory scaling ...")

# Make plots of memory scaling
command = pts_path + " plot_scaling memory total --debug --output '" + memory_path + "' --not_use_task_parallel --plot/not_connect_points --plot/figsize " + figsize + " --x_quantity processes"
#call(split_string_except_within_single_quotes(command)) # DOESN'T WORK??
call(command, shell=True, stdout=FNULL, stderr=STDOUT)

# -----------------------------------------------------------------
