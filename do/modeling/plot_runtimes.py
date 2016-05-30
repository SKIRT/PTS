#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_runtimes Make plots of the distribution of runtimes for the simulations that are part of
#  SKIRT radiative transfer modeling

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.distribution import Distribution
from pts.core.launch.timing import TimingTable

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Logging options
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path and the log path
arguments.path = fs.cwd()
log_path = fs.join(arguments.path, "log")

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(log_path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting plot_runtimes ...")

# -----------------------------------------------------------------

# Determine the path to the timing table
timing_table_path = fs.join(arguments.path, "fit", "timing.dat")

# Load the timing table
timing_table = TimingTable.read(timing_table_path)

# Keep a list of all the runtimes recorded for a certain remote host
runtimes_for_hosts = defaultdict(lambda: defaultdict(list))

# Loop over the entries in the timing table
for i in range(len(timing_table)):

    # Get the ID of the host and the cluster name for this particular simulation
    host_id = timing_table["Host id"][i]
    cluster_name = timing_table["Cluster name"][i]

    # Get the parallelization properties for this particular simulation
    cores = timing_table["Cores"][i]
    threads_per_core = timing_table["Threads per core"][i]
    processes = timing_table["Processes"][i]

    # Get the number of photon packages (per wavelength) used for this simulation
    packages = timing_table["Packages"][i]

    # Get the total runtime
    runtime = timing_table["Total runtime"][i]

    parallelization = (cores, threads_per_core, processes)

    runtimes_for_hosts[host_id][packages, parallelization].append(runtime)

# -----------------------------------------------------------------

bins = 25

# Loop over the different remote hosts
for host_id in runtimes_for_hosts:

    # Loop over the different configurations (packages, parallelization)
    for packages, parallelization in runtimes_for_hosts[host_id]:

        runtimes = runtimes_for_hosts[host_id][packages, parallelization]

        distribution = Distribution.from_values(runtimes, 15)

        # Plot the distribution
        distribution.plot()

# -----------------------------------------------------------------
