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
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time, filesystem, tables
from pts.core.basics.distribution import Distribution

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()

# Basic options
parser.add_argument("path", type=str, nargs='?', help="the modeling path")

# Logging options
parser.add_argument("--debug", action="store_true", help="enable debug logging mode")
parser.add_argument("--report", action='store_true', help='write a report file')

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# Set the modeling path
if arguments.path is None: arguments.path = filesystem.cwd()

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = filesystem.join(arguments.path, time.unique_name("log") + ".txt") if arguments.report else None

# Determine the log level
level = "DEBUG" if arguments.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting plot_runtimes ...")

# -----------------------------------------------------------------

# Determine the path to the runtime table
runtime_table_path = filesystem.join(arguments.path, "fit", "runtimes.dat")

# Load the runtime table
runtimes_table = tables.from_file(runtime_table_path, format="ascii.ecsv")

# Keep a list of all the runtimes recorded for a certain remote host
runtimes_for_hosts = defaultdict(lambda: defaultdict(list))

# Loop over the entries in the runtime table
# "Simulation name", "Host id", "Cluster name", "Cores", "Hyperthreads per core", "Processes", "Packages", "Runtime"
for i in range(len(runtimes_table)):

    # Get the ID of the host and the cluster name for this particular simulation
    host_id = runtimes_table["Host id"][i]
    cluster_name = runtimes_table["Cluster name"][i]

    # Get the parallelization properties for this particular simulation
    cores = runtimes_table["Cores"][i]
    threads_per_core = runtimes_table["Hyperthreads per core"][i]
    processes = runtimes_table["Processes"][i]

    # Get the number of photon packages (per wavelength) used for this simulation
    packages = runtimes_table["Packages"][i]

    # Get the total runtime
    runtime = runtimes_table["Runtime"][i]

    #parallelization = Parallelization(cores, threads_per_core, processes)

    parallelization = (cores, threads_per_core, processes)

    runtimes_for_hosts[host_id][packages, parallelization].append(runtime)

# -----------------------------------------------------------------

bins = 25

# Loop over the different remote hosts
for host_id in runtimes_for_hosts:

    # Loop over the different configurations (packages, parallelization)
    for packages, parallelization in runtimes_for_hosts[host_id]:

        print("host ID:", host_id)
        print("packages:", packages)
        print("parallelization:", parallelization)

        runtimes = runtimes_for_hosts[host_id][packages, parallelization]

        distribution = Distribution(runtimes, 15)

        #histogram = np.histogram(runtimes, bins=bins)

        #print(histogram)

        # the histogram of the data
        #n, bins, patches = plt.hist(runtimes, bins, normed=1, facecolor='green', alpha=0.75)

        print(distribution.counts)
        print(distribution.edges)

        #plt.hist(runtimes, distribution.edges)

        distribution.plot()

        #from scipy.signal import argrelextrema
        #from matplotlib.pyplot import *

        #np.random.seed()
        #x = np.random.random(50)
        #m = argrelextrema(x, np.greater)  # array of indexes of the locals maxima
        #y = [x[m] for i in m]
        #plot(x)
        #plot(m, y, 'rs')
        #show()


        #plt.show()

# -----------------------------------------------------------------
