#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.run Run a SKIRT or FitSKIRT simulation with the best performance, based on the current load of the
#  system.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import argparse
import psutil
import multiprocessing
import numpy as np

# Import the relevant PTS modules
from pts.log import Log
from pts.skirtexec import SkirtExec

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument('skifile', type=str, help='the name of the ski file')
parser.add_argument('inpath', type=str, help='the simulation input path')
parser.add_argument('outpath', type=str, help='the simulation output path')
parser.add_argument('--brief', action='store_true', help='add this option to enable brief console logging')
parser.add_argument('--verbose', action='store_true', help='add this option to enable verbose logging mode for SKIRT')

# Parse the command line arguments
args = parser.parse_args()
skifile = args.skifile
inpath = args.inpath
outpath = args.outpath
brief = args.brief
verbose = args.verbose

# -----------------------------------------------------------------

# Create a logger
log = Log()

# Determine the full path to the ski file
skifilepath = os.path.abspath(skifile)

# Determine the full path to the input and output directories
inputpath = os.path.abspath(inpath)
outputpath = os.path.abspath(outpath)

# Get the total number of processors on this system
total = multiprocessing.cpu_count()

# Get the load of the different processors
load = np.array(psutil.cpu_percent(percpu=True))/100.0

# Calculate the the number of full processors (where 2 processors with loads x and y contribute as a processor with load x+y)
full = load.sum()

# Get the number of free processors
free = total - full

# Get the currently available virtual memory (in gigabytes)
memory = psutil.virtual_memory().available / 1e9

# Inform the user
log.info("The number of currently available processors on this system is " + str(free))
log.info("The amount of currently available memory on this system is " + str(memory) + " gigabytes")

# Find the number of wavelengths and the number of dust cells from the ski file
wavelengths = 0
xpoints = 0
ypoints = 0
zpoints = 0
with open(skifilepath) as skifile:

    for line in skifile:

        if "minWavelength" in line and "maxWavelength" in line and "points" in line:

            wavelengths = int(line.split('points="')[1].split('"')[0])

        if "DustGridStructure" in line and "writeGrid" in line:

            xpoints = int(line.split('pointsX="')[1].split('"')[0])
            ypoints = int(line.split('pointsY="')[1].split('"')[0])
            zpoints = int(line.split('pointsZ="')[1].split('"')[0])

dustcells = xpoints * ypoints * zpoints

# Calculate the amount of required memory for this simulation (in gigabytes)
required = float(wavelengths) * float(dustcells) * 8 / 1e9

# Inform the user
log.info("The number of wavelengths for this simulation is " + str(wavelengths))
log.info("The number of dust cells for this simulation is " + str(dustcells))
log.info("The estimated memory requirement for this simulation is " + str(required) + " gigabytes")

# Calculate the maximum number of MPI processes
processes = int(memory/required)

if processes < 1:

    log.error("Not enough memory available to run this simulation")
    exit()

# Calculate the number of threads per process
threads = int(free / processes)

# Inform the user
log.info("The number of processes that will be used is " + str(processes))
log.info("The number of threads per process that will be used is " + str(threads))

# Create the SKIRT execution context
skirt = SkirtExec(log=log)

# Run the simulation
simulation = skirt.execute(skipattern=skifilepath, inpath=inputpath, outpath=outputpath, threads=threads, processes=processes, brief=brief, verbose=verbose)[0]

# Check whether the simulation finished
if simulation.status() != "Finished": raise ValueError("Simulation " + simulation.status())

# -----------------------------------------------------------------
