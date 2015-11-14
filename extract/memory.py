#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.memory Investigate the memory requirements of a particular SKIRT simulation

import subprocess
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np

from pts.skirtsimulation import SkirtSimulation
from pts.skifile import SkiFile

# This private helper function returns the datetime object corresponding to the time stamp in a line
def timestamp(line):

    date, time, dummy = line.split(None, 2)
    try:
        day, month, year = date.split('/')
    except ValueError:
        return None
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')
    
    return datetime(year=int(year), month=int(month), day=int(day),
                    hour=int(hour), minute=int(minute), second=int(second), microsecond=int(microsecond))

# Execute the skirtmem command and get the output
try:
    output = subprocess.check_output(["skirtmem", "PanSpiral/PanSpiral.ski", "-o", "PanSpiral"], stderr=subprocess.STDOUT)
except subprocess.CalledProcessError as e:
    output = e.output

# Create a simulation object from the simulation's output path
simulation = SkirtSimulation(prefix="PanSpiral.ski", outpath="PanSpiral")

# Create a ski file object
skifile = SkiFile("PanSpiral/PanSpiral.ski")

# Calculate the number of wavelengths and dust cells
Nlambda = skifile.nwavelengths()
Ncells = skifile.ncells()
Ndoubles = Nlambda * Ncells

# Get the list of logfiles for this simulation
logfilepaths = simulation.logfilepaths()

# Create a table to contain the timeline data from the log files
nprocs = len(logfilepaths)
data = [[] for _ in range(nprocs)]

# Keep track of the earliest time that is found in any of the log files
T0 = datetime.now()

time_list = []
memory_list = []

# Loop over all log files
for logfile in logfilepaths:

    # Get the rank of the process that created this log file (for the process with rank zero finding this rank will
    # not succeed, therefore use the try-except statements
    processrank = 0
    try: processrank = int(logfile[-7:-4])
    except ValueError: pass
    
    # Set the default value for the number of processes
    nprocs = 1
    
    # Iterate over each line in this log file
    for line in open(logfile):

        # Get the date and time information from this line in the log file
        time = timestamp(line)

        # Replace the start time with the time at which this simulation started, if earlier
        if "Starting simulation" in line:

            if time < T0: T0 = time

            # Get the number of processes with which the simulation was run
            if "with" in line: nprocs = int(line.split(' with ')[1].split()[0])

        memory = float(line.split(" (")[1].split(" GB)")[0])
        
        time_list.append(time)
        memory_list.append(memory)

seconds_list = []
for time in time_list:
    
    seconds = (time - T0).total_seconds()
    seconds_list.append(seconds)

plt.plot(seconds_list, memory_list)

times = []
deltas = []

# Iterate over each line in the output
for line in output.splitlines():
    
    # Get the time
    time = timestamp(line)
    
    # Skip the line
    if time is None: continue
    
    # Calculate the difference in number of seconds
    seconds = (time - T0).total_seconds()
    
    if line[26] == "+":
        
        memory = float(line.split("+")[1].split("GB")[0])
        times.append(seconds)
        deltas.append(memory)
        
    elif line[26] == "-":
        
        memory = -float(line.split("-", 1)[1].split("GB")[0])
        times.append(seconds)
        deltas.append(memory)

totals = np.cumsum(deltas)

#for i in range(len(times)):
#    print times[i], totals[i]

#plt.step(times, totals)
plt.fill_between(times, totals, color='green')
#plt.bar(times, totals, color='r')
plt.show()
