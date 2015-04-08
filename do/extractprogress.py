#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import standard modules
import argparse
from datetime import datetime

# Import relevant PTS modules
from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------

## This function extracts the progress from the simulation log files and writes them to file.
#  It takes the following arguments:
#
#  - simulationpath: the path to the simulation's log file
#  - progressfilepath: the path to a file that to where the runtimes should be written
#
def extract(simulationpath, progressfilepath):

    # Temporary
    phase = "stellar"

    # Open the progress file
    progressfile = open(progressfilepath, 'a')

    # Create a simulation object from the simulation's output path
    simulation = SkirtSimulation(outpath=simulationpath)

    # For the log file of each process
    for logfilepath in simulation.logfilepaths():

        # Initialize empty lists to contain the runtimes and corresponding percentages
        times = []
        percentages = []

        triggered = False

        # For each line in this log file
        for line in open(logfilepath):

            if not triggered and "Starting the {} emission phase".format(phase) in line:

                triggered = True
                process = line.split()[2][1:5] if "[P" in line else "P000"

            if triggered and " photon packages for " in line:

                startline = line

            if triggered and " photon packages: " in line and "Launched" in line:

                percentages.append( float(line.split()[-1][:-1]) )
                times.append(_timelapse(startline,line))

            if triggered and " Finished " in line: break

        # Add a line to the progress file for each (runtime, percentage) pair
        for runtime, percentage in zip(times, percentages):

            # Write the runtime and percentage to the progress file, with the process rank in the first column
            progressfile.write(str(process) + ' ' + str(runtime) + ' ' + str(percentage) + '\n')

    # Close the progress file
    progressfile.close()

# -----------------------------------------------------------------

# This private helper function returns the datetime object corresponding to the time stamp in a line
def _timestamp(line):

    date, time, dummy = line.split(None, 2)
    day, month, year = date.split('/')
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')
    return datetime(year=int(year), month=int(month), day=int(day),
                    hour=int(hour), minute=int(minute), second=int(second), microsecond=int(microsecond))

# This private helper function returns the difference in seconds between the time stamps in the two lines
def _timelapse(line1, line2):

    return (_timestamp(line2)-_timestamp(line1)).total_seconds()

# -----------------------------------------------------------------

# Execute the statements below if this script is run from the command line
if __name__ == "__main__":

    # Create the command-line parser and a set of subparsers
    parser = argparse.ArgumentParser()
    parser.add_argument('simulationpath', type=str, help='the path to the simulation output directory')
    parser.add_argument('progressfilepath', type=str, help='the path to the progress file')

    # Parse the command line arguments
    args = parser.parse_args()

    # Set the path of the simulation output directory
    simulationpath = args.simulationpath

    # Set the path to the results file
    progressfilepath = args.progressfilepath

    # Extract the timings
    extract(simulationpath, progressfilepath)

