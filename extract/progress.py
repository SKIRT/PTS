#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse
from datetime import datetime

# Import relevant PTS modules
from ..core.simulation import SkirtSimulation
from ..core import archive as arch

# *****************************************************************

class ProgressExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        pass

    # *****************************************************************

    def run(self, simulation, output_path):

        """
        This function ...
        :return:
        """

        pass

# *****************************************************************

phaseindices = {'stellar': 0, 'spectra': 1, 'dust': 2}

# *****************************************************************

## This function extracts the progress from the simulation log files and writes them to file.
#  It takes the following arguments:
#
#  - simulationpath: the path to the simulation's log file
#  - progressfilepath: the path to a file that to where the runtimes should be written
#
def extract(skifilename, outputpath, progressfilepath):

    # Open the progress file
    progressfile = open(progressfilepath, 'a')

    # Create a simulation object from the simulation's output path
    simulation = SkirtSimulation(prefix=skifilename, outpath=outputpath)

    # Get the list of logfiles for this simulation
    logfilepaths = simulation.logfilepaths()

    # Extract the progress of shooting stellar emission photons
    _extractphotonprogress(logfilepaths, progressfile, 'stellar')

    # Extract the progress of calculating the dust emission spectra
    _extractspectraprogress(logfilepaths, progressfile, simulation.staggered())

    # Extract the progress of shooting dust emission photons
    _extractphotonprogress(logfilepaths, progressfile, 'dust')

    # Close the progress file
    progressfile.close()

## This function
def _extractphotonprogress(logfiles, progressfile, phase):

    # For the log file of each process
    for logfile in logfiles:

        # Initialize empty lists to contain the runtimes and corresponding percentages
        times = []
        percentages = []

        # Get the rank of the process that created this log file (for the process with rank zero finding this rank will
        # not succeed, therefore use the try-except statements
        processrank = 0
        try: processrank = int(logfile[-7:-4])
        except ValueError: pass

        # We use the 'triggered' variable to indicate whether we are currently in the part of the log file
        # that is of interest for the progress of launching stellar/dust emission photons
        triggered = False

        # For each line in this log file
        for line in arch.opentext(logfile):

            # Check if this line signals the start of the stellar/dust emission phase
            if not triggered and "Starting the {} emission phase".format(phase) in line: triggered = True

            # If the emission phase has been triggered, check whether the shooting of photon packages has started
            if triggered and " photon packages for " in line: startline = line

            # If the emission phase has been triggered, acquire the progress of shooting photon packages
            if triggered and " photon packages: " in line and "Launched" in line:

                percentages.append( float(line.split()[-1][:-1]) )
                times.append(_timelapse(startline,line))

            # If the emission phase has been triggered, look for the line that signals the end of this phase
            if triggered and " Finished the {} emission phase".format(phase) in line: break

        # Add a line to the progress file for each (runtime, percentage) pair
        for runtime, percentage in zip(times, percentages):

            # Write the runtime and percentage to the progress file, with the process rank in the first column and
            # the simulation phase in the second column
            progressfile.write(str(phaseindices[phase]) + ' ' + str(processrank) + ' ' + str(runtime) + ' ' + str(percentage) + '\n')

## This function
def _extractspectraprogress(logfiles, progressfile, staggered=False):

    # Determine the number of processes used for the simulation
    numprocesses = len(logfiles)

    # For the log file of each process
    for logfile in logfiles:

        # Initialize empty lists to contain the runtimes and corresponding percentages
        times = []
        percentages = []

        # Get the rank of the process that created this log file (for the process with rank zero finding this rank will
        # not succeed, therefore use the try-except statements
        processrank = 0
        try: processrank = int(logfile[-7:-4])
        except ValueError: pass

        # We use the 'triggered' variable to indicate whether we are currently in the part of the log file
        # that is of interest for the progress of the dust emission spectra calculation
        triggered = False

        # For each line in this log file
        for line in arch.opentext(logfile):

            # Check if this line signals the start of the dust emission phase
            if not triggered and "Starting the dust emission phase" in line: triggered = True

            # If the dust emission phase has been triggered, check whether we are now at the calculation of the spectra
            if triggered and " Calculating dust emission spectra" in line: startline = line

            # If the dust emission phase has been triggered, look for the total number of library entries
            if triggered and " Library entries in use: " in line:

                totalentries = float(line.split()[-1])
                entriesperprocess = totalentries/numprocesses

            # If the dust emission phase has been triggered, acquire the progress of the spectra calculation
            if triggered and " Calculating emission for " in line:

                entry = float(line.split()[-1][:-3])

                if staggered: fraction = entry/totalentries
                else: fraction = (entry - processrank*entriesperprocess)/entriesperprocess

                percentages.append(100*fraction)
                times.append(_timelapse(startline,line))

            # If the dust emission phase has been triggered, look for the line that signals the end of the dust emission phase
            if triggered and " Dust emission spectra calculated" in line: break

         # Add a line to the progress file for each (runtime, percentage) pair
        for runtime, percentage in zip(times, percentages):

            # Write the runtime and percentage to the progress file, with the process rank in the first column and
            # the simulation phase in the second column
            progressfile.write(str(phaseindices['spectra']) + ' ' + str(processrank) + ' ' + str(runtime) + ' ' + str(percentage) + '\n')

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
    parser.add_argument('skifilepath', type=str, help='the path to the ski file of the simulation')
    parser.add_argument('outputpath', type=str, help='the path to the simulation output directory')
    parser.add_argument('progressfilepath', type=str, help='the path to the progress file')

    # Parse the command line arguments
    args = parser.parse_args()

    # Set the path to the simulation's ski file
    skifilepath = args.skifilepath

    # Set the path of the simulation output directory
    outputpath = args.outputpath

    # Set the path to the results file
    progressfilepath = args.progressfilepath

    # Extract the timings
    extract(skifilepath, outputpath, progressfilepath)
