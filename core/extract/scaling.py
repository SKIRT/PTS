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

# -----------------------------------------------------------------

class ScalingExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def run(self, simulation, output_path):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

## This function extracts the timings from the simulation log file and either returns them or writes them to file.
#  It takes the following arguments:
#
#  - logfilepath: the path to the simulation's log file
#  - processes: the number of processes with which the simulation was run
#  - threasd: the number of threads per process with which the simulation was run
#  - resultsfilepath: the path to a file that to where the runtimes should be written
#
def extract(logfilepath, processes=1, threads=1, resultsfilepath=""):

    # Create a dictionary with the runtimes for each simulation stage set to zero. The runtime will stay zero
    # for stages that are not encountered in the specified log file.
    keys = ('setup', 'stellar', 'dustselfabs', 'dustem', 'writing', 'total')
    runtimes = dict.fromkeys(keys, 0.0)

    # Read each line of the specified log file to search for the end of simulation stages
    for line in open(logfilepath):

        if 'Finished setup in' in line:

            # Add the runtime of the setup to the dictionary
            runtimes['setup'] = float(line.split(' in ')[1].split()[0])

        elif 'Finished the stellar emission phase in' in line:

            # Add the runtime of the stellar emission phase to the dictionary
            runtimes['stellar'] = float(line.split(' in ')[1].split()[0])

        elif 'Finished the dust self-absorption phase in' in line:

            # Add the runtime of the dust self-absorption phase to the dictionary
            runtimes['dustselfabs'] = float(line.split(' in ')[1].split()[0])

        elif 'Finished the dust emission phase in' in line:

            # Add the runtime of the dust emission phase to the dictionary
            runtimes['dustem'] = float(line.split(' in ')[1].split()[0])

        elif 'Finished writing results in' in line:

            # Add the runtime of the writing phase to the dictionary
            runtimes['writing'] = float(line.split(' in ')[1].split()[0])

        elif 'Finished simulation' in line:

            # Add the total runtime of the simulation to the dictionary
            runtimes['total'] = float(line.split(' in ')[1].split()[0])

    # Check whether a results file is specified
    if resultsfilepath:

        # Open the results file
        resultfile = open(resultsfilepath, 'a')

        # Add a line containing the runtimes for this number of processes and threads per process
        resultfile.write(str(processes) + ' ' + str(threads) + ' ' + str(processes*threads) + ' '
                         + str(runtimes['setup']) + ' ' + str(runtimes['stellar']) + ' ' + str(runtimes['dustselfabs'])
                         + ' ' + str(runtimes['dustem']) + ' ' + str(runtimes['writing']) + ' ' + str(runtimes['total'])
                         + '\n')

    else:

        # If no results file is given, return the dictionary of the extracted runtimes
        return runtimes

# -----------------------------------------------------------------

# Execute the statements below if this script is run from the command line
if __name__ == "__main__":

    # Create the command-line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('logfilepath', type=str, help='the path to the log file')
    parser.add_argument('processes', type=int, help='the number of processes used for the simulation')
    parser.add_argument('threads', type=int, help='the number of threads per process used for the simulation')
    parser.add_argument('resultsfilepath', type=str, help='the path to the results file')

    # Parse the command line arguments
    args = parser.parse_args()

    # Set the path of the simulation log file
    logfilepath = args.logfilepath

    # The number of processes and threads for this run
    processes = args.processes
    threads = args.threads

    # Set the path to the results file
    resultsfilepath = args.resultsfilepath

    # Extract the timings
    extract(logfilepath, processes, threads, resultsfilepath)

