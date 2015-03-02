#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import standard modules
import argparse

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

    # Set all runtimes to zero; for stages that are not encountered in the log file the runtime will stay zero.
    setuptime = 0.0
    stellartime = 0.0
    dustselfabstime = 0.0
    dustemissiontime = 0.0
    writingtime = 0.0
    simulationtime = 0.0

    # Read each line of the specified log file to search for the end of simulation stages
    for line in open(logfilepath):

        if 'Finished setup in' in line:
            setuptime = float(line.split(' in ')[1].split()[0])

        elif 'Finished the stellar emission phase in' in line:
            stellartime = float(line.split(' in ')[1].split()[0])

        elif 'Finished the dust self-absorption phase in' in line:
            dustselfabstime = float(line.split(' in ')[1].split()[0])

        elif 'Finished the dust emission phase in' in line:
            dustemissiontime = float(line.split(' in ')[1].split()[0])

        elif 'Finished writing results in' in line:
            writingtime = float(line.split(' in ')[1].split()[0])

        elif 'Finished simulation' in line:
            simulationtime = float(line.split(' in ')[1].split()[0])

    # Check whether a results file is specified
    if resultsfilepath:

        # Open the results file
        resultfile = open(resultsfilepath, 'a')

        # Add a line containing the runtimes for this number of processes and threads per process
        resultfile.write(str(processes) + ' ' + str(threads) + ' ' + str(processes*threads) + ' ' + str(setuptime)
                         + ' ' + str(stellartime) + ' ' + str(dustselfabstime) + ' ' + str(dustemissiontime)
                         + ' ' + str(writingtime) + ' ' + str(simulationtime) + '\n')

    else:

        # If no results file is given, return the extracted runtimes
        return [setuptime, stellartime, dustselfabstime, dustemissiontime, writingtime, simulationtime]

# -----------------------------------------------------------------

# Execute the statements below if this script is run from the command line
if __name__ == "__main__":

    # Create the command-line parser and a set of subparsers
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

