#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import standard modules
import argparse

# -----------------------------------------------------------------

## Extract the timings from the simulation log file
def extract(logfilepath, processes=1, threads=1, resultsfilepath=""):

    # Open the results file for appending the timings for this simulation
    if resultsfilepath:
        resultfile = open(resultsfilepath, 'a')

    for line in open(logfilepath):

        if 'Finished setup in' in line:
            setuptime = float(line.split(' in ')[1].split()[0])

        elif 'Finished the stellar emission phase in' in line:
            stellartime = float(line.split(' in ')[1].split()[0])

        elif 'Finished writing results in' in line:
            writingtime = float(line.split(' in ')[1].split()[0])

        elif 'Finished simulation' in line:
            simulationtime = float(line.split(' in ')[1].split()[0])

    if resultsfilepath:
        resultfile.write(str(processes) + ' ' + str(threads) + ' ' + str(processes*threads) + ' ' + str(setuptime) + ' ' + str(stellartime) + ' ' + str(writingtime) + ' ' + str(simulationtime) + '\n')
    else:
        return [setuptime, stellartime, writingtime, simulationtime]

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

