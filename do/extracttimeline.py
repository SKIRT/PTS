#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import standard modules
import argparse
from datetime import datetime

# Import relevant PTS modules
from pts.skirtsimulation import SkirtSimulation

# -----------------------------------------------------------------

# Define the indices used to identify the different simulation phases
phaseindices = {'setup': 0, 'stellar': 1, 'comm': 2, 'spectra': 3, 'dust': 4, 'writing': 5, 'waiting': 6}

# -----------------------------------------------------------------

## This function extracts the timeline of the different stages in a SKIRT simulation for the different processes
#  from the log files created by this simulation.
#
def extract(skifilename, outputpath, timelinefilepath=""):

    # Create a simulation object from the simulation's output path
    simulation = SkirtSimulation(prefix=skifilename, outpath=outputpath)

    # Get the list of logfiles for this simulation
    logfilepaths = simulation.logfilepaths()

    # Create a table to contain the timeline data from the log files
    nprocs = len(logfilepaths)
    data = [[] for _ in range(nprocs)]

    # Keep track of the earliest time that is found in any of the log files
    T0 = datetime.now()

    # Extract the profile
    for logfile in logfilepaths:

        # Get the rank of the process that created this log file (for the process with rank zero finding this rank will
        # not succeed, therefore use the try-except statements
        processrank = 0
        try: processrank = int(logfile[-7:-4])
        except ValueError: pass

        # Iterate over each line in this log file
        for line in open(logfile):

            # Get the date and time information from this line in the log file
            time = _timestamp(line)

            # Replace the start time with the time at which this simulation started, if earlier
            if "Starting simulation" in line:

                if time < T0: T0 = time

            # Check whether this line indicates the start of the setup
            elif "Starting setup" in line:

                data[processrank].append({'phase': phaseindices['setup'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the first part of the setup, and the start of the
            # subsequent waiting phase
            elif "Waiting for other processes to finish the calculation of the dust cell densities" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the waiting phase after the first part of the setup,
            # and the start of the communication of the dust densities
            elif "Starting communication of the dust densities" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['comm'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the communication of the dust densities, and the start
            # of the second part of the setup
            elif "Finished communication of the dust densities" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['setup'], 'start': time, 'end': None})

            # Check whether this line indicates the end of (the second part of) the setup, and the start of the
            # waiting phase after the setup
            elif "Waiting for other processes to finish the setup" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the waiting phase after the setup
            elif "Finished setup" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the start of the stellar emission phase
            elif "Starting the stellar emission phase" in line:

                data[processrank].append({'phase': phaseindices['stellar'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the stellar emission phase, and the start of the
            # waiting phase after the stellar emission phase
            elif "Waiting for other processes to finish the stellar emission phase" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the waiting phase after the stellar emission phase
            elif "Finished the stellar emission phase" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the start of the communication of the absorbed luminosities
            elif "Starting communication of the absorbed luminosities" in line:

                data[processrank].append({'phase': phaseindices['comm'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the communication of the absorbed luminosities
            elif "Finished communication of the absorbed luminosities" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the start of the dust emission spectra calculation
            elif "Library entries in use" in line:

                data[processrank].append({'phase': phaseindices['spectra'], 'start': time, 'end': None})


            # Check whether this line indicates the end of the dust emission spectra calculation, and the start
            # of the waiting phase after the dust emission spectra calculation
            elif "Waiting for other processes to finish the emission spectra calculation" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the waiting phase after the dust emission spectra
            # calculation, and the start of the communication of the emission luminosities
            elif "Starting communication of the dust emission spectra" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['comm'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the communication of the emission luminosities
            elif "Finished communication of the dust emission spectra" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the beginning of the dust emission phase or a dust self-absorption cycle
            elif "Dust emission spectra calculated." in line:

                data[processrank].append({'phase': phaseindices['dust'], 'start': time, 'end': None})

            # Check whether this line indicates the end of a dust self-absorption cycle, and the start of the
            # subsequent waiting phase
            elif "Waiting for other processes to finish this self-absorption cycle" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the waiting phase after a dust self-absorption cycle
            elif "Finished the" in line and "-stage dust self-absorption cycle" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the end of the dust emission phase, and the start of the waiting phase
            # after the dust emission phase
            elif "Waiting for other processes to finish the dust emission phase" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

                data[processrank].append({'phase': phaseindices['waiting'], 'start': time, 'end': None})

            # Check whether this line indicates the end of waiting phase after the dust emission phase
            elif "Finished the dust emission phase" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

            # Check whether this line indicates the start of the writing phase
            elif "Starting writing results" in line:

                data[processrank].append({'phase': phaseindices['writing'], 'start': time, 'end': None})

            # Check whether this line indicates the end of the writing phase
            elif "Finished writing results" in line:

                entries = len(data[processrank])
                data[processrank][entries-1]['end'] = time

    # Check whether a timeline file is specified
    if timelinefilepath:

        # Open the timeline file
        timelinefile = open(timelinefilepath, 'a')

        # Write a header
        timelinefile.write("# Times are in seconds, relative to T0 = " + T0.strftime("%Y-%m-%d %H:%M:%S.%f") + "\n")
        timelinefile.write("# Column 1: Process rank\n")
        timelinefile.write("# Column 2: Simulation phase (0=setup, 1=stellar emission, 2=communication, 3=spectra, 4=dust emission, 5=writing)\n")
        timelinefile.write("# Column 3: Start time (s)\n")
        timelinefile.write("# Column 4: End time (s)\n")

        # Loop over all the different ranks
        for rank, entries in enumerate(data):

            # Loop over all the different entries (different Timespans)
            for entry in entries:

                # If the end of this entry was not recorded, skip it
                if entry['end'] is None: continue

                phase_id = str(entry['phase'])
                phase_start = str((entry['start'] - T0).total_seconds())
                phase_end = str((entry['end'] - T0).total_seconds())
                timelinefile.write(str(rank) + " " + phase_id + " " + phase_start + " " + phase_end + "\n")

        # Close the progress file
        timelinefile.close()

    else:

        # If no timeline file is given, return the table of extracted timeline information
        return data

# -----------------------------------------------------------------

# This private helper function returns the datetime object corresponding to the time stamp in a line
def _timestamp(line):

    date, time, dummy = line.split(None, 2)
    day, month, year = date.split('/')
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')
    return datetime(year=int(year), month=int(month), day=int(day),
                    hour=int(hour), minute=int(minute), second=int(second), microsecond=int(microsecond))

# -----------------------------------------------------------------

# Execute the statements below if this script is run from the command line
if __name__ == "__main__":

    # Create the command-line parser and a set of subparsers
    parser = argparse.ArgumentParser()
    parser.add_argument('skifilepath', type=str, help='the path to the ski file of the simulation')
    parser.add_argument('outputpath', type=str, help='the path to the simulation output directory')
    parser.add_argument('timelinefilepath', type=str, help='the path to the timeline file')

    # Parse the command line arguments
    args = parser.parse_args()

    # Set the path to the simulation's ski file
    skifilepath = args.skifilepath

    # Set the path of the simulation output directory
    outputpath = args.outputpath

    # Set the path to the timeline file
    timelinefilepath = args.timelinefilepath

    # Extract the timeline
    extract(skifilepath, outputpath, timelinefilepath)

