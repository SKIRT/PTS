#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package extract.timeline

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import argparse
from datetime import datetime

# Import astronomical modules
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..simulation.simulation import SkirtSimulation

# -----------------------------------------------------------------

class TimeLineExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Attributes

        self.table = None

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None):

        """
        This function ...
        :return:
        """

        # Obtain the log files created by the simulation
        self.log_files = simulation.logfiles()

        # Perform the extraction
        self.extract()

        # Write the results
        if output_path is not None: self.write(output_path)

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Initialize lists for the columns
        process_list = []
        phase_list = []
        start_list = []
        end_list = []

        # Loop over all log files to determine the earliest recorded time
        t0 = datetime.now()
        for log_file in self.log_files:
            if log_file.t0 < t0: t0 = log_file.t0

        # Loop over the log files again and fill the column lists
        for log_file in self.log_files:

            # Get the process rank associated with this log file
            process = log_file.process

            # Keep track of the current phase while looping over the log file entries
            current_phase = log_file.contents["Phase"][0]
            process_list.append(process)
            phase_list.append(current_phase)
            start_list.append(log_file.t0)

            # Loop over all log file entries
            for j in range(len(log_file.contents)):

                # Get the description of the current simulation phase
                phase = log_file.contents["Phase"][j]

                # If a new phase is entered
                if phase != current_phase:

                    # Determine the current time
                    seconds = log_file.contents["Time"][j]

                    # Mark the end of the previous phase
                    end_list.append(seconds)

                    # Mark the start of the current phase
                    process_list.append(process)
                    phase_list.append(phase)
                    start_list.append(seconds)

                    # Update the current phase
                    current_phase = phase

        # Create the table data structures
        names = ['Process rank', 'Simulation phase', 'Start time', 'End time']
        data = [process_list, phase_list, start_list, end_list]

        # Create the table
        self.table = Table(data, names=names)
        self.table["Start time"].unit = "s"
        self.table["End time"].unit = "GB"

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :param output_path:
        :return:
        """

        # Write the table to file
        self.table.write(output_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def duration(self, phase, single=False):

        """
        This function ...
        :param phase:
        :return:
        """

        # Keep track of the total amount of time spent in the specified phase
        total = 0.0

        assert self.table["Process"][0] == 0

        # Loop over the table rows
        for i in range(len(self.table)):

            # Only add the contributions from the root process
            if self.table["Process"][i] > 0: break

            # Check whether the current entry corresponds to the desired phase
            if self.table["Simulation phase"][i] == phase:

                # Get the start and end time for the phase
                start = self.table["Start time"][i]
                end = self.table["End time"][i]

                # Calculate the time duration for this phase, returning it if single=True, otherwise add it to the total
                if single: return end - start
                else: total += end - start

        # Return the total amount of time spent in the specified phase
        return total

    # -----------------------------------------------------------------

    def duration_without(self, phases):

        """
        This function ...
        :param phase:
        :return:
        """

        # Create a list of phases is only one is given
        if isinstance(phases, basestring): phases = [phases]

        # Keep track of the total amount of time spent in phases other than the specified phase
        total = 0.0

        assert self.table["Process"][0] == 0

        # Loop over the table rows
        for i in range(len(self.table)):

            # Only add the contributions from the root process
            if self.table["Process"][i] > 0: break

            # Check whether the current entry corresponds to a phase different from the specified phase
            if self.table["Simulation phase"][i] not in phases:

                # Get the start and end time for this phase
                start = self.table["Start time"][i]
                end = self.table["End time"][i]

                # Add the duration to the total
                total += end - start

        # Return the total amount of time spent in phases other than the specified phase
        return total

    # -----------------------------------------------------------------

    @property
    def setup(self):

        """
        This function ...
        :return:
        """

        return self.duration("setup")

    # -----------------------------------------------------------------

    @property
    def stellar(self):

        """
        This function ...
        :return:
        """

        return self.duration("stellar", single=True)

    # -----------------------------------------------------------------

    @property
    def spectra(self):

        """
        This function ...
        :return:
        """

        return self.duration("spectra")

    # -----------------------------------------------------------------

    @property
    def dust(self):

        """
        This function ...
        :return:
        """

        return self.duration("dust")

    # -----------------------------------------------------------------

    @property
    def dustem(self):

        """
        This function ...
        :return:
        """

        # Loop over the table rows in opposite direction so that we know the first dust photon shooting phase is the
        # final dust emission phase
        for i in reversed(range(len(self.table))):

            # Only add the contributions from the root process
            if self.table["Process"][i] > 0: break

            # Check whether the current entry corresponds to the desired phase
            if self.table["Simulation phase"][i] == "dust":

                # Get the start and end time for the phase
                start = self.table["Start time"][i]
                end = self.table["End time"][i]

                # Calculate the time duration for the dust emission phase (only the shooting part)
                return end - start

    # -----------------------------------------------------------------

    @property
    def writing(self):

        """
        This function ....
        :return:
        """

        return self.duration("write")

    # -----------------------------------------------------------------

    @property
    def communication(self):

        """
        This function ...
        :return:
        """

        return self.duration("comm")

    # -----------------------------------------------------------------

    @property
    def waiting(self):

        """
        This function ...
        :return:
        """

        return self.duration("wait")

    # -----------------------------------------------------------------

    @property
    def other(self):

        """
        This function ...
        :return:
        """

        return self.duration(None)

    # -----------------------------------------------------------------

    @property
    def serial(self):

        """
        This function ...
        :return:
        """

        return self.setup + self.writing + self.other

    # -----------------------------------------------------------------

    @property
    def parallel(self):

        """
        This function ...
        :return:
        """

        return self.stellar + self.spectra + self.dust

    # -----------------------------------------------------------------

    @property
    def overhead(self):

        """
        This function ...
        :return:
        """

        return self.communication + self.waiting

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

        # Set the default value for the number of processes
        nprocs = 1

        # Iterate over each line in this log file
        for line in open(logfile):

            # Get the date and time information from this line in the log file
            time = _timestamp(line)

            # Replace the start time with the time at which this simulation started, if earlier
            if "Starting simulation" in line:

                if time < T0: T0 = time

                # Get the number of processes with which the simulation was run
                if "with" in line: nprocs = int(line.split(' with ')[1].split()[0])

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

                # If the number of processes is one, this line also marks the end of the calculation of the dust
                # emission spectra
                if nprocs == 1:

                    entries = len(data[processrank])
                    data[processrank][entries-1]['end'] = time

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

                # If the end of this entry was not recorded, take the start of the next entry
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

# *****************************************************************

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

# *****************************************************************
