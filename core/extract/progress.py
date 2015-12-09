#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module is used to extract simulation progress information from the simulation's log file
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# -----------------------------------------------------------------

class ProgressExtractor(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Attributes

        self.log_files = None
        self.table = None

    # -----------------------------------------------------------------

    @classmethod
    def open_table(cls, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Create a new ProgressExtractor instance
        extractor = cls()

        # Set the table attribute
        extractor.table = ascii.read(filepath)

        # Return the new ProgressExtractor instance
        return extractor

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None):

        """
        This function ...
        :return:
        """

        # Obtain the log files created by the simulation
        self.log_files = simulation.logfiles()

        # Determine whether the emission spectra calculation was performed using a staggered assignment scheme
        self.staggered = simulation.parameters().staggered()

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

        number_of_processes = None

        # Initialize lists for the columns
        process_list = []
        phase_list = []
        seconds_list = []
        progress_list = []

        # Loop over the log files again and fill the column lists
        for log_file in self.log_files:

            # Get the total number of processes
            if number_of_processes is None: number_of_processes = log_file.processes
            else: assert number_of_processes == log_file.processes

            # Get the process rank associated with this log file
            process = log_file.process

            stellar_start = None
            spectra_start = None
            dust_start = None

            first_spectra_phase = True
            last_dust_phase = False

            total_entries = None
            entries_per_process = None

            # Loop over the entries in the log file
            for i in range(len(log_file.contents)):

                # Get the description of the current simulation phase
                phase = log_file.contents["Phase"][i]

                # The log file entries corresponding to the stellar emission phase
                if phase == "stellar":

                    # The current log message
                    message = log_file.contents["Message"][i]

                    # If this is the log message that marks the very start of the stellar emission phase, record the associated time
                    if "photon packages for" in message: stellar_start = log_file.contents["Time"][i]

                    # If this is one of the log messages that log stellar emission progress
                    elif "Launched stellar emission photon packages" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - stellar_start).total_seconds()
                        seconds_list.append(seconds)

                        # Get the progress and add it to the list
                        progress = float(message.split("packages: ")[1].split("%"))
                        progress_list.append(progress)

                    elif "Finished the stellar emission phase" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - stellar_start).total_seconds()
                        seconds_list.append(seconds)

                        # Add 100% progress to the list
                        progress_list.append(100.0)

                # The log file entries corresponding to the stellar emission phase
                elif phase == "spectra" and first_spectra_phase:

                    # The current log message
                    message = log_file.contents["Message"][i]

                    # If this is the log message that marks the very start of the spectra calculation, record the associated time
                    if "Calculating dust emission spectra" in message: spectra_start = log_file.contents["Time"][i]

                    # If this log message states the total number of library entries that are used, record this number
                    elif "Library entries in use" in message:

                        # Get the total number of library entries in use and the number of entries per process
                        total_entries = int(message.split("use: ")[1].split(" out of")[0])
                        entries_per_process = total_entries / number_of_processes

                    elif "Calculating emission for" in message:

                        entry = float(message.split()[-1][:-3])

                        # Determine the progress
                        if self.staggered: fraction = entry / total_entries
                        else: fraction = (entry - process * entries_per_process) / entries_per_process

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - spectra_start).total_seconds()
                        seconds_list.append(seconds)

                        # Get the progress and add it to the list
                        progress = float(fraction*100.0)
                        progress_list.append(progress)

                    elif "Dust emission spectra calculated" in message:

                        # Add the process rank and phase entries
                        process_list.append(process)
                        phase_list.append(phase)

                        # Add the seconds entry
                        seconds = (log_file.contents["Time"][i] - spectra_start).total_seconds()
                        seconds_list.append(seconds)

                        # Add 100% progress to the list
                        progress_list.append(100.0)

                        # Indicate that the first spectra phase has already been processed (subsequent spectra phases can be ignored)
                        first_spectra_phase = False

                # The log file entries corresponding to the dust emission phase
                elif phase == "dust":

                    # The current log message
                    message = log_file.contents["Message"][i]

                    # Look for messages indicating whether this dust photon shooting phase corresponds to
                    # one of the dust self-absorption cycles or the actual dust emission phase
                    if "dust self-absorption cycle" in message: last_dust_phase = False
                    elif "Starting the dust emission phase" in message: last_dust_phase = True

                    # We only want to record the progress of the 'last' dust emission phase
                    if last_dust_phase:

                        # If this is the log message that marks the very start of the dust emission phase, record the associated time
                        if "photon packages for" in message: dust_start = log_file.contents["Time"][i]

                        # If this is one of the log messages that log dust emission progress
                        elif "Launched dust emission photon packages" in message:

                            # Add the process rank and phase entries
                            process_list.append(process)
                            phase_list.append(phase)

                            # Add the seconds entry
                            seconds = (log_file.contents["Time"][i] - dust_start).total_seconds()
                            seconds_list.append(seconds)

                            # Get the progress and add it to the list
                            progress = float(message.split("packages: ")[1].split("%"))
                            progress_list.append(progress)

                        elif "Finished the dust emission phase" in message:

                            # Add the process rank and phase entries
                            process_list.append(process)
                            phase_list.append(phase)

                            # Add the seconds entry
                            seconds = (log_file.contents["Time"][i] - dust_start).total_seconds()
                            seconds_list.append(seconds)

                            # Add 100% progress to the list
                            progress_list.append(100.0)

        # Create the table data structures
        names = ["Process rank", "Simulation phase", "Time", "Progress"]
        data = [process_list, phase_list, seconds_list, progress_list]

        # Create the table
        self.table = Table(data, names=names)
        self.table["Time"].unit = "s"
        self.table["Progress"].unit = "%"

    # -----------------------------------------------------------------

    def write(self, output_path):

        """
        This function ...
        :return:
        """

        # Write the table to file
        self.table.write(output_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the table to None
        self.table = None

# -----------------------------------------------------------------
