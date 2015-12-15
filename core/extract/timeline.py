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
from datetime import datetime

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

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

        # The list of log files created by the simulation
        self.log_files = None

        # The table containing the timeline information
        self.table = None

    # -----------------------------------------------------------------

    @classmethod
    def open_table(cls, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Create a new TimeLineExtractor instance
        extractor = cls()

        # Set the table attribute
        fill_values = ('--', '0', 'Simulation phase')
        extractor.table = ascii.read(filepath, fill_values=fill_values)

        # Return the new TimeLineExtractor instance
        return extractor

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
        t_0 = datetime.now()
        for log_file in self.log_files:
            if log_file.t_0 < t_0: t_0 = log_file.t_0

        # Loop over the log files again and fill the column lists
        for log_file in self.log_files:

            # Get the process rank associated with this log file
            process = log_file.process

            # Keep track of the current phase while looping over the log file entries
            current_phase = log_file.contents["Phase"][0]
            process_list.append(process)
            phase_list.append(current_phase)
            start_list.append((log_file.t_0 - t_0).total_seconds())

            # Loop over all log file entries
            for j in range(len(log_file.contents)):

                # Get the description of the current simulation phase
                phase = log_file.contents["Phase"][j]

                # If a new phase is entered
                if phase != current_phase:

                    # Determine the current time
                    seconds = (log_file.contents["Time"][j] - t_0).total_seconds()

                    # Mark the end of the previous phase
                    end_list.append(seconds)

                    # Mark the start of the current phase
                    process_list.append(process)
                    phase_list.append(phase)
                    start_list.append(seconds)

                    # Update the current phase
                    current_phase = phase

            # Add the last recorded time in the log file as the end of the last simulation phase
            end_list.append((log_file.t_last - t_0).total_seconds())

        # Create the table data structures
        names = ["Process rank", "Simulation phase", "Start time", "End time"]
        data = [process_list, phase_list, start_list, end_list]

        # Create the table
        self.table = Table(data, names=names, masked=True)
        self.table["Start time"].unit = "s"
        self.table["End time"].unit = "s"

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

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the table to None
        self.table = None

    # -----------------------------------------------------------------

    @property
    def processes(self):

        """
        This function ...
        :return:
        """

        return max(self.table["Process rank"]) + 1

    # -----------------------------------------------------------------

    def duration(self, phase, single=False):

        """
        This function ...
        :param phase:
        :return:
        """

        # Keep track of the total amount of time spent in the specified phase
        total = 0.0

        assert self.table["Process rank"][0] == 0

        # Loop over the table rows
        for i in range(len(self.table)):

            # Only add the contributions from the root process
            if self.table["Process rank"][i] > 0: break

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

        # If no phases are given, set an empty list
        if phases is None: phases = []

        # Create a list of phases is only one is given
        if isinstance(phases, basestring): phases = [phases]

        # Keep track of the total amount of time spent in phases other than the specified phase
        total = 0.0

        assert self.table["Process rank"][0] == 0

        # Loop over the table rows
        for i in range(len(self.table)):

            # Only add the contributions from the root process
            if self.table["Process rank"][i] > 0: break

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
    def total(self):

        """
        This function ...
        :return:
        """

        return self.duration_without(None)

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
            if self.table["Process rank"][i] > 0: break

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
