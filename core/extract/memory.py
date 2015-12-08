#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package extract.memory Extract memory information from SKIRT or SkirtMemory output

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from datetime import datetime

# Import astronomical modules
from astropy.table import Table
from astropy.io import ascii

# -----------------------------------------------------------------

class MemoryExtractor(object):

    """
    This function ...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        ## Attributes

        self.table = None

    # -----------------------------------------------------------------

    @classmethod
    def open_table(cls, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        # Create a new MemoryExtractor instance
        extractor = cls()

        # Set the table attribute
        fill_values = [('--', '0', 'Simulation phase'), ('--', '0', 'Array (de)allocation'), ('--', '0', 'Array ID')]
        extractor.table = ascii.read(filepath, fill_values=fill_values)

        # Return the new MemoryExtractor instance
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
        seconds_list = []
        memory_list = []
        delta_list = []
        id_list = []

        # Check whether allocation logging had been enabled for the simulation by looping over the entries of the first log file
        allocation_logging = False
        for i in range(len(self.log_files[0].contents)):

            # Check whether memory (de)allocation is reported in this entry
            if not allocation_logging and "GB for" in self.log_files[0].contents["Message"][i]:
                allocation_logging = True
                break

        # Loop over all log files to determine the earliest recorded time
        t_0 = datetime.now()
        for log_file in self.log_files:
            if log_file.t_0 < t_0: t_0 = log_file.t_0

        # Loop over the log files again and fill the column lists
        for log_file in self.log_files:

            unique_ids = 0
            address_to_id = dict()

            # Get the process rank associated with this log file
            process = log_file.process

            # Loop over all log file entries
            for j in range(len(log_file.contents)):

                # Calculate the number of seconds that have passed since the earliest recorded log time
                seconds = (log_file.contents["Time"][j] - t_0).total_seconds()

                # Fill in the column lists
                process_list.append(process)
                phase_list.append(log_file.contents["Phase"][j])
                seconds_list.append(seconds)

                memory_list.append(log_file.contents["Memory"][j])

                if allocation_logging:

                    # Test whether this log entry contains information about memory (de)allocation
                    message = log_file.contents["Message"][j]
                    if "GB for" in message:

                        # Get the amount of memory (de)allocated (in GB)
                        if log_file.contents["Message"][j][0] == "+": delta = float(message.split("+")[1].split(" GB")[0])
                        elif log_file.contents["Message"][j][0] == "-": delta = -float(message.split("-", 1)[1].split(" GB")[0])
                        else: raise ValueError("Cannot determine the amount of memory (de)allocation")

                        # Get the address of the associated Array
                        address = log_file.contents["Message"][j].split("for ")[1]

                        # If the address was not yet encountered, assign a new unique id to it (an simple integer)
                        if address not in address_to_id:
                            address_to_id[address] = unique_ids
                            unique_ids += 1

                        id = address_to_id[address]

                    else:

                        delta = None
                        id = None

                    delta_list.append(delta)
                    id_list.append(id)

        # Create the table data structures
        names = ['Process rank', 'Simulation phase', 'Simulation time', 'Memory usage']
        data = [process_list, phase_list, seconds_list, memory_list]

        # If memory (de)allocation logging was enabled, add the 2 additional columns
        if allocation_logging:

            data += [delta_list, id_list]
            names += ["Array (de)allocation", "Array ID"]

        # Create the table
        self.table = Table(data, names=names)
        self.table["Simulation time"].unit = "s"
        self.table["Memory usage"].unit = "GB"
        if allocation_logging: self.table["Array (de)allocation"].unit = "GB"

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

    @property
    def processes(self):

        """
        This function ...
        :return:
        """

        return max(self.table["Process rank"])+1

    # -----------------------------------------------------------------

    @property
    def peak(self):

        """
        This function ...
        :return:
        """

        return self.peak_per_process * self.processes

    # -----------------------------------------------------------------

    @property
    def peak_per_process(self):

        """
        This function ...
        :return:
        """

        return(max(self.table["Memory usage"]))

# -----------------------------------------------------------------
