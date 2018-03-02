#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.extract.memory Contains the MemoryUsageTable and MemoryExtractor classes. The latter is used for
# extracting memory information from SKIRT simulation log files into a MemoryUsageTable object.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from datetime import datetime

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.table import SmartTable

# -----------------------------------------------------------------

class NoMemoryData(Exception):

    """
    This class ...
    """

    def __init__(self, message, simulation_name=None):

        """
        Thisf unction ...
        :param message:
        :param simulation_name:
        """

        # Call the base class constructor with the parameters it needs
        super(NoMemoryData, self).__init__(message)

        # The simulation name
        self.simulation_name = simulation_name

# -----------------------------------------------------------------

class MemoryUsageTable(SmartTable):

    """
    This function ...
    """

    @classmethod
    def from_columns(cls, process_list, phase_list, seconds_list, memory_list, delta_list=None, id_list=None):

        """
        This function ...
        :param process_list:
        :param phase_list:
        :param seconds_list:
        :param memory_list:
        :param delta_list:
        :param id_list:
        :return:
        """

        names = ['Process rank', 'Phase', 'Simulation time', 'Memory usage']
        columns = [process_list, phase_list, seconds_list, memory_list]
        units = [None, None, "s", "GB"]

        # Add column of allocated memory
        if delta_list is not None:
            names.append("Array (de)allocation")
            columns.append(delta_list)
            units.append("GB")

        # Add column of allocation memory adresses
        if id_list is not None:
            names.append("Array ID")
            columns.append(id_list)
            units.append(None)

        # NEW
        return super(MemoryUsageTable, cls).from_columns(*columns, names=names, units=units)

    # -----------------------------------------------------------------

    @property
    def has_allocation_info(self):

        """
        This property ...
        :return:
        """

        return "Array (de)allocation" in self.colnames

    # -----------------------------------------------------------------

    @property
    def processes(self):

        """
        This function ...
        :return:
        """

        return max(self["Process rank"]) + 1

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

        return max(self["Memory usage"])

# -----------------------------------------------------------------

def extract_memory(simulation):

    """
    This function ...
    :param simulation:
    :return:
    """

    # Create a new MemoryExtractor instance
    extractor = MemoryExtractor()

    # Run the extractor and get the memory table
    extractor.run(simulation)

    # Return the table
    return extractor.table

# -----------------------------------------------------------------

def extract_memory_cwd():

    """
    This function ...
    :return:
    """

    from pts.core.simulation.simulation import createsimulations

    # Create a SkirtSimulation object based on a log file present in the current working directory
    simulation = createsimulations(single=True)

    # Extract
    return extract_memory(simulation)

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

        # -- Attributes --

        self.log_files = None
        self.table = None

        self.output_path = None

    # -----------------------------------------------------------------

    def run(self, simulation, output_path=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation, output_path=output_path)

        # 2. Perform the extraction
        self.extract()

        # 3. Write the results
        if self.output_path is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, simulation, output_path=None):

        """
        This function ...
        :param simulation:
        :param output_path:
        :return:
        """

        # Obtain the log files created by the simulation
        self.log_files = simulation.logfiles()

        # Check whether there are log files
        if not self.has_logfiles: raise ValueError("No log files found in simulation output")

        # Check whether the log files contain memory information
        if not self.log_files[0].has_memory: raise NoMemoryData("The log files don't contain memory information", simulation_name=simulation.name)

        # Set the output path
        self.output_path = output_path

    # -----------------------------------------------------------------

    @property
    def nlogfiles(self):

        """
        This function ...
        :return:
        """

        return len(self.log_files)

    # -----------------------------------------------------------------

    @property
    def has_logfiles(self):

        """
        This function ...
        :return:
        """

        return self.nlogfiles > 0

    # -----------------------------------------------------------------

    def extract(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting ...")

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

        # Set delta_list and id_list to None if not allocation logging
        if not allocation_logging:
            delta_list = None
            id_list = None

        # Create the memory usage table
        self.table = MemoryUsageTable.from_columns(process_list, phase_list, seconds_list, memory_list, delta_list, id_list)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the table to file
        self.table.saveto(self.output_path)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Set the table to None
        self.table = None

# -----------------------------------------------------------------
