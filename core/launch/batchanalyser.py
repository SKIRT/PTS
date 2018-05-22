#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.batchanalyser Contains the BatchAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ...core.basics.log import log
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable
from ..simulation.remote import get_simulation_for_host

# -----------------------------------------------------------------

class BatchAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(BatchAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The paths to the timing and memory table files
        self.timing_table_path = None
        self.memory_table_path = None

        # The timeline and memory usage tables
        self.timeline = None
        self.memory = None

        # The log file of the simulation
        self.log_file = None

        # The ski file corresponding to the simulation
        self.ski = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the log file of the simulation
        self.load_log_file()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the batch analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.timing_table_path = None
        self.memory_table_path = None
        self.timeline = None
        self.memory = None
        self.log_file = None
        self.ski = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchAnalyser, self).setup(**kwargs)

        # Make a local reference to the simulation object
        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        elif self.config.remote is not None and self.config.id is not None: self.load_simulation()
        else: raise ValueError("No simulation is specified")

        # Also make references to the timing and memory table files (for shortness of notation)
        self.timing_table_path = self.simulation.analysis.timing_table_path
        self.memory_table_path = self.simulation.analysis.memory_table_path

        # Make a reference to the timeline and memory usage tables
        self.timeline = kwargs.pop("timeline")
        self.memory = kwargs.pop("memory")

        # Load the ski file
        self.ski = self.simulation.parameters()

    # -----------------------------------------------------------------

    def load_simulation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the simulation ...")

        # Load simulation
        self.simulation = get_simulation_for_host(self.config.remote, self.config.id)

    # -----------------------------------------------------------------

    @property
    def has_timeline(self):

        """
        This function ...
        :return:
        """

        return self.timeline is not None

    # -----------------------------------------------------------------

    @property
    def has_memory(self):

        """
        This function ...
        :return:
        """

        return self.memory is not None

    # -----------------------------------------------------------------

    def load_log_file(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulation log file ...")

        # Get the log file produced by the simulation
        self.log_file = self.simulation.log_file

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the timing information
        if self.timing_table_path is not None: self.write_timing()

        # Write the memory information
        if self.memory_table_path is not None: self.write_memory()

    # -----------------------------------------------------------------

    def write_timing(self):

        """
        This function ...
        :return:
        """

        # Check
        if not self.has_timeline:
            if self.config.ignore_missing_data:
                log.warning("No timing information for this simulation: skipping ...")
                return
            else: raise RuntimeError("No timing information was found")

        # Inform the user
        log.info("Writing the timing information of the simulation ...")

        # Open the timing table
        timing_table = TimingTable.from_file(self.timing_table_path)

        # Already simulation with this name in the table?
        if timing_table.has_simulation(self.simulation.name):
            if self.config.replace: timing_table.remove_simulation(self.simulation.name)
            else: raise ValueError("Simulation with name '" + self.simulation.name + "' is already in the timing table")

        # Add an entry
        unique_name = timing_table.add_from_simulation(self.simulation, self.ski, self.log_file, self.timeline)

        # Check
        if unique_name != self.simulation.name: raise RuntimeError("The simulation did not have a unique name: this shouldn't happen")

        # Save the table
        timing_table.save()

    # -----------------------------------------------------------------

    def write_memory(self):

        """
        This function ...
        :return:
        """

        # Check
        if not self.has_memory:
            if self.config.ignore_missing_data:
                log.warning("No memory usage information for this simulation: skipping ...")
                return
            else: raise RuntimeError("No memory information was found")

        # Inform the user
        log.info("Writing the memory usage information of the simulation ...")

        # Open the memory table
        memory_table = MemoryTable.from_file(self.memory_table_path)

        # Already a simulation with this name in the table
        if memory_table.has_simulation(self.simulation.name):
            if self.config.replace: memory_table.remove_simulation(self.simulation.name)
            else: raise ValueError("Simulation with name '" + self.simulation.name + "' is already in the memory table")

        # Add an entry
        unique_name = memory_table.add_from_simulation(self.simulation, self.ski, self.log_file)

        # Check
        if unique_name != self.simulation.name: raise RuntimeError("The simulation did not have a unique name: this shouldn't happen")

        # Save the table
        memory_table.save()

# -----------------------------------------------------------------
