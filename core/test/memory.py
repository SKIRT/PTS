#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.memory Contains the MemoryTester class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ..basics.configurable import Configurable
from ..simulation.execute import SkirtExec
from ..simulation.skifile import SkiFile
from ..simulation.logfile import LogFile
from ..advanced.memoryestimator import MemoryEstimator
from ..launch.batchlauncher import BatchLauncher
from ..simulation.definition import SingleSimulationDefinition

# -----------------------------------------------------------------

class MemoryTester(Configurable):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(MemoryTester, self).__init__(config)

        # Local SKIRT execution environment
        self.skirt = SkirtExec()

        # Remote SKIRT batch execution environment
        self.launcher = BatchLauncher()

        # The memory estimator
        self.estimator = MemoryEstimator()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function
        self.setup(**kwargs)

        # Launch the simulations
        self.launch()

        # Load the data
        self.load_data()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MemoryTester, self).setup(**kwargs)

        # Setup the batch launcher
        self.setup_launcher()

    # -----------------------------------------------------------------

    def setup_launcher(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the batch launcher ...")

        from ..config.launch_batch import definition
        from ..basics.configuration import InteractiveConfigurationSetter

        #definition = ConfigurationDefinition(write_config=False)
        setter = InteractiveConfigurationSetter(self.class_name, add_logging=False)

        # Create new config
        self.launcher.config = setter.run(definition, prompt_optional=False)

        # Set working directory
        self.launcher.config.path = self.config.path

        # Set options for the batch launcher: basic options
        self.launcher.config.shared_input = True
        #self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = [self.config.remote]  # the remote hosts on which to run the simulations

        #self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        #self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file
        #self.launcher.config.retrieve_types = ["log"]  # only log files should be retrieved from the simulation output
        #self.launcher.config.keep = self.config.keep

        # Run in attached mode
        self.launcher.config.attached = True

        self.launcher.config.cores_per_process = 4

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        if self.config.ski is None:

            # Loop over all files in the current working directory
            #for ski_path, prefix in fs.files_in_path(self.config.path, extension="ski", returns=["path", "name"]):

                # Create simulation definition
                #definition = SingleSimulationDefinition(ski_path, output_path, input_path)

                # Add to the queue
                #self.launcher.add_to_queue(definition)

            # Run the batch launcher
            self.launcher.run()

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        if self.config.ski is None:

            # Loop over all xml files
            #for path, name in fs.files_in_path(self.config.path, extension="xml", returns=["path", "name"]):

                #original_name = name + ".xml"
                #new_name = name.replace("_parameters", "") + ".ski"
                #fs.rename_file(self.config.path, original_name, new_name)

            # Loop over all files in the current working directory
            for path, name in fs.files_in_path(self.config.path, extension="ski", returns=["path", "name"]):

                # Log path
                log_path = fs.join(self.config.path, name + "_log.txt")

                # Check if the log file is present
                if not fs.is_file(log_path): continue

                # Open the log file
                log_file = LogFile(log_path)

                # Get the peak memory usage
                memory = log_file.peak_memory

                # Load the ski file
                ski = SkiFile(path)

                # Configure the memory estimator tool
                self.estimator.config.ski = ski
                self.estimator.config.input = self.config.path
                self.estimator.config.ncells = None
                self.estimator.config.probe = False

                try:
                    # Run the estimator
                    self.estimator.run()
                except ValueError:
                    print(name)
                    exit()

                # Get the parallel and serial part of the memory
                parallel = self.estimator.parallel_memory
                serial = self.estimator.serial_memory

                print(memory, parallel, serial, parallel + serial)

        else: pass #

            # Calculate the procentual difference
            #diff = (consumed_memory - estimated_memory) / consumed_memory * 100.0

# -----------------------------------------------------------------
