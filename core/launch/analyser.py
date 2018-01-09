#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.launch.analyser Contains the SimulationAnalyser class, used for analysing simulation output.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable
from ..launch.basicanalyser import BasicAnalyser
from ..launch.batchanalyser import BatchAnalyser
from ..test.scalinganalyser import ScalingAnalyser
from ..basics.log import log
from ..simulation.simulation import RemoteSimulation
from ..tools import filesystem as fs
from ..simulation.remote import get_simulation_for_host

# -----------------------------------------------------------------

class SimulationAnalyser(Configurable):

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
        super(SimulationAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = BasicAnalyser()
        self.batch_analyser = BatchAnalyser()
        self.scaling_analyser = ScalingAnalyser()

    # -----------------------------------------------------------------

    @property
    def analysed_batch(self):

        """
        This function ...
        :return:
        """

        return self.simulation.analysed_batch

    # -----------------------------------------------------------------

    @property
    def analysed_scaling(self):

        """
        This function ...
        :return:
        """

        return self.simulation.analysed_scaling

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. If the simulation has no analysis options, finish the procedure right away
        if self.simulation.analysis is None: return

        # 3. Run the basic analysis
        if self.config.basic: self.analyse_basic()

        # 4. Run the batch analysis
        if self.simulation.from_batch and not self.analysed_batch: self.analyse_batch()

        # 5. Analyse the scaling, if the simulation is part of a scaling test
        if self.simulation.from_scaling_test and not self.analysed_scaling: self.analyse_scaling()

        # 6. Perform extra analysis
        if self.config.extra: self.analyse_extra()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationAnalyser, self).setup(**kwargs)

        # Make a local reference to the simulation object
        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        elif self.config.remote is not None and self.config.id is not None: self.load_simulation()
        else: raise ValueError("No simulation is specified")

        # Set flags
        self.basic_analyser.config.ignore_missing_data = self.config.ignore_missing_data
        self.batch_analyser.config.ignore_missing_data = self.config.ignore_missing_data

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

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the simulation analyser ...")

        # Set the simulation to None
        self.simulation = None

        # Clear the analysers
        self.basic_analyser.clear()
        self.scaling_analyser.clear()

    # -----------------------------------------------------------------

    def analyse_basic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the simulation output ...")

        # Run the analyser on the simulation
        self.basic_analyser.run(simulation=self.simulation)

    # -----------------------------------------------------------------

    def analyse_batch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the properties relevant for the batch of simulations ...")

        # Run the batch analyser on the simulation
        self.batch_analyser.run(simulation=self.simulation, timeline=self.basic_analyser.timeline, memory=self.basic_analyser.memory)

        # Set flag
        self.simulation.analysed_batch = True
        self.simulation.save()

    # -----------------------------------------------------------------

    def analyse_scaling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the scaling results ...")

        # Run the scaling analyser
        self.scaling_analyser.run(simulation=self.simulation, timeline=self.basic_analyser.timeline, memory=self.basic_analyser.memory)

        # Set flag
        self.simulation.analysed_scaling = True
        self.simulation.save()

    # -----------------------------------------------------------------

    def analyse_extra(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing extra analysis on the simulation output ...")

        # Loop over the 'extra' analyser classes that are defined for this simulation
        for analyser_class in self.simulation.analyser_classes:

            # Get name
            class_name = analyser_class.__name__

            # Check
            if class_name in self.simulation.analysed_extra:
                log.debug("Analysis with the " + class_name + " class has already been performed")
                continue

            # Debugging
            log.debug("Running the " + class_name + " on the simulation ...")

            # Create an instance of the analyser class
            analyser = analyser_class.for_simulation(self.simulation)

            # Run the analyser, giving this simulation analyser instance as an argument
            analyser.run(simulation_analyser=self)

            # Add name to analysed_extra
            self.simulation.analysed_extra.append(class_name)
            self.simulation.save()

        # Indicate that this simulation has been analysed
        self.simulation.analysed = True
        self.simulation.save()

        # If requested, remove the local output directory
        if isinstance(self.simulation, RemoteSimulation) and self.simulation.remove_local_output: fs.remove_directory(self.simulation.output_path)

# -----------------------------------------------------------------
