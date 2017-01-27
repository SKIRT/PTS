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
from ..tools.logging import log
from ..simulation.simulation import RemoteSimulation
from ..tools import filesystem as fs

# -----------------------------------------------------------------

class SimulationAnalyser(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(SimulationAnalyser, self).__init__(config)

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = BasicAnalyser()
        self.batch_analyser = BatchAnalyser()
        self.scaling_analyser = ScalingAnalyser()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # If the simulation has no analysis options, finish the procedure right away
        if self.simulation.analysis is None: return

        # 2. Run the basic analysis
        self.analyse_basic()

        # 3. Run the batch analysis
        if self.simulation.from_batch: self.analyse_batch()

        # 3. Analyse the scaling, if the simulation is part of a scaling test
        if self.simulation.from_scaling_test: self.analyse_scaling()

        # 4. Perform extra analysis
        self.analyse_extra()

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
        self.simulation = kwargs.pop("simulation")

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

            # Debugging
            log.debug("Running the " + analyser_class.__name__ + " on the simulation ...")

            # Create an instance of the analyser class
            analyser = analyser_class.for_simulation(self.simulation)

            # Run the analyser, giving this simulation analyser instance as an argument
            analyser.run(self)

        # Indicate that this simulation has been analysed
        self.simulation.analysed = True
        self.simulation.save()

        # If requested, remove the local output directory
        if isinstance(self.simulation, RemoteSimulation) and self.simulation.remove_local_output: fs.remove_directory(self.simulation.output_path)

# -----------------------------------------------------------------
