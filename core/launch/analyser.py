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
from ..basics.configurable import OldConfigurable
from ..launch.basicanalyser import BasicAnalyser
from ..launch.batchanalyser import BatchAnalyser
from ..test.scalinganalyser import ScalingAnalyser
from ..tools.logging import log

# -----------------------------------------------------------------

class SimulationAnalyser(OldConfigurable):

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
        super(SimulationAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = BasicAnalyser()
        self.batch_analyser = BatchAnalyser()
        self.scaling_analyser = ScalingAnalyser()

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :param simulation
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

        # If the simulation has no analysis options, finish the procedure right away
        if self.simulation.analysis is None:
            self.finish()
            return

        # 2. Run the basic analysis
        self.analyse_basic()

        # 3. Run the batch analysis
        if self.simulation.from_batch: self.analyse_batch()

        # 3. Analyse the scaling, if the simulation is part of a scaling test
        if self.simulation.from_scaling_test: self.analyse_scaling()

        # 5. Finish the analysis for the current simulation
        self.finish()

    # -----------------------------------------------------------------

    def setup(self, simulation):

        """
        This function ...
        :param simulation:
        :return:
        """

        # Call the setup function of the base class
        super(SimulationAnalyser, self).setup()

        # Make a local reference to the simulation object
        self.simulation = simulation

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
        self.basic_analyser.run(self.simulation)

    # -----------------------------------------------------------------

    def analyse_batch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the properties relevant for the batch of simulations ...")

        # Run the batch analyser on the simulation
        self.batch_analyser.run(self.simulation, self.basic_analyser.timeline, self.basic_analyser.memory)

    # -----------------------------------------------------------------

    def analyse_scaling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the scaling results ...")

        # Run the scaling analyser
        self.scaling_analyser.run(self.simulation, self.basic_analyser.timeline, self.basic_analyser.memory)

# -----------------------------------------------------------------
