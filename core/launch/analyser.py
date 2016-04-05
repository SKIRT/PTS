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
from ..test.scalinganalyser import ScalingAnalyser
from ...modeling.fitting.modelanalyser import ModelAnalyser
from ..tools.logging import log

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
        super(SimulationAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = BasicAnalyser()
        self.scaling_analyser = ScalingAnalyser()
        self.model_analyser = ModelAnalyser()

    # -----------------------------------------------------------------

    def run(self, simulation):

        """
        This function ...
        :param simulation
        :return:
        """

        # 1. Call the setup function
        self.setup(simulation)

        # 2. Run the basic analysis
        self.analyse_basic()

        # 3. Analyse the scaling, if the simulation is part of a scaling test
        if self.simulation.from_scaling_test: self.analyse_scaling()

        # 4. Analyse the goodness of fit of the radiative transfer model, if the simulation is part of a modeling run
        if self.simulation.from_modeling: self.analyse_model()

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
        self.model_analyser.clear()

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

    def analyse_scaling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the scaling results ...")

        # Run the scaling analyser
        self.scaling_analyser.run(self.simulation, self.basic_analyser.timeline_extractor, self.basic_analyser.memory_extractor)

    # -----------------------------------------------------------------

    def analyse_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Analysing the radiative transfer model ...")

        # Run the modeling analyser
        self.model_analyser.config.path = self.simulation.modeling_path
        self.model_analyser.run(self.simulation)

# -----------------------------------------------------------------
