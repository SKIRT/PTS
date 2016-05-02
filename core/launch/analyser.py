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
from ...modeling.fitting.modelanalyser import FitModelAnalyser
from ...modeling.analysis.bestmodelanalyser import BestModelAnalyser
from ..tools.logging import log
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
        super(SimulationAnalyser, self).__init__(config, "core")

        # -- Attributes --

        # Set the simulation object to None initially
        self.simulation = None

        # The lower-level analysers
        self.basic_analyser = BasicAnalyser()
        self.batch_analyser = BatchAnalyser()
        self.scaling_analyser = ScalingAnalyser()
        self.fit_model_analyser = FitModelAnalyser()
        self.best_model_analyser = BestModelAnalyser()

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

        # 3. Run the batch analysis
        if self.simulation.from_batch: self.analyse_batch()

        # 3. Analyse the scaling, if the simulation is part of a scaling test
        if self.simulation.from_scaling_test: self.analyse_scaling()

        # 4. Analyse the goodness of fit of the radiative transfer model, if the simulation is part of a modeling run
        if self.simulation.from_modeling: self.analyse_model()

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
        self.fit_model_analyser.clear()
        self.best_model_analyser.clear()

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
        self.batch_analyser.run(self.simulation, self.basic_analyser.timeline_extractor, self.basic_analyser.memory_extractor)

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

        # Determine the path to the output directory for the analysis of the best model within the radiative transfer
        # modeling environment
        analysis_out_path = fs.join(self.simulation.analysis.modeling_path, "analysis", "out")

        # If the output directory of the simulation corresponds to the analysis output directory, the simulation
        # corresponds to the best model in the current state of the radiative transfer modeling
        if self.simulation.output_path == analysis_out_path: self.analyse_best_model()

        # Else, the simulation is just one of the many simulations launched during the fitting step of the modeling
        else: self.analyse_fit_model()

    # -----------------------------------------------------------------

    def analyse_fit_model(self):

        """
        This function ...
        :return:
        """

        # Run the fit model analyser
        self.fit_model_analyser.config.path = self.simulation.analysis.modeling_path
        self.fit_model_analyser.run(self.simulation, self.basic_analyser.flux_calculator)

    # -----------------------------------------------------------------

    def analyse_best_model(self):

        """
        This function ...
        :return:
        """

        # Run the best model analyser
        self.best_model_analyser.config.path = self.simulation.analysis.modeling_path
        self.best_model_analyser.run(self.simulation)

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        # Indicate that this simulation has been analysed
        self.simulation.analysed = True
        self.simulation.save()

        # If requested, remove the local output directory
        if self.simulation.remove_local_output: fs.remove_directory(self.simulation.output_path)

# -----------------------------------------------------------------
