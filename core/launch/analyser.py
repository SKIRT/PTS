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
from ..launch.basicanalyser import BasicAnalyser, steps, extraction, plotting, misc
from ..launch.batchanalyser import BatchAnalyser
from ..test.scalinganalyser import ScalingAnalyser
from ..basics.log import log
from ..simulation.simulation import RemoteSimulation
from ..tools import filesystem as fs
from ..simulation.remote import get_simulation_for_host
from ..tools import sequences
from .options import extraction_names, plotting_names, misc_names

# -----------------------------------------------------------------

batch = "batch"
scaling = "scaling"

# -----------------------------------------------------------------

all_steps = steps + [batch, scaling]

# -----------------------------------------------------------------

def reanalyse_simulation(simulation, steps, features=None, not_steps=None, not_features=None):

    """
    This function ...
    :param simulation:
    :param steps:
    :param features:
    :param not_steps:
    :param not_features:
    :return:
    """

    # Only one step defined
    if len(steps) == 1:

        # Flag indicating whether this simulation has been analysed or not
        simulation.analysed = False

        # Get the step
        step = steps[0]

        # Extraction
        if step == extraction:

            # All features
            if features is None: simulation.analysed_extraction = []

            # Select features
            else: simulation.analysed_extraction = sequences.elements_not_in_other(extraction_names, features, check_existing=True)

        # Plotting
        elif step == plotting:

            # All features
            if features is None: simulation.analysed_plotting = []

            # Select features
            else: simulation.analysed_plotting = sequences.elements_not_in_other(plotting_names, features, check_existing=True)

        # Misc
        elif step == misc:

            # All features
            if features is None: simulation.analysed_misc = []

            # Select features
            else: simulation.analysed_misc = sequences.elements_not_in_other(misc_names, features, check_existing=True)

        # Batch
        elif step == batch:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

            # Set analysed_batch flag to False
            simulation.analysed_batch = False

        # Scaling
        elif step == scaling:

            # Check whether features are not defined
            if features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

            # Set analysed_scaling flag to False
            simulation.analysed_scaling = False

        # Invalid
        else: raise ValueError("Invalid step: '" + step + "'")

    # Multiple steps defined
    else:

        # Check whether features are not defined
        if features is not None: raise ValueError("Features cannot be specified with multiple steps")

        # Flag indicating whether this simulation has been analysed or not
        simulation.analysed = False

        # Reset extraction
        if extraction in steps: simulation.analysed_extraction = []

        # Reset plotting
        if plotting in steps: simulation.analysed_plotting = []

        # Reset misc
        if misc in steps: simulation.analysed_misc = []

        # Reset batch
        if batch in steps: simulation.analysed_batch = False

        # Reset scaling
        if scaling in steps: simulation.analysed_scaling = False

    # Now analyse the simulation
    analyse_simulation(simulation, not_steps=not_steps, not_features=not_features)

# -----------------------------------------------------------------

def analyse_simulation(simulation, not_steps=None, not_features=None):

    """
    This function ...
    :param simulation:
    :param not_steps:
    :param not_features:
    :return:
    """

    # Steps are defined
    if not_steps is not None:

        # One step is defined
        if len(not_steps) == 1:

            # Get the step
            not_step = not_steps[0]

            # Extraction
            if not_step == extraction:

                # Features are not specified
                if not_features is None: simulation.analysed_extraction = extraction_names

                # Features are defined
                else: sequences.extend_unique(simulation.analysed_extraction, not_features)

            # Plotting
            elif not_step == plotting:

                # Features are not specified
                if not_features is None: simulation.analysed_plotting = plotting_names

                # Features are defined
                else: sequences.extend_unique(simulation.analysed_plotting, not_features)

            # Misc
            elif not_step == misc:

                # Features are not specified
                if not_features is None: simulation.analysed_misc = misc_names

                # Features are defined
                else: sequences.extend_unique(simulation.analysed_misc, not_features)

            # Batch
            elif not_step == batch:

                # Check whether features are not defined
                if not_features is not None: raise ValueError("Cannot define features for the batch simulation analysis")

                # Set flag
                simulation.analysed_batch = True

            # Scaling
            elif not_step == scaling:

                # Check whether features are not defined
                if not_features is not None: raise ValueError("Cannot define features for scaling simulation analysis")

                # Set flag
                simulation.analysed_scaling = True

            # Invalid
            else: raise ValueError("Invalid step '" + not_step + "'")

        # Multiple steps are defined
        else:

            # Check whether features are not defined
            if not_features is not None: raise ValueError("Features cannot be specified with multiple steps")

            # Set all extracted
            if extraction in not_steps: simulation.analysed_extraction = extraction_names

            # Set all plotted
            if plotting in not_steps: simulation.analysed_plotting = plotting_names

            # Set all misc
            if misc in not_steps: simulation.analysed_misc = misc_names

            # Set batch
            if batch in not_steps: simulation.analysed_batch = True

            # Set scaling
            if scaling in not_steps: simulation.analysed_scaling = True

    # Steps are not defined
    elif not_features is not None: raise ValueError("Cannot specify features when step is not defined")

    # Create simulation analyser
    analyser = SimulationAnalyser()

    # Run the analyser on the simulation
    analyser.run(simulation=simulation)

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
