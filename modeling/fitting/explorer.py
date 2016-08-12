#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.explorer Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.launch.batchlauncher import BatchLauncher
from .modelgenerators.grid import GridModelGenerator
from .modelgenerators.initial import InitialModelGenerator
from .modelgenerators.genetic import GeneticModelGenerator
from .modelgenerators.instinctive import InstinctiveModelGenerator
from ...core.basics.filter import Filter
from ...core.tools import time
from ...core.basics.range import IntegerRange, RealRange, QuantityRange
from ...core.simulation.definition import SingleSimulationDefinition

# -----------------------------------------------------------------

class ParameterExplorer(FittingComponent):
    
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
        super(ParameterExplorer, self).__init__(config)

        # -- Attributes --

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The parameter ranges
        self.ranges = dict()

        # The generation index and name
        self.generation = None
        self.generation_name = None

        # The model generator
        self.generator = None

        # The ski file
        self.ski = None

        # The paths to the simulation input files
        self.input_paths = None

        # Initial generation is present
        self.has_initial = False

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 3. Set the parameter ranges
        self.set_ranges()

        # 4. Generate the model parameters
        self.generate_models()

        # 5. Set the paths to the input files
        self.set_input()

        # 5. Set the parallelization schemes for the different remote hosts
        self.set_parallelization()

        # 6. Estimate the runtimes for the different remote hosts
        self.estimate_runtimes()

        # 7. Launch the simulations for different parameter values
        self.launch()

        # 8. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ParameterExplorer, self).setup()

        # Check whether the initial generation has been created
        initial_generation_path = fs.join(self.fit_generations_path, "initial")
        self.has_initial = fs.is_directory(initial_generation_path)

        # Set options for the batch launcher
        self.set_launcher_options()

        # Set the model generator
        self.set_generator()

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Set options for the BatchLauncher: basic options
        self.launcher.config.shared_input = True  # The input directories (or files) for the different simulations are shared
        self.launcher.config.group_simulations = True  # group multiple simulations into a single job (because a very large number of simulations will be scheduled)
        self.launcher.config.remotes = self.config.remotes  # the remote hosts on which to run the simulations
        self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file

        # Set options for the BatchLauncher: simulation analysis options
        self.launcher.config.analysis.relative = True
        self.launcher.config.analysis.extraction.path = "res"
        self.launcher.config.analysis.misc.path = "res"  # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        self.launcher.config.analysis.plotting.path = "plot"  # The base directory where all of the simulations will have a seperate directory with the plotting analysis output
        self.launcher.config.analysis.extraction.timeline = True  # extract the simulation timeline
        self.launcher.config.analysis.plotting.seds = True  # Plot the output SEDs
        self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_path  # the path to the reference SED (for plotting the simulated SED against the reference points)
        # self.launcher.config.analysis.plotting.reference_sed = self.observed_sed_dustpedia_path # the path to the DustPedia SED
        # self.launcher.config.analysis.misc.fluxes = True  # Calculate observed fluxes
        # self.launcher.config.analysis.misc.images = True  # Make observed images
        self.launcher.config.analysis.misc.observation_filters = self.observed_filter_names  # The filters for which to create the observations
        self.launcher.config.analysis.plotting.format = "png"  # plot in PNG format so that an animation can be made from the fit SEDs

        # Set remote for the 'extra' simulations
        self.launcher.config.extra_remote = "nancy"

    # -----------------------------------------------------------------

    def set_generator(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the model generator ...")

        # Generate new models based on a simple grid (linear or logarithmic) of parameter values
        if self.config.generation_method == "grid":

            # Set a name for the generation
            self.generation_name = time.unique_name()

            # Create the model generator
            self.generator = GridModelGenerator()

        # Generate new models instinctively based on the current probability distribution of the parameters
        elif self.config.generation_method == "instinctive":

            # Set a name for the generation
            self.generation_name = time.unique_name()

            # Create the model generator
            self.generator = InstinctiveModelGenerator()

        # Generate new models using genetic algorithms
        elif self.config.generation_method == "genetic":

            if self.has_initial:

                self.generation = None # Determine from index of last generation
                self.generation_name = str()

                # Create the model generator
                self.generator = GeneticModelGenerator()

            else:

                self.generation_name = "initial"

                # Create the model generator
                self.generator = InitialModelGenerator()

    # -----------------------------------------------------------------

    def set_ranges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parameter ranges ...")

        # The given ranges are relative to the best or initial value
        if self.config.relative:

            # TODO: get parameter values of best fitting model

            # Get the current values in the ski file prepared by InputInitializer
            young_luminosity_guess, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
            ionizing_luminosity_guess, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
            dust_mass_guess = self.ski.get_dust_component_mass(0)

            # Set the range of the FUV luminosity of the young stellar population
            min_value = self.config.young[0] * young_luminosity_guess
            max_value = self.config.young[1] * young_luminosity_guess
            self.ranges["FUV young"] = RealRange(min_value, max_value)

            # Set the range of the FUV luminosity of the ionizing stellar population
            min_value = self.config.ionizing[0] * ionizing_luminosity_guess
            max_value = self.config.ionizing[1] * ionizing_luminosity_guess
            self.ranges["FUV ionizing"] = RealRange(min_value, max_value)

            # Set the range of the dust mass
            min_value = self.config.dust[0] * dust_mass_guess
            max_value = self.config.dust[1] * dust_mass_guess
            self.ranges["Dust mass"] = QuantityRange(min_value, max_value)

        else:

            # Set the ranges directly from the command-line arguments
            self.ranges["FUV young"] = self.config.young
            self.ranges["FUV ionizing"] = self.config.ionizing
            self.ranges["Dust mass"] = self.config.dust

    # -----------------------------------------------------------------

    def generate_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the model parameters ...")

        # Run the model generator
        self.generator.run()

    # -----------------------------------------------------------------

    def set_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the input paths ...")

        # Set the paths to the input maps and appropriate wavelength grid file
        self.input_paths = self.input_map_paths
        self.input_paths.append() # add the path to the wavelength grid

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Set the paths to the directories to contain the launch scripts (job scripts) for the different remote hosts
        for host_id in self.launcher.host_ids:
            script_dir_path = fs.join(self.fit_scripts_path, host_id)
            if not fs.is_directory(script_dir_path): fs.create_directory(script_dir_path)
            self.launcher.set_script_path(host_id, script_dir_path)

        # Set the paths to the screen output directories (for debugging) for remotes without a scheduling system for jobs
        for host_id in self.launcher.no_scheduler_host_ids: self.launcher.enable_screen_output(host_id)

        # Create a FUV filter object
        fuv = Filter.from_string("FUV")

        # Loop over the different parameter combinations
        for i in range(self.number_of_models):

            # Get the parameter values
            young_luminosity = self.parameters["FUV young"][i]
            ionizing_luminosity = self.parameters["FUV ionizing"][i]
            dust_mass = self.parameters["Dust mass"][i]

            # Create a unique name for this combination of parameter values
            simulation_name = time.unique_name()

            # Change the parameter values in the ski file
            self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv.centerwavelength() * Unit("micron"))
            self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv.centerwavelength() * Unit("micron"))
            self.ski.set_dust_component_mass(0, dust_mass)

            # Determine the directory for this simulation
            simulation_path = fs.join(self.fit_out_path, simulation_name)

            # Create the simulation directory
            fs.create_directory(simulation_path)

            # Create an 'out' directory within the simulation directory
            output_path = fs.join(simulation_path, "out")
            fs.create_directory(output_path)

            # Put the ski file with adjusted parameters into the simulation directory
            ski_path = fs.join(simulation_path, self.galaxy_name + ".ski")
            self.ski.saveto(ski_path)

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, self.input_paths, output_path)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug(" - ski path: " + definition.ski_path)
            log.debug(" - output path: " + definition.output_path)

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name)

            # Set scheduling options (for the different remote hosts with a scheduling system)
            for host_id in self.scheduling_options: self.launcher.set_scheduling_options(host_id, simulation_name, self.scheduling_options[host_id])

            # Add an entry to the parameter table
            self.table.add_row([simulation_name, young_luminosity, ionizing_luminosity, dust_mass])

        # Add simulations to the extra queue to calculate the contribution of the various stellar components
        self.add_contribution_simulations()

        # Add simulation to the extra queue to create simulated images
        self.add_images_simulation()

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in simulations:

            # Skip simulations for which analysis options are not defined (the contribution simulations and the images simulation)
            if simulation.analysis is None or simulation.name == "images": continue

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Save the simulation object
            simulation.save()

    # -----------------------------------------------------------------

    @property
    def number_of_models(self):

        """
        This function ...
        :return:
        """

        return self.generator.number_of_models

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------
