#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.parameterexploration Contains the ParameterExplorer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem
from ...core.simulation.arguments import SkirtArguments
from ...core.basics.filter import Filter
from ...core.simulation.skifile import SkiFile
from ...core.launch.batchlauncher import BatchLauncher

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

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski file
        self.ski = None

        # The parameter ranges
        self.young_luminosities = None
        self.ionizing_luminosities = None
        self.dust_masses = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new ParameterExplorer instance
        explorer = cls(arguments.config)

        # Set the output path
        explorer.config.path = arguments.path
        explorer.config.output_path = filesystem.join(arguments.path, "fit")

        # Set options for the young stellar population
        if arguments.young_nvalues is not None: explorer.config.young_stars.nvalues = arguments.young_nvalues
        if arguments.young_range is not None:
            explorer.config.young_stars.rel_min = arguments.young_range[0]
            explorer.config.young_stars.rel_max = arguments.young_range[1]
        if arguments.young_log: explorer.config.young_stars.scale = "log"
        else: explorer.config.young_stars.scale = "linear"

        # Set options for the ionizing stellar population
        if arguments.ionizing_nvalues is not None: explorer.config.ionizing_stars.nvalues = arguments.ionizing_nvalues
        if arguments.ionizing_range is not None:
            explorer.config.ionizing_stars.rel_min = arguments.ionizing_range[0]
            explorer.config.ionizing_stars.rel_max = arguments.ionizing_range[1]
        if arguments.ionizing_log: explorer.config.ionizing_stars = "log"
        else: explorer.config.ionizing_stars.scale = "linear"

        # Set options for the dust component
        if arguments.dust_nvalues is not None: explorer.config.dust.nvalues = arguments.dust_nvalues
        if arguments.dust_range is not None:
            explorer.config.dust.rel_min = arguments.dust_range[0]
            explorer.config.dust.rel_max = arguments.dust_range[1]
        if arguments.dust_log: explorer.config.dust.scale = "log"
        else: explorer.config.dust.scale = "linear"

        # Return the new instance
        return explorer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file
        self.load_ski()

        # 3. Set the ranges of the different fit parameters
        self.set_parameter_ranges()

        # 3. Launch the simulations for different parameter values
        self.simulate()

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Open the ski file (created by InputInitializer)
        self.ski = SkiFile(self.fit_ski_path)

    # -----------------------------------------------------------------

    def set_parameter_ranges(self):

        """
        This function ...
        :return:
        """

        # Get the current values in the ski file prepared by InputInitializer
        young_luminosity, young_filter = self.ski.get_stellar_component_luminosity("Young stars")
        ionizing_luminosity, ionizing_filter = self.ski.get_stellar_component_luminosity("Ionizing stars")
        dust_mass = self.ski.get_dust_component_mass(0)

        # Set the parameter ranges
        self.set_young_luminosity_range(young_luminosity)
        self.set_ionizing_luminosity_range(ionizing_luminosity)
        self.set_dust_mass_range(dust_mass)

    # -----------------------------------------------------------------

    def set_young_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Set the range of the FUV luminosity of the young stellar population
        min = self.config.young_stars.rel_min * luminosity
        max = self.config.young_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.young_stars.scale == "linear":
            self.young_luminosities = np.linspace(min, max, num=self.config.young_stars.nvalues, endpoint=True)
        elif self.config.young_stars.scale == "logarithmic":
            self.young_luminosities = np.logspace(min, max, num=self.config.young_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the young stellar luminosity values")

    # -----------------------------------------------------------------

    def set_ionizing_luminosity_range(self, luminosity):

        """
        This function ...
        :param luminosity:
        :return:
        """

        # Determine the minimum and maximum luminosity
        min = self.config.ionizing_stars.rel_min * luminosity
        max = self.config.ionizing_stars.rel_max * luminosity

        # Create a linear or logarithmic range of luminosities
        if self.config.ionizing_stars.scale == "linear":
            self.ionizing_luminosities = np.linspace(min, max, num=self.config.ionizing_stars.nvalues, endpoint=True)
        elif self.config.ionizing_stars.scale == "log":
            self.ionizing_luminosities = np.logspace(min, max, num=self.config.ionizing_stars.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the ionizing stellar luminosity values")

    # -----------------------------------------------------------------

    def set_dust_mass_range(self, mass):

        """
        This function ...
        :param mass:
        :return:
        """

        # Set the dust mass range
        min = self.config.dust.rel_min * mass
        max = self.config.dust.rel_max * mass

        # Create a linear or logarithmic range of dust masses
        if self.config.dust.scale == "linear":
            self.dust_masses = np.linspace(min, max, num=self.config.dust.nvalues, endpoint=True)
        elif self.config.dust.scale == "log":
            self.dust_masses = np.logspace(min, max, num=self.config.dust.nvalues, endpoint=True)
        else: raise ValueError("Invalid scale for the dust mass values")

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Scheduling options
        scheduling_options = {}
        scheduling_options["walltime"] = None

        # Create a FUV filter object
        fuv_filter = Filter.from_string("FUV")

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # The ski file pattern
        arguments.ski_pattern = None
        arguments.recursive = False
        arguments.relative = False

        # Input and output
        arguments.input_path = self.fit_in_path
        arguments.output_path = None

        # Parallelization settings
        arguments.parallel.threads = threads
        arguments.parallel.processes = processes

        # Loop over the different values of the young stellar luminosity
        for young_luminosity in self.young_luminosities:

            # Loop over the different values of the ionizing stellar luminosity
            for ionizing_luminosity in self.ionizing_luminosities:

                # Loop over the different values of the dust mass
                for dust_mass in self.dust_masses:

                    # Create a unique name for this combination of parameter values
                    simulation_name = str(young_luminosity) + "_" + str(ionizing_luminosity) + "_" + str(dust_mass)

                    # Change the parameter values in the ski file
                    self.ski.set_stellar_component_luminosity("Young stars", young_luminosity, fuv_filter) # TODO: units!
                    self.ski.set_stellar_component_luminosity("Ionizing stars", ionizing_luminosity, fuv_filter) # TODO: units!
                    self.ski.set_dust_component_mass(0, dust_mass) # TODO: units !

                    # Determine the directory for this simulation
                    simulation_path = filesystem.join(self.fit_out_path, simulation_name)

                    # Create the simulation directory
                    filesystem.create_directory(simulation_path)

                    # Put the ski file with adjusted parameters into the simulation directory
                    ski_path = filesystem.join(simulation_path, self.galaxy_name + ".ski")
                    self.ski.saveto(ski_path)

                    # Create a directory 'out' in the simulation directory to contain the SKIRT output
                    ouput_path = filesystem.join(simulation_path, "out")
                    filesystem.create_directory(output_path)

                    # Adjust the SKIRT arguments
                    arguments.ski_pattern = ski_path
                    arguments.output_path = output_path

                    # Put the parameters in the queue and get the simulation object
                    self.launcher.add_to_queue(arguments, scheduling_options)

        # Run the launcher, schedules the simulations
        simulations = self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.modeling_path = self.config.path

            # Save the simulation object
            simulation.save()

# -----------------------------------------------------------------
