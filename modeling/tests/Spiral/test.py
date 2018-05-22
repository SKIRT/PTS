#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.units.parsing import parse_unit as u
from pts.core.basics.configuration import Configuration
from pts.core.simulation.skifile import SkiFile
from pts.core.basics.range import QuantityRange, RealRange
from pts.core.basics.map import Map
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "determining parameters based on mock observations of a simple spiral galaxy model"

# -----------------------------------------------------------------

# Determine the ski path
ski_path = fs.join(this_dir_path, "spiral.ski")

# Get the initial dust mass of the exponential disk with spiral structure
ski = SkiFile(ski_path)
dust_mass = ski.get_labeled_value("exp_dustmass")

# -----------------------------------------------------------------

class SpiralTest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SpiralTest, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        super(SpiralTest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Determine the simulation output path
        simulation_output_path = "./ref"

        # Settings
        settings_launch = dict()
        settings_launch["ski"] = ski_path
        settings_launch["output"] = simulation_output_path
        settings_launch["create_output"] = True

        # Input
        input_launch = dict()

        # Construct the command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")

        # Add the command
        #commands.append(launch)

        launcher = self.run_command(launch)

    # -----------------------------------------------------------------

    def create_sed(self):

        """
        This function ...
        :return:
        """

        # Create simulation
        #prefix = name = "spiral"
        output_path = simulation_output_path
        #simulation = SkirtSimulation(prefix, outpath=output_path, ski_path=ski_path, name=name)

        # Settings
        settings_sed = dict()
        settings_sed["spectral_convolution"] = False

        # Input
        input_sed = dict()
        input_sed["simulation_output_path"] = simulation_output_path
        input_sed["output_path"] = "."

        # Construct the command
        create_sed = Command("observed_fluxes", "create the mock SED", settings_sed, input_sed, cwd=".")

        # Add the command
        #commands.append(create_sed)

        calculator = self.run_command(create_sed)

        # Determine the path to the mock SED
        mock_sed_path = "spiral_earth_fluxes.dat"

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # -----------------------------------------------------------------
        # SETUP THE MODELLING
        # -----------------------------------------------------------------

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "sed"
        settings_setup["name"] = "Spiral"
        settings_setup["fitting_host_ids"] = None

        # Create object config
        object_config = dict()
        object_config["ski"] = ski_path

        # Create input dict for setup
        input_setup = dict()
        input_setup["object_config"] = object_config
        input_setup["sed"] = mock_sed_path

        # Construct the command
        stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

        # Add the command
        commands.append(stp)

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Settings
        settings_model = dict()
        settings_model["ngenerations"] = 4
        settings_model["nsimulations"] = 20
        settings_model["fitting_settings"] = {"spectral_convolution": False}

        # Input

        # Get free parameter names
        ski = SkiFile(ski_path)
        free_parameter_names = ski.labels

        # Get fitting filter names
        #filter_names = sed.filter_names()

        # Set descriptions
        descriptions = Map()
        descriptions["exp_dustmass"] = "dust mass of the exponential disk with spiral structure"

        # Set types
        types = Map()
        types["exp_dustmass"] = "dust mass"

        # Set units
        units = Map()
        units["exp_dustmass"] = u("Msun")

        # Set the range of the dust mass
        dustmass_range = QuantityRange(0.1*dust_mass, 100*dust_mass)

        # Create input dict for model
        input_model = dict()
        input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
        input_model["descriptions_config"] = Configuration(descriptions=descriptions)
        input_model["types_config"] = Configuration(types=types)
        input_model["units_config"] = Configuration(units=units)
        input_model["ranges_config"] = Configuration(exp_dustmass_range=dustmass_range)
        #input_model["filters_config"] = Configuration(filters=filter_names)

        # Fitting initializer config
        input_model["initialize_config"] = Configuration(npackages=1e4)

        # Add dict of input for 'model' command to the list
        #input_dicts.append(input_model)

        # Construct the command
        command = Command("model", "perform the modelling", settings_model, input_model, "./Spiral")

        # Add the command
        #commands.append(command)

        modeler = self.run_command(command)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    """

    return

# -----------------------------------------------------------------
