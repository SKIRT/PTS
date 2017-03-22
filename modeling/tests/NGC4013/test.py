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
from pts.core.tools import network
from pts.core.tools.logging import log
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "Fitting the galaxy NGC4013 using a flattened Sersic profile for the central bulge and a double exponential for the stellar disk and dust disk"

# -----------------------------------------------------------------

input_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013.tar.gz"
all_url = "http://www.skirt.ugent.be/downloads/tutorial_NGC4013_complete.tar.gz"

# -----------------------------------------------------------------

ski_path = fs.join(this_dir_path, "ngc4013.ski")
fski_path = fs.join(this_dir_path, "ngc4013.fski")

# -----------------------------------------------------------------

class NGC4013Test(TestImplementation):

    """
    This class ...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(NGC4013Test, self).__init__(config, interactive)

        # The path with the test data
        self.data_path = None

        # The path for the fitskirt run
        self.reference_path = None
        self.reference_output_path = None
        self.reference_ski_path = None
        self.reference_fski_path = None

        # The FitSKIRT launcher
        self.launcher = None

        # The modeler
        self.modeler = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the data
        self.get_data()

        # Write
        self.write()

        # 3. Launch with FitSKIRT
        self.launch_fitskirt()

        # Setup the modeling
        self.setup_modelling()

        # Model
        self.model()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # CAll the setup function of the base class
        super(NGC4013Test, self).setup(**kwargs)

        # Create the data directory
        self.data_path = fs.create_directory_in(self.path, "data")

        # Create the reference directory and subdirectories
        self.reference_path = fs.create_directory_in(self.path, "ref")
        self.reference_output_path = fs.create_directory_in(self.reference_path, "out")
        self.reference_ski_path = fs.join(self.data_path, "NGC4013.ski")
        self.reference_fski_path = fs.join(self.reference_path, "NGC4013.fski")

    # -----------------------------------------------------------------

    def get_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Downloading the input ...")

        # Download the input
        network.download_and_decompress_directory(all_url, self.data_path, progress_bar=True, into_root=True)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write ski
        self.write_ski()

        # Write fski
        self.write_fski()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        fs.copy_file(ski_path, self.reference_ski_path)

    # -----------------------------------------------------------------

    def write_fski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fski file ...")

        fs.copy_file(fski_path, self.reference_fski_path)

    # -----------------------------------------------------------------

    def launch_fitskirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching FitSKIRT ...")

        # Create configuration
        config = dict()

        config["ski"] = self.reference_ski_path
        config["fski"] = self.reference_fski_path

        config["input"] = self.data_path
        config["output"] = self.reference_output_path

        # Input
        input_dict = dict()

        command = Command("fitskirt", "run the reference fitting with FitSKIRT", config, input_dict, cwd=".")

        # Run the command
        self.launcher = self.run_command(command)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # -----------------------------------------------------------------
        # SETUP MODELLING
        # -----------------------------------------------------------------

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "images"
        settings_setup["name"] = "NGC4013"
        settings_setup["fitting_host_ids"] = None

        # -----------------------------------------------------------------

        # Input

        # Create object config
        object_config = dict()
        ski_path = fs.join(this_dir_path, "NGC4013.ski")
        object_config["ski"] = ski_path

        # Create input dict for setup
        input_setup = dict()
        input_setup["object_config"] = object_config
        #input_setup["sed"] = sed

        # Construct the command
        setup_command = Command("setup", "setup the modelling", settings_setup, input_setup, ".")

        # Add the command
        #commands.append(setup_command)

        tool = self.run_command(setup_command)

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Settings for 'model_sed'
        settings_model = dict()
        settings_model["ngenerations"] = "4"

        # -----------------------------------------------------------------

        # Create input dict for model
        input_model = dict()
        #input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
        #input_model["descriptions_config"] = Configuration(descriptions=descriptions)
        #input_model["types_config"] = Configuration(types=types)
        #input_model["units_config"] = Configuration(units=units)
        #input_model["ranges_config"] = Configuration(luminosity_range=luminosity_range, dustmass_range=dustmass_range, grainsize_range=grainsize_range, fsil_range=fsil_range)
        #input_model["filters_config"] = Configuration(filters=filter_names)

        # Fitting initializer config
        #input_model["initialize_config"] = Configuration(npackages=1e4)

        # Construct the command
        model_command = Command("model", "perform the modelling", settings_model, input_model, cwd="./NGC4013")

        # Add the command
        #commands.append(model_command)

        self.modeler = self.run_command(model_command)

# -----------------------------------------------------------------
# TEST FUNCTION
# -----------------------------------------------------------------

def test(temp_path):

    """
    This function ...
    :param temp_path:
    """

    return

# -----------------------------------------------------------------
