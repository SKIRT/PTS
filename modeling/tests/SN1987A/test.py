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
from pts.core.data.sed import ObservedSED
from pts.core.units.parsing import parse_unit as u
from pts.core.basics.configuration import Configuration
from pts.core.simulation.skifile import SkiFile
from pts.core.basics.range import QuantityRange, RealRange
from pts.core.basics.map import Map
from pts.do.commandline import Command
from pts.core.test.implementation import TestImplementation
from pts.core.basics.log import log

# -----------------------------------------------------------------

this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "modeling of the Supernova 1987A with genetic algorithms and 4 free parameters"

# -----------------------------------------------------------------

class SN1987ATest(TestImplementation):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SN1987ATest, self).__init__(*args, **kwargs)

        # The modeler
        self.modeler = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup the modelling
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

        # Call the setup function of the base class
        super(SN1987ATest, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting up the modeling ...")

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "sed"
        settings_setup["name"] = "SN1987A"
        settings_setup["fitting_host_ids"] = None

        # -----------------------------------------------------------------

        # Input

        # Construct the observed SED
        sed = ObservedSED(photometry_unit="Jy")
        sed.add_point("Pacs 70", 0.0455 * u("Jy"), 0.0034 * u("Jy"))
        sed.add_point("Pacs 100", 0.0824 * u("Jy"), 0.0045 * u("Jy"))
        sed.add_point("Pacs 160", 0.1530 * u("Jy"), 0.0090 * u("Jy"))
        sed.add_point("SPIRE 250", 0.1107 * u("Jy"), 0.0252 * u("Jy"))
        sed.add_point("SPIRE 350", 0.0693 * u("Jy"), 0.0228 * u("Jy"))
        sed.add_point("ALMA 440um", 0.0500 * u("Jy"), 0.0150 * u("Jy"))
        sed.add_point("ALMA 870um", 0.0050 * u("Jy"), 0.0010 * u("Jy"))

        # Create object config
        object_config = dict()
        ski_path = fs.join(this_dir_path, "SN1987A.ski")
        object_config["ski"] = ski_path

        # Create input dict for setup
        input_setup = dict()
        input_setup["object_config"] = object_config
        input_setup["sed"] = sed

        # Construct the command
        setup_command = Command("setup", "setup the modelling", settings_setup, input_setup, cwd=".")

        # Add the command
        #commands.append(setup_command)

        tool = self.run_command(setup_command)

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the modelling ...")

        # Settings
        settings_model = dict()
        settings_model["ngenerations"] = 4
        settings_model["nsimulations"] = 20
        settings_model["fitting_settings"] = {"spectral_convolution": False}

        # -----------------------------------------------------------------

        # Input

        # Get free parameter names
        ski = SkiFile(ski_path)
        free_parameter_names = ski.labels

        # Get fitting filter names
        filter_names = sed.filter_names()

        # Set descriptions
        descriptions = Map()
        descriptions["luminosity"] = "total luminosity of the SN"
        descriptions["dustmass"] = "total dust mass"
        descriptions["grainsize"] = "dust grain size"
        descriptions["fsil"] = "dust silicate fraction"

        # Set types
        types = Map()
        types["luminosity"] = "luminosity"
        types["dustmas"] = "mass"
        types["grainsize"] = "grainsize"
        types["fsil"] = "dimless"

        # Set units
        units = Map()
        units["luminosity"] = u("Lsun")
        units["dustmass"] = u("Msun")
        units["grainsize"] = u("micron")
        units["fsil"] = None

        # Set ranges
        luminosity_range = QuantityRange(100, 1000, "Lsun")
        dustmass_range = QuantityRange(0.3, 5, "Msun")
        grainsize_range = QuantityRange(0.1, 5, "micron")
        fsil_range = RealRange(0.1, 100)

        # Create input dict for model
        input_model = dict()
        input_model["parameters_config"] = Configuration(free_parameters=free_parameter_names)
        input_model["descriptions_config"] = Configuration(descriptions=descriptions)
        input_model["types_config"] = Configuration(types=types)
        input_model["units_config"] = Configuration(units=units)
        input_model["ranges_config"] = Configuration(luminosity_range=luminosity_range, dustmass_range=dustmass_range, grainsize_range=grainsize_range, fsil_range=fsil_range)
        input_model["filters_config"] = Configuration(filters=filter_names)

        # Fitting initializer config
        input_model["initialize_config"] = Configuration(npackages=1e4)

        # Construct the command
        model_command = Command("model", "perform the modelling", settings_model, input_model, "./SN1987A")

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
