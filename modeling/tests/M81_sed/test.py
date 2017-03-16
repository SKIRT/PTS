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
from pts.core.test.implementation import TestImplementation
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.core.tools.logging import log
from pts.modeling.basics.instruments import SEDInstrument
from pts.modeling.tests.base import M81TestBase, fitting_filter_names

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting the parameters of a model of M81 based on only a mock observed SED"

# -----------------------------------------------------------------

class M81SEDTest(M81TestBase):

    """
    This class runs the test on M81, but by only adjusting the normalizations (not by creating a model),
    and fitting to a mock observed SED
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(M81SEDTest, self).__init__(config, interactive)

        # The galaxy properties
        self.properties = None

        # Bulge and disk model
        self.bulge = None
        self.disk = None

        # The input maps
        self.maps = dict()

        # The instrument
        self.instrument = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the properties
        self.load_properties()

        # 3. Load the components
        self.load_components()

        # Create instrument
        self.create_instrument()

        self.create_deprojections()

        self.write()

        self.plot()

        self.launch_reference()

        self.setup_modelling()

        self.model()

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the SED instrument ...")

        azimuth = 0.

        # Create the SED instrument
        self.instrument = SEDInstrument(distance=self.galaxy_distance, inclination=self.galaxy_inclination, azimuth=azimuth, position_angle=self.galaxy_position_angle)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski()

        # Write the input
        self.write_input()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the wavelengths
        self.plot_wavelengths()

        # Plot the filters
        self.plot_filters()

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Inform the user

# -----------------------------------------------------------------
