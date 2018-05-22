#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.solve.fitter Contains the Fitter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from .extremizer import Extremizer

# -----------------------------------------------------------------

def fit_function(function):

    """
    This function ...
    :param function: 
    :return: 
    """

    pass

# -----------------------------------------------------------------

def chi_squared(function, x, y, parameters):

    """
    This function ...
    :param function: 
    :param x: 
    :param y:
    :param parameters:
    :return: 
    """

    pass

# -----------------------------------------------------------------

class Fitter(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param config: 
        :param interactive: 
        """

        # Call the constructor of the base class
        super(Fitter, self).__init__(*args, **kwargs)

        # The function to be minimized or maximized
        self.function = None

        # The initial values and ranges
        self.initial_values = None
        self.ranges = None

        # The data
        self.data = None

        # The minimizer
        self.minimizer = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 2. Find minimum chi squared
        self.minimize()

        # 3. Plot
        self.plot()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(Fitter, self).setup(**kwargs)

        # Get input
        self.function = kwargs.pop("function")
        self.initial_values = kwargs.pop("initial_values")
        self.ranges = kwargs.pop("ranges")
        self.data = kwargs.pop("data")

    # -----------------------------------------------------------------

    def minimize(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Minimizing the chi-squared ...")

        # Settings
        settings = dict()
        settings["min_or_max"] = "min"

        # Input
        input_dict = dict()

        # Create the extremizer
        self.minimizer = Extremizer(settings)

        # Run the extremizer
        self.minimizer.run(**input_dict)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Plotting ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------
