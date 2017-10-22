#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.plotter Contains the MapsPlotter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent

# -----------------------------------------------------------------

class MapsPlotter(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(MapsPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        self.plot_colours()

        self.plot_ssfr()

        self.plot_tir()

        self.plot_attenuation()

        self.plot_old()

        self.plot_dust()

        self.plot_young()

        self.plot_ionizing()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsPlotter, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def colours_scale(self):

        """
        This function ...
        :return:
        """

        return "squared"

    # -----------------------------------------------------------------

    def plot_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps(scale=self.colours_scale, share_limits=False)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

    # -----------------------------------------------------------------

    def plot_ssfr(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_tir(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_attenuation(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_old(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_dust(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_young(self):

        """
        Thins function ...
        :return:
        """

    # -----------------------------------------------------------------

    def plot_ionizing(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------
