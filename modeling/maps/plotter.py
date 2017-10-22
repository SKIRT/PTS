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

        # 2. Colour maps
        if self.has_colour_maps: self.load_colours()
        if self.has_colour_maps: self.plot_colours()

        # 3. sSFR maps
        if self.has_ssfr_maps: self.load_ssfr()
        if self.has_ssfr_maps: self.plot_ssfr()

        # 4. TIR maps
        if self.has_tir_maps: self.load_tir()
        if self.has_tir_maps: self.plot_tir()

        # 5. Attenuation maps
        if self.has_attenuation_maps: self.load_attenuation()
        if self.has_attenuation_maps: self.plot_attenuation()

        # 6. Old stellar maps
        if self.has_old_maps: self.load_old()
        if self.has_old_maps: self.plot_old()

        # 7. Dust maps
        if self.has_dust_maps: self.load_dust()
        if self.has_dust_maps: self.plot_dust()

        # 8. Young stellar maps
        if self.has_young_maps: self.load_young()
        if self.has_young_maps: self.plot_young()

        # 9. Ionizing stellar maps
        if self.has_ionizing_maps: self.load_ionizing()
        if self.has_ionizing_maps: self.plot_ionizing()

    # -----------------------------------------------------------------

    @property
    def has_colour_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_colour_maps

    # -----------------------------------------------------------------

    @property
    def has_ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_ssfr_maps

    # -----------------------------------------------------------------

    @property
    def has_tir_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_tir_maps

    # -----------------------------------------------------------------

    @property
    def has_attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_attenuation_maps

    # -----------------------------------------------------------------

    @property
    def has_old_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_old_maps

    # -----------------------------------------------------------------

    @property
    def has_young_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_young_maps

    # -----------------------------------------------------------------

    @property
    def has_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.has_ionizing_maps

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

    def clear(self):

        """
        This function ...
        :return:
        """

        self.maps = None
        self.origins = None
        self.methods = None

    # -----------------------------------------------------------------

    def load_colours(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the colour maps ...")

        # Load the maps
        self.maps = self.static_collection.colour_maps

        # Load the origins
        self.origins = self.get_colours_origins()

        # Load the methods
        self.origins = self.get_colours_methods()

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
        log.info("Plotting the colour maps ...")

        # Plot the maps
        self.plot_maps(scale=self.colours_scale, share_limits=False)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the sSFR maps ...")

        # Load the maps
        self.maps = self.static_collection.ssfr_maps

        # Load the origins
        self.origins = self.get_ssfr_origins()

        # Load the methods
        self.origins = self.get_ssfr_methods()

    # -----------------------------------------------------------------

    @property
    def ssfr_scale(self):

        """
        This function ...
        :return:
        """

        return "squared"

    # -----------------------------------------------------------------

    def plot_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR maps ...")

        # Plot the maps
        self.plot_maps(scale=self.ssfr_scale, share_limits=False)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the TIR maps ...")

        # Load the maps
        self.maps = self.static_collection.tir_maps

        # Load the origins
        self.origins = self.get_tir_origins()

        # Load the methods
        self.origins = self.get_tir_methods()

    # -----------------------------------------------------------------

    @property
    def tir_scale(self):

        """
        This function ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    def plot_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the TIR maps ...")

        # Plot the maps
        self.plot_maps(scale=self.tir_scale)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the NaNs masks
        self.plot_nans()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the attenuation maps ...")

        # Load the maps
        self.maps = self.static_collection.attenuation_maps

        # Load the origins
        self.origins = self.get_attenuation_origins()

        # Load the methods
        self.origins = self.get_attenuation_methods()

    # -----------------------------------------------------------------

    @property
    def attenuation_scale(self):

        """
        This function ...
        :return:
        """

        return "linear"

    # -----------------------------------------------------------------

    def plot_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the attenuation maps ...")

        # Plot the maps
        self.plot_maps(scale=self.attenuation_scale)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the extra maps
        self.plot_extra_maps(scale="sqrt")

        # Plot the NaNs masks
        self.plot_nans()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_old(self):

        """
        This function ...
        :return:
        """

        # Informt he user
        log.info("Loading the old stellar maps ...")

        # Load the maps
        self.maps = self.static_collection.old_maps

        # Load the origins
        self.origins = self.get_old_origins()

        # Load the methods
        self.origins = self.get_old_methods()

    # -----------------------------------------------------------------

    @property
    def old_scale(self):

        """
        Thisf unction ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    def plot_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the old stellar maps ...")

        # Plot the maps
        self.plot_maps(scale=self.old_scale)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the NaNs masks
        self.plot_nans()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust maps ...")

        # Load the maps
        self.maps = self.static_collection.dust_maps

        # Load the origins
        self.origins = self.get_dust_origins()

        # Load the methods
        self.origins = self.get_dust_methods()

    # -----------------------------------------------------------------

    @property
    def dust_scales(self):

        """
        Thisf unction ...
        :return:
        """

        from .dust import blackbody, emission, attenuation, hot
        scales = dict()
        scales[blackbody] = "linear"
        scales[emission] = "linear"
        scales[attenuation] = "linear"
        scales[hot] = "log"
        return scales

    # -----------------------------------------------------------------

    def plot_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust maps ...")

        # Plot the maps
        self.plot_maps(scales=self.dust_scales)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the NaNs masks
        self.plot_nans()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the young stellar maps ...")

        # Load the maps
        self.maps = self.static_collection.young_maps

        # Load the origins
        self.origins = self.get_young_origins()

        # Load the methods
        self.origins = self.get_young_methods()

    # -----------------------------------------------------------------

    @property
    def young_scale(self):

        """
        This function ...
        :return:
        """

        return "sqrt"

    # -----------------------------------------------------------------

    def plot_young(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the young stellar maps ...")

        # Plot the maps
        self.plot_maps(scale=self.young_scale, mask_negatives=True)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the negative pixel masks
        self.plot_negatives()

        # Plot the NaN pixel masks
        self.plot_nans()

        # Clear the maps
        self.clear()

    # -----------------------------------------------------------------

    def load_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ionizing stellar maps ...")

        # Load the maps
        self.maps = self.static_collection.ionizing_maps

        # Load the origins
        self.origins = self.get_ionizing_origins()

        # Load the methods
        self.origins = self.get_ionizing_methods()

    # -----------------------------------------------------------------

    @property
    def ionizing_scale(self):

        """
        Thisfunction ...
        :return:
        """

        return "log"

    # -----------------------------------------------------------------

    def plot_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the ionizing stellar maps ...")

        # Plot the maps
        self.plot_maps(scale=self.ionizing_scale)

        # Plot the contours
        self.plot_contours(filled=True)

        # Plot the radial profiles
        self.plot_profiles()

        # Plot the NaNs masks
        self.plot_nans()

        # Clear the maps
        self.clear()

# -----------------------------------------------------------------
