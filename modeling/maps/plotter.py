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

        # Maps sub path
        self._maps_sub_path = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self._maps_sub_path

    # -----------------------------------------------------------------

    @property
    def colours(self):

        """
        This function ...
        :return:
        """

        return "colours" in self.config.types and self.has_colour_maps

    # -----------------------------------------------------------------

    @property
    def ssfr(self):

        """
        Thinsfunction ...
        :return:
        """

        return "ssfr" in self.config.types and self.has_ssfr_maps

    # -----------------------------------------------------------------

    @property
    def tir(self):

        """
        This function ...
        :return:
        """

        return "tir" in self.config.types and self.has_tir_maps

    # -----------------------------------------------------------------

    @property
    def attenuation(self):

        """
        This function ...
        :return:
        """

        return "attenuation" in self.config.types and self.has_attenuation_maps

    # -----------------------------------------------------------------

    @property
    def old(self):

        """
        Thisf ucntion ...
        :return:
        """

        return "old" in self.config.types and self.has_old_maps

    # -----------------------------------------------------------------

    @property
    def dust(self):

        """
        Thisf unction ...
        :return:
        """

        return "dust" in self.config.types and self.has_dust_maps

    # -----------------------------------------------------------------

    @property
    def young(self):

        """
        Thisf unction ...
        :return:
        """

        return "young" in self.config.types and self.has_young_maps

    # -----------------------------------------------------------------

    @property
    def ionizing(self):

        """
        Thisfunctio n...
        :return:
        """

        return "ionizing" in self.config.types and self.has_ionizing_maps

    # -----------------------------------------------------------------

    @property
    def maps_plotting(self):

        """
        This function ...
        :return:
        """

        return "maps" in self.config.features

    # -----------------------------------------------------------------

    @property
    def contours_plotting(self):

        """
        This function ...
        :return:
        """

        return "contours" in self.config.features

    # -----------------------------------------------------------------

    @property
    def profiles_plotting(self):

        """
        Thisf unction ...
        :return:
        """

        return "profiles" in self.config.features

    # -----------------------------------------------------------------

    @property
    def nans_plotting(self):

        """
        This function ...
        :return:
        """

        return "nans" in self.config.features

    # -----------------------------------------------------------------

    @property
    def negatives_plotting(self):

        """
        This funtion ...
        :return:
        """

        return "negatives" in self.config.features

    # -----------------------------------------------------------------

    @property
    def extra_maps_plotting(self):

        """
        This function ...
        :return:
        """

        return "extra" in self.config.features

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Colour maps
        if self.colours: self.process_colours()

        # 3. sSFR maps
        if self.ssfr: self.process_ssfr()

        # 4. TIR maps
        if self.tir: self.process_tir()

        # 5. Attenuation maps
        if self.attenuation: self.process_attenuation()

        # 6. Old stellar maps
        if self.old: self.process_old()

        # 7. Dust maps
        if self.dust: self.process_dust()

        # 8. Young stellar maps
        if self.young: self.process_young()

        # 9. Ionizing stellar maps
        if self.ionizing: self.process_ionizing()

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

        # Reset
        self._maps_sub_path = None

        # Clear
        self.maps = None
        self.origins = None
        self.methods = None
        self.extra_maps = None

    # -----------------------------------------------------------------

    def process_colours(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Processing colour maps ...")

        # Load the maps
        self.load_colours()

        # Plot the maps
        self.plot_colours()

        # Clear the maps
        self.clear()

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

        # Set path
        self._maps_sub_path = self.maps_colours_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.colours_scale, share_limits=False, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing sSFR maps ...")

        # Load maps
        self.load_ssfr()

        # Plot the maps
        self.plot_ssfr()

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

    @property
    def ssfr_cmap(self):

        """
        Thisn function ...
        :return:
        """

        return "jet"

    # -----------------------------------------------------------------

    @property
    def ssfr_interval(self):

        """
        This function ...
        :return:
        """

        return [0, 10]

    # -----------------------------------------------------------------

    @property
    def ssfr_min(self):

        """
        This function ...
        :return:
        """

        return 0

    # -----------------------------------------------------------------

    @property
    def ssfr_soft_max(self):

        """
        This function ...
        :return:
        """

        return 10

    # -----------------------------------------------------------------

    def plot_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR maps ...")

        # Set path
        self._maps_sub_path = self.maps_ssfr_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.ssfr_scale, cmap=self.ssfr_cmap, strict_vmin=self.ssfr_min, soft_vmax=self.ssfr_soft_max, share_limits=False, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing TIR maps ...")

        # Load the TIR maps
        self.load_tir()

        # Plot the TIR maps
        self.plot_tir()

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

    @property
    def tir_cmap(self):

        """
        Thisfunction ...
        :return:
        """

        return "viridis"

    # -----------------------------------------------------------------

    def plot_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the TIR maps ...")

        # Set path
        self._maps_sub_path = self.maps_tir_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.tir_scale, cmap=self.tir_cmap, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the NaNs masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the attenuation maps ...")

        # Load maps
        self.load_attenuation()

        # Plot maps
        self.plot_attenuation()

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

        # Get current extra maps
        self.extra_maps = self.get_attenuation_extra_maps()

    # -----------------------------------------------------------------

    @property
    def attenuation_scale(self):

        """
        This function ...
        :return:
        """

        return "linear"

    # -----------------------------------------------------------------

    @property
    def attenuation_cmap(self):

        """
        This function ...
        :return:
        """

        return "cool"

    # -----------------------------------------------------------------

    @property
    def tir_to_fuv_cmap(self):

        """
        This function ...
        :return:
        """

        return "ds9heat"

    # -----------------------------------------------------------------

    def plot_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the attenuation maps ...")

        # Set sub path
        self._maps_sub_path = self.maps_attenuation_path

        # Define name for extra maps
        self.extra_maps_name = "TIRtoFUV"

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.attenuation_scale, cmap=self.attenuation_cmap, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the extra maps
        if self.extra_maps_plotting: self.plot_extra_maps(scale="sqrt", format=self.config.format, clear_other_formats=True, cmap=self.tir_to_fuv_cmap)

        # Plot the NaNs masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing old stellar maps ...")

        # Load
        self.load_old()

        # Plot
        self.plot_old()

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

    @property
    def old_cmap(self):

        """
        This function ...
        :return:
        """

        #return "Wistia"
        return "afmhot"

    # -----------------------------------------------------------------

    def plot_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the old stellar maps ...")

        # Set path
        self._maps_sub_path = self.maps_old_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.old_scale, cmap=self.old_cmap, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the NaNs masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_dust(self):

        """
        This function ...
        :return:
        """

        # Inorm the user
        log.info("Processing dust maps ...")

        # Load dust maps
        self.load_dust()

        # Plot dust maps
        self.plot_dust()

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

    @property
    def dust_cmaps(self):

        """
        Thisf unction ...
        :return:
        """

        from .dust import blackbody, emission, attenuation, hot
        cmaps = dict()
        cmaps[blackbody] = "jet"
        cmaps[emission] = "jet"
        cmaps[attenuation] = "gist_ncar"
        cmaps[hot] = "jet"
        return cmaps

    # -----------------------------------------------------------------

    def plot_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust maps ...")

        # Set path
        self._maps_sub_path = self.maps_dust_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scales=self.dust_scales, cmaps=self.dust_cmaps, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the NaNs masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_young(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Processing young stellar maps ...")

        # Load
        self.load_young()

        # Plot
        self.plot_young()

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

    @property
    def young_cmap(self):

        """
        This function ...
        :return:
        """

        return "plasma"

    # -----------------------------------------------------------------

    def plot_young(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the young stellar maps ...")

        # Set path
        self._maps_sub_path = self.maps_young_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.young_scale, cmap=self.young_cmap, mask_negatives=True, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the negative pixel masks
        if self.negatives_plotting: self.plot_negatives(format=self.config.format, clear_other_formats=True)

        # Plot the NaN pixel masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

    # -----------------------------------------------------------------

    def process_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing ionizing stellar maps ...")

        # Load
        self.load_ionizing()

        # Plot
        self.plot_ionizing()

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

    @property
    def ionizing_cmap(self):

        """
        This function ...
        :return:
        """

        #return "ds9he"
        return "summer"

    # -----------------------------------------------------------------

    def plot_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the ionizing stellar maps ...")

        # Set sub path
        self._maps_sub_path = self.maps_ionizing_path

        # Plot the maps
        if self.maps_plotting: self.plot_maps(scale=self.ionizing_scale, cmap=self.ionizing_cmap, format=self.config.format, clear_other_formats=True)

        # Plot the contours
        if self.contours_plotting: self.plot_contours(filled=True, format=self.config.format, clear_other_formats=True)

        # Plot the radial profiles
        if self.profiles_plotting: self.plot_profiles(format=self.config.format, clear_other_formats=True)

        # Plot the NaNs masks
        if self.nans_plotting: self.plot_nans(format=self.config.format, clear_other_formats=True)

# -----------------------------------------------------------------
