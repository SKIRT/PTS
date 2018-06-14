#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.viewer Contains the MapsViewer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapMakingComponent
from ...magic.view.multi import MultiImageViewer

# -----------------------------------------------------------------

ncolumns = 3
image_width = 300
image_height = 300
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

class MapsViewer(MapMakingComponent):

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
        super(MapsViewer, self).__init__(*args, **kwargs)

        # Maps sub path
        self._maps_sub_path = None

        # The map paths
        self.paths = None

        # The viewer
        self.viewer = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self._maps_sub_path

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the map paths
        self.load_map_paths()

        # 3. View
        self.view()

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
        super(MapsViewer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Reset
        self._maps_sub_path = None

    # -----------------------------------------------------------------

    @property
    def colours(self):

        """
        This function ...
        :return:
        """

        return self.config.which == "colours"

    # -----------------------------------------------------------------

    @property
    def ssfr(self):

        """
        Thsifunction ...
        :return:
        """

        return self.config.which == "ssfr"

    # -----------------------------------------------------------------

    @property
    def tir(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.which == "tir"

    # -----------------------------------------------------------------

    @property
    def attenuation(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.which == "attenuation"

    # -----------------------------------------------------------------

    @property
    def dust(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.config.which == "dust"

    # -----------------------------------------------------------------

    @property
    def old(self):

        """
        This function ...
        :return:
        """

        return self.config.which == "old"

    # -----------------------------------------------------------------

    @property
    def young(self):

        """
        Thisfunction ...
        :return:
        """

        return self.config.which == "young"

    # -----------------------------------------------------------------

    @property
    def ionizing(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.which == "ionizing"

    # -----------------------------------------------------------------

    @property
    def components(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.which == "components"

    # -----------------------------------------------------------------

    def load_map_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map paths ...")

        # Colour maps
        if self.colours: self.paths = self.get_colour_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # sSFR maps
        elif self.ssfr: self.paths = self.get_ssfr_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # TIR maps
        elif self.tir: self.paths = self.get_tir_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Attenuation maps
        elif self.attenuation: self.paths = self.get_attenuation_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Old maps
        elif self.old: self.paths = self.get_old_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Dust maps
        elif self.dust: self.paths = self.get_dust_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Young maps
        elif self.young: self.paths = self.get_young_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Ionizing stellar maps
        elif self.ionizing: self.paths = self.get_ionizing_map_paths(flatten=True, method=self.config.method, startswith=self.config.startswith, factors=self.config.factors)

        # Component maps
        elif self.components: self.paths = self.get_component_map_paths(flatten=True)

        # Invalid
        else: raise ValueError("Invalid value for 'which': '" + self.config.which + "'")

    # -----------------------------------------------------------------

    @property
    def scale(self):

        """
        This function ...
        :return:
        """

        if self.config.scale is not None: return self.config.scale
        else:

            if self.colours: return self.colours_scale
            elif self.ssfr: return self.ssfr_scale
            elif self.tir: return self.tir_scale
            elif self.attenuation: return self.attenuation_scale
            elif self.old: return self.old_scale
            elif self.dust: return self.dust_scales["attenuation"]
            elif self.young: return self.young_scale
            elif self.ionizing: return self.ionizing_scale
            elif self.components: return

    # -----------------------------------------------------------------

    @property
    def colormap(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.colormap is not None: return self.config.colormap
        else:

            if self.colours: return self.colours_js9_cmap
            elif self.ssfr: return self.ssfr_js9_cmap
            elif self.tir: return self.tir_js9_cmap
            elif self.attenuation: return self.attenuation_js9_cmap
            elif self.old: return self.old_js9_cmap
            elif self.dust: return self.dust_js9_cmap
            elif self.young: return self.young_js9_cmap
            elif self.ionizing: return self.ionizing_js9_cmap
            elif self.components: return

    # -----------------------------------------------------------------

    @property
    def zoom(self):

        """
        This function ...
        :return:
        """

        return self.config.zoom

    # -----------------------------------------------------------------

    def view(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Viewing ...")

        # Initialize the viewer
        self.viewer = MultiImageViewer()

        # Set settings
        self.viewer.config.scale = self.scale
        self.viewer.config.colormap = self.colormap
        self.viewer.config.zoom = self.zoom

        # Set advanced settings
        self.viewer.config.preload_all = self.config.preload_all
        self.viewer.config.preload = self.config.preload
        self.viewer.config.dynamic = self.config.dynamic

        # Region?
        if self.config.truncation_ellipse: region = self.truncation_ellipse
        else: region = None

        # Run the viewer
        self.viewer.run(paths=self.paths, region=region)

# -----------------------------------------------------------------
