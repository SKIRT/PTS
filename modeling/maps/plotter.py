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
from .component import MapMakingComponent
from .selectioncomponent import MapsSelectionComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class MapsPlotter(MapMakingComponent):

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

        # Set path
        self._maps_sub_path = self.maps_colours_path

        # Load the maps
        self.load_colours()

        # Plot the maps
        self.plot_colours(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                          format=self.config.format)

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

    def process_ssfr(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing sSFR maps ...")

        # Set path
        self._maps_sub_path = self.maps_ssfr_path

        # Load maps
        self.load_ssfr()

        # Plot the maps
        self.plot_ssfr(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                       format=self.config.format)

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

    def process_tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing TIR maps ...")

        # Set path
        self._maps_sub_path = self.maps_tir_path

        # Load the TIR maps
        self.load_tir()

        # Plot the TIR maps
        self.plot_tir(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                      nans=self.nans_plotting, format=self.config.format)

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

    def process_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the attenuation maps ...")

        # Set sub path
        self._maps_sub_path = self.maps_attenuation_path

        # Define name for extra maps
        self.extra_maps_name = "TIRtoFUV"

        # Load maps
        self.load_attenuation()

        # Plot maps
        self.plot_attenuation(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                              extra=self.extra_maps_plotting, nans=self.nans_plotting, format=self.config.format)

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

    def process_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing old stellar maps ...")

        # Set path
        self._maps_sub_path = self.maps_old_path

        # Load
        self.load_old()

        # Plot
        self.plot_old(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                      negatives=self.negatives_plotting, nans=self.nans_plotting, format=self.config.format)

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

    def process_dust(self):

        """
        This function ...
        :return:
        """

        # Inorm the user
        log.info("Processing dust maps ...")

        # Set path
        self._maps_sub_path = self.maps_dust_path

        # Load dust maps
        self.load_dust()

        # Plot dust maps
        self.plot_dust(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                       negatives=self.negatives_plotting, nans=self.nans_plotting, format=self.config.format)

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

    def process_young(self):

        """
        This function ...
        :return:
        """

        # Inform the suer
        log.info("Processing young stellar maps ...")

        # Set path
        self._maps_sub_path = self.maps_young_path

        # Load
        self.load_young()

        # Plot
        self.plot_young(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                        negatives=self.negatives_plotting, nans=self.nans_plotting, format=self.config.format)

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

    def process_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing ionizing stellar maps ...")

        # Set sub path
        self._maps_sub_path = self.maps_ionizing_path

        # Define name for extra maps
        self.extra_maps_name = "HalphaToHot"

        # Load
        self.load_ionizing()

        # Plot
        self.plot_ionizing(maps=self.maps_plotting, contours=self.contours_plotting, profiles=self.profiles_plotting,
                           extra=self.extra_maps_plotting, negatives=self.negatives_plotting, nans=self.nans_plotting,
                           format=self.config.format)

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

class ComponentMapsPlotter(MapsSelectionComponent):

    """
    This class ..
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ComponentMapsPlotter, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load
        self.load()

        # 3. Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsPlotter, self).setup(**kwargs)

        # Replot
        if self.config.replot:

            from .selectioncomponent import map_plot_filename, edgeon_plot_filename, deprojected_plot_filename, deprojected_skirt_plot_filename

            if self.plot_old:

                if self.maps_plotting: fs.remove_files_in_path(self.old_plots_path, exact_name=map_plot_filename, recursive=True)
                if self.deprojected_plotting: fs.remove_files_in_path(self.old_plots_path, exact_name=deprojected_plot_filename, recursive=True)
                if self.deprojected_skirt_plotting: fs.remove_files_in_path(self.old_plots_path, exact_name=deprojected_skirt_plot_filename, recursive=True)
                if self.edgeon_plotting: fs.remove_files_in_path(self.old_plots_path, exact_name=edgeon_plot_filename, recursive=True)

            if self.plot_young:

                if self.maps_plotting: fs.remove_files_in_path(self.young_plots_path, exact_name=map_plot_filename, recursive=True)
                if self.deprojected_plotting: fs.remove_files_in_path(self.young_plots_path, exact_name=deprojected_plot_filename, recursive=True)
                if self.deprojected_skirt_plotting: fs.remove_files_in_path(self.young_plots_path, exact_name=deprojected_skirt_plot_filename, recursive=True)
                if self.edgeon_plotting: fs.remove_files_in_path(self.young_plots_path, exact_name=edgeon_plot_filename, recursive=True)

            if self.plot_ionizing:

                if self.maps_plotting: fs.remove_files_in_path(self.ionizing_plots_path, exact_name=map_plot_filename, recursive=True)
                if self.deprojected_plotting: fs.remove_files_in_path(self.ionizing_plots_path, exact_name=deprojected_plot_filename, recursive=True)
                if self.deprojected_skirt_plotting: fs.remove_files_in_path(self.ionizing_plots_path, exact_name=deprojected_skirt_plot_filename, recursive=True)
                if self.edgeon_plotting: fs.remove_files_in_path(self.ionizing_plots_path, exact_name=edgeon_plot_filename, recursive=True)

            if self.plot_dust:

                if self.maps_plotting: fs.remove_files_in_path(self.dust_plots_path, exact_name=map_plot_filename, recursive=True)
                if self.deprojected_plotting: fs.remove_files_in_path(self.dust_plots_path, exact_name=deprojected_plot_filename, recursive=True)
                if self.deprojected_skirt_plotting: fs.remove_files_in_path(self.dust_plots_path, exact_name=deprojected_skirt_plot_filename, recursive=True)
                if self.edgeon_plotting: fs.remove_files_in_path(self.dust_plots_path, exact_name=edgeon_plot_filename, recursive=True)

    # -----------------------------------------------------------------

    @property
    def plot_old(self):

        """
        Thisfunction ...
        :return:
        """

        return "old" in self.config.types

    # -----------------------------------------------------------------

    @property
    def plot_young(self):

        """
        Thisf unction ...
        :return:
        """

        return "young" in self.config.types

    # -----------------------------------------------------------------

    @property
    def plot_ionizing(self):

        """
        This function ...
        :return:
        """

        return "ionizing" in self.config.types

    # -----------------------------------------------------------------

    @property
    def plot_dust(self):

        """
        Thisf unction ...
        :return:
        """

        return "dust" in self.config.types

    # -----------------------------------------------------------------

    def load(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading input ...")

        # Load the maps
        if self.maps_plotting: self.load_maps()

        # Load the deprojected maps
        if self.deprojected_plotting: self.load_deprojected()

        # Deprojected with SKIRT
        if self.deprojected_skirt_plotting: self.load_deprojected_skirt()

        # Edgeon
        if self.edgeon_plotting: self.load_edgeon()

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Old
        if self.plot_old: self.load_old_maps()

        # Young
        if self.plot_young: self.load_young_maps()

        # Ionizing
        if self.plot_ionizing: self.load_ionizing_maps()

        # Dust
        if self.plot_dust: self.load_dust_maps()

    # -----------------------------------------------------------------

    def load_old_maps(self):

        """
        Thisfunction ...
        :return:
        """

        # Load
        self.old_maps = self.get_component_old_maps()

    # -----------------------------------------------------------------

    def load_young_maps(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_maps = self.get_component_young_maps()

    # -----------------------------------------------------------------

    def load_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_maps = self.get_component_ionizing_maps()

    # -----------------------------------------------------------------

    def load_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Load
        self.dust_maps = self.get_component_dust_maps()

    # -----------------------------------------------------------------

    def load_deprojected(self):

        """
        This function ...
        :return:
        """

        # Old
        if self.plot_old: self.load_old_deprojected()

        # Young
        if self.plot_young: self.load_young_deprojected()

        # Ionizing
        if self.plot_ionizing: self.load_ionizing_deprojected()

        # Dust
        if self.plot_dust: self.load_dust_deprojected()

    # -----------------------------------------------------------------

    def load_old_deprojected(self):

        """
        This function ...
        :return:
        """

        # Load
        self.old_deprojected = self.get_component_old_deprojected()

    # -----------------------------------------------------------------

    def load_young_deprojected(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_deprojected = self.get_component_young_deprojected()

    # -----------------------------------------------------------------

    def load_ionizing_deprojected(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_deprojected = self.get_component_ionizing_deprojected()

    # -----------------------------------------------------------------

    def load_dust_deprojected(self):

        """
        This function ...
        :return:
        """

        # Load
        self.dust_deprojected = self.get_component_dust_deprojected()

    # -----------------------------------------------------------------

    def load_deprojected_skirt(self):

        """
        Thisf unction ...
        :return:
        """

        # Old
        if self.plot_old: self.load_old_deprojected_skirt()

        # Young
        if self.plot_young: self.load_young_deprojected_skirt()

        # Ionizing
        if self.plot_ionizing: self.load_ionizing_deprojected_skirt()

        # Dust
        if self.plot_dust: self.load_dust_deprojected_skirt()

    # -----------------------------------------------------------------

    def load_old_deprojected_skirt(self):

        """
        Thisf unction ...
        :return:
        """

        # Load
        self.old_deprojected_skirt = self.get_component_old_deprojected_skirt()

    # -----------------------------------------------------------------

    def load_young_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_deprojected_skirt = self.get_component_young_deprojected_skirt()

    # -----------------------------------------------------------------

    def load_ionizing_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_deprojected_skirt = self.get_component_ionizing_deprojected_skirt()

    # -----------------------------------------------------------------

    def load_dust_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Load
        self.dust_deprojected_skirt = self.get_component_dust_deprojected_skirt()

    # -----------------------------------------------------------------

    def load_edgeon(self):

        """
        This function ...
        :return:
        """

        # Old
        if self.plot_old: self.load_old_edgeon()

        # Young
        if self.plot_young: self.load_young_edgeon()

        # Ionizing
        if self.plot_ionizing: self.load_ionizing_edgeon()

        # Dust
        if self.plot_dust: self.load_dust_edgeon()

    # -----------------------------------------------------------------

    def load_old_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.old_edgeon_skirt = self.get_component_old_edgeon()

    # -----------------------------------------------------------------

    def load_young_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.young_edgeon_skirt = self.get_component_young_edgeon()

    # -----------------------------------------------------------------

    def load_ionizing_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.ionizing_edgeon_skirt = self.get_component_ionizing_edgeon()

    # -----------------------------------------------------------------

    def load_dust_edgeon(self):

        """
        This function ...
        :return:
        """

        # Load
        self.dust_edgeon_skirt = self.get_component_dust_edgeon()

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
    def masks_plotting(self):

        """
        This function ...
        :return:
        """

        return "masks" in self.config.features

    # -----------------------------------------------------------------

    @property
    def deprojected_plotting(self):

        """
        Thisf unction ...
        :return:
        """

        return "deprojected" in self.config.features

    # -----------------------------------------------------------------

    @property
    def deprojected_skirt_plotting(self):

        """
        Thisf unction ...
        :return:
        """

        return "deprojected_skirt" in self.config.features

    # -----------------------------------------------------------------

    @property
    def edgeon_plotting(self):

        """
        This function ...
        :return:
        """

        return "edgeon" in self.config.features

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps
        if self.maps_plotting: self.plot_components_maps(format=self.config.format)

        # Masks
        if self.masks_plotting: self.plot_components_masks(format=self.config.format)

        # Deprojected
        if self.deprojected_plotting: self.plot_components_deprojected(format=self.config.format)

        # Deprojected with SKIRT
        if self.deprojected_skirt_plotting: self.plot_components_deprojected_skirt(format=self.config.format)

        # Edgeon
        if self.edgeon_plotting: self.plot_components_edgeon(format=self.config.format)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        Thisf unction ...
        :return:
        """

        return None

# -----------------------------------------------------------------
