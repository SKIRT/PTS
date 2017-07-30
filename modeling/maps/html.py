#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html Contains the AllMapsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import MapsComponent
from ..html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size
from ..core.environment import map_sub_names, colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ...core.tools import filesystem as fs
from ...core.tools.html import HTMLPage, SimpleTable

# -----------------------------------------------------------------

ncolumns = 2
colour_map = "jet"

# -----------------------------------------------------------------

class AllMapsPageGenerator(MapsComponent):

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
        super(AllMapsPageGenerator, self).__init__(*args, **kwargs)

        # Paths
        self.plots_path = None
        self.colour_plots_path = None
        self.ssfr_plots_path = None
        self.tir_plots_path = None
        self.attenuation_plots_path = None
        self.old_plots_path = None
        self.young_plots_path = None
        self.ionizing_plots_path = None
        self.dust_plots_path = None

        # The plots
        self.colour_plots = dict()
        self.ssfr_plots = dict()
        self.tir_plots = dict()
        self.attenuation_plots = dict()
        self.old_plots = dict()
        self.young_plots = dict()
        self.ionizing_plots = dict()
        self.dust_plots = dict()

        # The tables
        self.colour_table = None
        self.ssfr_table = None
        self.tir_table = None
        self.attenuation_table = None
        self.old_table = None
        self.young_table = None
        self.ionizing_table = None
        self.dust_table = None

        # The page
        self.page = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Make plots
        self.make_plots()

        # Make the tables
        self.make_tables()

        # Generate the page
        self.generate_page()

        # 5. Writing
        self.write()

        # Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AllMapsPageGenerator, self).setup(**kwargs)

        # Make directory to contain the plots
        self.plots_path = fs.join(self.maps_html_path, "plots")
        if fs.is_directory(self.plots_path): fs.clear_directory(self.plots_path)
        else: fs.create_directory(self.plots_path)

        # Create directories for each type of map
        self.colour_plots_path = fs.create_directory_in(self.plots_path, colours_name)
        self.ssfr_plots_path = fs.create_directory_in(self.plots_path, ssfr_name)
        self.tir_plots_path = fs.create_directory_in(self.plots_path, tir_name)
        self.attenuation_plots_path = fs.create_directory_in(self.plots_path, attenuation_name)
        self.old_plots_path = fs.create_directory_in(self.plots_path, old_name)
        self.young_plots_path = fs.create_directory_in(self.plots_path, young_name)
        self.ionizing_plots_path = fs.create_directory_in(self.plots_path, ionizing_name)
        self.dust_plots_path = fs.create_directory_in(self.plots_path, dust_name)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Maps"

    # -----------------------------------------------------------------

    @property
    def colour_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.colour_maps_flat

    # -----------------------------------------------------------------

    @property
    def ssfr_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.ssfr_maps_flat

    # -----------------------------------------------------------------

    @property
    def tir_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.tir_maps_flat

    # -----------------------------------------------------------------

    @property
    def attenuation_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.attenuation_maps_flat

    # -----------------------------------------------------------------

    @property
    def old_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.old_maps_flat

    # -----------------------------------------------------------------

    @property
    def young_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.young_maps_flat

    # -----------------------------------------------------------------

    @property
    def ionizing_maps(self):

        """
        Thisn function ...
        :return:
        """

        return self.static_collection.ionizing_maps_flat

    # -----------------------------------------------------------------

    @property
    def dust_maps(self):

        """
        This function ...
        :return:
        """

        return self.static_collection.dust_maps_flat

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

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make colours plots
        if self.has_colour_maps: self.make_colour_plots()

        # Make sSFR plots
        if self.has_ssfr_maps: self.make_ssfr_plots()

        # TIR
        if self.has_tir_maps: self.make_tir_plots()

        # Attenuation
        if self.has_attenuation_maps: self.make_attenuation_plots()

        # Old stellar maps
        if self.has_old_maps: self.make_old_plots()

        # Young stellar maps
        if self.has_young_maps: self.make_young_plots()

        # Ionizing stellar maps
        if self.has_ionizing_maps: self.make_ionizing_plots()

        # Dust maps
        if self.has_dust_maps: self.make_dust_plots()

    # -----------------------------------------------------------------

    def make_colour_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Determine path
            filepath = fs.join(self.colour_plots_path, name + ".png")

            # Save as PNG image
            self.colour_maps[name].saveto_png(filepath, colours=self.config.colours, absolute_alpha=True)

    # -----------------------------------------------------------------

    def make_ssfr_plots(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            # Determine path
            filepath = fs.join(self.ssfr_plots_path, name + ".png")

            # Save as PNG image
            self.ssfr_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_tir_plots(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Making plots of the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Determine path
            filepath = fs.join(self.tir_plots_path, name + ".png")

            # Save as PNG image
            self.tir_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_attenuation_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Determine path
            filepath = fs.join(self.attenuation_plots_path, name + ".png")

            # Save as PNG image
            self.attenuation_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_old_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Determine path
            filepath = fs.join(self.old_plots_path, name + ".png")

            # Save as PNG image
            self.old_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_young_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Determine path
            filepath = fs.join(self.young_plots_path, name + ".png")

            # Save as PNG image
            self.young_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_ionizing_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Determine path
            filepath = fs.join(self.ionizing_plots_path, name + ".png")

            # Save as PNG image
            self.ionizing_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_dust_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plot of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Determine path
            filepath = fs.join(self.dust_plots_path, name + ".png")

            # Save as PNG image
            self.dust_maps[name].saveto_png(filepath, colours=self.config.colours)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Colours
        if self.has_colour_maps: self.make_colour_table()

        # sSFR
        if self.has_ssfr_maps: self.make_ssfr_table()

        # TIR
        if self.has_tir_maps: self.make_tir_table()

        # Attenuation
        if self.has_attenuation_maps: self.make_attenuation_table()

        # Old stars
        if self.has_old_maps: self.make_old_table()

        # Young stars
        if self.has_young_maps: self.make_young_table()

        # Ionizing stars
        if self.has_ionizing_maps: self.make_ionizing_table()

        # Dust
        if self.has_dust_maps: self.make_dust_table()

    # -----------------------------------------------------------------

    def make_colour_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of colour maps ...")

        cells = []

        # Loop over the maps
        #for name, map in self.get_colour_maps(flatten=True).items():

        # Make
        self.colour_table = SimpleTable.rasterize(cells, ncolumns=ncolumns)

    # -----------------------------------------------------------------

    def make_ssfr_table(self):

        """

        :return:
        """

        # Inform the user
        log.info("Making the table of sSFR maps ...")

        # Make
        self.ssfr_table = SimpleTable.rasterize(cells, ncolumns=ncolumns)

    # -----------------------------------------------------------------

    def make_tir_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of TIR maps ...")

        # Make
        #self.tir_table =

    # -----------------------------------------------------------------

    def make_attenuation_table(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def make_old_table(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def make_young_table(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def make_ionizing_table(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def make_dust_table(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Create the page
        self.page = HTMLPage(self.title, style=page_style, css_path=stylesheet_url, footing=self.footing)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

# -----------------------------------------------------------------
