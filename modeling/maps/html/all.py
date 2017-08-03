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
from ....core.tools.logging import log
from ..component import MapsComponent
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size
from ...core.environment import map_sub_names, colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ....core.tools import filesystem as fs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts
from ....core.tools import browser
from ....core.tools.stringify import tostr
from ....core.tools.utils import lazyproperty
from ....core.tools import numbers
from ....core.basics.range import RealRange

# -----------------------------------------------------------------

plots_name = "plots"
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

        # The image info
        self.colour_info = dict()
        self.ssfr_info = dict()
        self.tir_info = dict()
        self.attenuation_info = dict()
        self.old_info = dict()
        self.young_info = dict()
        self.ionizing_info = dict()
        self.dust_info = dict()

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

        # Get image info
        if self.config.info: self.get_info()

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
        self.plots_path = fs.join(self.maps_html_path, plots_name)
        if fs.is_directory(self.plots_path):
            if self.config.replot: fs.clear_directory(self.plots_path)
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
    def image_width(self):

        """
        This fucntion ...
        :return:
        """

        #return 150
        return None

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        This function ...
        :return:
        """

        return 300

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

    def get_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting image info ...")

        # Colour
        self.get_colour_info()

        # SSFR
        self.get_ssfr_info()

        # TIR
        self.get_tir_info()

        # Attenuation
        self.get_attenuation_info()

        # Old
        self.get_old_info()

        # Young
        self.get_young_info()

        # Ionizing
        self.get_ionizing_info()

        # Dut
        self.get_dust_info()

    # -----------------------------------------------------------------

    def get_colour_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Get info
            info = get_image_info(self.colour_maps[name])

            # Make list
            code = html.unordered_list(info)

            # Add info
            self.colour_info[name] = code

    # -----------------------------------------------------------------

    def get_ssfr_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the sSFR maps ...")

    # -----------------------------------------------------------------

    def get_tir_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the TIR maps ...")

    # -----------------------------------------------------------------

    def get_attenuation_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the attenuation maps ...")

    # -----------------------------------------------------------------

    def get_old_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the old stellar maps ...")

    # -----------------------------------------------------------------

    def get_young_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the young stellar maps ...")

    # -----------------------------------------------------------------

    def get_ionizing_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the ionizing stellar maps ...")

    # -----------------------------------------------------------------

    def get_dust_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the dust maps ...")

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.softening_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_radius(self):

        """
        This function ...
        :return:
        """

        return numbers.geometric_mean(self.config.softening_start, 1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(self.config.softening_start / self.softening_radius, 1. / self.softening_radius)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.colour_maps[name].wcs, self.colour_maps[name].xsize, self.colour_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.colour_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.colour_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.colour_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.ssfr_maps[name].wcs, self.ssfr_maps[name].xsize, self.ssfr_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.ssfr_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.ssfr_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.ssfr_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.tir_maps[name].wcs, self.tir_maps[name].xsize, self.tir_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.tir_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.tir_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.tir_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.attenuation_maps[name].wcs, self.attenuation_maps[name].xsize, self.attenuation_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.attenuation_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.attenuation_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.attenuation_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.old_maps[name].wcs, self.old_maps[name].xsize, self.old_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.old_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.old_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.old_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.young_maps[name].wcs, self.young_maps[name].xsize, self.young_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.young_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.young_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.young_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.ionizing_maps[name].wcs, self.ionizing_maps[name].xsize, self.ionizing_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.ionizing_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.ionizing_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.ionizing_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

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

            # Check if plot is already made
            if fs.is_file(filepath):

                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Get the truncation mask and mask out the pixel beyond the truncation limit
            wcs, xsize, ysize = self.dust_maps[name].wcs, self.dust_maps[name].xsize, self.dust_maps[name].ysize
            ellipse = self.truncation_ellipse.to_pixel(wcs)
            mask = ellipse.to_mask(xsize, ysize).inverse()
            self.dust_maps[name][mask] = 0.0

            # Make RGBA image
            rgba = self.dust_maps[name].to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
            rgba.soften_edges(self.softening_ellipse.to_pixel(self.dust_maps[name].wcs), self.softening_range)

            # Save
            rgba.saveto(filepath)

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "realtable"

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
        for name in self.colour_maps:

            # Determine the relative path
            path = fs.join(plots_name, colours_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.colour_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_ssfr_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of sSFR maps ...")

        cells = []

        # Loop over the maps
        for name in self.ssfr_maps:

            # Determine the relative path
            path = fs.join(plots_name, ssfr_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.ssfr_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_tir_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of TIR maps ...")

        cells = []

        # Loop over the maps
        for name in self.tir_maps:

            # Determine the relative path
            path = fs.join(plots_name, tir_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.tir_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_attenuation_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of attenuation maps ...")

        cells = []

        # Loop over the maps
        for name in self.attenuation_maps:

            # Determine the relative path
            path = fs.join(plots_name, attenuation_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.attenuation_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_old_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of old stellar maps ...")

        cells = []

        # Loop over the maps
        for name in self.old_maps:

            # Determine the relative path
            path = fs.join(plots_name, old_name, name + ".png")

            # Make iamge
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.old_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_young_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of young stellar maps ...")

        cells = []

        # Loop over the maps
        for name in self.young_maps:

            # Determine the relative path
            path = fs.join(plots_name, young_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.young_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_ionizing_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of ionizing stellar maps ...")

        cells = []

        # Loop over the maps
        for name in self.ionizing_maps:

            # Determine the relative path
            path = fs.join(plots_name, ionizing_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.ionizing_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_dust_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of dust maps ...")

        cells = []

        # Loop over the maps
        for name in self.dust_maps:

            # Determine the relative path
            path = fs.join(plots_name, dust_name, name + ".png")

            # Make image
            image = html.image(path, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Make cell
            cell = ""
            cell += html.center(name)
            cell += html.newline
            cell += image

            # Add
            cells.append(cell)

        # Make
        self.dust_table = SimpleTable.rasterize(cells, ncolumns=ncolumns, css_class=self.table_class)

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create the page
        self.page = HTMLPage(self.title, style=page_style, css_path=css_paths, javascript_path=javascripts, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes))

        self.page += html.newline

        # Add the tables
        #self.page += self.table

        # Add the colours table
        self.page += "COLOURS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.colour_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the sSFR table
        self.page += "SSFR"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.ssfr_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the TIR table
        self.page += "TIR"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.tir_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the attenuation table
        self.page += "ATTENUATION"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.attenuation_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the old table
        self.page += "OLD STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.old_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the young table
        self.page += "YOUNG STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.young_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the ionizing table
        self.page += "IONIZING STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.ionizing_table
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the dust table
        self.page += "DUST"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.dust_table
        self.page += html.make_line("heavy")
        self.page += html.newline

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

        # Write the page
        self.write_page()

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Save
        self.page.saveto(self.all_maps_html_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.all_maps_html_page_path)

# -----------------------------------------------------------------

def get_image_info(frame):

    """
    This function ...
    :param frame:
    :return:
    """

    info = []

    fltr = frame.filter
    wavelength = fltr.wavelength if fltr is not None else None
    pixelscale = frame.average_pixelscale
    fwhm = frame.fwhm

    # Get filesize
    filesize = fs.file_size(frame.path).to("MB")

    info.append("Filter: " + tostr(fltr))
    info.append("Wavelength: " + tostr(wavelength))
    info.append("Unit: " + tostr(frame.unit))
    info.append("Pixelscale: " + tostr(pixelscale))
    info.append("PSF filter: " + frame.psf_filter_name)
    info.append("FWHM: " + tostr(fwhm))
    info.append("Dimensions: " + str((frame.xsize, frame.ysize)))
    info.append("File size: " + tostr(filesize))

    # Return the info
    return info

# -----------------------------------------------------------------
