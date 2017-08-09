#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html.summary Contains the MapsSummaryPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import MapsComponent
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size
from ...core.environment import map_sub_names, colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ....core.tools import filesystem as fs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts, JS9Spawner, make_replace_nans_infs
from ....core.tools import browser
from ....core.tools.utils import lazyproperty
from ....core.tools import numbers
from ....core.basics.range import RealRange
from ....magic.tools.info import get_image_info_strings, get_image_info

# -----------------------------------------------------------------

background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

class MapsSummaryPageGenerator(MapsComponent):

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
        super(MapsSummaryPageGenerator, self).__init__(*args, **kwargs)

        # The image info
        self.colour_info = dict()
        self.ssfr_info = dict()
        self.tir_info = dict()
        self.attenuation_info = dict()
        self.old_info = dict()
        self.young_info = dict()
        self.ionizing_info = dict()
        self.dust_info = dict()

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
        self.get_info()

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
        super(MapsSummaryPageGenerator, self).setup(**kwargs)

        # Set the number of allowed open file handles
        #fs.set_nallowed_open_files(self.config.nopen_files)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Maps summary"

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
        if self.has_colour_maps: self.get_colour_info()

        # SSFR
        if self.has_ssfr_maps: self.get_ssfr_info()

        # TIR
        if self.has_tir_maps: self.get_tir_info()

        # Attenuation
        if self.has_attenuation_maps: self.get_attenuation_info()

        # Old
        if self.has_old_maps: self.get_old_info()

        # Young
        if self.has_young_maps: self.get_young_info()

        # Ionizing
        if self.has_ionizing_maps: self.get_ionizing_info()

        # Dust
        if self.has_dust_maps: self.get_dust_info()

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
            #info = get_image_info_strings(name, self.colour_maps[name])
            info = get_image_info(name, self.colour_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

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

        # Loop over the maps
        for name in self.ssfr_maps:

            # Get info
            #info = get_image_info_strings(name, self.ssfr_maps[name])
            info = get_image_info(name, self.ssfr_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.ssfr_info[name] = code

    # -----------------------------------------------------------------

    def get_tir_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Get info
            #info = get_image_info_strings(name, self.tir_maps[name])
            info = get_image_info(name, self.tir_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.tir_info[name] = code

    # -----------------------------------------------------------------

    def get_attenuation_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Get info
            #info = get_image_info_strings(name, self.attenuation_maps[name])
            info = get_image_info(name, self.attenuation_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.attenuation_info[name] = code

    # -----------------------------------------------------------------

    def get_old_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Get info
            #info = get_image_info_strings(name, self.old_maps[name])
            info = get_image_info(name, self.old_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.old_info[name] = code

    # -----------------------------------------------------------------

    def get_young_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get info
            #info = get_image_info_strings(name, self.young_maps[name])
            info = get_image_info(name, self.young_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.young_info[name] = code

    # -----------------------------------------------------------------

    def get_ionizing_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Get info
            #info = get_image_info_strings(name, self.ionizing_maps[name])
            info = get_image_info(name, self.ionizing_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.ionizing_info[name] = code

    # -----------------------------------------------------------------

    def get_dust_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting info about the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Get info
            #info = get_image_info_strings(name, self.dust_maps[name])
            info = get_image_info(name, self.dust_maps[name], path=False)

            # Make list
            #code = html.unordered_list(info)
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.dust_info[name] = code

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

    def make_rgba_plot(self, name, frame, filepath):

        """
        This function ...
        :param frame:
        :param filepath:
        :return:
        """

        # Debugging
        log.debug("Making an RGBA plot from the '" + name + "' map at '" + filepath + "' ...")

        # Crop the frame
        frame = frame.cropped_to(self.truncation_box, factor=self.config.cropping_factor)

        # Get the truncation mask and mask out the pixel beyond the truncation limit
        wcs, xsize, ysize = frame.wcs, frame.xsize, frame.ysize
        ellipse = self.truncation_ellipse.to_pixel(wcs)
        mask = ellipse.to_mask(xsize, ysize).inverse()

        #from ....magic.tools import plotting
        #plotting.plot_mask(mask)
        #plotting.plot_mask(self.truncation_box.to_pixel(wcs).to_mask(xsize, ysize).inverse())

        frame[mask] = 0.0

        # Make RGBA image
        rgba = frame.to_rgba(scale=self.config.scale, colours=self.config.colours, absolute_alpha=True)
        rgba.soften_edges(self.softening_ellipse.to_pixel(wcs), self.softening_range)

        # Save
        rgba.saveto(filepath)

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

            # Set the filepath
            self.colour_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.colour_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.colour_maps[name], filepath)

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

            # Set the filepath
            self.ssfr_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.ssfr_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.ssfr_maps[name], filepath)

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

            # Set the filepath
            self.tir_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.tir_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.tir_maps[name], filepath)

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

            # Set the filepath
            self.attenuation_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.attenuation_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.attenuation_maps[name], filepath)

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

            # Set the filepath
            self.old_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.old_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.old_maps[name], filepath)

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

            # Set the filepath
            self.young_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.young_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.young_maps[name], filepath)

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

            # Set the filepath
            self.ionizing_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.ionizing_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.ionizing_maps[name], filepath)

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

            # Set the filepath
            self.dust_plots_paths[name] = filepath

            # Determine the relative path
            relpath = fs.relative_to(filepath, self.maps_html_path)

            # Make image plot
            self.dust_plots[name] = html.image(relpath, alttext=name, height=self.image_height, width=self.image_width, hover=None)

            # Check if plot is already made
            if fs.is_file(filepath):
                if self.config.replot: fs.remove_file(filepath)
                else: continue

            # Make the plot
            self.make_rgba_plot(name, self.dust_maps[name], filepath)

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "realtable"

    # -----------------------------------------------------------------

    def make_view(self, name, filepath, plot):

        """
        This function ...
        :param name:
        :param filepath:
        :param plot
        :return:
        """

        # Debugging
        log.debug("Making a view of the '" + name + "' map ...")

        settings = dict()
        settings["scale"] = self.config.scale
        settings["colormap"] = self.config.colormap
        settings["zoom"] = self.config.zoom

        # Get region in image coordinates
        #region = self.disk_ellipse.to_pixel(self.coordinate_systems[name])
        #regions_for_loader = region if self.config.load_regions else None

        # Add the region
        #self.ellipses[name] = region

        regions_for_loader = None

        # Set text
        text = plot

        # Create the loader
        loader = JS9Spawner.from_path(text, name, filepath, settings=settings, button=False,
                                      menubar=self.config.menubar, colorbar=self.config.colorbar,
                                      regions=regions_for_loader, background_color=background_color,
                                      replace=True, width=self.view_height, height=self.view_height,
                                      center=True, replace_nans=True, replace_infs=True)

        #display_id = self.loaders[name].display_id
        #self.windows[name] = self.loaders[name].placeholder

        # Set load info
        #load_info[display_id] = (name, path, regions_for_loader)
        #images[display_id] = self.loaders[name].image
        #placeholders[display_id] = self.loaders[name].spawn_div_name

        # Return the loader
        return loader

    # -----------------------------------------------------------------

    def relative_path(self, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        return fs.relative_to(filepath, self.maps_html_path)

    # -----------------------------------------------------------------

    def make_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views ...")

        # Colours
        if self.has_colour_maps: self.make_colour_views()

        # sSFR
        if self.has_ssfr_maps: self.make_ssfr_views()

        # TIR
        if self.has_tir_maps: self.make_tir_views()

        # Attenuation
        if self.has_attenuation_maps: self.make_attenuation_views()

        # Old
        if self.has_old_maps: self.make_old_views()

        # Young
        if self.has_young_maps: self.make_young_views()

        # Ionizing
        if self.has_ionizing_maps: self.make_ionizing_views()

        # Dust
        if self.has_dust_maps: self.make_dust_views()

    # -----------------------------------------------------------------

    def make_colour_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.colour_plots_paths[name])
            else: path = self.colour_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.colour_plots[name])

            # Add
            self.colour_views[name] = view

    # -----------------------------------------------------------------

    def make_ssfr_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.ssfr_plots_paths[name])
            else: path = self.ssfr_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.ssfr_plots[name])

            # Add
            self.ssfr_views[name] = view

    # -----------------------------------------------------------------

    def make_tir_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.tir_plots_paths[name])
            else: path = self.tir_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.tir_plots[name])

            # Add
            self.tir_views[name] = view

    # -----------------------------------------------------------------

    def make_attenuation_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.attenuation_plots_paths[name])
            else: path = self.attenuation_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.attenuation_plots[name])

            # Add
            self.attenuation_views[name] = view

    # -----------------------------------------------------------------

    def make_old_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.old_plots_paths[name])
            else: path = self.old_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.old_plots[name])

            # Add
            self.old_views[name] = view

    # -----------------------------------------------------------------

    def make_young_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.young_plots_paths[name])
            else: path = self.young_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.young_plots[name])

            # Add
            self.young_views[name] = view

    # -----------------------------------------------------------------

    def make_ionizing_views(self):

        """
        This functino ...
        :return:
        """

        # Inform the user
        log.info("Making views of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.ionizing_plots_paths[name])
            else: path = self.ionizing_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.ionizing_plots[name])

            # Add
            self.ionizing_views[name] = view

    # -----------------------------------------------------------------

    def make_dust_views(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making views of the dust maps ...")

        # Loopover the maps
        for name in self.dust_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.dust_plots_paths[name])
            else: path = self.dust_maps[name].path

            # Make the view
            view = self.make_view(name, path, self.dust_plots[name])

            # Add
            self.dust_views[name] = view

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

        # Make
        self.colour_table = SimpleTable(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_ssfr_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of sSFR maps ...")

        # Make
        self.ssfr_table = SimpleTable(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_tir_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of TIR maps ...")

        # Make
        self.tir_table = SimpleTable(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_attenuation_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of attenuation maps ...")

        # Make
        self.attenuation_table = SimpleTable.rasterize(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_old_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of old stellar maps ...")

        # Make
        self.old_table = SimpleTable.rasterize(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_young_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of young stellar maps ...")

        # Make
        self.young_table = SimpleTable.rasterize(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_ionizing_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of ionizing stellar maps ...")

        # Make
        self.ionizing_table = SimpleTable.rasterize(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def make_dust_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of dust maps ...")

        # Make
        self.dust_table = SimpleTable.rasterize(cells, css_class=self.table_class)

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Generating the page ...")

        # Create list of css scripts
        css_paths = css_scripts[:]
        css_paths.append(stylesheet_url)

        # Create CSS for the page width
        css = make_page_width(page_width)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths, javascript_path=javascripts, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes))
        self.page += html.newline

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
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the TIR table
        self.page += "TIR"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.tir_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the attenuation table
        self.page += "ATTENUATION"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.attenuation_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the old table
        self.page += "OLD STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.old_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the young table
        self.page += "YOUNG STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.young_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the ionizing table
        self.page += "IONIZING STARS"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.ionizing_table
        self.page += html.newline
        self.page += html.newline
        self.page += html.make_line("heavy")
        self.page += html.newline

        # Add the dust table
        self.page += "DUST"
        self.page += html.newline
        self.page += html.line
        self.page += html.newline
        self.page += self.dust_table
        self.page += html.newline
        self.page += html.newline
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
