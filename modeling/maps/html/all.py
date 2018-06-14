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

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import MapMakingComponent
from ...html.component import stylesheet_url, page_style
from ...core.environment import colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name
from ....core.tools import filesystem as fs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts, JS9Spawner, make_replace_nans_infs, make_load_region, make_replace_infs_by_nans
from ....core.tools import browser
from ....core.tools.utils import lazyproperty
from ....core.tools import numbers
from ....core.basics.range import RealRange
from ....magic.tools.info import get_image_info
from ....magic.core.frame import Frame
from ....magic.core.image import Image
from ....core.tools.serialization import load_dict, write_dict

# -----------------------------------------------------------------

plots_name = "plots"
ncolumns = 2
background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

intervals_filename = "intervals.dat"

# -----------------------------------------------------------------

class AllMapsPageGenerator(MapMakingComponent):

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

        # The regions
        self.colour_regions = dict()
        self.ssfr_regions = dict()
        self.tir_regions = dict()
        self.attenuation_regions = dict()
        self.old_regions = dict()
        self.young_regions = dict()
        self.ionizing_regions = dict()
        self.dust_regions = dict()

        # The plots
        self.colour_plots = dict()
        self.ssfr_plots = dict()
        self.tir_plots = dict()
        self.attenuation_plots = dict()
        self.old_plots = dict()
        self.young_plots = dict()
        self.ionizing_plots = dict()
        self.dust_plots = dict()

        # The views
        self.colour_views = dict()
        self.ssfr_views = dict()
        self.tir_views = dict()
        self.attenuation_views = dict()
        self.old_views = dict()
        self.young_views = dict()
        self.ionizing_views = dict()
        self.dust_views = dict()

        # Buttons
        self.colour_buttons = defaultdict(list)
        self.ssfr_buttons = defaultdict(list)
        self.tir_buttons = defaultdict(list)
        self.attenuation_buttons = defaultdict(list)
        self.old_buttons = defaultdict(list)
        self.young_buttons = defaultdict(list)
        self.ionizing_buttons = defaultdict(list)
        self.dust_buttons = defaultdict(list)

        # The tables
        self.colour_table = None
        self.ssfr_table = None
        self.tir_table = None
        self.attenuation_table = None
        self.old_table = None
        self.young_table = None
        self.ionizing_table = None
        self.dust_table = None

        # Plot paths
        self.colour_plots_paths = dict()
        self.ssfr_plots_paths = dict()
        self.tir_plots_paths = dict()
        self.attenuation_plots_paths = dict()
        self.old_plots_paths = dict()
        self.young_plots_paths = dict()
        self.ionizing_plots_paths = dict()
        self.dust_plots_paths = dict()

        # The page
        self.page = None

        # The intervals
        self.colour_plot_intervals = dict()
        self.ssfr_plot_intervals = dict()
        self.tir_plot_intervals = dict()
        self.attenuation_plot_intervals = dict()
        self.old_plot_intervals = dict()
        self.young_plot_intervals = dict()
        self.ionizing_plot_intervals = dict()
        self.dust_plot_intervals = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get image info
        if self.config.info: self.get_info()

        # 3. Load the regions
        self.get_regions()

        # 4. Load the intervals
        self.get_intervals()

        # 5. Make plots
        self.make_plots()

        # 6. Make the views
        self.make_views()

        # 7. Make buttons for extra functionality
        self.make_buttons()

        # 8. Make the tables
        self.make_tables()

        # 9. Generate the page
        self.generate_page()

        # 10. Writing
        self.write()

        # 11. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AllMapsPageGenerator, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

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
        return 0.4 * page_width

    # -----------------------------------------------------------------

    @property
    def view_width(self):

        """
        This function ...
        :return:
        """

        return 0.5 * page_width

    # -----------------------------------------------------------------

    @property
    def image_height(self):

        """
        This function ...
        :return:
        """

        #return 300
        return None

    # -----------------------------------------------------------------

    @property
    def view_height(self):

        """
        This function ...
        :return:
        """

        return None

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

    @property
    def has_dust_maps(self):

        """
        Thisn function ...
        :return:
        """

        return self.static_collection.has_dust_maps

    # -----------------------------------------------------------------

    def get_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting image info ...")

        # Colour
        if self.has_colour_maps:
            if self.has_colour_info: self.load_colour_info()
            else: self.get_colour_info()

        # SSFR
        if self.has_ssfr_maps:
            if self.has_ssfr_info: self.load_ssfr_info()
            else: self.get_ssfr_info()

        # TIR
        if self.has_tir_maps:
            if self.has_tir_info: self.load_tir_info()
            else: self.get_tir_info()

        # Attenuation
        if self.has_attenuation_maps:
            if self.has_attenuation_info: self.load_attenuation_info()
            else: self.get_attenuation_info()

        # Old
        if self.has_old_maps:
            if self.has_old_info: self.load_old_info()
            else: self.get_old_info()

        # Young
        if self.has_young_maps:
            if self.has_young_info: self.load_young_info()
            else: self.get_young_info()

        # Ionizing
        if self.has_ionizing_maps:
            if self.has_ionizing_info: self.load_ionizing_info()
            else: self.get_ionizing_info()

        # Dust
        if self.has_dust_maps:
            if self.has_dust_info: self.load_dust_info()
            else: self.get_dust_info()

    # -----------------------------------------------------------------

    @property
    def colour_info_path(self):

        """
        This function ....
        :return:
        """

        return fs.join(self.colour_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_colour_info(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.colour_info_path)

    # -----------------------------------------------------------------

    def load_colour_info(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Loading colour maps info ...")

        #return load_dict(self.colour_info_path)
        self.colour_info = load_dict(self.colour_info_path)

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
            info = get_image_info(name, self.colour_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.colour_info[name] = code

        # Write the info
        write_dict(self.colour_info, self.colour_info_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.ssfr_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_info(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ssfr_info_path)

    # -----------------------------------------------------------------

    def load_ssfr_info(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Loading sSFR maps info ...")

        #return load_dict(self.ssfr_info_path)
        self.ssfr_info = load_dict(self.ssfr_info_path)

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
            info = get_image_info(name, self.ssfr_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.ssfr_info[name] = code

        # Write the info
        write_dict(self.ssfr_info, self.ssfr_info_path)

    # -----------------------------------------------------------------

    @property
    def tir_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.tir_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_tir_info(self):

        """
        Thisfnuction ...
        :return:
        """

        return fs.is_file(self.tir_info_path)

    # -----------------------------------------------------------------

    def load_tir_info(self):

        """
        This funciton ...
        :return:
        """

        # Inform the user
        log.info("Loading TIR maps info ...")

        #return load_dict(self.tir_info_path)
        self.tir_info = load_dict(self.tir_info_path)

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
            info = get_image_info(name, self.tir_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.tir_info[name] = code

        # Write the info
        write_dict(self.tir_info, self.tir_info_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_info_path(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.join(self.attenuation_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_attenuation_info(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.attenuation_info_path)

    # -----------------------------------------------------------------

    def load_attenuation_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading attenuation maps info ...")

        #return load_dict(self.attenuation_info_path)
        self.attenuation_info = load_dict(self.attenuation_info_path)

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
            info = get_image_info(name, self.attenuation_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.attenuation_info[name] = code

        # Write the info
        write_dict(self.attenuation_info, self.attenuation_info_path)

    # -----------------------------------------------------------------

    @property
    def old_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.old_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_info(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_info_path)

    # -----------------------------------------------------------------

    def load_old_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading old stellar maps info ...")

        #return load_dict(self.old_info_path)
        self.old_info = load_dict(self.old_info_path)

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
            info = get_image_info(name, self.old_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.old_info[name] = code

        # Write the info
        write_dict(self.old_info, self.old_info_path)

    # -----------------------------------------------------------------

    @property
    def young_info_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.young_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_info(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.young_info_path)

    # -----------------------------------------------------------------

    def load_young_info(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Loading young stellar maps info ...")

        #return load_dict(self.young_info_path)
        self.young_info = load_dict(self.young_info_path)

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
            info = get_image_info(name, self.young_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.young_info[name] = code

        # Write the info
        write_dict(self.young_info, self.young_info_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_info_path(self):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.ionizing_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_ionizing_info(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_info_path)

    # -----------------------------------------------------------------

    def load_ionizing_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading ionizing stellar maps info ...")

        #return load_dict(self.ionizing_info_path)
        self.ionizing_info = load_dict(self.ionizing_info_path)

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
            info = get_image_info(name, self.ionizing_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.ionizing_info[name] = code

        # Write the info
        write_dict(self.ionizing_info, self.ionizing_info_path)

    # -----------------------------------------------------------------

    @property
    def dust_info_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.dust_plots_path, "info.dat")

    # -----------------------------------------------------------------

    @property
    def has_dust_info(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_info_path)

    # -----------------------------------------------------------------

    def load_dust_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading dust map info ...")

        #return load_dict(self.dust_info_path)
        self.dust_info = load_dict(self.dust_info_path)

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
            info = get_image_info(name, self.dust_maps[name], path=False)

            # Make list
            code = html.dictionary(info, key_color=key_color)

            # Add info
            self.dust_info[name] = code

        # Write the info
        write_dict(self.dust_info, self.dust_info_path)

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

    @property
    def alpha(self):

        """
        This function ...
        :return:
        """

        if self.config.alpha: return "combined"
        else: return None

    # -----------------------------------------------------------------

    @property
    def peak_alpha(self):

        """
        This function ...
        :return:
        """

        return 2.5

    # -----------------------------------------------------------------

    def make_rgba_plot(self, name, frame, filepath, around_zero=False, scale=None, interval="pts", colours="jet",
                       alpha="absolute", peak_alpha=1.):

        """
        This function ...
        :param name:
        :param frame:
        :param filepath:
        :param around_zero:
        :param scale: if None, configured value is used
        :param interval:
        :param colours:
        :param alpha:
        :param peak_alpha:
        :return:
        """

        # Debugging
        log.debug("Making an RGBA plot from the '" + name + "' map at '" + filepath + "' ...")

        # Get primary frame
        if isinstance(frame, Image): frame = frame.primary
        elif isinstance(frame, Frame): pass
        else: raise ValueError("Invalid argument of type '" + str(type(frame)) + "'")

        # Crop the frame
        frame = frame.cropped_to(self.truncation_box, factor=self.config.cropping_factor, out_of_bounds="expand")

        # Get the truncation mask and mask out the pixel beyond the truncation limit
        wcs, xsize, ysize = frame.wcs, frame.xsize, frame.ysize
        ellipse = self.truncation_ellipse.to_pixel(wcs)
        mask = ellipse.to_mask(xsize, ysize).inverse()
        frame[mask] = 0.0

        # Around zero?
        if around_zero: symmetric = True
        else: symmetric = False

        # Make RGBA image
        rgba, vmin, vmax = frame.to_rgba(scale=scale, colours=colours, around_zero=around_zero, symmetric=symmetric, alpha=alpha, interval=interval, return_minmax=True, peak_alpha=peak_alpha)
        rgba.soften_edges(self.softening_ellipse.to_pixel(wcs), self.softening_range)

        # Save
        rgba.saveto(filepath)

        # Return vmin and vmax
        return vmin, vmax

    # -----------------------------------------------------------------

    def get_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the regions ...")

        # Colour
        if self.has_colour_maps: self.get_colour_regions()

        # sSFR
        if self.has_ssfr_maps: self.get_ssfr_regions()

        # TIR
        if self.has_tir_maps: self.get_tir_regions()

        # Attenuation
        if self.has_attenuation_maps: self.get_attenuation_regions()

        # Old stars
        if self.has_old_maps: self.get_old_regions()

        # Young stars
        if self.has_young_maps: self.get_young_regions()

        # Ionizing stars
        if self.has_ionizing_maps: self.get_ionizing_regions()

        # Dust
        if self.has_dust_maps: self.get_dust_regions()

    # -----------------------------------------------------------------

    def get_colour_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.colour_maps[name].wcs)

            # Add
            self.colour_regions[name] = region

    # -----------------------------------------------------------------

    def get_ssfr_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.ssfr_maps[name].wcs)

            # Add
            self.ssfr_regions[name] = region

    # -----------------------------------------------------------------

    def get_tir_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.tir_maps[name].wcs)

            # Add
            self.tir_regions[name] = region

    # -----------------------------------------------------------------

    def get_attenuation_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.attenuation_maps[name].wcs)

            # Add
            self.attenuation_regions[name] = region

    # -----------------------------------------------------------------

    def get_old_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.old_maps[name].wcs)

            # Add
            self.old_regions[name] = region

    # -----------------------------------------------------------------

    def get_young_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.young_maps[name].wcs)

            # Add
            self.young_regions[name] = region

    # -----------------------------------------------------------------

    def get_ionizing_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.ionizing_maps[name].wcs)

            # Add
            self.ionizing_regions[name] = region

    # -----------------------------------------------------------------

    def get_dust_regions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the regions for the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Get truncation region in image coordinates
            region = self.truncation_ellipse.to_pixel(self.dust_maps[name].wcs)

            # Add
            self.dust_regions[name] = region

    # -----------------------------------------------------------------

    def get_intervals(self):

        """
        This function ...
        :return:
        """

        # Inform the uer
        log.info("Loading the intervals ...")

        # Colours
        if self.has_colour_maps and self.has_colour_intervals: self.load_colour_intervals()

        # sSFR
        if self.has_ssfr_maps and self.has_ssfr_intervals: self.load_ssfr_intervals()

        # TIR
        if self.has_tir_maps and self.has_tir_intervals: self.load_tir_intervals()

        # Attenuation
        if self.has_attenuation_maps and self.has_attenuation_intervals: self.load_attenuation_intervals()

        # Old
        if self.has_old_maps and self.has_old_intervals: self.load_old_intervals()

        # Young
        if self.has_young_maps and self.has_young_intervals: self.load_young_intervals()

        # Ionizing
        if self.has_ionizing_maps and self.has_ionizing_intervals: self.load_ionizing_intervals()

        # Dust
        if self.has_dust_maps and self.has_dust_intervals: self.load_dust_intervals()

    # -----------------------------------------------------------------

    @property
    def colour_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.colour_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_colour_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.colour_intervals_path) and fs.has_lines(self.colour_intervals_path)

    # -----------------------------------------------------------------

    def load_colour_intervals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading colour map intervals ...")

        #return load_dict(self.colour_intervals_path)
        self.colour_plot_intervals = load_dict(self.colour_intervals_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.ssfr_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_intervals(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.ssfr_intervals_path) and fs.has_lines(self.ssfr_intervals_path)

    # -----------------------------------------------------------------

    def load_ssfr_intervals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading sSFR map intervals ...")

        #return load_dict(self.ssfr_intervals_path)
        self.ssfr_plot_intervals = load_dict(self.ssfr_intervals_path)

    # -----------------------------------------------------------------

    @property
    def tir_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.tir_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_tir_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.tir_intervals_path) and fs.has_lines(self.tir_intervals_path)

    # -----------------------------------------------------------------

    def load_tir_intervals(self):

        """
        This ufnciton ..
        :return:
        """

        # Inform the user
        log.info("Loading TIR map intervals ...")

        #return load_dict(self.tir_intervals_path)
        self.tir_plot_intervals = load_dict(self.tir_intervals_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.attenuation_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_attenuation_intervals(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.is_file(self.attenuation_intervals_path) and fs.has_lines(self.attenuation_intervals_path)

    # -----------------------------------------------------------------

    def load_attenuation_intervals(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Loading attenuation map intervals ...")

        #return load_dict(self.attenuation_intervals_path)
        self.attenuation_plot_intervals = load_dict(self.attenuation_intervals_path)

    # -----------------------------------------------------------------

    @property
    def old_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.old_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_old_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.old_intervals_path) and fs.has_lines(self.old_intervals_path)

    # -----------------------------------------------------------------

    def load_old_intervals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading old stellar map intervals ...")

        #return load_dict(self.old_intervals_path)
        self.old_plot_intervals = load_dict(self.old_intervals_path)

    # -----------------------------------------------------------------

    @property
    def young_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.young_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_young_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.young_intervals_path) and fs.has_lines(self.young_intervals_path)

    # -----------------------------------------------------------------

    def load_young_intervals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading young stellar map intervals ...")

        #return load_dict(self.young_intervals_path)
        self.young_plot_intervals = load_dict(self.young_intervals_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.ionizing_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_ionizing_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ionizing_intervals_path) and fs.has_lines(self.ionizing_intervals_path)

    # -----------------------------------------------------------------

    def load_ionizing_intervals(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Loading ionizing stellar map intervals ...")

        #return load_dict(self.ionizing_intervals_path)
        self.ionizing_plot_intervals = load_dict(self.ionizing_intervals_path)

    # -----------------------------------------------------------------

    @property
    def dust_intervals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.dust_plots_path, intervals_filename)

    # -----------------------------------------------------------------

    @property
    def has_dust_intervals(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_intervals_path) and fs.has_lines(self.dust_intervals_path)

    # -----------------------------------------------------------------

    def load_dust_intervals(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Loading dust map intervals ...")

        #return load_dict(self.dust_intervals_path)
        self.dust_plot_intervals = load_dict(self.dust_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.colour_maps[name], filepath, colours=self.colours_cmap,
                                             scale=self.colours_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.colour_plot_intervals[name] = (vmin, vmax)

        # Save the intervals
        write_dict(self.colour_plot_intervals, self.colour_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.ssfr_maps[name], filepath, colours=self.ssfr_cmap,
                                             scale=self.ssfr_scale, interval=self.ssfr_interval, alpha=self.alpha,
                                             peak_alpha=self.peak_alpha)

            # Set the interval
            self.ssfr_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.ssfr_plot_intervals, self.ssfr_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.tir_maps[name], filepath, colours=self.tir_cmap,
                                             scale=self.tir_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.tir_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.tir_plot_intervals, self.tir_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.attenuation_maps[name], filepath, colours=self.attenuation_cmap,
                                             scale=self.attenuation_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.attenuation_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.attenuation_plot_intervals, self.attenuation_intervals_path)

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
                if self.config.replot or name not in self.old_plot_intervals: fs.remove_file(filepath)
                else: continue

            # Make the plot
            vmin, vmax = self.make_rgba_plot(name, self.old_maps[name], filepath, colours=self.old_cmap,
                                             scale=self.old_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.old_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.old_plot_intervals, self.old_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.young_maps[name], filepath, colours=self.young_cmap,
                                             scale=self.young_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.young_plot_intervals[name] = (vmin, vmax)

        # Write the interrvals
        write_dict(self.young_plot_intervals, self.young_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.ionizing_maps[name], filepath, colours=self.ionizing_cmap,
                                             scale=self.ionizing_scale, alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.ionizing_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.ionizing_plot_intervals, self.ionizing_intervals_path)

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
            vmin, vmax = self.make_rgba_plot(name, self.dust_maps[name], filepath, colours=self.dust_cmaps["attenuation"],
                                             scale=self.dust_scales["attenuation"], alpha=self.alpha, peak_alpha=self.peak_alpha)

            # Set the interval
            self.dust_plot_intervals[name] = (vmin, vmax)

        # Write the intervals
        write_dict(self.dust_plot_intervals, self.dust_intervals_path)

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "realtable"

    # -----------------------------------------------------------------

    def make_view(self, name, filepath, plot, scale, colormap, zoom, next_zoom=None, interval=None):

        """
        This function ...
        :param name:
        :param filepath:
        :param plot:
        :param scale:
        :param colormap:
        :param zoom:
        :param next_zoom:
        :param interval:
        :return:
        """

        # Debugging
        log.debug("Making a view of the '" + name + "' map ...")

        settings = dict()
        settings["scale"] = scale
        settings["colormap"] = colormap
        settings["zoom"] = zoom

        regions_for_loader = None

        # Set text
        text = plot

        # Determine filepath
        filepath = fs.in_localhost_path(filepath, self.config.path)

        # Create the loader
        loader = JS9Spawner.from_path(text, name, filepath, settings=settings, button=False,
                                      menubar=self.config.menubar, colorbar=self.config.colorbar,
                                      regions=regions_for_loader, background_color=background_color,
                                      replace=True, width=self.view_height, height=self.view_height,
                                      center=True, replace_nans=True, replace_infs=True, zoom=next_zoom,
                                      interval=interval)

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

            # Determine image name
            image_name = "colours___" + name

            # Get the interval
            interval = self.colour_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.colour_plots[name], self.colours_scale,
                                  self.colours_js9_cmap, self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "ssfr___" + name

            # Get the interval
            interval = self.ssfr_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.ssfr_plots[name], self.ssfr_scale, self.ssfr_js9_cmap,
                                  self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "tir___" + name

            # Get the interval
            interval = self.tir_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.tir_plots[name], self.tir_scale, self.tir_js9_cmap,
                                  self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "attenuation___" + name

            # Get the interval
            interval = self.attenuation_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.attenuation_plots[name], self.attenuation_scale,
                                  self.attenuation_js9_cmap, self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "old___" + name

            # Get the interval
            interval = self.old_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.old_plots[name], self.old_scale, self.old_js9_cmap,
                                  self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "young___" + name

            # Get the interval
            interval = self.young_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.young_plots[name], self.young_scale, self.young_js9_cmap,
                                  self.config.zoom, interval=interval)

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

            # Determine image name
            image_name = "ionizing___" + name

            # Get the interval
            interval = self.ionizing_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.ionizing_plots[name], self.ionizing_scale,
                                  self.ionizing_js9_cmap, self.config.zoom, interval=interval)

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

        # Loop over the maps
        for name in self.dust_maps:

            # Get path
            if self.config.view_png: path = self.relative_path(self.dust_plots_paths[name])
            else: path = self.dust_maps[name].path

            # Determine image name
            image_name = "dust___" + name

            # Get the interval
            interval = self.dust_plot_intervals[name]

            # Make the view
            view = self.make_view(image_name, path, self.dust_plots[name], self.dust_scales["attenuation"],
                                  self.dust_js9_cmap, self.config.zoom, interval=interval)

            # Add
            self.dust_views[name] = view

    # -----------------------------------------------------------------

    def make_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons ...")

        # Colours
        if self.has_colour_maps: self.make_colour_buttons()

        # sSFR
        if self.has_ssfr_maps: self.make_ssfr_buttons()

        # TIR
        if self.has_tir_maps: self.make_tir_buttons()

        # Attenuation
        if self.has_attenuation_maps: self.make_attenuation_buttons()

        # Old stellar maps
        if self.has_old_maps: self.make_old_buttons()

        # Young stellar maps
        if self.has_young_maps: self.make_young_buttons()

        # Ionizing stellar maps
        if self.has_ionizing_maps: self.make_ionizing_buttons()

        # Dust maps
        if self.has_dust_maps: self.make_dust_buttons()

    # -----------------------------------------------------------------

    def make_colour_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            ## NANS/INFS

            # Create nan/infs replacer button
            button_id = self.colour_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.colour_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.colour_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.colour_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.colour_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.colour_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.colour_views[name].image_name)
            load_region = make_load_region(self.colour_regions[name], display=self.colour_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.colour_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_ssfr_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            ## NANS/INFS

            # Create nan/infs replacer button
            button_id = self.ssfr_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.ssfr_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.ssfr_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.ssfr_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.ssfr_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.ssfr_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.ssfr_views[name].image_name)
            load_region = make_load_region(self.ssfr_regions[name], display=self.ssfr_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.ssfr_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_tir_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.tir_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.tir_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.tir_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.tir_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.tir_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.tir_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.tir_views[name].image_name)
            load_region = make_load_region(self.tir_regions[name], display=self.tir_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.tir_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_attenuation_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.attenuation_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.attenuation_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.attenuation_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.attenuation_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.attenuation_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.attenuation_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.attenuation_views[name].image_name)
            load_region = make_load_region(self.attenuation_regions[name], display=self.attenuation_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.attenuation_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_old_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.old_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.old_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.old_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.old_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.old_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.old_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.old_views[name].image_name)
            load_region = make_load_region(self.old_regions[name], display=self.old_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.old_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_young_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.young_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.young_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.young_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.young_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.young_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.young_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.young_views[name].image_name)
            load_region = make_load_region(self.young_regions[name], display=self.young_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.young_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_ionizing_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.ionizing_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.ionizing_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.ionizing_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.ionizing_views[name].display_id)

            # Create the button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.ionizing_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.ionizing_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.ionizing_views[name].image_name)
            load_region = make_load_region(self.ionizing_regions[name], display=self.ionizing_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.ionizing_buttons[name].append(region_button)

    # -----------------------------------------------------------------

    def make_dust_buttons(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making buttons for dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            ## NANS/INFS

            # Create nans/infs replacer button
            button_id = self.dust_views[name].image_name + "nansinfs"
            replace_function_name = "replace_infs_nans_" + html.make_usable(self.dust_views[name].image_name)
            if self.config.replace_nans: replace_nans_infs = make_replace_nans_infs(self.dust_views[name].display_id)
            else: replace_nans_infs = make_replace_infs_by_nans(self.dust_views[name].display_id)

            # Create button
            button = html.make_script_button(button_id, "Replace infs/nans", replace_nans_infs, replace_function_name)
            self.dust_buttons[name].append(button)

            ## REGIONS

            # Regions button
            region_button_id = self.dust_views[name].image_name + "regionsbutton"
            load_region_function_name = "load_regions_" + html.make_usable(self.dust_views[name].image_name)
            load_region = make_load_region(self.dust_regions[name], display=self.dust_views[name].display_id, movable=False,
                                           rotatable=False, removable=False, resizable=False, quote_character="'")

            # Create region loader
            region_button = html.make_script_button(region_button_id, "Load regions", load_region, load_region_function_name)
            self.dust_buttons[name].append(region_button)

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.colour_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.colour_views[name]))

            # Add info
            cell += html.newline
            cell += self.colour_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.ssfr_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.ssfr_views[name]))

            # Add info
            cell += html.newline
            cell += self.ssfr_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.tir_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.tir_views[name]))

            # Add info
            cell += html.newline
            cell += self.tir_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.attenuation_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.attenuation_views[name]))

            # Add info
            cell += html.newline
            cell += self.attenuation_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.old_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.old_views[name]))

            # Add info
            cell += html.newline
            cell += self.old_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.young_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.young_views[name]))

            # Add info
            cell += html.newline
            cell += self.young_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.ionizing_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.ionizing_views[name]))

            # Add info
            cell += html.newline
            cell += self.ionizing_info[name]

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

            # Make cell
            cell = ""

            # Add buttons
            for button in self.dust_buttons[name]: cell += html.center(str(button))

            # Add the view
            cell += html.center(str(self.dust_views[name]))

            # Add info
            cell += html.newline
            cell += self.dust_info[name]

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

        test_path = fs.join(self.config.path, "maps", "raw", "attenuation", "buat", "NUV__single_SPIRE_PSW.fits")
        localhost_test_path = fs.in_localhost_path(test_path, self.config.path)
        self.page += html.hyperlink(localhost_test_path, "get test")

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

        # Serve and show
        #with browser.serve_local_host(): browser.open_path(self.all_maps_html_page_path)

        thread = browser.start_localhost()
        browser.open_path(self.all_maps_html_page_path)

        raw_input("ENTER...")

# -----------------------------------------------------------------
