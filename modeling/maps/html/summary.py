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

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ..component import MapsComponent
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size, sortable_url, preview_url
from ....core.tools import filesystem as fs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts, JS9Spawner, make_replace_nans_infs
from ....core.tools import browser
from ....magic.tools.info import get_image_info_strings, get_image_info
from ....core.basics.table import SmartTable
from ....core.basics.composite import SimplePropertyComposite
from ....core.tools.utils import lazyproperty
from .all import plots_name, colours_name, ssfr_name, tir_name, attenuation_name, old_name, young_name, ionizing_name, dust_name

# -----------------------------------------------------------------

background_color = "white"
key_color = "#4180d3"

# -----------------------------------------------------------------

page_width = 600

# -----------------------------------------------------------------

thumbnail_title = "Thumbnail"

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

        # The map thumbnails
        self.colour_thumbnails = dict()
        self.ssfr_thumbnails = dict()
        self.tir_thumbnails = dict()
        self.attenuation_thumbnails = dict()
        self.old_thumbnails = dict()
        self.young_thumbnails = dict()
        self.ionizing_thumbnails = dict()
        self.dust_thumbnails = dict()

        # The map previews
        self.colour_previews = dict()
        self.ssfr_previews = dict()
        self.tir_previews = dict()
        self.attenuation_previews = dict()
        self.old_previews = dict()
        self.young_previews = dict()
        self.ionizing_previews = dict()
        self.dust_previews = dict()

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

        # Make the thumbnails
        if self.config.thumbnails: self.make_thumbnails()

        # Make the previews
        if self.config.previews: self.make_previews()

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
        fs.set_nallowed_open_files(self.config.nopen_files)

    # -----------------------------------------------------------------

    @property
    def plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, plots_name)

    # -----------------------------------------------------------------

    @property
    def colour_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, colours_name)

    # -----------------------------------------------------------------

    @property
    def has_colour_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.colour_plots_path) and not fs.is_empty(self.colour_plots_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, ssfr_name)

    # -----------------------------------------------------------------

    @property
    def has_ssfr_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.ssfr_plots_path) and not fs.is_empty(self.ssfr_plots_path)

    # -----------------------------------------------------------------

    @property
    def tir_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, tir_name)

    # -----------------------------------------------------------------

    @property
    def has_tir_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.tir_plots_path) and not fs.is_empty(self.tir_plots_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, attenuation_name)

    # -----------------------------------------------------------------

    @property
    def has_attenuation_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.attenuation_plots_path) and not fs.is_empty(self.attenuation_plots_path)

    # -----------------------------------------------------------------

    @property
    def old_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, old_name)

    # -----------------------------------------------------------------

    @property
    def has_old_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.old_plots_path) and not fs.is_empty(self.old_plots_path)

    # -----------------------------------------------------------------

    @property
    def young_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, young_name)

    # -----------------------------------------------------------------

    @property
    def has_young_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.young_plots_path) and not fs.is_empty(self.young_plots_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_plots_path(self):

        """
        Thisfunction ...
        :return:
        """

        return fs.join(self.plots_path, ionizing_name)

    # -----------------------------------------------------------------

    @property
    def has_ionizing_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.ionizing_plots_path) and not fs.is_empty(self.ionizing_plots_path)

    # -----------------------------------------------------------------

    @property
    def dust_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, dust_name)

    # -----------------------------------------------------------------

    @property
    def has_dust_plots(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.dust_plots_path) and not fs.is_empty(self.dust_plots_path)

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

    @lazyproperty
    def info_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["path"] = False
        kwargs["name"] = False
        kwargs["xsize"] = False
        kwargs["ysize"] = False
        kwargs["psf_filter"] = False
        kwargs["filesize"] = False
        kwargs["filter"] = False
        kwargs["wavelength"] = False
        return kwargs

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
            info = get_image_info(name, self.colour_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.colour_info[name] = info

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
            info = get_image_info(name, self.ssfr_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.ssfr_info[name] = info

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
            info = get_image_info(name, self.tir_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.tir_info[name] = info

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
            info = get_image_info(name, self.attenuation_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.attenuation_info[name] = info

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
            info = get_image_info(name, self.old_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.old_info[name] = info

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
            info = get_image_info(name, self.young_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.young_info[name] = info

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
            info = get_image_info(name, self.ionizing_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.ionizing_info[name] = info

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
            info = get_image_info(name, self.dust_maps[name], **self.info_kwargs)

            # Make list
            #code = html.unordered_list(info)
            #code = html.dictionary(info, key_color=key_color)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.dust_info[name] = info

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "sortable"

    # -----------------------------------------------------------------

    def relative_path(self, filepath):

        """
        This function ...
        :param filepath:
        :return:
        """

        return fs.relative_to(filepath, self.maps_html_path)

    # -----------------------------------------------------------------

    def make_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails ...")

        # Colours
        if self.has_colour_maps and self.has_colour_plots: self.make_colour_thumbnails()

        # sSFR
        if self.has_ssfr_maps and self.has_ssfr_plots: self.make_ssfr_thumbnails()

        # TIR
        if self.has_tir_maps and self.has_tir_plots: self.make_tir_thumbnails()

        # Attenuation
        if self.has_attenuation_maps and self.has_attenuation_plots: self.make_attenuation_thumbnails()

        # Old
        if self.has_old_maps and self.has_old_plots: self.make_old_thumbnails()

        # Young
        if self.has_young_maps and self.has_young_plots: self.make_young_thumbnails()

        # Ionizing
        if self.has_ionizing_maps and self.has_ionizing_plots: self.make_ionizing_thumbnails()

        # Dust
        if self.has_dust_maps and self.has_dust_plots: self.make_dust_thumbnails()

    # -----------------------------------------------------------------

    def make_colour_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Find the plot file
            path = fs.join(self.colour_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.colour_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_ssfr_thumbnails(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            # Find the plot file
            path = fs.join(self.ssfr_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.ssfr_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_tir_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Find the plot file
            path = fs.join(self.tir_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.tir_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_attenuation_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Find the plot file
            path = fs.join(self.attenuation_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.attenuation_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_old_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Find the plot file
            path = fs.join(self.old_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.old_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_young_thumbnails(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Find the plot file
            path = fs.join(self.young_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.young_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_ionizing_thumbnails(self):

        """
        This function
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Find the plot file
            path = fs.join(self.ionizing_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.ionizing_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_dust_thumbnails(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Find the plot file
            path = fs.join(self.dust_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.dust_thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews ...")

        # Colour maps
        if self.has_colour_maps: self.make_colour_previews()

        # sSFR maps
        if self.has_ssfr_maps: self.make_ssfr_previews()

        # TIR maps
        if self.has_tir_maps: self.make_tir_previews()

        # Attenuation
        if self.has_attenuation_maps: self.make_attenuation_previews()

        # Old stars
        if self.has_old_maps: self.make_old_previews()

        # Young stars
        if self.has_young_maps: self.make_young_previews()

        # Ionizing stars
        if self.has_ionizing_maps: self.make_ionizing_previews()

        # Dust
        if self.has_dust_maps: self.make_dust_previews()

    # -----------------------------------------------------------------

    def make_colour_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the colour maps ...")

        # Loop over the maps
        for name in self.colour_maps:

            # Get the path
            path = html.get_image_url(self.colour_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.colour_thumbnails[name])

            # Add
            self.colour_previews[name] = preview

    # -----------------------------------------------------------------

    def make_ssfr_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the sSFR maps ...")

        # Loop over the maps
        for name in self.ssfr_maps:

            # Get the path
            path = html.get_image_url(self.ssfr_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.ssfr_thumbnails[name])

            # Add
            self.ssfr_previews[name] = preview

    # -----------------------------------------------------------------

    def make_tir_previews(self):

        """
        This function ..
        """

        # Inform the user
        log.info("Making previews for the TIR maps ...")

        # Loop over the maps
        for name in self.tir_maps:

            # Get the path
            path = html.get_image_url(self.tir_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.tir_thumbnails[name])

            # Add
            self.tir_previews[name] = preview

    # -----------------------------------------------------------------

    def make_attenuation_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the attenuation maps ...")

        # Loop over the maps
        for name in self.attenuation_maps:

            # Get the path
            path = html.get_image_url(self.attenuation_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.attenuation_thumbnails[name])

            # Add
            self.attenuation_previews[name] = preview

    # -----------------------------------------------------------------

    def make_old_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Get the path
            path = html.get_image_url(self.old_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.old_thumbnails[name])

            # Add
            self.old_previews[name] = preview

    # -----------------------------------------------------------------

    def make_young_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Get the path
            path = html.get_image_url(self.young_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.young_thumbnails[name])

            # Add
            self.young_previews[name] = preview

    # -----------------------------------------------------------------

    def make_ionizing_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Get the path
            path = html.get_image_url(self.ionizing_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.ionizing_thumbnails[name])

            # Add
            self.ionizing_previews[name] = preview

    # -----------------------------------------------------------------

    def make_dust_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews for the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Get the path
            path = html.get_image_url(self.dust_thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.dust_thumbnails[name])

            # Add
            self.dust_previews[name] = preview

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

        # Get info for each map and the map labels
        infos = self.colour_info.values()
        labels = self.colour_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.colour_previews: thumbnails.append(html.center(self.colour_previews[label]))
            elif label in self.colour_thumbnails: thumbnails.append(html.center(self.colour_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Colour map"
        self.colour_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                        extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_ssfr_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of sSFR maps ...")

        # Get info for each map and the map labels
        infos = self.ssfr_info.values()
        labels = self.ssfr_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.ssfr_previews: thumbnails.append(html.center(self.ssfr_previews[label]))
            elif label in self.ssfr_thumbnails: thumbnails.append(html.center(self.ssfr_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "sSFR map"
        self.ssfr_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                      extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_tir_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of TIR maps ...")

        # Get info for each map and the map labels
        infos = self.tir_info.values()
        labels = self.tir_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.tir_previews: thumbnails.append(html.center(self.tir_previews[label]))
            elif label in self.tir_thumbnails: thumbnails.append(html.center(self.tir_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "TIR map"
        self.tir_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                     extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_attenuation_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of attenuation maps ...")

        # Get info for each map and the map labels
        infos = self.attenuation_info.values()
        labels = self.attenuation_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.attenuation_previews: thumbnails.append(html.center(self.attenuation_previews[label]))
            elif label in self.attenuation_thumbnails: thumbnails.append(html.center(self.attenuation_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Attenuation map"
        self.attenuation_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                             extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_old_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of old stellar maps ...")

        # Get info for each map and the map labels
        infos = self.old_info.values()
        labels = self.old_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.old_previews: thumbnails.append(html.center(self.old_previews[label]))
            elif label in self.old_thumbnails: thumbnails.append(html.center(self.old_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Old stellar map"
        self.old_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                     extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_young_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of young stellar maps ...")

        # Get info for each map and the map labels
        infos = self.young_info.values()
        labels = self.young_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.young_previews: thumbnails.append(html.center(self.young_previews[label]))
            elif label in self.young_thumbnails: thumbnails.append(html.center(self.young_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Young stellar map"
        self.young_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                       extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_ionizing_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of ionizing stellar maps ...")

        # Get info for each map and the map labels
        infos = self.ionizing_info.values()
        labels = self.ionizing_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.ionizing_previews: thumbnails.append(html.center(self.ionizing_previews[label]))
            elif label in self.ionizing_thumbnails: thumbnails.append(html.center(self.ionizing_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Ionizing stellar map"
        self.ionizing_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                          extra_column=thumbnails, extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    def make_dust_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the table of dust maps ...")

        # Get info for each map and the map labels
        infos = self.dust_info.values()
        labels = self.dust_info.keys()

        # Set the thumbnails
        thumbnails = []
        for label in labels:
            if label in self.dust_previews: thumbnails.append(html.center(self.dust_previews[label]))
            elif label in self.dust_thumbnails: thumbnails.append(html.center(self.dust_thumbnails[label]))
            else: thumbnails.append("")

        # Make the table
        label = "Dust map"
        self.dust_table = SimpleTable.from_composites(infos, css_class=self.table_class, labels=labels, label=label,
                                                    extra_column=thumbnails, extra_column_label=thumbnail_title)

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

        # Make javascripts urls
        javascript_paths = javascripts[:]
        javascript_paths.append(sortable_url)
        javascript_paths.append(preview_url)

        # Create the page
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths, javascript_path=javascript_paths, footing=updated_footing())

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
        self.page.saveto(self.maps_summary_html_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.maps_summary_html_page_path)

# -----------------------------------------------------------------
