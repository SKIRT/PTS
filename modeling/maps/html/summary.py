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
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size, sortable_url
from ....core.tools import filesystem as fs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import html
from ....magic.view.html import javascripts, css_scripts, JS9Spawner, make_replace_nans_infs
from ....core.tools import browser
from ....magic.tools.info import get_image_info_strings, get_image_info
from ....core.basics.table import SmartTable
from ....core.basics.composite import SimplePropertyComposite

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
        fs.set_nallowed_open_files(self.config.nopen_files)

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
            info = get_image_info(name, self.ssfr_maps[name], path=False)

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
            info = get_image_info(name, self.tir_maps[name], path=False)

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
            info = get_image_info(name, self.attenuation_maps[name], path=False)

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
            info = get_image_info(name, self.old_maps[name], path=False)

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
            info = get_image_info(name, self.young_maps[name], path=False)

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
            info = get_image_info(name, self.ionizing_maps[name], path=False)

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
            info = get_image_info(name, self.dust_maps[name], path=False)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Colour map")

        # Make
        self.colour_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="sSFR map")

        # Make
        self.ssfr_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="TIR map")

        # Make
        self.tir_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Attenuation map")

        # Make
        self.attenuation_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Old stellar map")

        # Make
        self.old_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Young stellar map")

        # Make
        self.young_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Ionizing stellar map")

        # Make
        self.ionizing_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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

        # Create the table
        table = SmartTable.from_composites(*infos, labels=labels, label="Dust map")

        # Make
        self.dust_table = SimpleTable(table.as_tuples(), header_row=table.column_names, css_class=self.table_class)

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
