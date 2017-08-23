#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html.dust Contains the DustMapsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ..selectioncomponent import MapsSelectionComponent
from .all import plots_name, dust_name
from ....core.tools import html
from ....magic.tools.info import get_image_info_from_header_file
from ....core.basics.composite import SimplePropertyComposite
from ....core.tools.utils import lazyproperty
from ....core.tools import types
from ...html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size, sortable_url, preview_url
from ....magic.view.html import javascripts, css_scripts, JS9Spawner, make_replace_nans_infs
from ....core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ....core.tools import browser
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

page_width = 600
thumbnail_title = "Thumbnail"

# -----------------------------------------------------------------

class DustMapsPageGenerator(MapsSelectionComponent):

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
        super(DustMapsPageGenerator, self).__init__(*args, **kwargs)

        # The info
        self.info = dict()

        # The thumbnails
        self.thumbnails = dict()

        # The previews
        self.previews = dict()

        # The categories
        self.categories = dict()

        # The tables
        self.tables = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Get the image info
        self.get_info()

        # Categorize
        self.categorize()

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
        super(DustMapsPageGenerator, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

    # -----------------------------------------------------------------

    @property
    def map_names(self):

        """
        Thisf unction ...
        :return:
        """

        return self.dust_map_names

    # -----------------------------------------------------------------

    @property
    def map_paths(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_paths

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
        log.info("Getting the map info ...")

        # Loop over the maps
        for name in self.map_names:

            # Get the path
            path = self.map_paths[name]

            # Get info
            info = get_image_info_from_header_file(name, path, **self.info_kwargs)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.info[name] = info

    # -----------------------------------------------------------------

    def categorize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing ...")

        # Loop over the maps
        for name in self.map_names:

            # Get the methods
            methods = self.dust_map_methods[name]
            reversed_methods = list(reversed(methods))

            # Loop over the methods
            parent_dict = self.categories
            for method in reversed_methods[:-1]:

                if method not in parent_dict: parent_dict[method] = dict()
                parent_dict = parent_dict[method]

            # Add map name to the list
            last_method = methods[0]
            if last_method not in parent_dict: parent_dict[last_method] = []
            parent_dict[last_method].append(name)

        # Print
        #print(self.categories)

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
    def dust_plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.plots_path, dust_name)

    # -----------------------------------------------------------------

    def make_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails ...")

        # Loop over the maps
        for name in self.map_names:

            # Find the plot file
            path = fs.join(self.dust_plots_path, name + ".png")

            # Check
            if not fs.is_file(path): continue

            # Make image
            image = html.image(path, alttext=name, height=self.config.thumbnail_height)

            # Add the image
            self.thumbnails[name] = image

    # -----------------------------------------------------------------

    def make_previews(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making previews ...")

        # Loop over the maps
        for name in self.map_names:

            # Get the path
            path = html.get_image_url(self.thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.thumbnails[name])

            # Add
            self.previews[name] = preview

    # -----------------------------------------------------------------

    @property
    def table_class(self):

        """
        This function ...
        :return:
        """

        return "sortable"

    # -----------------------------------------------------------------

    @lazyproperty
    def category_lists(self):

        """
        This function ...
        :return:
        """

        return get_paths_through_dictionary(self.categories)

    # -----------------------------------------------------------------

    @lazyproperty
    def category_tuples(self):

        """
        This function ...
        :return:
        """

        return get_paths_through_dictionary(self.categories, as_tuples=True)

    # -----------------------------------------------------------------

    def get_map_paths_for_categories(self, categories):

        """
        This function ...
        :param categories:
        :return:
        """

        return get_hierarchical_dictionary_element(self.categories, categories)

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the tables ...")

        # Loop over the different category lists
        for categories in self.category_tuples:

            # Get the map paths
            map_names = self.get_map_paths_for_categories(categories)

            # Fill the info and labels
            infos = []
            labels = []
            for name in map_names:
                infos.append(self.info[name])
                labels.append(name)

            # Set thumbnails
            thumbnails = []
            for label in labels:
                if label in self.previews: thumbnails.append(html.center(self.previews[label]))
                elif label in self.thumbnails: thumbnails.append(html.center(self.thumbnails[label]))
                else: thumbnails.append("")

            # Make the table
            label = "Dust map"
            self.tables[categories] = html.SimpleTable.from_composites(infos, css_class=self.table_class,
                                                                   labels=labels, label=label,
                                                                   extra_column=thumbnails,
                                                                   extra_column_label=thumbnail_title)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Dust maps"

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        This function ...
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
        self.page = HTMLPage(self.title, css=css, style=page_style, css_path=css_paths,
                             javascript_path=javascript_paths, footing=updated_footing())

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(classes=classes, images=False))
        self.page += html.newline

        for categories in self.tables:

            self.page += tostr(categories)
            self.page += html.newline
            self.page += self.tables[categories]
            self.page += html.newline

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
        self.page.saveto(self.dust_maps_html_page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        browser.open_path(self.dust_maps_html_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return:
        """

        return None

# -----------------------------------------------------------------

def get_hierarchical_dictionary_element(dictionary, path):

    """
    Thisf unction ...
    :param dictionary:
    :param path:
    :return:
    """

    #print("PATH", path)
    parent_dict = dictionary
    for element in path[:-1]: parent_dict = parent_dict[element]
    return parent_dict[path[-1]]

# -----------------------------------------------------------------

def get_paths_through_dictionary(dictionary, as_tuples=False):

    """
    This function ...
    :param dictionary:
    :param as_tuples:
    :return:
    """

    paths = []

    for key in dictionary:

        path = [key]

        if types.is_dictionary(dictionary[key]):

            key_paths = get_paths_through_dictionary(dictionary[key])
            complete_paths = [path + key_path for key_path in key_paths]
            paths.extend(complete_paths)

        else: paths.append(path)

    # Return the paths
    if as_tuples: return [tuple(path) for path in paths]
    else: return paths

# -----------------------------------------------------------------
