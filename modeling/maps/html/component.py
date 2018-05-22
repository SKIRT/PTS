#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.html.component Contains the abstract ComponentMapsPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import gc
from abc import ABCMeta, abstractproperty
import numpy as np

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ..selectioncomponent import MapsSelectionComponent
from .all import plots_name
from ....core.tools import html
from ....magic.tools.info import get_image_info_from_header_file
from ....core.basics.composite import SimplePropertyComposite
from ....core.tools.utils import lazyproperty
from ....core.tools import types
from ...html.component import stylesheet_url, page_style, sortable_url, preview_url
from ....magic.view.html import javascripts, css_scripts
from ....core.tools.html import HTMLPage, updated_footing, make_page_width
from ....core.tools import browser
from ....core.tools.stringify import tostr
from ....magic.core.image import Image
from ....core.tools import strings
from ....core.tools import numbers
from ....core.tools.parsing import real
from ....core.basics import containers
from ....core.tools import sequences

# -----------------------------------------------------------------

page_width = 600
thumbnail_title = "Thumbnail"
thumbnails_title = "Thumbnails"

# -----------------------------------------------------------------

class ComponentMapsPageGenerator(MapsSelectionComponent):

    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ComponentMapsPageGenerator, self).__init__(*args, **kwargs)

        # Set of invalid maps
        self.invalid = set()

        # Set of maps with too many zeroes
        self.zero = set()

        # Set of maps that are constant
        self.constant = set()

        # Set of maps with too many negative values
        self.negative = set()

        # Number of negative pixels per map
        self.nnegatives = dict()

        # The FWHMs and pixelscales
        self.fwhms = dict()
        self.pixelscales = dict()

        # The filtered map names
        self.filtered_map_names = None

        # Hidden
        self.hidden = None

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Filter the images
        if self.config.filter: self.filter_maps()

        # 3. Get the image info
        self.get_info()

        # 4. Categorize
        self.categorize()

        # 5. Make the thumbnails
        if self.config.thumbnails: self.make_thumbnails()

        # 6. Make the previews
        if self.config.previews: self.make_previews()

        # 7. Make the tables
        self.make_tables()

        # 8. Generate the page
        self.generate_page()

        # 9. Writing
        self.write()

        # 10. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    @property
    def nmaps(self):

        """
        This function ...
        :return:
        """

        return len(self.map_names)

    # -----------------------------------------------------------------

    @property
    def no_maps(self):

        """
        This function ...
        :return:
        """

        return self.nmaps == 0

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsPageGenerator, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

        # CHECK WHETHER THERE ARE MAPS
        if self.no_maps: raise RuntimeError("There are no maps yet. Run the appropriate command first to create the maps.")

        # Set filtered map names
        self.filtered_map_names = self.map_names

    # -----------------------------------------------------------------

    @abstractproperty
    def map_names(self):

        """
        Thisf unction ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def map_paths(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def filter_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filtering the maps ...")

        # Filter invalid, zero and negative pixels
        self.filter_invalid()

        # Filter resolution
        self.filter_resolution()

        # Set hidden maps
        self.set_hidden()

    # -----------------------------------------------------------------

    @lazyproperty
    def central_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.central_ellipse_factor

    # -----------------------------------------------------------------

    def filter_invalid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filtering invalid maps ...")

        # Loop over the maps
        for name in self.map_names:

            # Get the path
            path = self.map_paths[name]

            # Load the map
            #the_map = Frame.from_file(path)
            image = Image.from_file(path)
            the_map = image.primary
            nans = image.masks["nans"] if "nans" in image.masks else None
            negatives = image.masks["negatives"] if "negatives" in image.masks else None

            # Get the FWHM and pixelscale
            fwhm = the_map.fwhm
            pixelscale = the_map.average_pixelscale

            # Set FWHM and pixelscale
            self.fwhms[name] = fwhm
            self.pixelscales[name] = pixelscale

            # Check complete map for NaNS
            if the_map.all_nans or (nans is not None and nans.all_masked): # Old, new way
                log.warning("The '" + name + "' map contains only NaN values")
                self.invalid.add(name)

            # INFS SHOULDN'T HAPPEN
            # # Check complete map for infs
            # if the_map.all_infs:
            #     log.warning("The '" + name + "' map contains only infinities")
            #     self.invalid.add(name)

            # Check complete map for zeros
            if the_map.all_zeroes:
                log.warning("The '" + name + "' map contains only zeros")
                self.zero.add(name)

            # Check whether complete map is constant
            if the_map.is_constant:
                log.warning("The '" + name + "' map is constant everywhere")
                self.constant.add(name)

            # Check whether complete map is negative
            if the_map.all_negatives or (negatives is not None and negatives.all_masked): # Old, new way
                log.warning("The '" + name + "' map contains only negative values")
                self.negative.add(name)

            # Get center ellipse in pixel coordinates
            center_ellipse = self.central_ellipse.to_pixel(the_map.wcs)
            npixels_in_ellipse = the_map.npixels_in(center_ellipse)

            # Check NaNs within center ellipse
            #if the_map.nnans_in(center_ellipse) / npixels_in_ellipse > self.config.ninvalid_pixels_tolerance: # old
            if nans is not None and nans.nmasked_in(center_ellipse) / npixels_in_ellipse > self.config.ninvalid_pixels_tolerance:
                log.warning("The '" + name + "' map contains too many NaN values within the central ellipse")
                self.invalid.add(name)

            # INFS SHOULDN'T HAPPEN
            # Check infinities within center ellipse
            # if the_map.ninfs_in(center_ellipse) / npixels_in_ellipse > self.config.ninvalid_pixels_tolerance:
            #     log.warning("The '" + name + "' map contains too many infinities within the central ellipse")
            #     self.invalid.add(name)

            # Check zeros within center ellipse
            if the_map.nzeroes_in(center_ellipse) / npixels_in_ellipse > self.config.nzero_pixels_tolerance:
                log.warning("The '" + name + "' map contains too many zero values within the central ellipse")
                self.zero.add(name)

            # Check whether constant within center ellipse
            if the_map.is_constant_in(center_ellipse):
                log.warning("The '" + name + "' map is constant within the central ellipse")
                self.constant.add(name)

            # Check negatives within central ellipse
            #relative_nnegatives = the_map.nnegatives_in(center_ellipse) / npixels_in_ellipse # old
            relative_nnegatives = negatives.nmasked_in(center_ellipse) / npixels_in_ellipse if negatives is not None else 0.0
            if relative_nnegatives > self.config.nnegative_pixels_tolerance:
                log.warning("The '" + name + "' map contains too many negative values within the central ellipse")
                self.negative.add(name)

            # NEW
            self.nnegatives[name] = relative_nnegatives

            # Check central pixel
            if np.isnan(the_map.center_value) or (nans is not None and nans.data[the_map.pixel_center.y, the_map.pixel_center.x]): # old, new way
                log.warning("The '" + name + "' map has NaN at the center pixel")
                self.invalid.add(name)

            # INFS SHOULDN'T HAPPEN
            # # Check central pixel
            # if np.isinf(the_map.center_value):
            #     log.warning("The '" + name + "' map has infinity at the center pixel")
            #     self.invalid.add(name)

            # Clean
            gc.collect()

    # -----------------------------------------------------------------

    def filter_resolution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Filter the maps based on resolution ...")

        # Clear filtered
        self.filtered_map_names = []

        # Loop over the maps
        for name in self.map_names:

            # Check pixelscale
            if self.config.max_pixelscale is not None and self.pixelscales[name] > self.config.max_pixelscale: continue

            # Check FWHM
            if self.config.max_fwhm is not None and self.fwhms[name] > self.config.max_fwhm: continue

            # Add
            self.filtered_map_names.append(name)

    # -----------------------------------------------------------------

    def set_hidden(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the hidden maps ...")

        # Initialize hidden
        self.hidden = set()

        # Loop over the maps
        for name in self.filtered_map_names:

            # Skip invalid?
            if self.config.hide_invalid and name in self.invalid: self.hidden.add(name)

            # Skip zero
            if self.config.hide_zero and name in self.zero: self.hidden.add(name)

            # Skip constant
            if self.config.hide_constant and name in self.constant: self.hidden.add(name)

            # Skip negative
            if self.config.hide_negative and name in self.negative: self.hidden.add(name)

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
        kwargs["unit"] = False
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
        for name in self.filtered_map_names:

            # Skip
            if self.hidden is not None and name in self.hidden: continue

            # Get the path
            path = self.map_paths[name]

            # Get info
            info = get_image_info_from_header_file(name, path, **self.info_kwargs)

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            self.info[name] = info

    # -----------------------------------------------------------------

    @abstractproperty
    def component_map_methods(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def component_map_origins(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def categorize(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing ...")

        # Loop over the maps
        for name in self.filtered_map_names:

            # Skip
            if self.hidden is not None and name in self.hidden: continue

            # Get the methods
            methods = self.component_map_methods[name]
            reversed_methods = list(reversed(methods))

            # # CHECK FOR ISSUE WHERE LAST TWO ARE THE SAME -> NOT NECESSARY ANYMORE?
            # if len(reversed_methods) >= 2 and reversed_methods[0] == reversed_methods[1]:
            #     log.warning("Last two methods are the same " + tostr(methods))
            #     reversed_methods = reversed_methods[1:]

            #print(reversed_methods)

            # Loop over the methods
            parent_dict = self.categories
            #print(self.categories)
            for method in reversed_methods[:-1]:
                #print("METHOD", method)
                if method not in parent_dict: parent_dict[method] = dict()
                parent_dict = parent_dict[method]
                #print("PARENT DICT", parent_dict)

            # Add map name to the list
            last_method = methods[0]
            #print("LAST METHOD", last_method)
            # try: print("LAST METHOD", parent_dict[last_method])
            # except: pass

            if last_method not in parent_dict: parent_dict[last_method] = []
            parent_dict[last_method].append(name)

            # FIXING FOR YOUNG:
            # if types.is_dictionary(parent_dict):
            #     if last_method not in parent_dict: parent_dict[last_method] = []
            #     parent_dict[last_method].append(name)
            # else:
            #     parent_dict.append(name)

    # -----------------------------------------------------------------

    @property
    def plots_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.maps_html_path, plots_name)

    # -----------------------------------------------------------------

    @abstractproperty
    def component_plots_path(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def sub_name(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def make_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails ...")

        # Loop over the maps
        for name in self.filtered_map_names:

            # Skip
            if self.hidden is not None and name in self.hidden: continue

            # Determine full plot name
            if self.has_methods_for_sub_name(self.sub_name):
                method = self.component_map_methods[name][-1]
                if name.startswith(method): full_name = name
                else: full_name = method + "_" + name # add the method just as in the AllMapsPageGenerator
            else: full_name = name

            # Find the plot file
            path = fs.join(self.component_plots_path, full_name + ".png")

            # Check
            if not fs.is_file(path):
                log.error("Didn't find '" + path + "'")
                raise RuntimeError("The plot for the '" + name + "' map is not present. Run generate_all_maps_page first.")

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
        for name in self.filtered_map_names:

            # Skip
            if self.hidden is not None and name in self.hidden: continue

            # Get the path
            path = html.get_image_url(self.thumbnails[name])

            # Make preview
            preview = html.image_preview(path, self.thumbnails[name], title=name)

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

    @property
    def tostr_kwargs(self):

        """
        This function ...
        :return:
        """

        kwargs = dict()
        kwargs["round"] = True
        # kwargs["scientific"] = True  NO: let tostr decide
        kwargs["ndigits"] = 3
        kwargs["delimiter"] = ", "
        kwargs["html"] = True
        return kwargs

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

            # Group maps with same name except for factor
            if self.config.group_factors:

                # Initialize dictinoary
                grouped = dict()

                # Loop over all the maps
                for name in map_names:

                    if "__" not in name:
                        grouped[name] = name
                        continue
                    base_name, last_part = strings.split_at_last(name, "__")

                    if not numbers.is_number(last_part):
                        grouped[name] = name
                        continue

                    factor = real(last_part)
                    if base_name not in grouped: grouped[base_name] = dict()
                    grouped[base_name][factor] = name

                # Sort for each base name based on the factors values
                for base_name in grouped:
                    if not types.is_dictionary(grouped[base_name]): continue
                    grouped[base_name] = containers.ordered_by_key(grouped[base_name])

            # No grouping
            else: grouped = {name: name for name in map_names}

            # Fill the info and labels
            infos = []
            names = []
            origins = []
            #for name in map_names:
            for name in grouped:
                # name is actually base_name

                # Get one of the, or the only map name
                if types.is_dictionary(grouped[name]):
                    one_name = grouped[name][grouped[name].keys()[0]]
                    grouped_names = grouped[name].values()
                else:
                    one_name = name
                    grouped_names = [name]

                infos.append(self.info[one_name])

                names.append(grouped_names)

                origins_string = tostr(self.component_map_origins[one_name], delimiter=", ")
                origins.append(origins_string)

            # Set thumbnails
            thumbnails = []
            for grouped_names in names:
                thumbnails_group = ""
                for name in grouped_names:
                    if name in self.previews: thumbnails_group += self.previews[name]
                    elif name in self.thumbnails: thumbnails_group += self.thumbnails[name]
                    elif len(grouped_names) > 1: raise RuntimeError("Don't have a thumbnail for the '" + name + "' map")
                thumbnails_group = html.center(thumbnails_group)
                thumbnails.append(thumbnails_group)

            # Set extra columns
            if self.config.group_factors:

                # Set factors
                factors = []
                for name in grouped:
                    # name is actually base name
                    if types.is_dictionary(grouped[name]): factors_group = tostr(grouped[name].keys(), **self.tostr_kwargs)
                    else: factors_group = ""
                    factors.append(factors_group)

                if sequences.all_false(factors):
                    extra_columns = [thumbnails]
                    extra_column_labels = [thumbnail_title]

                else:
                    # Set fractions of negatives
                    nnegatives = []
                    for name in grouped:
                        # name is actually base_name
                        if types.is_dictionary(grouped[name]): nnegatives_group = tostr([self.nnegatives[single_name] * 100. for single_name in grouped[name].values()], **self.tostr_kwargs)
                        else: nnegatives_group = ""
                        nnegatives.append(nnegatives_group)

                    extra_columns = [factors, thumbnails, nnegatives]
                    extra_column_labels = ["Factors", thumbnails_title, "Negative values (%)"]
            else:
                extra_columns = [thumbnails]
                extra_column_labels = [thumbnail_title]

            # Make the table
            label = "Origins"
            if not self.config.group_factors:
                # WE KNOW THAT EACH ENTRY IN NAMES IS A LIST OF ONLY ONE
                strike_rows = []
                for group_names in names:
                    name = sequences.get_singleton(group_names)
                    if name in self.invalid: strike = True
                    elif name in self.constant: strike = True
                    elif name in self.zero: strike = True
                    elif name in self.negative: strike = True
                    else: strike = False
                    strike_rows.append(strike)
            else: strike_rows = None

            # Make the table
            self.tables[categories] = html.SimpleTable.from_composites(infos, css_class=self.table_class,
                                                                       labels=origins, label=label,
                                                                       extra_columns=extra_columns,
                                                                       extra_column_labels=extra_column_labels,
                                                                       strike_rows=strike_rows, tostr_kwargs=self.tostr_kwargs,
                                                                       wrap=False)

    # -----------------------------------------------------------------

    @abstractproperty
    def title(self):

        """
        This function ...
        :return:
        """

        pass

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

        # Loop over the different category lists
        for categories in self.tables:

            # Add table with title
            self.page += html.big(html.bold(tostr(categories, delimiter=" :: ")))
            self.page += html.newline
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

        # Write the ellipses
        self.write_ellipse()

        # Write the page
        self.write_page()

    # -----------------------------------------------------------------

    def write_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the central ellipse ...")

        # Determine filepath
        filename = "central_ellipse_" + self.sub_name + ".reg"
        filepath = fs.join(self.maps_html_path, filename)

        # Write the region
        self.central_ellipse.saveto(filepath)

    # -----------------------------------------------------------------

    @abstractproperty
    def page_path(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_page(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the page ...")

        # Save
        self.page.saveto(self.page_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open in browser
        with browser.serve_local_host(): browser.open_path(self.page_path)

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
