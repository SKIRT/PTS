#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.data Contains the DataPageGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import HTMLPageComponent, sortable_table_class
from ...core.tools import html
from ...dustpedia.core.properties import DustPediaProperties
from ...dustpedia.core.database import DustPediaDatabase
from ...dustpedia.core.database import get_account
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...magic.core.remote import RemoteFrame
from ...core.launch.pts import execute_pts_remote
from ...core.filter.filter import parse_filter
from ..html.component import stylesheet_url, page_style, sortable_url, preview_url
from ...magic.view.html import javascripts, css_scripts
from ...core.tools.html import HTMLPage, SimpleTable, updated_footing, make_page_width
from ...core.filter.broad import BroadBandFilter
from ...core.tools.utils import lazyproperty
from ...magic.tools.info import get_image_info_from_header_file, get_image_info_from_remote_header_file
from ..data.component import galex, sdss, twomass, spitzer, wise, herschel, planck, other, halpha, data_origins
from ...core.tools import strings
from ...magic.tools import headers
from ...core.basics.composite import SimplePropertyComposite

# -----------------------------------------------------------------

page_width = 600
thumbnail_title = "Thumbnail"

# -----------------------------------------------------------------

username, password = get_account()

# -----------------------------------------------------------------

class DataPageGenerator(HTMLPageComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DataPageGenerator, self).__init__(*args, **kwargs)

        # The image paths
        self.paths = None
        self.additional_paths = None

        # The properties table
        self.properties_table = None

        # The image paths
        self.image_paths = dict()

        # The additional image paths
        self.additional_image_paths = dict()

        # Thumbnails
        self.thumbnails = dict()

        # Previews
        self.previews = dict()

        # The table or tables
        self.table = None
        self.tables = OrderedDict()
        self.additional_table = None

    # -----------------------------------------------------------------

    @lazyproperty
    def image_names(self):

        """
        This function ...
        :return:
        """

        # Sorted on wavelength
        return sorted(self.image_paths.keys(), key=lambda name: parse_filter(name).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    @lazyproperty
    def all_image_paths(self):

        """
        This function ...
        :return:
        """

        paths = dict()
        paths.update(**self.image_paths)
        paths.update(**self.additional_image_paths)
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def all_image_names(self):

        """
        This function ...
        :return:
        """

        return sorted(self.all_image_paths.keys(), key=lambda name: parse_filter(name).wavelength.to("micron").value)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the image paths
        self.get_image_paths()

        # 3. Make plots
        self.make_plots()

        # 4. Make previews
        if self.config.thumbnails: self.make_thumbnails()

        # 5. Make previews
        if self.config.previews: self.make_previews()

        # 6. Make tables
        self.make_tables()

        # 7. Generaet the html
        self.generate()

        # 8. Write
        self.write()

        # Show the page
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DataPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def get_image_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the image paths ...")

        # Get paths and additional paths
        self.paths = self.get_data_image_paths_with_cached(lazy=True, not_origins=[other, halpha], attached=True)
        # paths = self.get_data_image_paths_with_cached_for_origin(, lazy=True)
        self.additional_paths = self.get_data_image_paths_with_cached(lazy=True, origins=[other, halpha], attached=True)

        # Present filenames
        filenames = {fs.name(filepath): name for name, filepath in self.paths.items()}

        # Get the image names and URLs from the DustPedia archive
        database = DustPediaDatabase()
        database.login(username, password)
        all_urls = database.get_image_names_and_urls(self.ngc_name_nospaces, error_maps=False)

        #print("filenames", filenames)

        # Check the paths
        # Loop over all DustPedia images
        for name in all_urls:

            #print(name)
            #continue

            # Skip error maps: normally shouldn't happen!
            if "_Error" in name: continue

            # Skip DSS
            if "DSS" in name and "SDSS" not in name: continue

            # Find name in paths dictionary
            if name in filenames: continue # OK

            # Parse the filter
            fltr = headers.get_filter(name)
            if fltr is None: raise ValueError("Could not find the filter for the '" + name + "' image")
            filter_name = str(fltr)

            # Error is missing
            if not self.config.fetch_missing: raise IOError("The '" + filter_name + "' image is missing from the data/images directory (and remote location)")

            # Warning
            log.warning("The '" + filter_name + "' image is missing from the data/images directory: fetching ...")

            # Not present

            # Add url to the dictionary
            if "pacs" in name.lower() or "spire" in name.lower(): origin = herschel

            # Not Herschel
            else: origin = strings.contains_which(name, data_origins)

            # Determine local directory for origin
            origin_path = fs.join(self.data_images_path, origin)

            # Determine the path to the image file
            path = fs.join(origin_path, name)

            # DOWNLOAD
            database.download_image_from_url(all_urls[name], path)

            # Set the path
            self.paths[filter_name] = path

    # -----------------------------------------------------------------

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make image plots
        self.make_image_plots()

        # Make additional image plots
        self.make_additional_image_plots()

    # -----------------------------------------------------------------

    @lazyproperty
    def session(self):

        """
        This function ...
        :return:
        """

        # Start python session on cache remote
        if self.has_cache_host:
            if self.config.use_session: return self.cache_remote.start_python_session(output_path=self.cache_remote.pts_temp_path, attached=True, new_connection_for_attached=True)
            else: return None
        else: return None

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        if self.has_cache_host:
            if self.config.use_session: return None
            else: return self.cache_remote
        else: return None

    # -----------------------------------------------------------------

    def make_image_plot(self, name, path, output_path):

        """
        This function ...
        :param name:
        :param path:
        :param output_path:
        :return:
        """

        # Debugging
        log.debug("Making a plot for the '" + name + "' image ...")

        # Local
        if fs.is_file(path):

            # Load the frame
            frame = Frame.from_file(path)

            # Downsample
            if frame.xsize > self.config.max_npixels or frame.ysize > self.config.max_npixels:
                downsample_factor = max(frame.xsize, frame.ysize) / float(self.config.max_npixels)
                frame.downsample(downsample_factor)

            # Plot
            frame.saveto_png(output_path, colours=self.config.colours, interval=self.config.interval, scale=self.config.scale,
                         alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Remote session
        elif self.session is not None:

            # Check if the file exists remotely
            if not self.session.is_file(path): raise ValueError("The file [" + path + "] does not exist locally or remotely")

            # Load the frame
            frame = RemoteFrame.from_remote_file(path, self.session)

            # Downsample
            if frame.xsize > self.config.max_npixels or frame.ysize > self.config.max_npixels:
                downsample_factor = max(frame.xsize, frame.ysize) / float(self.config.max_npixels)
                frame.downsample(downsample_factor)

            # Plot
            frame.saveto_png(output_path, colours=self.config.colours, interval=self.config.interval, scale=self.config.scale,
                         alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha)

        # Remote
        elif self.remote is not None:

            # Check if the file exists remotely
            if not self.remote.is_file(path): raise ValueError("The file [" + path + "] does not exist locally or remotely")

            # Temporary remote path
            temp_output_path = fs.join(self.remote.pts_temp_path, name + ".png")

            # Run the PTS command to create the PNG
            execute_pts_remote(self.remote, "fits_to_png", path, output=temp_output_path, show=False, show_output=True,
                               colours=self.config.colours, interval=self.config.interval, scale=self.config.scale,
                               alpha=self.config.alpha_method, peak_alpha=self.config.peak_alpha, max_npixels=self.config.max_npixels)

            # Check whether the remote file exists
            if not self.remote.is_file(temp_output_path): raise RuntimeError("Remote file does not exist")

            # Retrieve the file
            directory_path = fs.directory_of(output_path)
            new_name = fs.name(output_path)
            output_path = self.remote.download_file_to(temp_output_path, directory_path, remove=True, new_name=new_name)

        # Not found!
        else: raise ValueError("File '" + path + "' not found")

    # -----------------------------------------------------------------

    def make_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the image plots ...")

        # Loop over the images
        for name in self.paths:

            # Get the image path
            path = self.paths[name]

            # Determine the path
            output_path = fs.join(self.images_data_path, name + ".png")

            # Set the path
            self.image_paths[name] = output_path

            # CHECK WHETHER THE FILE ALREADY EXISTS
            if fs.is_file(output_path) and not self.config.replot:
                log.success("Plot of the '" + name + "' image is already present: not replotting")
                continue

            # Make image plot
            self.make_image_plot(name, path, output_path)

    # -----------------------------------------------------------------

    def make_additional_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the additional image plots ...")

        # Loop over the additional images
        for name in self.additional_paths:

            # Get image path
            path = self.additional_paths[name]

            # Determine the path
            output_path = fs.join(self.images_data_path, name + ".png")

            # Set the path
            self.additional_image_paths[name] = output_path

            # CHECK WHETHER THE FILE ALREADY EXISTS
            if fs.is_file(output_path) and not self.config.replot:
                log.success("Plot of the '" + name + "' image is already present: not replotting")
                continue

            # Make image plot
            self.make_image_plot(name, path, output_path)

    # -----------------------------------------------------------------

    def make_thumbnails(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making thumbnails ...")

        # Loop over the maps
        for name in self.all_image_names:

            # Get the path
            path = self.all_image_paths[name]

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
        for name in self.all_image_names:

            # Get the path
            path = self.all_image_paths[name]

            # Make preview
            preview = html.image_preview(path, self.thumbnails[name], title=name)

            # Add
            self.previews[name] = preview

    # -----------------------------------------------------------------

    @property
    def nadditional_images(self):

        """
        This function ...
        :return:
        """

        return len(self.additional_image_paths)

    # -----------------------------------------------------------------

    @property
    def has_additional(self):

        """
        This function ...
        :return:
        """

        return self.nadditional_images > 0

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the tables ...")

        # Group
        if self.config.group_observatories: self.make_tables_observatories()

        # Don't group
        else: self.make_table()

        # Additional images
        if self.has_additional: self.make_additional_table()

    # -----------------------------------------------------------------

    def make_tables_observatories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables for each observatory ...")

        # Get DustPedia data properties
        properties = DustPediaProperties()

        # Make tables
        tables = properties.create_tables(self.ngc_name, username, password)

        # Create HTML tables
        for observatory in tables:

            # # Set the thumbnails
            thumbnails = []

            # for label in labels[method]:
            for band in tables[observatory]["Band"]:

                # Determine filter name / image name
                image_name = str(BroadBandFilter.from_instrument_and_band(observatory, band))

                if image_name in self.previews: thumbnails.append(html.center(self.previews[image_name]))
                elif image_name in self.thumbnails: thumbnails.append(html.center(self.thumbnails[image_name]))
                else: thumbnails.append("")

            # Set extra columns
            extra_columns = [thumbnails]
            extra_column_labels = [thumbnail_title]

            # Create the table
            table = html.SimpleTable.from_table(tables[observatory], css_class=sortable_table_class,
                                                tostr_kwargs=self.tostr_kwargs, extra_columns=extra_columns,
                                                extra_column_labels=extra_column_labels)

            # Add the talbe
            self.tables[observatory] = table

    # -----------------------------------------------------------------

    def make_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making table for all observatories ...")

        # Get DustPedia data properties
        properties = DustPediaProperties()

        # Get table
        table = properties.create_table(self.ngc_name, username, password)

        # # Set the thumbnails
        thumbnails = []

        # for label in labels[method]:
        for instrument, band in zip(table["Instrument"], table["Band"]):

            # Determine filter name / image name
            image_name = str(BroadBandFilter.from_instrument_and_band(instrument, band))

            if image_name in self.previews: thumbnails.append(html.center(self.previews[image_name]))
            elif image_name in self.thumbnails: thumbnails.append(html.center(self.thumbnails[image_name]))
            else: thumbnails.append("")

        # Set extra columns
        extra_columns = [thumbnails]
        extra_column_labels = [thumbnail_title]

        # Create the table
        self.table = html.SimpleTable.from_table(table, css_class=sortable_table_class, tostr_kwargs=self.tostr_kwargs,
                                                 extra_columns=extra_columns, extra_column_labels=extra_column_labels)

    # -----------------------------------------------------------------

    def make_additional_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making table of additional images ...")

        infos = []
        thumbnails = []

        #print(self.additional_paths)

        # Loop over the additional images
        for name in self.additional_paths:

            # Get image info
            path = self.additional_paths[name]
            if fs.is_file(path): info = get_image_info_from_header_file(name, path, path=False)
            elif self.session is not None: info = get_image_info_from_remote_header_file(name, path, self.session, path=False)
            else: info = get_image_info_from_remote_header_file(name, path, self.remote, path=False)

            if name in self.previews: thumbnails.append(html.center(self.previews[name]))
            elif name in self.thumbnails: thumbnails.append(html.center(self.thumbnails[name]))
            else: thumbnails.append("")

            # Make property composite
            info = SimplePropertyComposite.from_dict(info)

            # Add info
            infos.append(info)

        # Set extra columns
        extra_columns = [thumbnails]
        extra_column_labels = [thumbnail_title]

        # Create the table
        self.additional_table = html.SimpleTable.from_composites(infos, css_class=sortable_table_class, tostr_kwargs=self.tostr_kwargs,
                                                 extra_columns=extra_columns, extra_column_labels=extra_column_labels)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        # Generate the page
        self.generate_page()

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

        #classes = dict()
        #classes["JS9Menubar"] = "data-backgroundColor"
        self.page += html.center(html.make_theme_button(images=False))
        self.page += html.newline

        # Add the table (if existing)
        if self.table is not None:

            self.page += "DUSTPEDIA DATASET"
            self.page += html.newline + html.newline
            self.page += self.table
            self.page += html.newline + html.newline

        # Add tables per observatory
        for observatory in self.tables:

            self.page += observatory.upper()
            self.page += html.newline + html.newline

            # Add the table
            self.page += self.tables[observatory]
            self.page += html.newline + html.newline

        # Add additional table
        if self.additional_table is not None:

            self.page += "ADDITIONAL DATA"
            self.page += html.newline + html.newline

            # Add the table
            self.page += self.additional_table

            # New lines
            self.page += html.newline + html.newline

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write status page
        self.write_page()

    # -----------------------------------------------------------------

    @property
    def page_path(self):

        """
        This function ...
        :return:
        """

        return self.data_page_path

# -----------------------------------------------------------------
