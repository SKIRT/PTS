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

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from .component import HTMLPageComponent, table_class
from ...core.tools import html
from ...dustpedia.core.properties import DustPediaProperties
from ...dustpedia.core.database import get_account
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...core.filter.filter import parse_filter
from ...magic.core.remote import RemoteFrame

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

        # The properties table
        self.properties_table = None

        # The images table
        self.images_table = None

        # The image paths
        self.image_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Setup
        self.setup(**kwargs)

        # Make plots
        self.make_plots()

        # Make tables
        self.make_tables()

        # Generaet the html
        self.generate()

        # Write
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

    def make_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making plots ...")

        # Make image plots
        self.make_image_plots()

    # -----------------------------------------------------------------

    def make_image_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the image plots ...")

        # Start python session on cache remote
        if self.has_cache_host:
            session = remote = None
            if self.config.use_session: session = self.cache_remote.start_python_session(output_path=self.cache_remote.pts_temp_path, attached=True, new_connection_for_attached=True)
            else: remote = self.cache_remote
        else: session = remote = None

        # Loop over the images
        paths = self.get_data_image_paths_with_cached(lazy=True)
        #print("paths", paths)
        for name in paths:

            path = paths[name]

            # Determine the path
            output_path = fs.join(self.images_path, name + ".png")

            # Local
            if fs.is_file(path): Frame.from_file(paths[name]).saveto_png(output_path)

            # Remote session
            elif session is not None: RemoteFrame.from_remote_file(path, session).saveto_png(output_path)

            # Remote
            elif remote is not None:

                # Temporary remote path
                temp_output_path = fs.join(remote.pts_temp_path, name + ".png")

                # Run the PTS command to create the PNG
                remote.execute_pts("fits_to_png", path, output=temp_output_path, show_output=True)

                # Check whether the remote file exists
                if not remote.is_file(temp_output_path): raise RuntimeError("Remote file does not exist")

                # Retrieve the file
                output_path = remote.download_file_to(temp_output_path, self.images_path, remove=True)

            # Not found!
            else: raise ValueError("File '" + path + "' not found")

            # Add the path
            self.image_paths[name] = output_path

    # -----------------------------------------------------------------

    def make_tables(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making tables ...")

        # Make properties table
        self.make_properties_table()

        # Make images table
        self.make_images_table()

    # -----------------------------------------------------------------

    def make_properties_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the data properties table ...")

        # Get DustPedia username and password
        username, password = get_account()
        properties = DustPediaProperties()

        # Get the table
        table = properties.create_table(self.ngc_name, username, password)

        # Create the html table
        self.properties_table = html.SimpleTable(table.as_tuples(), table.column_names, css_class=table_class, tostr_kwargs=self.tostr_kwargs)

    # -----------------------------------------------------------------

    @property
    def image_names(self):

        """
        This function ...
        :return:
        """

        return sorted(self.image_paths.keys(), key=lambda filter_name: parse_filter(filter_name).wavelength)

    # -----------------------------------------------------------------

    def make_images_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the images table ...")

        # Create the code for each image
        cells = []

        for name in self.image_names:

            # Get the relative image path
            image_path = self.image_paths[name]
            relative_path = fs.relative_to(image_path, self.environment.html_path)

            # Set title
            text = ""
            text += html.center_template.format(text=html.bold_template.format(text=name))
            text += html.newline + html.newline
            text += html.image(relative_path, width=400, height=400)

            cells.append(text)

        # Create table from cells
        self.images_table = html.SimpleTable.rasterize(cells, 4, tostr_kwargs=self.tostr_kwargs)

    # -----------------------------------------------------------------

    def generate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating the HTML ...")

        body = self.heading + html.newline
        body += str(self.properties_table)

        body += html.line + html.newline

        body += str(self.images_table)

        body += self.footing

        # Make page
        self.make_page(body)

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
