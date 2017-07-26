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
from ...core.tools import html
from ..html.component import stylesheet_url, page_style, table_class, hover_table_class, top_title_size, title_size

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

        # Generate the page
        self.generate_page()

        # 5. Writing
        self.write()

        # Show
        self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AllMapsPageGenerator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Maps"

    # -----------------------------------------------------------------

    def generate_page(self):

        """
        Thisfunction ...
        :return:
        """

        body = ""
        #body += html.line + html.newline
        #body += str(self.images_table)


        #body += self.footing

        # Make page
        # Create contents
        contents = dict()
        contents["title"] = self.title
        contents["head"] = html.link_stylesheet_header_template.format(url=stylesheet_url)
        contents["body"] = body
        contents["style"] = page_style

        # Make page
        self.page = html.page_template.format(**contents)

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
