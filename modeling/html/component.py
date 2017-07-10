#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.html.component Contains the HTMLPageComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import webbrowser
from abc import ABCMeta, abstractmethod, abstractproperty

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import html
from ...core.tools import filesystem as fs
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import time

# -----------------------------------------------------------------

stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"

# -----------------------------------------------------------------

table_class = "realtable"

# -----------------------------------------------------------------

page_style = "ugentstyle"

# -----------------------------------------------------------------

top_title_size = 24
title_size = 20

# -----------------------------------------------------------------

class HTMLPageComponent(GalaxyModelingComponent):

    """
    This function ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HTMLPageComponent, self).__init__(*args, **kwargs)

        # Page
        self.page = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(HTMLPageComponent, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @abstractmethod
    def make_tables(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractmethod
    def make_plots(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def images_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_images_path

    # -----------------------------------------------------------------

    @property
    def status_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_status_path

    # -----------------------------------------------------------------

    @property
    def status_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.status_page_path)

    # -----------------------------------------------------------------

    @property
    def data_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, "data.html")

    # -----------------------------------------------------------------

    @property
    def data_page_name(self):

        """
        THisn function ...
        :return:
        """

        return fs.name(self.data_page_path)

    # -----------------------------------------------------------------

    @property
    def model_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, "model.html")

    # -----------------------------------------------------------------

    @property
    def model_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.model_page_path)

    # -----------------------------------------------------------------

    @property
    def maps_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, "maps.html")

    # -----------------------------------------------------------------

    @property
    def maps_page_name(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.maps_page_path)

    # -----------------------------------------------------------------

    @property
    def top_title(self):

        """
        This function ...
        :return:
        """

        # Create title
        return "High resolution 3D radiative transfer modeling of DustPedia galaxies"

    # -----------------------------------------------------------------

    @property
    def html_top_title(self):

        """
        This function ...
        :return:
        """

        return html.fontsize_template.format(size=top_title_size, text=html.underline_template.format(text=self.top_title))

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        #return "Modeling of " + self.galaxy_name
        return self.galaxy_name + " (" + self.ngc_name_nospaces + ")"

    # -----------------------------------------------------------------

    @property
    def html_title(self):

        """
        This function ...
        :return:
        """

        # Create title
        return html.fontsize_template.format(size=title_size, text=html.underline_template.format(text=self.title))

    # -----------------------------------------------------------------

    @property
    def heading(self):

        """
        This function ...
        :return:
        """

        heading = ""
        heading += self.html_top_title + html.newline + html.newline + html.line + html.newline
        heading += self.html_title + html.newline
        return heading

    # -----------------------------------------------------------------

    @property
    def footing(self):

        """
        This function ...
        :return:
        """

        text = "Last updated " + time.pretty_date().lower()

        footing = ""
        footing += html.newline + html.line + html.newline + html.center_template.format(text=text)
        return footing

    # -----------------------------------------------------------------

    @property
    def style(self):

        """
        This function ...
        :return:
        """

        return page_style

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

    def make_page(self, body):

        """
        This function ...
        :param body:
        :return:
        """

        # Create contents
        contents = dict()
        contents["title"] = self.title
        contents["head"] = html.link_stylesheet_header_template.format(url=stylesheet_url)
        contents["body"] = body
        contents["style"] = self.style

        # Make page
        self.page = html.page_template.format(**contents)

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

        # Write
        fs.write_text(self.page_path, self.page)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the page ...")

        # Open
        webbrowser._tryorder = ["safari"]
        webbrowser.open(self.page_path, new=2)

# -----------------------------------------------------------------
