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

# -----------------------------------------------------------------

stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"

# -----------------------------------------------------------------

table_class = "realtable"

# -----------------------------------------------------------------

page_style = "ugentstyle"

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
    def status_page_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.html_status_path

    # -----------------------------------------------------------------

    @property
    def models_page_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.environment.html_path, "models.html")

    # -----------------------------------------------------------------

    @property
    def title(self):

        """
        This function ...
        :return:
        """

        return "Modeling of " + self.galaxy_name

    # -----------------------------------------------------------------

    @property
    def html_title(self):

        """
        This function ...
        :return:
        """

        # Create title
        return html.fontsize_template.format(size=20, text=html.underline_template.format(text="Modeling of " + self.galaxy_name))

    # -----------------------------------------------------------------

    @property
    def style(self):

        """
        This function ...
        :return:
        """

        return page_style

    # -----------------------------------------------------------------

    @abstractproperty
    def page_path(self):

        """
        This function ...
        :return:
        """

        pass

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
