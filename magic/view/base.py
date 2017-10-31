#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.view.base Contains the ImageViewer base class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty, abstractmethod

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from pts.magic.view.html import make_replace_infs_by_nans, make_replace_negatives_by_nans
from pts.core.tools import html
from pts.core.tools import browser

# -----------------------------------------------------------------

stylesheet_url = "http://users.ugent.be/~sjversto/stylesheet.css"
background_color = "white"
page_style = "ugentstyle"

# -----------------------------------------------------------------

class ImageViewer(Configurable):

    """
    This class ...
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
        super(ImageViewer, self).__init__(*args, **kwargs)

        # The preloader
        self.preloader = None

        # The regions button
        self.regions_button = None

        # The page
        self.page = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImageViewer, self).setup(**kwargs)

        # Check
        if self.config.show and not self.config.page: raise ValueError("Cannot show the page when generating the page is disabled")

    # -----------------------------------------------------------------

    @abstractproperty
    def has_regions(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @property
    def theme_button(self):

        """
        This function ...
        :return:
        """

        classes = dict()
        classes["JS9Menubar"] = "data-backgroundColor"
        classes["JS9Menubar JS9Plugin"] = "data-backgroundColor"
        return html.make_theme_button(classes=classes)

    # -----------------------------------------------------------------

    def make_infs_button(self, image_name, display_name):

        """
        This function ...
        :param image_name:
        :param display_name:
        :return:
        """

        # Create nan/infs replacer button
        button_id = image_name + "nansinfs"
        replace_function_name = "replace_infs_nans_" + html.make_usable(image_name)
        replace_nans_infs = make_replace_infs_by_nans(display_name)

        # Create the button
        return html.make_script_button(button_id, "Replace infs", replace_nans_infs, replace_function_name)

    # -----------------------------------------------------------------

    def make_negatives_button(self, image_name, display_name):

        """
        This function ...
        :param image_name:
        :param display_name:
        :return:
        """

        # Create nan/infs replacer button
        button_id = image_name + "negatives"
        replace_function_name = "replace_negatives_nans_" + html.make_usable(image_name)
        replace_nans_negatives = make_replace_negatives_by_nans(display_name)

        # Create the button
        return html.make_script_button(button_id, "Replace negatives", replace_nans_negatives, replace_function_name)

    # -----------------------------------------------------------------

    @abstractmethod
    def _initialize_page(self):

        """
        Thisfunction ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def show(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Open the page
        browser.open_page(self.page)

# -----------------------------------------------------------------
