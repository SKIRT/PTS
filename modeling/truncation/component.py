#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.truncation.component Contains the TruncationComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

analytics_name = "analytics"
html_name = "html"

# -----------------------------------------------------------------

ellipse_page_filename = "ellipse.html"
significance_page_filename = "significance.html"

# -----------------------------------------------------------------

class TruncationComponent(GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(TruncationComponent, self).__init__(*args, **kwargs)

        # The path to the truncation/analytics directory
        self.truncation_analytics_path = None

        # HTML directory
        self.truncation_html_path = None

        # Ellipse page
        self.ellipse_page_path = None
        self.significance_page_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TruncationComponent, self).setup(**kwargs)

        # Set the path to the truncation/analytics directory
        self.truncation_analytics_path = fs.create_directory_in(self.truncation_path, analytics_name)

        # Set the path to the HTML directory
        self.truncation_html_path = fs.create_directory_in(self.truncation_path, html_name)

        # Page paths
        self.ellipse_page_path = fs.join(self.truncation_html_path, ellipse_page_filename)
        self.significance_page_path = fs.join(self.truncation_html_path, significance_page_filename)

# -----------------------------------------------------------------
