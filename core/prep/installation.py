#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.installation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import urllib

# Import the relevant PTS classes and modules
from ..basics.configurable import Configurable

# -----------------------------------------------------------------

class SkirtInstaller(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SkirtInstaller, self).__init__()

        ## Attributes

        self.has_qt = False
        
    # -----------------------------------------------------------------
    
    def run(self):
        
        """
        This function ...
        """

        # 1. Call the setup function
        self.setup()

        # 2.
        self.check_qt()

        # 3.
        if not self.has_qt: self.install_qt()

        # 4.
        self.get()

        # 5.
        self.install()

    # -----------------------------------------------------------------

    def check_qt(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def install_qt(self):

        """
        This function ...
        :return:
        """

        link = "http://download.qt.io/official_releases/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.gz"

        # Determine the path for the Qt source code
        path = os.path.join(None, "qt.tar.gz")

        # Create a
        urllib.urlretrieve (self.config.qt_link, path)

    # -----------------------------------------------------------------

    def get(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def install(self):

        """
        This function ...
        :return:
        """


# -----------------------------------------------------------------
