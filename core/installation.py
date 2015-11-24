#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.installation
#

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import urllib

# -----------------------------------------------------------------

class SkirtInstallation(object):
    
    """
    This class ...
    """
    
    def __init__(self):
        
        """
        The constructor ...
        """
        
        self.has_qt = False
        
    # -----------------------------------------------------------------
    
    def run(self):
        
        """
        This function ...
        """

        self.check_qt()

        if not self.has_qt: self.install_qt()

        self.get()

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
