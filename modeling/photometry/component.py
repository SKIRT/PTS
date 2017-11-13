#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.component Contains the PhotometryComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class PhotometryComponent(GalaxyModelingComponent):
    
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
        super(PhotometryComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The path to the flux differences table
        self.phot_differences_path = None

        # The path to the flux errors table
        self.phot_errors_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PhotometryComponent, self).setup(**kwargs)

        # Set ...
        self.phot_differences_path = fs.join(self.phot_path, "differences.dat")

        # Set ...
        self.phot_errors_path = fs.join(self.phot_path, "errors.dat")

# -----------------------------------------------------------------
