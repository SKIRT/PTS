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
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class PhotometryComponent(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PhotometryComponent, self).__init__(config)

        # -- Attributes --

        # Phot/temp directory
        self.phot_temp_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PhotometryComponent, self).setup()

        # Set ...
        self.phot_temp_path = fs.create_directory_in(self.phot_path, "temp")

# -----------------------------------------------------------------
