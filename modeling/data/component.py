#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.component Contains the DataComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class DataComponent(ModelingComponent):
    
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
        super(DataComponent, self).__init__(config)

        # -- Attributes --

        # The path to the galaxy info file
        self.galaxy_info_path = None

        # Different origins
        self.data_origins = ["GALEX", "SDSS", "Halpha", "2MASS", "Spitzer", "WISE", "Herschel", "Planck"]

        # The paths for the images from different origins
        self.data_galex_path = None
        self.data_sdss_path = None
        self.data_halpha_path = None
        self.data_2mass_path = None
        self.data_spitzer_path = None
        self.data_wise_path = None
        self.data_herschel_path = None
        self.data_planck_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataComponent, self).setup()

        # Set the path to the galaxy info file
        self.galaxy_info_path = fs.join(self.data_path, "info.dat")

        # The paths for the images from different origins
        self.data_galex_path = fs.create_directory_in(self.data_path, "GALEX")
        self.data_sdss_path = fs.create_directory_in(self.data_path, "SDSS")
        self.data_ha_path = fs.create_directory_in(self.data_path, "Halpha")
        self.data_2mass_path = fs.create_directory_in(self.data_path, "2MASS")
        self.data_spitzer_path = fs.create_directory_in(self.data_path, "Spitzer")
        self.data_wise_path = fs.create_directory_in(self.data_path, "WISE")
        self.data_herschel_path = fs.create_directory_in(self.data_path, "Herschel")

# -----------------------------------------------------------------
