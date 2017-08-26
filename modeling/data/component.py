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
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from ...core.tools import strings

# -----------------------------------------------------------------

galex = "GALEX"
sdss = "SDSS"
twomass = "2MASS"
spitzer = "Spitzer"
wise = "WISE"
herschel = "Herschel"
planck = "Planck"
other = "Other"

halpha = "Halpha"

# -----------------------------------------------------------------

data_origins = [galex, sdss, halpha, twomass, spitzer, wise, herschel, planck, other]

# -----------------------------------------------------------------

def instrument_to_origin(instrument):

    """
    This function ...
    :param instrument:
    :return:
    """

    if instrument.lower() == "pacs": return "Herschel"
    elif instrument.lower() == "spire": return "Herschel"
    elif instrument.lower() == "lfi": return "Planck"
    elif instrument.lower() == "hfi": return "Planck"
    else: return strings.find_any_case(instrument, data_origins)

# -----------------------------------------------------------------

class DataComponent(GalaxyModelingComponent):
    
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
        super(DataComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Different origins
        self.data_origins = data_origins

        # The paths to the data/images/ directories for the different origins
        self.data_images_paths = dict()

        # Determine the path
        self.urls_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DataComponent, self).setup(**kwargs)

        # Set ...
        for origin in self.data_origins: self.data_images_paths[origin] = fs.create_directory_in(self.data_images_path, origin)

        # Set the urls path
        self.urls_path = fs.join(self.data_images_path, "urls.dat")

# -----------------------------------------------------------------
