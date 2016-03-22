#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.component Contains the MapsComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem

# -----------------------------------------------------------------

class MapsComponent(ModelingComponent):
    
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
        super(MapsComponent, self).__init__(config)

        # -- Attributes --

        # The path to the maps/solar directory
        self.maps_solar_path = None

        # The path to the maps/intermediate directory
        self.maps_intermediate_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MapsComponent, self).setup()

        # Set the output path (= the maps directory)
        self.config.output_path = self.maps_path

        # Set the path to the maps/solar directory
        self.maps_solar_path = filesystem.join(self.maps_path, "solar")

        # Set the path to the maps/intermediate directory
        self.maps_intermediate_path = filesystem.join(self.maps_path, "intermediate")

        # Create the maps/solar and maps/intermediate directories
        filesystem.create_directories([self.maps_solar_path, self.maps_intermediate_path])

# -----------------------------------------------------------------
