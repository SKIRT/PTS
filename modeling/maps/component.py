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

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class MapsComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(MapsComponent, self).__init__(config)

        # -- Attributes --

        # The path to the maps/dust directory
        self.maps_dust_path = None

        # The path to the maps/old directory
        self.maps_old_path = None

        # The path to the maps/young directory
        self.maps_young_path = None

        # The path to the maps/ionizing directory
        self.maps_ionizing_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MapsComponent, self).setup()

        # Set the path to the maps/dust directory
        self.maps_dust_path = fs.create_directory_in(self.maps_path, "dust")

        # Set the path to the maps/old directory
        self.maps_old_path = fs.create_directory_in(self.maps_path, "old")

        # Set the path to the maps/young directory
        self.maps_young_path = fs.create_directory_in(self.maps_path, "young")

        # Set the path to the maps/ionizing directory
        self.maps_ionizing_path = fs.create_directory_in(self.maps_path, "ionizing")

# -----------------------------------------------------------------
