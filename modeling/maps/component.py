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
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class MapsComponent(GalaxyModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(MapsComponent, self).__init__(config, interactive)

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

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsComponent, self).setup(**kwargs)

        # Set the path to the maps/dust directory
        self.maps_dust_path = fs.create_directory_in(self.maps_path, "dust")

        # Set the path to the maps/old directory
        self.maps_old_path = fs.create_directory_in(self.maps_path, "old")

        # Set the path to the maps/young directory
        self.maps_young_path = fs.create_directory_in(self.maps_path, "young")

        # Set the path to the maps/ionizing directory
        self.maps_ionizing_path = fs.create_directory_in(self.maps_path, "ionizing")

# -----------------------------------------------------------------

def get_dust_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "dust")

# -----------------------------------------------------------------

def get_dust_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_dust_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_old_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "old")

# -----------------------------------------------------------------

def get_old_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_old_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_young_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "young")

# -----------------------------------------------------------------

def get_young_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_young_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------

def get_ionizing_stars_maps_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "maps", "ionizing")

# -----------------------------------------------------------------

def get_ionizing_stellar_map_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.files_in_path(get_ionizing_stars_maps_path(modeling_path), extension="fits", returns="name")

# -----------------------------------------------------------------
