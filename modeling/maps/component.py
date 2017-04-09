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

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from ...magic.convolution.aniano import AnianoKernels
from ...core.tools.logging import log

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

        # The path to the maps/colours directory
        self.maps_colours_path = None

        # The path to the maps/tir directory
        self.maps_tir_path = None

        # The path to the maps/attenuation directory
        self.maps_attenuation_path = None

        # The path to the maps/old directory
        self.maps_old_path = None

        # The path to the maps/young directory
        self.maps_young_path = None

        # The path to the maps/ionizing directory
        self.maps_ionizing_path = None

        # The path to the maps/dust directory
        self.maps_dust_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MapsComponent, self).setup(**kwargs)

        # The path to the maps/colours directory
        self.maps_colours_path = fs.create_directory_in(self.maps_path, "colours")

        # The path to the maps/TIR directory
        self.maps_tir_path = fs.create_directory_in(self.maps_path, "tir")

        # The path to the maps/attenuation directory
        self.maps_attenuation_path = fs.create_directory_in(self.maps_path, "attenuation")

        # Set the path to the maps/old directory
        self.maps_old_path = fs.create_directory_in(self.maps_path, "old")

        # Set the path to the maps/young directory
        self.maps_young_path = fs.create_directory_in(self.maps_path, "young")

        # Set the path to the maps/ionizing directory
        self.maps_ionizing_path = fs.create_directory_in(self.maps_path, "ionizing")

        # Set the path to the maps/dust directory
        self.maps_dust_path = fs.create_directory_in(self.maps_path, "dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def aniano(self):

        """
        This function ...
        :return: 
        """

        return AnianoKernels()

    # -----------------------------------------------------------------

    @staticmethod
    def convert_to_same_unit(*frames, **kwargs):

        """
        This function ...
        :param frames:
        :param kwargs:
        :return: 
        """

        # Inform the user
        log.info("Converting frames to the same unit ...")

        # Get the unit of the first frame
        unit = frames[0].unit

        # Initialize list for converted frames
        new_frames = []

        # Convert all
        for frame in frames: new_frames.append(frame.converted_to(unit, **kwargs))

        # Return the new set of frames
        return new_frames

    # -----------------------------------------------------------------

    @staticmethod
    def rebin_to_highest_pixelscale(*frames):

        """
        This function ...
        :param frames: 
        :return: 
        """

        # Inform the user
        log.info("Rebinning frames to the coordinate system with the highest pixelscale ...")

        highest_pixelscale = None
        highest_pixelscale_wcs = None

        # Loop over the frames
        for frame in frames:

            wcs = frame.wcs
            if highest_pixelscale is None or wcs.average_pixelscale > highest_pixelscale:

                highest_pixelscale = wcs.average_pixelscale
                highest_pixelscale_wcs = wcs

        # Initialize list for rebinned frames
        new_frames = []

        # Rebin
        for frame in frames:

            if frame.wcs == highest_pixelscale_wcs: new_frames.append(frame.copy())
            else:
                if frame.unit.is_per_angular_area: rebinned = frame.rebinned(highest_pixelscale_wcs)
                else:
                    rebinned = frame.converted_to_corresponding_angular_area_unit()
                    rebinned.rebin(highest_pixelscale_wcs)
                new_frames.append(rebinned)

        # Return the rebinned frames
        return new_frames

    # -----------------------------------------------------------------

    def convolve_to_highest_fwhm(self, *frames):

        """
        This function ...
        :param frames: 
        :return: 
        """

        # Inform the user
        log.info("Convolving frames to the resolution of the frame with the highest FWHM ...")

        highest_fwhm = None
        highest_fwhm_filter = None

        # Loop over the frames
        for frame in frames:

            if highest_fwhm is None or frame.fwhm > highest_fwhm:

                highest_fwhm = frame.fwhm
                highest_fwhm_filter = frame.filter

        # Initialize list for convolved frames
        new_frames = []

        # Convolve
        for frame in frames:

            if frame.filter == highest_fwhm_filter: new_frames.append(frame.copy())
            else: new_frames.append(frame.convolved(self.aniano.get_kernel(frame.filter, highest_fwhm_filter)))

        # Return the convolved frames
        return new_frames

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
