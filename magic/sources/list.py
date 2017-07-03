#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.list Contains the GalaxyList and StarList classes.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.mask import Mask

# -----------------------------------------------------------------

class GalaxyList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        self.principal = None
        self.companions = []
        self.other = []

    # -----------------------------------------------------------------

    def append(self, galaxy):

        """
        This function ...
        :param galaxy:
        :return:
        """

        if galaxy.principal:
            if self.principal is not None: raise ValueError("Principal galaxy is already set")
            self.principal = galaxy
        elif galaxy.companion: self.companions.append(galaxy)
        else: self.other.append(galaxy)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        self.principal = None
        self.companions = []
        self.other = []

    # -----------------------------------------------------------------

    def get_positions(self, wcs):

        """
        This function ...
        :return:
        """

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the galaxies
        for galaxy in self:

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(galaxy.pixel_position(wcs))

        # Return the list
        return positions

    # -----------------------------------------------------------------

    def get_principal_mask(self, frame_or_wcs):

        """
        This function ...
        :param frame_or_wcs:
        :return:
        """

        # Create a new mask with the dimensions of the frame
        mask = Mask.empty_like(frame_or_wcs)

        # Add the principal galaxy's mask to the total mask
        mask[self.principal.detection.cutout.y_slice, self.principal.detection.cutout.x_slice] = self.principal.detection.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def get_companion_mask(self, frame_or_wcs):

        """
        This function ...
        :return:
        """

        # Create a new mask with the dimension of the frame
        mask = Mask.empty_like(frame_or_wcs)

        # Loop over all companion galaxies
        for galaxy in self.companions:

            # Check if the galaxy has a source and add its mask to the total mask
            if galaxy.has_detection: mask[galaxy.detection.cutout.y_slice, galaxy.detection.cutout.x_slice] = galaxy.detection.mask

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        count = 0

        if self.principal is not None: count += 1
        count += len(self.companions)
        count += len(self.other)

        return count

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        if self.principal is not None: yield self.principal
        for companion in self.companions: yield companion
        for other in self.other: yield other

# -----------------------------------------------------------------

class StarList(object):

    """
    This function ...
    """

    def __init__(self):

        """
        This function ...
        """

        self.stars = []

    # -----------------------------------------------------------------

    def append(self, star):

        """
        This function ...
        :param star:
        :return:
        """

        self.stars.append(star)

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        self.stars = []

    # -----------------------------------------------------------------

    def get_positions(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Initialize a list to contain the object positions
        positions = []

        # Loop over the galaxies
        for skyobject in self.stars:

            # Calculate the pixel coordinate in the frame and add it to the list
            positions.append(skyobject.pixel_position(wcs))

        # Return the list
        return positions

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.stars)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for star in self.stars: yield star

# -----------------------------------------------------------------
