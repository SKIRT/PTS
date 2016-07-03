#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.source Contains the Source class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .image import Image
from ...modeling.core.sed import SED

# -----------------------------------------------------------------

class DataCube(Image):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(DataCube, self).__init__()

        # The wavelength grid
        self.wavelength_grid = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, image_path, wavelength_grid):

        """
        This function ...
        :param image_path:
        :param wavelength_grid:
        :return:
        """

        # Call the corresponding base class function
        datacube = super(DataCube, cls).from_file(image_path, always_call_first_primary=False)

        # Check wavelength grid size
        assert len(wavelength_grid) == datacube.nframes

        # Set the wavelength grid
        datacube.wavelength_grid = wavelength_grid

        # Loop over the frames
        for i in range(datacube.nframes):

            # Frame name
            frame_name = "frame" + str(i)

            # Set the wavelength of the frame
            datacube.frames[frame_name].wavelength = datacube.wavelength_grid[i]

    # -----------------------------------------------------------------

    def wavelengths(self, unit=None, asarray=False, add_unit=True):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid.wavelengths(unit, asarray, add_unit)

    # -----------------------------------------------------------------

    def asarray(self, axis=3):

        """
        This function ...
        :return:
        """

        # Get a list that contains the frames
        frame_list = self.frames.as_list()

        # Stack the frames into a 3D numpy array
        if axis == 3: return np.dstack(frame_list)
        elif axis == 2: return np.hstack(frame_list)
        elif axis == 1: return np.vstack(frame_list)
        elif axis == 0: return np.stack(frame_list)
        else: raise ValueError("'axis' parameter should be integer 0-3")

    # -----------------------------------------------------------------

    def convert_to_fluxdensity(self, new_unit):

        """
        This function ...
        :return:
        """

        # Loop over the wavelengths
        index = 0
        for wavelength in self.wavelength_grid.wavelengths():

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Divide this frame by the wavelength in micron
            self.frames[frame_name] /= wavelength.to("micron")

            # Set the new unit
            self.frames[frame_name].unit = new_unit

            # Increment the index
            index += 1

    # -----------------------------------------------------------------

    def to_sed(self):

        """
        This function ...
        :return:
        """

        sed = SED()

        # Loop over the wavelengths
        index = 0
        for wavelength in self.wavelength_grid.wavelengths():

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            total_flux = self.frames[frame_name].sum() * self.frames[frame_name].unit

            sed.add_entry(wavelength, total_flux)

        # Return the SED
        return sed

# -----------------------------------------------------------------
