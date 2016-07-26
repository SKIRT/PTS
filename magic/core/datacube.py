#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.datacube Contains the DataCube class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .image import Image
from .frame import Frame
from ...modeling.core.sed import SED, ObservedSED
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.tools.logging import log
from ..basics.mask import Mask
from ...core.basics.errorbar import ErrorBar

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

    @classmethod
    def from_files(cls, paths):

        """
        This function ...
        :param paths: paths of frames
        :return:
        """

        frames = [Frame.from_file(path) for path in paths]
        return cls.from_frames(frames)

    # -----------------------------------------------------------------

    @classmethod
    def from_frames(cls, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Create a datacube instance
        datacube = cls()

        # The indices of the frames, sorted on wavelength
        sorted_indices = sorted(range(len(frames)), key=lambda i: frames[i].filter.pivotwavelength())

        # The list of wavelengths
        wavelengths = []

        # Add the frames
        for index in sorted_indices:

            # Add the frame
            frame_name = "frame" + str(index)
            datacube.add_frame(frames[index], frame_name)

            # Add the wavelength
            wavelengths.append(frames[index].filter.pivotwavelength())

        # Create the wavelength grid
        datacube.wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths, unit="micron")

        # Return the datacube
        return datacube

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

    def local_sed(self, region):

        """
        This function ...
        :param region:
        :return:
        """

        # Initialize the SED
        sed = SED()

        # Create a mask from the region (or shape)
        mask = region.to_mask(self.xsize, self.ysize)

        # Loop over the wavelengths
        index = 0
        for wavelength in self.wavelengths():

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixels that belong to the region
            flux = np.sum(self.frames[frame_name][mask]) * self.unit

            # Add an entry to the SED
            sed.add_entry(wavelength, flux)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def pixel_sed(self, x, y):

        """
        This function ...
        :param x:
        :param y:
        :return:
        """

        # Initialize the SED
        #sed = SED()
        sed = ObservedSED()

        # Loop over the wavelengths
        index = 0
        for wavelength in self.wavelengths():

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Get the flux in the pixel
            flux = self.frames[frame_name][y, x] * self.unit

            # Add an entry to the SED
            #sed.add_entry(wavelength, flux)
            errorbar = ErrorBar(0.0, 0.0)
            sed.add_entry(self.frames[frame_name].filter, flux, errorbar)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def global_sed(self, mask=None):

        """
        This function ...
        :return:
        """

        # Determine the mask
        if isinstance(mask, basestring): inverse_mask = self.masks[mask].inverse()
        elif isinstance(mask, Mask): inverse_mask = mask.inverse()
        elif mask is None: inverse_mask = None
        else: raise ValueError("Mask must be string or Mask (or None) instead of " + str(type(mask)))

        # Initialize the SED
        sed = SED()

        # Loop over the wavelengths
        index = 0
        for wavelength in self.wavelengths():

            # Determine the name of the frame in the datacube
            frame_name = "frame" + str(index)

            # Calculate the total flux
            if mask is not None: total_flux = self.frames[frame_name].sum() * self.unit
            else: total_flux = np.sum(self.frames[frame_name][inverse_mask]) * self.unit

            # Add an entry to the SED
            sed.add_entry(wavelength, total_flux)

            # Increment the index
            index += 1

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def convolve_with_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with the " + str(fltr) + " filter ...")

        # Convert the datacube to a numpy array where wavelength is the third dimension
        array = self.asarray()

        # Calculate the observed image frame
        data = fltr.convolve(self.wavelengths, array)
        frame = Frame(data)

        # Set the unit of the frame
        frame.unit = self.unit

        # Return the resulting frame
        return frame

    # -----------------------------------------------------------------

    def convolve_with_filters(self, filters):

        """
        This function ...
        :param filters:
        :return:
        """

        # Inform the user
        log.info("Convolving the datacube with " + str(len(filters)) + " different filter ...")

        # Initialize list to contain the output frames
        frames = []

        # Convert the datacube to a numpy array where wavelength is the third dimension
        array = self.asarray()

        # Loop over the filters
        for fltr in filters:

            # Debugging
            log.debug("Convolving the datacube with the " + str(fltr) + " filter ...")

            # Calculate the observed image frame
            data = fltr.convolve(self.wavelengths, array)
            frame = Frame(data)

            # Set the unit of the frame
            frame.unit = self.unit

            # Set the filter
            frame.filter = fltr

            # Set the wcs
            frame.wcs = self.wcs

            # Add the frame to the list
            frames.append(frame)

        # Return the list of resulting frames
        return frames

    # -----------------------------------------------------------------

    def _end_as_array(self):

        """
        This function ...
        :return:
        """

        self._asarray = None

# -----------------------------------------------------------------
