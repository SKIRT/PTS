#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.wavelengthgrids Contains the WavelengthGridGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import bisect
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools.logging import log
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..core.emissionlines import EmissionLines

# -----------------------------------------------------------------

class WavelengthGridGenerator(FittingComponent):
    
    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(WavelengthGridGenerator, self).__init__()

        # -- Attributes --

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # ...

        # 17. Writing
        #self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Create the low-resolution wavelength grid
        self.create_low_res_wavelength_grid()

        # Create the high-resolution wavelength grid
        self.create_high_res_wavelength_grid()

    # -----------------------------------------------------------------

    def create_low_res_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the low-resolution wavelength grid ...")

        # Verify the grid parameters
        if self.config.wavelengths.npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
        if self.config.wavelengths.npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
        if self.config.wavelengths.min <= 0: raise ValueError("the shortest wavelength should be positive")
        if (self.config.wavelengths.min_zoom <= self.config.wavelengths.min
            or self.config.wavelengths.max_zoom <= self.config.wavelengths.min_zoom
            or self.config.wavelengths.max <= self.config.wavelengths.max_zoom):
                raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

        logmin = np.log10(float(self.config.wavelengths.min))
        logmax = np.log10(float(self.config.wavelengths.max))
        logmin_zoom = np.log10(float(self.config.wavelengths.min_zoom))
        logmax_zoom = np.log10(float(self.config.wavelengths.max_zoom))

        # Build the high- and low-resolution grids independently
        base_grid = np.logspace(logmin, logmax, num=self.config.wavelengths.npoints, endpoint=True, base=10)
        zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=self.config.wavelengths.npoints_zoom, endpoint=True, base=10)

        # Merge the two grids
        total_grid = []

        # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength < self.config.wavelengths.min_zoom: total_grid.append(wavelength)

        # Add the wavelengths of the high-resolution grid
        for wavelength in zoom_grid: total_grid.append(wavelength)

        # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength > self.config.wavelengths.max_zoom: total_grid.append(wavelength)

        # Add the central wavelengths of the filters used for normalizing the stellar components
        bisect.insort(total_grid, self.i1.centerwavelength())
        bisect.insort(total_grid, self.fuv.centerwavelength())

        # Create table for the low-resolution wavelength grid
        self.lowres_wavelength_grid = WavelengthGrid.from_wavelengths(total_grid)

    # -----------------------------------------------------------------

    def create_high_res_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the high-resolution wavelength grid ...")

        # Emission lines instance
        lines = EmissionLines()

        # Create the high-resolution wavelength grid
        self.highres_wavelength_grid = lines.create_wavelength_grid()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------

    