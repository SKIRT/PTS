#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization Contains the InputInitializer

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import the relevant PTS classes and modules
from ..core import ModelingComponent
from ...core.tools import inspection
from ...core.simulation.skifile import SkiFile

# -----------------------------------------------------------------

template_ski_path = os.path.join(inspection.pts_dat_dir("modeling"), "ski", "template.ski")

# -----------------------------------------------------------------

class InputInitializer(ModelingComponent):
    
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
        super(InputInitializer, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski_file = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the template ski file
        self.load_template()

        # 3. Create the wavelength grid
        self.create_wavelength_grid()

        # 4. Copy the input maps
        self.copy_maps()

        # 5. Place ski file
        self.place_ski_file()

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Open the template ski file
        self.ski_file = SkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Verify the grid parameters
        if self.config.wavelengths.npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
        if self.config.wavelengths.npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
        if self.config.wavelengths.min <= 0: raise ValueError("the shortest wavelength should be positive")
        if (self.config.wavelengths.min_zoom <= self.config.wavelengths.min
            or self.config.wavelengths.max_zoom <= self.config.wavelengths.min_zoom
            or self.config.wavelengths.max <= self.config.wavelengths.max_zoom):
                raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

        # Build the high- and low-resolution grids independently
        base_grid = np.logspace(float(self.config.wavelengths.min), float(self.config.wavelengths.max), num=self.config.wavelengts.npoints, endpoint=True, base=10.)
        zoom_grid = np.logspace(float(self.config.wavelengths.min_zoom), float(self.config.wavelengths.max_zoom), num=self.config.wavelengths.npoints_zoom, endpoint=True, base=10.)

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

        # Add the wavelengths of the bands of interest

    # -----------------------------------------------------------------

    def copy_maps(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def place_ski_file(self):

        """
        This function ...
        :return:
        """



# -----------------------------------------------------------------
