#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launch Contains the BestModelLauncher class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem, tables
from ...core.simulation.skifile import SkiFile
from ...core.tools.logging import log
from ..basics.instruments import FullInstrument
from ...magic.basics.vector import Position

# -----------------------------------------------------------------

class BestModelLauncher(AnalysisComponent):
    
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
        super(BestModelLauncher, self).__init__(config)

        # -- Attributes --

        # The path to the directory with the best model parameters
        self.best_path = None

        # The ski file for the best model
        self.ski = None

        # The wavelength grid
        self.wavelength_grid = None

        # The instruments
        self.instruments = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new BestModelLauncher instance
        launcher = cls()

        # Set the modeling path
        launcher.config.path = arguments.path

        # Return the new instance
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file describing the best model
        self.load_ski()

        # 3. Create the wavelength grid
        self.create_wavelength_grid()

        # 4. Create the instruments
        self.create_instruments()

        # 5. Adjust the ski file
        self.adjust_ski()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BestModelLauncher, self).setup()

        # The path to the directory with the best model parameters
        self.best_path = filesystem.join(self.fit_path, "best")

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file for the best fitting model ...")

        # Determine the path to the best model ski file
        path = filesystem.join(self.best_path, self.galaxy_name + ".ski")

        # Load the ski file
        self.ski = SkiFile(path)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Verify the grid parameters
        if self.config.wavelengths.npoints < 2: raise ValueError(
            "the number of points in the low-resolution grid should be at least 2")
        if self.config.wavelengths.npoints_zoom < 2: raise ValueError(
            "the number of points in the high-resolution subgrid should be at least 2")
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
        zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=self.config.wavelengths.npoints_zoom, endpoint=True,
                                base=10)

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

        # Create table for the wavelength grid
        self.wavelength_grid = tables.new([total_grid], names=["Wavelength"])

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # SKIRT:  incl.  azimuth PA
        # XY-plane	0	 0	    90
        # XZ-plane	90	 -90	0
        # YZ-plane	90	 0	    0

        # Determine the instrument properties
        distance = self.parameters.distance
        inclination = self.parameters.inclination
        azimuth = 0.0
        position_angle = self.parameters.disk.PA  # SAME PA AS THE DISK, BUT TILT THE BULGE W.R.T. THE DISK
        pixels_x = self.reference_wcs.xsize
        pixels_y = self.reference_wcs.ysize
        pixel_center = self.parameters.center.to_pixel(self.reference_wcs)
        # center = Position(0.5*pixels_x - pixel_center.x - 0.5, 0.5*pixels_y - pixel_center.y - 0.5) # when not convolved ...
        center = Position(0.5 * pixels_x - pixel_center.x - 1,
                          0.5 * pixels_y - pixel_center.y - 1)  # when convolved ...
        center_x = center.x * Unit("pix")
        center_y = center.y * Unit("pix")
        center_x = (center_x * self.reference_wcs.pixelscale.x.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * self.reference_wcs.pixelscale.y.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        field_x_angular = self.reference_wcs.pixelscale.x.to("deg/pix") * pixels_x * Unit("pix")
        field_y_angular = self.reference_wcs.pixelscale.y.to("deg/pix") * pixels_y * Unit("pix")
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Create the 'earth' instrument
        self.instruments["earth"] = FullInstrument(distance, inclination, azimuth, position_angle, field_x_physical,
                                                   field_y_physical, pixels_x, pixels_y, center_x, center_y)

        # Create the face-on instrument
        position_angle = Angle(90., "deg")
        self.instruments["faceon"] = FullInstrument(distance, 0.0, 0.0, position_angle, field_x_physical,
                                                    field_y_physical, pixels_x, pixels_y, center_x, center_y)

        # Create the edge-on instrument
        # azimuth = Angle(-90., "deg")
        self.instruments["edgeon"] = FullInstrument(distance, 90.0, 0.0, 0.0, field_x_physical, field_y_physical,
                                                    pixels_x, pixels_y, center_x, center_y)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Enable all writing options for analysis
        self.ski.enable_all_writing_options()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski()

        # Write the wavelength grid
        self.write_wavelength_grid()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save the ski file
        self.ski.saveto(self.analysis_ski_path)

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        self.wavelength_grid.rename_column("Wavelength", str(len(self.wavelength_grid)))  # Trick to have the number of wavelengths in the first line (required for SKIRT)
        tables.write(self.wavelength_grid, self.analysis_wavelengths_path, format="ascii")

# -----------------------------------------------------------------
