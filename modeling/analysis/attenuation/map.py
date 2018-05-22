#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.attenuation.map Contains the AttenuationMapAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AttenuationAnalysisComponent
from ....core.basics.log import log
from ....core.tools import filesystem as fs
from ....core.simulation.wavelengthgrid import WavelengthGrid
from ....magic.core.datacube import DataCube
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

class AttenuationMapAnalyser(AttenuationAnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AttenuationMapAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The wavelength grid
        self.wavelength_grid = None

        # The total FUV frame
        self.total = None

        # The transparent FUV frame
        self.transparent = None

        # The attenuation map
        self.attenuation = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the wavelength grid
        self.load_wavelength_grid()

        # Load the frames
        self.load_frames()

        # Make the map
        self.make_map()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AttenuationMapAnalyser, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid file produced by SKIRT ...")

        # Determine the path to the wavelength grid file in the analysis simulation output path
        wavelengths_path = fs.join(self.analysis_run.out_path, self.galaxy_name + "_wavelengths.dat")

        # Load the wavelength grid as a table
        self.wavelength_grid = WavelengthGrid.from_skirt_output(wavelengths_path)

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the frames ...")

        # Load the total FUV frame
        self.load_total()

        # Load the transparent FUV frame
        self.load_transparent()

    # -----------------------------------------------------------------

    def load_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading ...")

        # Detemrine the path to the datacube
        path = fs.join(self.analysis_run.out_path, self.galaxy_name + "_earth_total.fits")

        # The total datacube
        datacube = DataCube.from_file(path, self.wavelength_grid)

        # Get the total frame at FUV wavelength
        self.total = datacube.get_frame_for_wavelength(self.fuv_filter.pivot)

    # -----------------------------------------------------------------

    def load_transparent(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the datacube
        path = fs.join(self.analysis_run.out_path, self.galaxy_name + "_earth_transparent.fits")

        # The transparent datacube
        datacube = DataCube.from_file(path, self.wavelength_grid)

        # Get the transparent frame at FUV wavelength
        self.transparent = datacube.get_frame_for_wavelength(self.fuv_filter.pivot)

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the FUV attenuation map ...")

        # Calculate the attenuations
        self.attenuation = Frame(-2.5 * np.log10(self.total.data / self.transparent.data))

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write attenuation map
        self.write_map()

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the modelled FUV attenuation map ...")

        # Determine the path and save
        path = fs.join(self.attenuation_map_path, "FUV attenuation.fits")
        self.attenuation.saveto(path)

# -----------------------------------------------------------------
