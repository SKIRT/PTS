#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.tir.single Contains the SingleBandTIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....magic.core.frame import Frame
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

possible_filters = ["IRAC I4", "MIPS 24mu", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250"]

# -----------------------------------------------------------------

class SingleBandTIRMapMaker(Configurable):

    """
    This class...
    """

    def __init__(self, config=None, ):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(SingleBandTIRMapMaker, self).__init__()

        # -- Attributes --

        # The frames
        self.frames = dict()

        # The error maps
        self.errors = dict()

        # The maps
        self.maps = dict()

        # The Galametz TIR calibration object
        self.galametz = GalametzTIRCalibration()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the data
        self.load_data()

        # 4. Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SingleBandTIRMapMaker, self).setup()

        # Get the distance
        self.distance = kwargs.pop("distance")

    # -----------------------------------------------------------------

    @lazyproperty
    def available_filters(self):

        """
        This function ...
        :return: 
        """

        filters = []

        # Loop over the colours
        for fltr in self.config.filters:

            # If no image is avilalbe for this filters, skip
            if not self.dataset.has_frame_for_filter(fltr): continue

            # otherwise, add to the list of filters
            filters.append(fltr)

        # Return the available filters
        return filters

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        # Loop over the filters
        for fltr in self.available_filters:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.dataset.get_frame_for_filter(fltr)
            self.frames[fltr] = frame

            # Load the error map
            errors = self.dataset.load_errormap_for_filter(fltr)
            self.errors[fltr] = errors

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Making the TIR maps ...")

        # Get the galaxy distance
        distance = self.galaxy_properties.distance

        # Loop over the frames
        for fltr in self.frames:

            # Debugging
            log.debug("Making TIR map from the '" + str(fltr) + "' frame ...")

            # Get the parameters
            a, b = self.galametz.get_parameters_single_brightness(fltr)

            # Convert to neutral intrinsic surface brightness
            frame = self.frames[fltr].convert_to("W/kpc2", density=True, distance=distance, brightness=True,
                                                 density_strict=True, brightness_strict=True)

            # Calculate the TIR map in W/kpc2 (intrinsic surface brightness)
            logtir = a * np.log(frame.data) + b
            tir = Frame(10**logtir)
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frame.wcs

            # Set the TIR map
            self.maps[fltr] = tir

# -----------------------------------------------------------------
