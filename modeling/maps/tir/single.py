#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.tir.single Contains the SingleBandTIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent
from ....core.basics.unit import parse_unit as u
from ....magic.misc.galametz import GalametzTIRCalibration
from ....magic.core.frame import Frame

# -----------------------------------------------------------------

possible_filters = ["IRAC I4", "MIPS 24mu", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250"]

# -----------------------------------------------------------------

class SingleBandTIRMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

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

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the data
        self.load_data()

        # 4. Make the maps
        self.make_maps()

        # 6. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SingleBandTIRMapMaker, self).setup()

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
            a, b = self.galametz.get_parameters_single(fltr)

            # Convert to neutral luminosity
            frame = self.frames[fltr].convert_to("W", density=True, distance=distance)

            # Calculate the TIR
            logtir = a * np.log(frame.data) + b
            tir = Frame(10**logtir)
            tir.unit = u("W", density=True)
            tir.wcs = frame.wcs

            # Set the TIR map
            self.maps[fltr] = tir

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the maps
        for fltr in self.maps:

            # Determine the path
            path = fs.join(self.maps_tir_path, str(fltr) + ".fits")

            # Save
            self.maps[fltr].saveto(path)

# -----------------------------------------------------------------
