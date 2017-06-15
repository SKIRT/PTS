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

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....magic.core.frame import Frame
from ....core.basics.configurable import Configurable
from ....magic.core.list import FrameList
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

def make_map(frame, errors=None):

    """
    This function ...
    :param frame: 
    :param errors:
    :return: 
    """

    # Create list of frames and error maps
    frames = FrameList(frame)
    if errors is not None: error_maps = FrameList(errors)
    else: error_maps = None

    # Create the maker
    maker = SingleBandTIRMapMaker()

    # Run the map maker
    maker.run(frames=frames, errors=error_maps)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

class SingleBandTIRMapMaker(Configurable):

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
        super(SingleBandTIRMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frames
        self.frames = None

        # The error maps
        self.errors = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The Galametz TIR calibration object
        self.galametz = GalametzTIRCalibration()

        # The distance
        self.distance = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

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

        # Get the input
        self.frames = kwargs.pop("frames")
        self.errors = kwargs.pop("errors", None)
        self.distance = kwargs.pop("distance")

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        THis function ...
        :return: 
        """

        # Inform the user
        log.info("Making the TIR maps ...")

        # Loop over the frames
        for fltr in self.filters:

            # Debugging
            log.debug("Making TIR map from the '" + str(fltr) + "' frame ...")

            # Get the parameters
            a, b = self.galametz.get_parameters_single_brightness(fltr)

            # Get the frame
            frame = self.frames[fltr]

            # Convert to neutral intrinsic surface brightness
            frame = frame.converted_to("W/kpc2", density=True, distance=self.distance, brightness=True,
                                                 density_strict=True, brightness_strict=True)

            # Debugging
            #log.debug("Conversion factor: " + str(factor))

            # Calculate the TIR map in W/kpc2 (intrinsic surface brightness)
            logtir = a * np.log(frame.data) + b
            tir = Frame(10**logtir)
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True) # TIR can only be bolometric, right??
            #tir.unit = u("W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frame.wcs

            # Set the name
            name = tostr(fltr, delimiter="_")

            # Set the TIR map
            self.maps[name] = tir

            # Set the origins
            self.origins[name] = [fltr]

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------
