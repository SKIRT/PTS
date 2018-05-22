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
from ....core.basics.log import log
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....magic.core.frame import Frame
from ....core.basics.configurable import Configurable
from ....magic.core.list import FrameList
from ....core.tools.stringify import tostr
from ....magic.core.image import Image

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

        # The method name
        self.method_name = None

        # The frames
        self.frames = None

        # The error maps
        self.errors = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # The Galametz TIR calibration object
        self.galametz = GalametzTIRCalibration()

        # Region of interest
        self.region_of_interest = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

        # Get the method name
        self.method_name = kwargs.pop("method_name", None)

        # Get maps that have already been created
        if "maps" in kwargs: self.maps = kwargs.pop("maps")

        # Get region of interest
        self.region_of_interest = kwargs.pop("region_of_interest", None)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    @property
    def has_method_name(self):

        """
        This function ...
        :return:
        """

        return self.method_name is not None

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

            # Set the name
            name = tostr(fltr, delimiter="_")

            # Set the origins, also when already created map
            self.origins[name] = [fltr]

            # Set the methods
            if self.has_method_name: self.methods[name] = [self.method_name]

            # Check whether a TIR map is already present
            if name in self.maps:
                log.success("The " + name + " TIR map is already created: not creating it again")
                continue

            # Debugging
            log.debug("Making TIR map from the '" + str(fltr) + "' frame ...")

            # Get the parameters
            a, b = self.galametz.get_parameters_single_brightness(fltr)
            
            # Get the frame
            frame = self.frames[fltr]

            # Convert to neutral intrinsic surface brightness
            frame = frame.converted_to("W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)

            # Debugging
            #log.debug("Conversion factor: " + str(factor))

            # Calculate the TIR map in W/kpc2 (intrinsic surface brightness)
            logtir = a * np.log10(frame.data) + b

            # Create TIR map in linear units
            tir = Frame(10**logtir)

            # Set the unit and wcs
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True) # TIR can only be bolometric, right??
            #tir.unit = u("W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frame.wcs

            # Set other properties
            tir.distance = frame.distance
            tir.pixelscale = frame.pixelscale
            tir.psf_filter = frame.psf_filter
            tir.fwhm = frame.fwhm

            # Interpolate NaNs
            #tir.interpolate_nans(max_iterations=None)
            nans = tir.interpolate_nans_if_below(min_max_in=self.region_of_interest)
            tir.replace_negatives(0.0)
            image = Image()
            image.add_frame(tir, "tir")
            if nans is not None: image.add_mask(nans, "nans")

            # Set the TIR map
            self.maps[name] = image

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
