#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.dust.tir.multi Contains the MultiBandTIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....magic.core.frame import linear_combination
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....core.tools import sequences
from ....core.basics.configurable import Configurable
from ....magic.core.list import FrameList

# -----------------------------------------------------------------

def make_map(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    # Create the
    maker = MultiBandTIRMapMaker()

    # Set input
    frames = FrameList(*args)
    errors = kwargs.pop("errors", None)
    if errors is not None: errors = FrameList(*errors)
    distance = kwargs.pop("distance")
    lengths = [len(args)]

    # Run the maker
    maker.run(frames=frames, errors=errors, lengths=lengths, distance=distance)

    # Return the single map
    return maker.single_map

# -----------------------------------------------------------------

def make_maps(*args, **kwargs):

    """
    This function ...
    :param args: 
    :param kwargs:
    :return: 
    """

    # Create the maker
    maker = MultiBandTIRMapMaker()

    # Set input
    frames = FrameList(*args)
    errors = kwargs.pop("errors", None)
    if errors is not None: errors = FrameList(**errors)
    distance = kwargs.pop("distance")

    # Run the maker
    maker.run(frames=frames, errors=errors, distance=distance)

    # Return the maps
    return maker.maps

# -----------------------------------------------------------------

class MultiBandTIRMapMaker(Configurable):

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
        super(MultiBandTIRMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Maps/dust/tirfuv path
        self.maps_tirfuv_path = None

        # Frames and error maps
        self.frames = None
        self.errors = None

        # The TIR maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The galametz calibration object
        self.galametz = GalametzTIRCalibration()

        # The distance
        self.distance = None

        # Lengths
        self.lengths = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Make the TIR maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(MultiBandTIRMapMaker, self).setup(**kwargs)

        # Get input
        self.frames = kwargs.pop("frames")
        self.errors = kwargs.pop("errors", None)

        # Get distance
        self.distance = kwargs.pop("distance")

        # The number of filters to consider as combinations
        self.lengths = kwargs.pop("lengths", [2,3])

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
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the TIR maps ...")

        # Loop over each combination of 2 or 3 filters
        for filters in sequences.combinations(self.filters, self.lengths):

            # Check if the combination if possible
            if not self.galametz.has_combination_multi_brightness(*filters): continue

            # Get the parameters
            coefficients = self.galametz.get_parameters_multi_brightness(*filters)

            # Get the frames
            frames = []
            for fltr in filters:

                # Convert the frame to neutral intrinsic surface brightness and add it to the list
                frame = self.frames[fltr].converted_to("W/kpc2", density=True, brightness=True, density_strict=True,
                                                       brightness_strict=True, distance=self.distance)
                frames.append(frame)

            # Calculate the TIR
            tir = linear_combination(frames, coefficients)
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frames[0].wcs

            # Determine keys
            combination = tuple([str(fltr) for fltr in filters])

            # Set the TIR map
            self.maps[combination] = tir

            # Set the origins
            self.maps[combination] = filters

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
