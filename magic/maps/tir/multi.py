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
from ....core.basics.log import log
from ....magic.core.frame import linear_combination
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....core.tools import sequences
from ....core.basics.configurable import Configurable
from ....magic.core.list import FrameList
from ....core.tools.stringify import tostr
from ....magic.core.image import Image

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
    lengths = [len(args)]

    # Run the maker
    maker.run(frames=frames, errors=errors, lengths=lengths)

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

    # Run the maker
    maker.run(frames=frames, errors=errors)

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

        # The method name
        self.method_name = None

        # Maps/dust/tirfuv path
        self.maps_tirfuv_path = None

        # Frames and error maps
        self.frames = None
        self.errors = None

        # The TIR maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # The galametz calibration object
        self.galametz = GalametzTIRCalibration()

        # Lengths
        self.lengths = None

        # Region of interest
        self.region_of_interest = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

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

        # Get maps that have already been created
        if "maps" in kwargs: self.maps = kwargs.pop("maps")

        # The number of filters to consider as combinations
        self.lengths = kwargs.pop("lengths", [2,3])

        # Get the method name
        self.method_name = kwargs.pop("method_name", None)

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
        This ufnction ...
        :return:
        """

        return self.method_name is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the TIR maps ...")

        # Loop over each combination of 2 or 3 filters
        for filters in sequences.combinations(self.filters, self.lengths): # returns lists, not tuples, now

            #print(filters, type(filters))

            # Check if the combination if possible
            if not self.galametz.has_combination_multi_brightness(*filters): continue

            # Convert filter combination to string
            key = tostr(filters, delimiter="__", value_delimiter="_")

            # Set the origins (also when the map is already present)
            self.origins[key] = filters  # cannot be tuple (for serialization reasons of the origins dict): OK: sequences.combinations now returns list instead of tuple

            # Set the methods
            if self.has_method_name: self.methods[key] = [self.method_name]

            # Check whether a TIR map is already present
            if key in self.maps:
                log.success("The " + key + " TIR map is already created: not creating it again")
                continue

            # Debugging
            log.debug("Making a TIR map with the filters " + tostr(filters) + " ...")

            # Get the parameters
            coefficients = self.galametz.get_parameters_multi_brightness(*filters)

            # Debugging
            log.debug("The calibration coefficients are " + tostr(coefficients) + " ...")

            # Get the frames
            #frames = []
            #for fltr in filters:
                # Convert the frame to neutral intrinsic surface brightness and add it to the list
                #frame = self.frames[fltr].converted_to("W/kpc2", density=True, brightness=True, density_strict=True,
                #                                       brightness_strict=True, distance=self.distance)
                #frames.append(frame)

            #print(self.frames)
            #print("filters", filters, type(filters))

            # Create frame list
            #frames = FrameList(self.frames[fltr_a], self.frames[fltr_b])
            frames = self.frames[filters] # get sub list

            #print(frames)

            # Convolve
            frames.convolve_to_highest_fwhm()

            # Rebin the frames to the same pixelgridt
            frames.rebin_to_highest_pixelscale()

            # Convert the frames to the same unit
            frames.convert_to_same_unit(unit="W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)

            # Get the frames
            #frame_a, frame_b = frames[0], frames[1]

            # Calculate the TIR
            tir = linear_combination(frames, coefficients)
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True) # TIR can only be bolometric, right??
            #tir.unit = u("W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frames[0].wcs

            #print("test", tir.distance, tir.pixelscale, tir.psf_filter, tir.fwhm)

            # Set other properties
            tir.distance = frames.distance
            tir.pixelscale = frames.pixelscale
            tir.psf_filter = frames.psf_filter
            tir.fwhm = frames.fwhm

            # Determine keys
            #combination = tuple([fltr for fltr in filters])

            # Interpolate NaNs (although normally there aren't any)
            #tir.interpolate_nans(max_iterations=None)
            nans = tir.interpolate_nans_if_below(min_max_in=self.region_of_interest)
            tir.replace_negatives(0.0)
            image = Image()
            image.add_frame(tir, "tir")
            if nans is not None: image.add_mask(nans, "nans")

            # Set the TIR map
            self.maps[key] = image

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
