#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.ionizingstars.ionizing Contains the IonizingStellarMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from copy import copy

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ....core.filter.filter import parse_filter
from ....core.tools import sequences
from ...core.list import NamedFrameList
from ...core.image import Image
from ...core.mask import union

# -----------------------------------------------------------------

def make_map(halpha, hot=None):

    """
    This function ...
    :param halpha:
    :param hot:
    :return: 
    """

    # Create the maker
    maker = IonizingStellarMapsMaker()

    # Set maps of hot dust
    if hot is not None: hots = {"mips24": hot}
    else: hots = None

    # Run the map maker
    maker.run(halpha=halpha, hots=hots)

    # Return
    if hot is not None: return maker.maps["mips24"]
    else: return maker.single_map

# -----------------------------------------------------------------

class IonizingStellarMapsMaker(Configurable):

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
        super(IonizingStellarMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Input
        self.halpha = None
        self.hots = None
        self.hots_origins = None
        self.hots_methods = None
        self.hots_negatives = None
        self.region_of_interest = None

        # Halpha negatives and nans masks
        self.halpha_nans = None
        self.halpha_negatives = None

        # THE MAPS OF IONIZING STARS
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # The H-alpha to hot dust ratio maps
        self.halphatohots = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 2. Process the H-alpha map
        if self.has_halpha: self.process_halpha()

        # 3. Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(IonizingStellarMapsMaker, self).setup()

        # Get the input
        self.halpha = kwargs.pop("halpha")
        self.hots = kwargs.pop("hots", None)
        self.hots_origins = kwargs.pop("hots_origins", None)
        self.hots_methods = kwargs.pop("hots_methods", None)
        self.hots_negatives = kwargs.pop("hots_negatives", None)
        self.region_of_interest = kwargs.pop("region_of_interest", None)

        # Get already existing maps
        self.maps = kwargs.pop("maps", dict())

        # Get already calculated H-alpha to hot dust ratio maps
        self.halphatohots = kwargs.pop("halpha_to_hots", dict())

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return: 
        """

        return self.hots_origins is not None

    # -----------------------------------------------------------------

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        return self.hots_methods is not None

    # -----------------------------------------------------------------

    @property
    def has_negatives(self):

        """
        This function ...
        :return:
        """

        return self.hots_negatives is not None

    # -----------------------------------------------------------------

    def has_negatives_for_name(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.has_negatives and name in self.hots_negatives and self.hots_negatives[name] is not None

    # -----------------------------------------------------------------

    def process_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the H-alpha image ...")

        # Make copy if necessary
        if self.config.interpolate_halpha or self.config.halpha_smoothing_factor: self.halpha = self.halpha.copy()

        # Interpolate?
        if self.config.interpolate_halpha:

            # Interpolate NaNs
            self.halpha_nans = self.halpha.interpolate_nans(max_iterations=None)
            if not self.halpha_nans.has_masked: self.halpha_nans = None

            # Interpolate negative values
            self.halpha_negatives = self.halpha.interpolate_negatives(max_iterations=None)
            if not self.halpha_negatives.has_masked: self.halpha_negatives = None

        # Smooth?
        if self.config.smooth_halpha: self.halpha.smooth(self.config.halpha_smoothing_factor)

    # -----------------------------------------------------------------

    @property
    def has_hots(self):

        """
        This function ...
        :return:
        """

        return self.hots is not None

    # -----------------------------------------------------------------

    @property
    def has_halpha(self):

        """
        This function ...
        :return:
        """

        return self.halpha is not None

    # -----------------------------------------------------------------

    @property
    def halpha_with_hots(self):

        """
        This function ...
        :return:
        """

        return self.has_hots and self.has_halpha

    # -----------------------------------------------------------------

    @property
    def only_halpha(self):

        """
        This function ...
        :return:
        """

        return self.has_halpha and not self.has_hots

    # -----------------------------------------------------------------

    @property
    def only_hots(self):

        """
        This function ...
        :return:
        """

        return self.has_hots and not self.has_halpha

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the map of ionizing stars ...")

        # H-alpha with hot dust contribution
        if self.halpha_with_hots: self.make_maps_halpha_hot()

        # Only halpha
        if self.has_halpha: self.make_map_halpha()

        # Only hot dust contribution, and no Halpha to combine it with
        if self.only_hots: self.make_maps_hot()

    # -----------------------------------------------------------------

    def make_maps_halpha_hot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps of ionizing stars based on H-alpha emission and hot dust maps ...")

        # Set method name
        method_name = "halpha-hot"

        # NEW: MAKE MAP OF IONIZING STARS FOR VARIOUS DIFFERENT maps of hot dust
        for name in self.hots:

            # Set the origin
            if self.has_origins:
                origins = copy(self.hots_origins[name])
                sequences.append_unique(origins, parse_filter("Halpha"))
                self.origins[name] = origins

            # Set the method
            if self.has_methods:
                methods = copy(self.hots_methods[name])
                methods.append(method_name)
                self.methods[name] = methods

            # Check whether a map is already present
            if name in self.maps:
                log.success("The " + name + " ionizing stellar map is already created: not creating it again")
                continue

            # Get the hot dust map
            hot = self.hots[name]

            # Uniformize the hot dust and H-alpha map
            frames = NamedFrameList(hot=hot, halpha=self.halpha)
            frames.convolve_and_rebin()  # convolve and rebin

            # Convert to erg/s
            frames["halpha"].convert_to("erg/s")
            # convert to neutral spectral (frequency/wavelength) density (nu * Lnu or lambda * Llambda)
            frames["hot"].convert_to("erg/s", density=True, density_strict=True)

            # Calculate the ratio of the halpha and the hot dust map
            ratio = frames["halpha"] / frames["hot"]
            ratio.replace_infs_by_nans()
            ratio.replace_by_nans_where_greater_than(1)
            ratio.cutoff_below_zero()

            # Add the ratio
            self.halphatohots[name] = ratio

            # Calculate ionizing stars map and ratio
            # CALZETTI et al., 2007
            ionizing = frames["halpha"] + 0.031 * frames["hot"]

            # Normalize the ionizing stellar map
            ionizing.normalize()

            # Create image
            image = Image()
            image.add_frame(ionizing, "ionizing")

            # Collect negative pixels masks
            negatives = []
            if self.halpha_negatives is not None: negatives.append(self.halpha_negatives.rebinned(ionizing.wcs))
            if self.has_negatives_for_name(name): negatives.append(self.hots_negatives[name].rebinned(ionizing.wcs))

            # Create and add negatives mask
            if len(negatives) > 0:
                negatives = union(*negatives)
                image.add_mask(negatives, "negatives")

            # Add hot dust negatives map
            if self.has_negatives_for_name(name): image.add_mask(self.hots_negatives[name].rebinned(ionizing.wcs), "hot_negatives")

            # Add NaNs mask
            if self.halpha_nans is not None:
                nans = self.halpha_nans.rebinned(ionizing.wcs)
                image.add_mask(nans, "nans")

            # Add the map of ionizing stars
            self.maps[name] = image

    # -----------------------------------------------------------------

    def make_map_halpha(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making map of ionizing stars based on H-alpha emission ...")

        # Set method name
        method_name = "halpha"

        # Take the processed halpha image
        ionizing = self.halpha.copy()

        # Normalize the map
        ionizing.normalize()

        # Create image
        image = Image()
        image.add_frame(ionizing, "ionizing")
        if self.halpha_nans is not None: image.add_mask(self.halpha_nans, "nans")
        if self.halpha_negatives is not None: image.add_mask(self.halpha_negatives, "negatives")

        # Add the map of ionizing stars
        self.maps[method_name] = image

        # Set origins
        if self.has_origins: self.origins[method_name] = [parse_filter("Halpha")]

        # Set method
        if self.has_methods: self.methods[method_name] = [method_name]

    # -----------------------------------------------------------------

    def make_maps_hot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making map of ionizing stars based on hot dust maps ...")

        # Set method name
        method_name = "hot"

        # MAKE MAP OF IONIZING STARS FOR VARIOUS DIFFERENT maps of hot dust
        for name in self.hots:

            # Set the origin
            if self.has_origins:
                origins = copy(self.hots_origins[name])
                self.origins[name] = origins

            # Set the method
            if self.has_methods:
                methods = copy(self.hots_methods[name])
                methods.append(method_name)
                self.methods[name] = methods

            # Check whether a map is already present
            if name in self.maps:
                log.warning("The " + name + " ionizing stellar map is already created: not creating it again")
                continue

            # Get the hot dust map
            ionizing = self.hots[name].copy()
            #ionizing[ionizing < 0.0] = 0.0

            # Normalize the dust map
            ionizing.normalize()

            # Create image
            image = Image()
            image.add_frame(ionizing, "ionizing")
            if self.has_negatives_for_name(name): image.add_mask(self.hots_negatives[name], "negatives")

            # Add the map of ionizing stars
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
