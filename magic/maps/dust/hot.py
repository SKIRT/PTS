#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.hot Contains the HotDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from copy import copy

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ...core.list import NamedFrameList
from ....core.tools import sequences
from ....core.units.parsing import parse_unit as u
from ...core.image import Image

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

def make_map(mips24, old, factor):

    """
    This function ...
    :return: 
    """

    # Create the map maker
    maker = HotDustMapsMaker()

    # Set input
    factors = [factor]

    # Run the map maker
    maker.run(mips24=mips24, old=old, factors=factors)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

def make_maps(mips24, old, factors):

    """
    This fucntion ...
    :return: 
    """

    # Create the map maker
    maker = HotDustMapsMaker()

    # Run the map maker
    maker.run(mips24=mips24, old=old, factors=factors)

    # Return the maps
    return maker.maps

# -----------------------------------------------------------------

class HotDustMapsMaker(Configurable):

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
        super(HotDustMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The mips 24 image
        self.mips24 = None

        # The maps of the old stellar disk
        self.old = None

        # The origins
        self.old_origins = None

        # The methods
        self.old_methods = None

        # Factors
        self.factors = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # Method name
        self.method_name = None

        # Region of interest
        self.region_of_interest = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # 2. Make the maps
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ....
        :param kwargs: 
        :return: 
        """

        # Call the setup function of the base class
        super(HotDustMapsMaker, self).setup(**kwargs)

        # Get the input
        self.mips24 = kwargs.pop("mips24")

        # Maps of old stars and their origins
        self.old = kwargs.pop("old")
        self.old_origins = kwargs.pop("old_origins", None)
        self.old_methods = kwargs.pop("old_methods", None)

        # The method name
        self.method_name = kwargs.pop("method_name", None)
        if self.has_methods and self.method_name is None: raise ValueError("Method name has to be specified when methods are given")

        # Set factors
        self.factors = kwargs.pop("factors")

        # Get region of interest
        self.region_of_interest = kwargs.pop("region_of_interest", None)

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return:
        """

        return self.old_origins is not None

    # -----------------------------------------------------------------

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        return self.old_methods is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of hot dust ...")

        # Loop over the different old stellar maps
        for old_name in self.old:

            # Debugging
            log.debug("Creating maps of hot dust based on the '" + old_name + "' old stellar map ...")

            # Normalize the old map
            normalized_old = self.old[old_name].normalized()

            # Uniformize the MIPS 24 micron image and the old disk map
            frames = NamedFrameList(old=normalized_old, mips24=self.mips24)
            frames.convolve_and_rebin(unitless="old")

            # CHECK IF OLD IS STILL NORMALIZED
            if not frames["old"].is_normalized():
                log.warning("Need to re-normalize the old stellar map")
                frames["old"].normalize()

            # Loop over the different factors
            for factor in self.factors:

                # Debugging
                log.debug("Creating maps of host dust with a MIPS 24mu correction factor of " + str(factor) + " ...")

                # Determine name
                name = old_name + "__" + str(factor)

                # Debugging
                log.debug("Name for the map: " + name)

                # Set origin
                if self.has_origins:

                    origins = [self.mips24.filter]
                    old_origins = copy(self.old_origins[old_name])

                    # Add old origins
                    sequences.extend_unique(origins, old_origins)

                    # Add the origins
                    self.origins[name] = origins

                # Set method
                if self.has_methods:

                    methods = copy(self.old_methods[old_name])
                    methods.append(self.method_name)
                    self.methods[name] = methods

                # Check whether a map is already present
                if name in self.maps:
                    log.success("The " + name + " hot dust map is already created: not creating it again")
                    continue

                # Calculate the corrected 24 micron image
                hot_dust = make_corrected_24mu_map(frames["mips24"], frames["old"], factor, normalize_in=self.region_of_interest)

                # Interpolate negatives
                negatives = hot_dust.interpolate_negatives_if_below(min_max_in=self.region_of_interest)
                hot_dust.replace_negatives(0.0)  # if any left
                # DON'T Normalize!!

                # Create image
                image = Image()
                image.add_frame(hot_dust, "hot")
                if negatives is not None: image.add_mask(negatives, "negatives")

                # Add the image to the dictionary
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

def make_corrected_24mu_map(mips24, disk, factor, normalize_in=None):

    """
    This function ...
    :param mips24:
    :param disk:
    :param factor:
    :param normalize_in:
    :return:
    """

    # Inform the user
    log.info("Subtracting the old stellar contribution from the 24 micron emission map with a factor of " + str(factor) + " ...")
    # Lu et al. 2014: 48% voor MIPS 24

    # Calculate sum of MIPS 24 in the specified region
    if normalize_in is not None: normalization_value = mips24.sum_in(normalize_in, add_unit=False)
    else: normalization_value = mips24.sum(add_unit=False)

    # Total contribution in solar units
    total_contribution = factor * normalization_value

    # Subtract the disk contribution to the 24 micron image
    new_mips = mips24 - total_contribution * disk # disk image is normalized

    # Return the new 24 micron frame
    return new_mips

# -----------------------------------------------------------------
