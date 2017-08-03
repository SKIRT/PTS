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

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.units.parsing import parse_unit as u
from ....core.basics.configurable import Configurable
from ....core.filter.filter import parse_filter
from ....core.tools import sequences
from ...core.list import NamedFrameList

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

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

        # THE MAPS OF IONIZING STARS
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 3. Make the map
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

        # Get already existing maps
        self.maps = kwargs.pop("maps", dict())

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

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # H-ALPHA HAS BEEN CONVERTED TO LSUN (ABOVE)

        # Inform the user
        log.info("Making the map of ionizing stars ...")

        # With hot dust contribution
        if self.hots is not None: self.make_maps_hot()

        # Only halpha
        self.make_map_halpha()

    # -----------------------------------------------------------------

    def make_maps_hot(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making maps of ionizing stars based on H-alpha emission and hot dust maps ...")

        method_name = "hot"

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
                log.warning("The " + name + " ionizing stellar map is already created: not creating it again")
                continue

            # Young ionizing stars = Ha + 0.031 x MIPS24_corrected
            #best_corrected_24mu_map = self.corrected_24mu_maps[factor]

            # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
            #best_corrected_24mu_map[best_corrected_24mu_map < 0.0] = 0.0

            hot = self.hots[name]
            hot[hot < 0.0] = 0.0

            # Uniformize the hot dust and H-alpha map
            frames = NamedFrameList(hot=hot, halpha=self.halpha)
            frames.convolve_and_rebin()  # convolve and rebin

            # SEE CALZETTI et al., 2007

            # Convert to erg/s
            frames["halpha"].convert_to("erg/s")
            # convert to neutral spectral (frequency/wavelength) density (nu * Lnu or lambda * Llambda)
            frames["hot"].convert_to("erg/s", density=True, density_strict=True)

            # Calculate ionizing stars map and ratio
            ionizing = frames["halpha"] + 0.031 * frames["hot"]

            # Normalize the dust map
            ionizing.normalize()

            # Add the map of ionizing stars
            self.maps[name] = ionizing

        # ionizing_ratio = self.ha / (0.031*mips_young_stars)

        # MASK NEGATIVE AND LOW SIGNAL-TO-NOISE PIXELS

        # Set pixels to zero with low signal-to-noise in the H Alpha image
        # ionizing[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0
        # ionizing_ratio[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0

        # Set pixels to zero with low signal-to-noise in the 24 micron image
        # ionizing[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0
        # ionizing_ratio[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        # ionizing[ionizing < 0.0] = 0.0
        # ionizing_ratio[ionizing < 0.0] = 0.0

        # New
        #ionizing[self.cutoff_masks["Halpha"]] = 0.0

        # Set the ionizing stars map
        #self.map = ionizing

    # -----------------------------------------------------------------

    def make_map_halpha(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making map of ionizing stars based on H-alpha emission ...")

        method_name = "halpha"

        ionizing = self.halpha.copy()
        ionizing[ionizing < 0.0] = 0.0

        # Normalize the dust map
        ionizing.normalize()

        # Add the map of ionizing stars
        self.maps[method_name] = ionizing

        # Set origins
        if self.has_origins: self.origins[method_name] = [parse_filter("Halpha")]

        # Set method
        if self.has_methods: self.methods[method_name] = [method_name]

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
