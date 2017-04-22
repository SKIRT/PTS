#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust.tir.multi Contains the MultiBandTIRMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ....core.tools import filesystem as fs
from ....core.tools.logging import log
from ..component import MapsComponent
from ....magic.core.frame import linear_combination
from ....core.units.parsing import parse_unit as u
from ....magic.calibrations.galametz import GalametzTIRCalibration
from ....core.tools import sequences

# -----------------------------------------------------------------

possible_filters = ["MIPS 24mu", "Pacs blue", "Pacs green", "Pacs red", "SPIRE 250"]

# -----------------------------------------------------------------

class MultiBandTIRMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self):

        """
        The constructor ...
        :return:
        """

        # Call the constructor of the base class
        super(MultiBandTIRMapMaker, self).__init__()

        # -- Attributes --

        # Maps/dust/tirfuv path
        self.maps_tirfuv_path = None

        # Frames and error maps
        self.frames = dict()
        self.errors = dict()

        # The TIR maps
        self.maps = dict()

        # The galametz calibration object
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

        # 4. Make the TIR maps
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
        super(MultiBandTIRMapMaker, self).setup()

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

            # Load the error map, if present
            if not self.dataset.has_errormap_for_for_filter(fltr): continue
            errors = self.dataset.get_errormap_for_filter(fltr)
            self.errors[fltr] = errors

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the TIR maps ...")

        # Get the galaxy distance
        distance = self.galaxy_properties.distance

        # Loop over each combination of 2 or 3 filters
        for filters in sequences.combinations(self.available_filters, lengths=[2,3]):

            # Check if the combination if possible
            if not self.galametz.has_combination_multi_brightness(*filters): continue

            # Get the parameters
            coefficients = self.galametz.get_parameters_multi_brightness(*filters)

            # Get the frames
            frames = []
            for fltr in filters:

                # Convert the frame to neutral intrinsic surface brightness and add it to the list
                frame = self.frames[fltr].converted_to("W/kpc2", density=True, brightness=True, density_strict=True, brightness_strict=True)
                frames.append(frame)

            # Calculate the TIR
            tir = linear_combination(frames, coefficients)
            tir.unit = u("W/kpc2", density=False, brightness=True, density_strict=True, brightness_strict=True)
            tir.wcs = frames[0].wcs

            # Determine keys
            combination = tuple([str(fltr) for fltr in filters])

            # Set the TIR map
            self.maps[combination] = tir

    # -----------------------------------------------------------------

    #def make_tir_old(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Inform the user
        #log.info("Creating the TIR map in W/m2 units ...")

        ## GET THE GALAMETZ PARAMETERS
        #a, b, c = self.galametz.get_parameters_multi("MIPS 24mu", "Pacs blue", "Pacs red")

        #assert a == 2.133
        #assert b == 0.681
        #assert c == 1.125

        # MIPS, PACS BLUE AND PACS RED CONVERTED TO LSUN (ABOVE)
        # Galametz (2013) formula for Lsun units
        #tir_map = a * self.frames["MIPS 24mu"] + b * self.frames["Pacs blue"] + c * self.frames["Pacs red"]

        ## Convert the TIR map from Lsun to W / m2

        #conversion_factor = 1.0

        # Conversion from Lsun to W

        #conversion_factor *= solar_luminosity.to("W").value

        # Conversion from W [LUMINOSITY] to W / m2 [FLUX]
        #distance = self.galaxy_properties.distance
        #conversion_factor /= (4. * np.pi * distance**2).to("m2").value

        ## CONVERT AND SET NEW UNIT

        #self.tir_si = Frame(tir_map * conversion_factor)
        #self.tir_si.unit = "W/m2"

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the TIR map
        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the TIR maps ...")

        # Loop over the combinations
        for combination in self.maps:

            # Determine name
            name = "-".join(combination)

            # Determine path
            path = fs.join(self.maps_tir_path, name + ".fits")

            # Save the map
            self.maps[combination].saveto(path)

# -----------------------------------------------------------------
