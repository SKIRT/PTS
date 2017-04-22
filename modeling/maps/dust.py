#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.dust Contains the DustMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import MapsComponent
from ...core.tools.logging import log
from ...magic.maps.dust.blackbody import BlackBodyDustMapMaker
from ...magic.maps.dust.emission import EmissionDustMapMaker
from ...magic.maps.dust.attenuation import AttenuationDustMapMaker
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class DustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(DustMapMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The TIR to FUV ratio (in log)
        #self.log_tir_to_fuv = None

        # The dust maps (and error maps)
        self.maps = dict()
        self.error_maps = dict()

        # The image of significance masks
        #self.significance = Image()

        # The cutoff mask
        #self.cutoff_mask = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Calculate the significance masks
        #self.calculate_significance()

        # 3.. Make a dust map based on black body pixel fitting
        if self.config.make_black_body: self.make_black_body()

        # 4. Make a dust map simply based on FIR / submm emission in a certain band
        if self.config.make_emission: self.make_emission()

        # 6. Make a dust map based on UV attenuation
        if self.config.make_attenuation: self.make_attenuation()

        # Make the cutoff mask
        #self.make_cutoff_mask()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustMapMaker, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance masks
        if self.config.fuv_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("GALEX FUV", self.config.fuv_significance), "GALEX_FUV")
        if self.config.mips24_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("MIPS 24mu", self.config.mips24_significance), "MIPS_24mu")
        if self.config.pacs70_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("Pacs blue", self.config.pacs70_significance), "Pacs_blue")
        if self.config.pacs160_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("Pacs red", self.config.pacs160_significance), "Pacs_red")
        if self.config.h_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("2MASS H", self.config.h_significance), "2MASS_H")

    # -----------------------------------------------------------------

    def make_black_body(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on black-body fitting to the FIR/submm SED ...")

        # Create the black body dust map maker
        maker = BlackBodyDustMapMaker(self.config.black_body)

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps["black-body"] = maker.maps
        self.error_maps["black-body"] = maker.error_maps

    # -----------------------------------------------------------------

    def make_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the emission in a certain FIR/submm band ...")

        # Create the emission dust map maker
        maker = EmissionDustMapMaker()

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps["emission"] = maker.maps

    # -----------------------------------------------------------------

    def make_attenuation(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on the UV attenuation ...")

        # Create the Attenuation dust map maker
        maker = AttenuationDustMapMaker(self.config.cortese)

        # Run the maker
        maker.run()

        # Add the dust map to the dictionary
        self.maps["attenuation"] = maker.maps

    # -----------------------------------------------------------------

    def make_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the cutoff mask ...")

        # Combine the significance masks
        high_significance = self.significance.intersect_masks()

        # Fill holes
        if self.config.remove_holes: high_significance.fill_holes()

        # Set
        self.cutoff_mask = high_significance.inverse()

    # -----------------------------------------------------------------

    def cutoff_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cutting-off the map at low significance of the data ...")

        # Set zero outside of significant pixels
        self.map[self.cutoff_mask] = 0.0

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

        # Write the error maps
        self.write_error_maps()

        # Write the significance mask
        #self.write_significance_masks()

        # Write the cutoff mask
        #self.write_cutoff_mask()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps ...")

        # Loop over the methods
        for method in self.maps:

            # Create a directory
            path = fs.create_directory_in(self.maps_dust_path, method)

            # Loop over the maps
            for name in self.maps[method]:

                # Determine path
                map_path = fs.join(path, name + ".fits")

                # Save the map
                self.maps[method][name].saveto(map_path)

    # -----------------------------------------------------------------

    def write_error_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the error maps (with different methods) ...")

        # Loop over the methods
        for method in self.maps:

            # Create a directory
            path = fs.create_directory_in(self.maps_dust_path, method)

            # Loop over the maps
            for name in self.error_maps[method]:

                # Determine path
                map_path = fs.join(path, name + "_error.fits")

                # Save the map
                self.maps[method][name].saveto(map_path)

    # -----------------------------------------------------------------

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.saveto(self.dust_significance_path)

    # -----------------------------------------------------------------

    def write_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cutoff mask ...")

        # Write
        self.cutoff_mask.saveto(self.dust_cutoff_path)

# -----------------------------------------------------------------
