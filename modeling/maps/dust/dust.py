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
from ..component import MapsComponent
from ....core.tools.logging import log
from .blackbody import BlackBodyDustMapMaker
from .emission import EmissionDustMapMaker
from .buat import BuatDustMapMaker
from .cortese import CorteseDustMapMaker
from ....core.tools import filesystem as fs
from .tirtofuv import TIRtoFUVMapMaker
from ....magic.core.image import Image

# -----------------------------------------------------------------

class DustMapMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(DustMapMaker, self).__init__(config)

        # -- Attributes --

        # The TIR to FUV ratio (in log)
        self.log_tir_to_fuv = None

        # The dust maps
        self.maps = dict()

        # The best dust map
        self.map = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Calculate the significance masks
        self.calculate_significance()

        # 2. Make the TIR to FUV map
        if self.config.make_buat or self.config.make_cortese: self.make_tir_to_fuv()

        # 3.. Make a dust map based on black body pixel fitting
        if self.config.make_black_body: self.make_black_body()

        # 4. Make a dust map simply based on FIR / submm emission in a certain band
        if self.config.make_emission: self.make_emission()

        # 5. Make a dust map based on Buat
        if self.config.make_buat: self.make_buat()

        # 6. Make a dust map based on Cortese
        if self.config.make_cortese: self.make_cortese()

        # Make the final map
        self.make_map()

        # Make the cutoff mask
        self.make_cutoff_mask()

        # Cut-off the map
        self.cutoff_map()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustMapMaker, self).setup()

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

    def make_tir_to_fuv(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the TIR to FUV map ...")

        # Create the maker
        maker = TIRtoFUVMapMaker()

        # Run the maker
        maker.run()

        # Set the TIR to FUV map
        self.log_tir_to_fuv = maker.log_tir_to_fuv

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
        self.maps["black-body"] = maker.map

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
        self.maps["emission"] = maker.map

    # -----------------------------------------------------------------

    def make_buat(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on Buat et al. ...")

        # Create the Buat dust map maker
        maker = BuatDustMapMaker(self.config.buat)

        # Run the maker
        maker.run(self.log_tir_to_fuv)

        # Add the dust map to the dictionary
        self.maps["buat"] = maker.map

    # -----------------------------------------------------------------

    def make_cortese(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a dust map based on Cortese et al. ...")

        # Create the Cortese dust map maker
        maker = CorteseDustMapMaker(self.config.cortese)

        # Run the maker
        maker.run(self.log_tir_to_fuv)

        # Add the dust map to the dictionary
        self.maps["cortese"] = maker.map

    # -----------------------------------------------------------------

    def make_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Selecting the best dust map ...")

        # Set the best dust map
        self.map = self.maps[self.config.best_method].copy()

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

        # Write the final dust map
        self.write_map()

        # Write the significance mask
        self.write_significance_masks()

        # Write the cutoff mask
        self.write_cutoff_mask()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps (with different methods) ...")

        # Loop over the maps
        for label in self.maps:

            # Determine the path
            path = fs.join(self.maps_dust_path, label + ".fits")

            # Save the dust map
            self.maps[label].save(path)

    # -----------------------------------------------------------------

    def write_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the final dust map ...")

        # Write the final dust map
        self.map.save(self.dust_map_path)

    # -----------------------------------------------------------------

    def write_significance_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the significance masks ...")

        # Write
        self.significance.save(self.dust_significance_path)

    # -----------------------------------------------------------------

    def write_cutoff_mask(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cutoff mask ...")

        # Write
        self.cutoff_mask.save(self.dust_cutoff_path)

# -----------------------------------------------------------------
