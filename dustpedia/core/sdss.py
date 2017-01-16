#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.sdss Contains the SDSSMosaicMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools.logging import log
from .dataprocessing import DustPediaDataProcessing
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

sdss_bands = ["u", "g", "r", "i", "z"]

# -----------------------------------------------------------------

class SDSSMosaicMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(SDSSMosaicMaker, self).__init__(config)

        # The DustPedia data processing instance
        self.dpdp = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Create directories
        self.create_directories()

        # 2. Do the mosaicing
        self.mosaic()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SDSSMosaicMaker, self).setup()

        # Create the DustPedia data processing instance
        self.dpdp = DustPediaDataProcessing()

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the temporary directory
        temp_path = fs.join(fs.home(), time.unique_name("SDSS_" + galaxy_name + "_" + band))

        # Create the temporary directory
        fs.create_directory(temp_path)

        # temp_path = fs.join(fs.home(), "SDSS_NGC3031_u_2016-08-10--16-48-33-775")

        ###

        # RAW PATH
        raw_path = fs.join(temp_path,
                           "raw")  # frames with HDU0, HDU1, HDU2 .. (primary, calib, sky, ...) IN NANOMAGGIES PER PIXEL
        fs.create_directory(raw_path)

        # FIELDS PATH
        fields_path = fs.join(temp_path, "fields")  # photoField files
        fs.create_directory(fields_path)

        # COUNTS PATH
        counts_path = fs.join(temp_path, "counts")  # frames IN COUNTS (DN)
        fs.create_directory(counts_path)

        # POISSON PATH
        poisson_path = fs.join(temp_path, "poisson_nmgy")  # error maps in NANOMAGGIES PER PIXEL
        fs.create_directory(poisson_path)

        # REBINNED PATH
        rebinned_path = fs.join(temp_path,
                                "rebinned")  # images with frame and corresponding error map in NANOMAGGIES, REBINNED
        fs.create_directory(rebinned_path)

        # FOOTPRINTS PATH
        footprints_path = fs.join(temp_path, "footprints")
        fs.create_directory(footprints_path)

        # MOSAIC PATH
        mosaics_path = fs.join(temp_path, "mosaics")
        fs.create_directory(mosaics_path)

        # RESULT PATH
        results_path = fs.join(temp_path, "results")
        fs.create_directory(results_path)

    def mosaic(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making SDSS mosaic(s) ...")

        # If the band is not specified, do all bands
        if self.config.band is None:

            # Loop over all bands and make the mosaics and Poisson frames
            for band in sdss_bands: self.mosaic_band(band)

        # Do just the specified band otherwise
        else: self.mosaic_band(self.config.band)

    # -----------------------------------------------------------------

    def mosaic_band(self, band):

        """
        This function ...
        :param band:
        :return:
        """

        # Inform the user
        log.info("Making the SDSS " + band + " mosaic ...")

        # Make the mosaic for the specified band
        # just place all the images (for the different bands) in the same output directory, but specifiy the band in the filename
        self.dpdp.make_sdss_mosaic_and_poisson_frame(self.config.galaxy_name, band, self.output_path)

# -----------------------------------------------------------------
