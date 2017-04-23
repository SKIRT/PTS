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
import numpy as np

# Import astronomical modules
from astropy import constants

# Import the relevant PTS classes and modules
from ....core.tools.logging import log
from ....core.tools import filesystem as fs
from ....core.basics.distribution import Distribution
from ....core.plot.distribution import DistributionPlotter
from ....magic.region.composite import PixelCompositeRegion
from ....magic.region.list import PixelRegionList
from ....magic.core.image import Image
from ....core.units.parsing import parse_unit as u
from ....core.basics.configurable import Configurable

# -----------------------------------------------------------------

speed_of_light = constants.c
solar_luminosity = 3.846e26 * u("W")

# -----------------------------------------------------------------

def make_map():

    """
    This function ...
    :return: 
    """

    maker = IonizingStellarMapsMaker()

    maker.run()

# -----------------------------------------------------------------

class IonizingStellarMapsMaker(Configurable):

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
        super(IonizingStellarMapsMaker, self).__init__(config, interactive)

        # -- Attributes --

        # The MIPS 24 micron image and error map IN SOLAR UNITS
        self.mips24 = None
        self.mips24_errors = None

        # The Halpha map IN SOLAR UNITS
        self.halpha = None

        # The NORMALIZED (to unity) DISK IMAGE
        self.disk = None

        # The maps of the corrected 24 micron emission
        self.corrected_24mu_maps = dict()

        # The distributions of the pixel values of the corrected 24mu maps
        self.corrected_24mu_distributions = dict()

        # The map of ionizing stars
        #self.map = None

        # NEW: THE MAPS OF IONIZING STARS
        self.maps = dict()

        # The path to the maps/ionizing/24mu directory
        self.maps_ionizing_24mu_path = None

        # The path to the maps/ionizing/solar directory
        self.maps_ionizing_solar_path = None

        # The region of the pixels used for plotting the distributions of pixel values
        self.distribution_region = None

        # The image of significance masks
        self.significance = Image()

        # The cutoff mask
        self.cutoff_mask = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the necessary frames
        self.load_frames()

        # Calculate the significance masks
        self.calculate_significance()

        # 3. Make the map
        self.make_maps()

        # ...
        #self.create_distribution_region()
        #self.make_distributions()

        # 4. Normalize the map
        #self.normalize_map()

        # Make the cutoff mask
        #self.make_cutoff_mask()

        # 5. Cut-off map
        #self.cutoff_map()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(IonizingStellarMapsMaker, self).setup()

        # Set paths
        #self.maps_ionizing_24mu_path = fs.create_directory_in(self.maps_ionizing_path, "24mu")
        #self.maps_ionizing_solar_path = fs.create_directory_in(self.maps_ionizing_path, "solar")

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary data ...")

        # Load the MIPS 24 micron image and convert to solar units
        self.load_mips()

        # Load the H alpha image and convert to solar units
        self.load_halpha()

        # Load the disk image and normalize to unity
        self.load_disk()

    # -----------------------------------------------------------------

    def calculate_significance(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the significance masks ...")

        # Get the significance mask
        if self.config.mips24_significance > 0: self.significance.add_mask(self.dataset.get_significance_mask("MIPS 24mu", self.config.mips24_significance), "MIPS_24mu")
        if self.config.halpha_significance > 0: self.significance.add_mask(self.get_halpha_significance_mask(self.config.halpha_significance), "Halpha")

    # -----------------------------------------------------------------

    def load_mips(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the MIPS 24 micron image and converting to solar units ...")

        # Get MIPS 24 micron frame and error map
        self.mips24 = self.dataset.get_frame("MIPS 24mu")  # in original MJy/sr units
        self.mips24_errors = self.dataset.get_errormap("MIPS 24mu")  # in original MJy/sr units

        ## CONVERT TO LSUN

        # Get the galaxy distance
        distance = self.galaxy_properties.distance

        # Get pixelscale and wavelength
        pixelscale = self.mips24.average_pixelscale
        wavelength = self.mips24.filter.pivot

        # Conversion from MJy / sr to Jy / sr
        conversion_factor = 1e6

        # Conversion from Jy / sr to Jy / pix(2)
        conversion_factor *= (pixelscale ** 2).to("sr/pix2").value

        # Conversion from Jy / pix to W / (m2 * Hz) (per pixel)
        conversion_factor *= 1e-26

        # Conversion from W / (m2 * Hz) (per pixel) to W / (m2 * m) (per pixel)
        conversion_factor *= (speed_of_light / wavelength ** 2).to("Hz/m").value

        # Conversion from W / (m2 * m) (per pixel) [SPECTRAL FLUX] to W / m [SPECTRAL LUMINOSITY]
        conversion_factor *= (4. * np.pi * distance ** 2).to("m2").value

        # Conversion from W / m [SPECTRAL LUMINOSITY] to W [LUMINOSITY]
        conversion_factor *= wavelength.to("m").value

        # Conversion from W to Lsun
        conversion_factor *= 1. / solar_luminosity.to("W").value

        ## DO THE CONVERSION

        self.mips24 *= conversion_factor
        self.mips24_errors *= conversion_factor

    # -----------------------------------------------------------------

    def load_halpha(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the H-alpha image and converting to solar units ...")

        # Get the H-alpha image
        self.halpha = self.masked_halpha_frame

        # Convert from erg/s to Lsun
        self.halpha.convert_to("Lsun")

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Get disk frame
        self.disk = self.masked_disk_frame

        # Normalize the disk image
        self.disk.normalize()
        self.disk.unit = None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # H-ALPHA HAS BEEN CONVERTED TO LSUN (ABOVE)

        # Inform the user
        log.info("Making the map of ionizing stars ...")

        # NEW: MAKE MAP OF IONIZING STARS FOR VARIOUS DIFFERENT FACTORS
        for factor in self.corrected_24mu_maps:

            # Young ionizing stars = Ha + 0.031 x MIPS24_corrected
            best_corrected_24mu_map = self.corrected_24mu_maps[factor]

            # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
            best_corrected_24mu_map[best_corrected_24mu_map < 0.0] = 0.0

            # Calculate ionizing stars map and ratio
            ionizing = self.halpha + 0.031 * best_corrected_24mu_map

            # Add the map of ionizing stars
            self.maps[factor] = ionizing


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

    def create_distribution_region(self):

        """
        This function ...
        :return:
        """

        disk_ellipse = self.disk_ellipse.to_pixel(self.mips24.wcs)
        inner_ellipse = disk_ellipse * self.config.histograms_annulus_range.min
        outer_ellipse = disk_ellipse * self.config.histograms_annulus_range.max
        composite = PixelCompositeRegion(outer_ellipse, inner_ellipse)
        region = PixelRegionList()
        region.append(composite)
        self.distribution_region = region

    # -----------------------------------------------------------------

    def normalize_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Normalizing the map of ionizing stars ...")

        # Normalize the dust map
        self.map.normalize()
        self.map.unit = None

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
