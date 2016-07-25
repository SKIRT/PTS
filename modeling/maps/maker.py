#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.mapmaking Contains the MapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import numpy as np

# Import astronomical modules
from astropy.units import Unit
from photutils import detect_sources

# Import the relevant PTS classes and modules
from ...magic.basics.mask import Mask
from ...magic.core.image import Image
from ...magic.core.frame import Frame
from .component import MapsComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from .dust.dust import DustMapMaker
from .stars.old import OldStellarMapMaker
from .stars.young import YoungStellarMapMaker
from .stars.ionizing import IonizingStellarMapMaker

# -----------------------------------------------------------------

class MapMaker(MapsComponent):

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
        super(MapMaker, self).__init__(config)

        # -- Attributes --

        # The distance to the galaxy
        self.distance_mpc = None

        # Input images
        self.images = dict()

        # Cutoff masks
        self.cutoff_masks = dict()

        # Bulge and disk
        self.disk = None
        self.bulge = None

        # Set the low signal-to-noise mask to None initially
        self.mask = None

        # The map of dust
        self.dust = None

        # The map of old stars
        self.old_stars = None

        # The map of young stars
        self.young_stars = None

        # The map of ionizing stars
        self.ionizing_stars = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new MapMaker instance
        maker = cls(arguments.config)

        # Set the input and output path
        maker.config.path = arguments.path

        # A single map name can be specified so the procedure is only run for that map
        maker.config.single_map = arguments.map

        # Return the new instance
        return maker

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the input images
        self.load_images()

        # Calculate the signal-to-noise
        self.create_cutoff_masks()
        self.write_cutoff_masks()

        # 3. Cut-off the low signal-to-noise pixels
        #self.cutoff_low_snr()

        # 4. If requested, save the maps with masked
        #if self.config.save_cutoff_maps: self.save_cutoff_maps()

        # 5. Convert maps to solar luminosity units
        self.convert_units()

        # Make the dust map
        self.make_dust_map()

        # Make the old stars map
        self.make_old_stars_map()

        # Make the young (non-ionizing) stars map
        self.make_young_stars_map()

        # Make the ionizing young stars map
        self.make_ionizing_stars_map()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(MapMaker, self).setup()

        # Get the galaxy distance
        self.distance_mpc = self.galaxy_parameters.distance.to("Mpc").value

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Load the GALEX FUV image
        self.load_image("GALEX FUV", "FUV")

        # Load the SDSS i image
        self.load_image("SDSS i", "i")

        # Load the 2MASS H image
        self.load_image("2MASS H", "H")

        # Load the H-alpha image
        self.load_image("Mosaic Halpha", "Halpha")

        # Load the 24 micron image
        self.load_image("MIPS 24mu", "24mu")

        # Load the IRAC 3.6 image
        self.load_image("IRAC I1", "3.6mu")

        # Load the PACS 70 image
        self.load_image("Pacs blue", "70mu")

        # load the PACS 160 image
        self.load_image("Pacs red", "160mu")

        # Load the disk image
        self.load_disk()

        # Load the bulge image
        self.load_bulge()

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the disk image ...")

        # Determine the path to the disk image
        path = fs.join(self.truncation_path, "disk.fits")

        # Load the disk image
        self.disk = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_bulge(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the bulge image ...")

        # Determine the path to the bulge image
        path = fs.join(self.truncation_path, "bulge.fits")

        # Load the bulge image
        self.bulge = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_image(self, image_name, image_id):

        """
        This function ...
        :param image_name:
        :param image_id:
        """

        # Inform the user
        log.info("Loading the " + image_name + " image ...")

        # Determine the full path to the image
        path = fs.join(self.truncation_path, image_name + ".fits")

        # Check whether the image is present
        if not fs.is_file(path): raise IOError("Could not find the " + image_name + " image")

        # Open the image
        image = Image.from_file(path)

        # Assert that the units are MJy/sr
        if not "Halpha" in image_name: assert image.unit == Unit("MJy/sr")

        # Add the image to the dictionary
        self.images[image_id] = image

    # -----------------------------------------------------------------

    def create_cutoff_masks(self):

        """
        This function ...
        :return:
        """

        sigma_level = 3.0

        # Loop over all images
        for name in self.images:

            # Calculate the signal-to-noise ratio in each pixel
            snr = self.images[name].frames.primary / self.images[name].frames.errors

            # Calculate the snr > sigma level mask and add it to the dictionary
            self.cutoff_masks[name] = Mask(snr < sigma_level)

    # -----------------------------------------------------------------

    def write_cutoff_masks(self):

        """
        This function ...
        :return:
        """

        # Loop over all cutoff masks
        for name in self.cutoff_masks:

            # Get the mask
            mask = self.cutoff_masks[name]

            # Determine the path to the FITS file
            path = fs.join(self.maps_cutoff_path, name + ".fits")

            # Save the mask as a FITS file
            Frame(mask.astype(float)).save(path)

    # -----------------------------------------------------------------

    def cutoff_low_snr(self):

        """
        This function ...
        :return:
        """

        # Check whether the reference image is present
        if os.path.isfile(self.config.cutoff.reference_path):

            # Open the reference image
            reference = Image.from_file(self.config.cutoff.reference_path)

            # Check whether the errors frame is present
            assert self.config.errors in reference.frames, "An error map could not be found for the reference image"

            # Create a mask for the pixels with a signal-to-noise ratio of 5 or less
            data = reference.frames[self.config.primary] < self.config.cutoff.level*reference.frames[self.config.errors]
            self.mask = Mask(data)

            # If requested, remove holes from the cut-off mask
            if self.config.cutoff.remove_holes:

                # Save the mask as a FITS file
                Frame(self.mask.astype(float)).save(self.config.saving.cutoff_mask_with_holes_path)

                # Perform the segmentation
                segments = detect_sources(self.mask.astype(float), 0.5, 1).data

                # Save segments
                Frame(segments.astype(float)).save(self.config.saving.cutoff_mask_segments_path)

                # Find the label of the largest segment (=the background)
                label_counts = np.bincount(segments.flatten())
                background_label = np.argmax(label_counts)

                # Create a mask for the holes identified as background
                holes = copy.deepcopy(self.mask)
                holes[segments == background_label] = False

                # Save holes mask
                Frame(holes.astype(float)).save(self.config.saving.cutoff_mask_holes_path)

                # Remove holes from the mask
                self.mask[holes] = False

            # Save the mask as a FITS file
            Frame(self.mask.astype(float)).save(self.config.saving.cutoff_mask_path)

        # If not, raise an error
        else: raise IOError("The prepared reference image could not be found")

        # Cut-off the input images at the same contour level
        for name in self.images: self.images[name].apply_mask(self.mask, 0.0)

        # Cut-off the bulge and disk map at the same contour level
        self.disk[self.mask] = 0.0
        self.bulge[self.mask] = 0.0

    # -----------------------------------------------------------------

    def save_cutoff_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving input images set to zero outside the low signal-to-noise contour level")

        # Save each of the frames
        #self.h.save(self.config.saving.h_cutoff_path)
        #self.fuv.save(self.config.saving.fuv_cutoff_path)
        #self.ha.save(self.config.saving.ha_cutoff_path)
        #self.irac.save(self.config.saving.irac_cutoff_path)
        #self.mips.save(self.config.saving.mips_cutoff_path)
        #self.pacsblue.save(self.config.saving.pacsblue_cutoff_path)
        #self.pacsred.save(self.config.saving.pacsred_cutoff_path)
        #self.disk.save(self.config.saving.disk_cutoff_path)
        #self.bulge.save(self.config.saving.bulge_cutoff_path)

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the H-alpha, 24mu, 70mu and 160mu images to solar luminosities ...")

        # FUV is not converted to Lsun

        # Convert the H-alpha image to solar luminosities
        self.convert_halpha_to_solar()

        # Convert the MIPS image to solar luminosities
        #self.convert_mips_to_solar()

        # Convert the PACS images to solar luminosities
        #self.convert_pacs_to_solar()

        # Save the images that are converted to solar units

        #halpa_path = fs.join(self.maps_solar_path, self.images["Halpha"].name + ".fits")
        #self.images["Halpha"].save(halpa_path)

        #mips_path = fs.join(self.maps_solar_path, self.images["24mu"].name + ".fits")
        #self.images["24mu"].save(mips_path)

        #pacs70_path = fs.join(self.maps_solar_path, self.images["70mu"].name + ".fits")
        #self.images["70mu"].save(pacs70_path)

        #pacs160_path = fs.join(self.maps_solar_path, self.images["160mu"].name + ".fits")
        #self.images["160mu"].save(pacs160_path)

    # -----------------------------------------------------------------

    def convert_halpha_to_solar(self):

        """
        This function converts the H alpha image from
        :return:
        """

        # Debugging
        log.debug("Converting the H-alpha image to solar luminosities ...")

        # Convert from erg/s to Lsun
        self.images["Halpha"].convert_to("Lsun")

    # -----------------------------------------------------------------

    def make_dust_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Create the dust map maker
        maker = DustMapMaker()

        # Run the the map maker
        self.dust = maker.run()

    # -----------------------------------------------------------------

    def make_old_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of old stars ...")

        # Create the old stellar map maker
        maker = OldStellarMapMaker()

        # Run the map maker
        self.old_stars = maker.run()

    # -----------------------------------------------------------------

    def make_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of young stars ...")

        # Create the young stellar map maker
        maker = YoungStellarMapMaker()

        # Run the map maker
        self.young_stars = maker.run()

    # -----------------------------------------------------------------

    def make_ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of ionizing young stars ...")

        # Create the ionizing stellar map maker
        maker = IonizingStellarMapMaker()

        # Run the map maker
        self.ionizing_stars = maker.run()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the output maps
        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the output maps ...")

        # Determine the path to the dust map
        dust_path = fs.join(self.maps_path, "dust.fits")

        # Save the dust map
        self.dust.save(dust_path)

        # Determine the path to the old stars map
        old_stars_path = fs.join(self.maps_path, "old_stars.fits")

        # Save the old stars map
        self.old_stars.save(old_stars_path)

        # Determine the path to the young stars map
        young_stars_path = fs.join(self.maps_path, "young_stars.fits")

        # Save the young stars map
        self.young_stars.save(young_stars_path)

        # Determine the path to the ionizing stars map
        ionizing_stars_path = fs.join(self.maps_path, "ionizing_stars.fits")

        # Save the ionizing stars map
        self.ionizing_stars.save(ionizing_stars_path)

# -----------------------------------------------------------------
