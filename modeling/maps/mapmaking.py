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
from ...core.tools import filesystem, inspection, tables
from ..decomposition.decomposition import load_parameters
from ...core.tools.logging import log

# -----------------------------------------------------------------

# The path to the table containing the parameters from Cortese et. al 2008
cortese_table_path = filesystem.join(inspection.pts_dat_dir("modeling"), "cortese.dat")

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

        # The structural galaxy parameters
        self.parameters = None
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

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

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

        # Determine the path to the parameters file
        path = filesystem.join(self.components_path, "parameters.dat")

        # Load the structural parameters
        self.parameters = load_parameters(path)

        # Get the galaxy distance
        self.distance_mpc = self.parameters.distance.to("Mpc").value

        # Load the Cortese et. al 2008 table
        self.cortese = tables.from_file(cortese_table_path)

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
        path = filesystem.join(self.truncation_path, "disk.fits")

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
        path = filesystem.join(self.truncation_path, "bulge.fits")

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
        path = filesystem.join(self.truncation_path, image_name + ".fits")

        # Check whether the image is present
        if not filesystem.is_file(path): raise IOError("Could not find the " + image_name + " image")

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
            path = filesystem.join(self.maps_cutoff_path, name + ".fits")

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
        self.h.save(self.config.saving.h_cutoff_path)
        self.fuv.save(self.config.saving.fuv_cutoff_path)
        self.ha.save(self.config.saving.ha_cutoff_path)
        self.irac.save(self.config.saving.irac_cutoff_path)
        self.mips.save(self.config.saving.mips_cutoff_path)
        self.pacsblue.save(self.config.saving.pacsblue_cutoff_path)
        self.pacsred.save(self.config.saving.pacsred_cutoff_path)
        self.disk.save(self.config.saving.disk_cutoff_path)
        self.bulge.save(self.config.saving.bulge_cutoff_path)

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
        self.convert_mips_to_solar()

        # Convert the PACS images to solar luminosities
        self.convert_pacs_to_solar()

        # Save the images that are converted to solar units

        halpa_path = filesystem.join(self.maps_solar_path, self.images["Halpha"].name + ".fits")
        self.images["Halpha"].save(halpa_path)

        mips_path = filesystem.join(self.maps_solar_path, self.images["24mu"].name + ".fits")
        self.images["24mu"].save(mips_path)

        pacs70_path = filesystem.join(self.maps_solar_path, self.images["70mu"].name + ".fits")
        self.images["70mu"].save(pacs70_path)

        pacs160_path = filesystem.join(self.maps_solar_path, self.images["160mu"].name + ".fits")
        self.images["160mu"].save(pacs160_path)

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

    def convert_mips_to_solar(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Converting the 24 micron image to solar luminosities ...")

        # Calculate conversion factors from MJy/sr to solar luminosities
        exponent = (2.0*np.log10(2.85/206264.806247)) - 20.0 + np.log10(3e8/24e-6) + np.log10(4.0*np.pi) + (2.0*np.log10(self.distance_mpc*3.08567758e22)) - np.log10(3.846e26)

        # Multiply the image
        self.images["24mu"] *= 10.**exponent

        # Set the new unit
        self.images["24mu"].unit = "Lsun"

    # -----------------------------------------------------------------

    def convert_pacs_to_solar(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Converting the 70 and 160 micron images to solar luminosities ...")

        # Calculate the conversion factors
        exponent_70 = (2.0*np.log10(2.85/206264.806247)) - 20.0 + np.log10(3e8/70e-6) + np.log10(4.0*np.pi) + (2.0*np.log10(self.distance_mpc*3.08567758e22)) - np.log10(3.846e26)
        exponent_160 = (2.0*np.log10(2.85/206264.806247)) - 20.0 + np.log10(3e8/160e-6) + np.log10(4.0*np.pi) + (2.0*np.log10(self.distance_mpc*3.08567758e22)) - np.log10(3.846e26)

        # Convert the units of the 70 micron and 160 micron images from MJy/sr to solar luminosities
        self.images["70mu"] *= 10.0**exponent_70
        self.images["160mu"] *= 10.0**exponent_160

        # Set the new unit
        self.images["70mu"].unit = "Lsun"
        self.images["160mu"].unit = "Lsun"

    # -----------------------------------------------------------------

    def make_dust_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust map ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # Creat the FUV attenuation map according to the calibration in Cortese et. al 2008
        a_fuv_cortese = self.create_afuv_cortese("FUV-i")

        # Set attenuation to zero where the original FUV map is smaller than zero
        a_fuv_cortese[self.images["FUV"].frames.primary <= 0.0] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        a_fuv_cortese[a_fuv_cortese < 0.0] = 0.0

        # Cutoff
        a_fuv_cortese[self.cutoff_masks["160mu"]] = 0.0

        # Set the A(FUV) map as the dust map
        self.dust = a_fuv_cortese

    # -----------------------------------------------------------------

    def make_old_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of old stars ...")

        # Old stars = IRAC3.6 - bulge
        # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

        # The relative contribution of the bulge to the 3.6mu emission
        #bulge_rel_contribution = self.parameters.bulge.f

        # Total flux of the IRAC 3.6mu image
        #total_flux = np.sum(self.images["3.6mu"].frames.primary)

        # Calculate factor
        #factor = bulge_rel_contribution * total_flux / np.sum(self.bulge)

        # Create the old stars map
        #old_stars = self.images["3.6mu"].frames.primary - factor * self.bulge

        # Convert the 3.6 micron image from MJy/sr to Jy/sr
        conversion_factor = 1.0
        conversion_factor *= 1e6

        # Convert the 3.6 micron image from Jy / sr to Jy / pixel
        pixelscale = self.images["3.6mu"].xy_average_pixelscale
        pixel_factor = (1.0/pixelscale**2).to("pix2/sr").value
        conversion_factor /= pixel_factor
        self.images["3.6mu"] *= conversion_factor
        self.images["3.6mu"].unit = "Jy"

        i1_jy_path = filesystem.join(self.maps_intermediate_path, "i1_jy.fits")
        self.images["3.6mu"].save(i1_jy_path)

        # Subtract bulge
        #old_stars = self.images["3.6mu"].frames.primary - (self.bulge * 1.5)
        old_stars = self.images["3.6mu"].frames.primary - self.bulge

        bulge_residual = self.images["3.6mu"].frames.primary - self.disk
        bulge_residual_path = filesystem.join(self.maps_intermediate_path, "bulge_residual.fits")
        bulge_residual.save(bulge_residual_path)

        # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
        #old_stars[self.irac < self.config.old_stars.irac_snr_level*self.irac_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        old_stars[old_stars < 0.0] = 0.0

        # Mask pixels outside of the low signal-to-noise contour
        #old_stars[self.mask] = 0.0

        # Set the old stars map
        self.old_stars = old_stars

    # -----------------------------------------------------------------

    def make_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of young stars ...")

        # Calculate the non ionizing young stars map from the FUV data
        non_ionizing_stars = self.get_fuv_young_stars_map()

        # Set the young stars map
        self.young_stars = non_ionizing_stars

    # -----------------------------------------------------------------

    def make_ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the map of ionizing young stars ...")

        # H-ALPHA HAS BEEN CONVERTED TO LSUN (ABOVE)

        # Young ionizing stars = Ha + 0.031 x MIPS24

        # CALCULATE THE IONIZING STARS MAP BASED ON THE CONVERTED H ALPHA AND THE DISK-SUBTRACTED 24 MICRON IMAGE

        # Calculate the young stellar contribution to the 24 micron image
        mips_young_stars = self.get_mips_young_stars_map()

        # Save the mips_young_stars map
        mips_young_path = filesystem.join(self.maps_intermediate_path, "24mu_young.fits")
        mips_young_stars.save(mips_young_path)

        # Calculate ionizing stars map and ratio
        ionizing = self.images["Halpha"].frames.primary + 0.031 * mips_young_stars

        #ionizing_ratio = self.ha / (0.031*mips_young_stars)

        # MASK NEGATIVE AND LOW SIGNAL-TO-NOISE PIXELS

        # Set pixels to zero with low signal-to-noise in the H Alpha image
        #ionizing[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0
        #ionizing_ratio[self.ha < self.config.ionizing_stars.ha_snr_level*self.ha_errors] = 0.0

        # Set pixels to zero with low signal-to-noise in the 24 micron image
        #ionizing[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0
        #ionizing_ratio[self.mips < self.config.ionizing_stars.mips_snr_level*self.mips_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        #ionizing[ionizing < 0.0] = 0.0
        #ionizing_ratio[ionizing < 0.0] = 0.0

        # New
        ionizing[self.cutoff_masks["Halpha"]] = 0.0

        # Set the ionizing stars map
        self.ionizing_stars = ionizing

    # -----------------------------------------------------------------

    def get_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the TIR map ...")

        # MIPS, PACS BLUE AND PACS RED CONVERTED TO LSUN (ABOVE)

        # Galametz (2013) formula for Lsun units
        tir_data = 2.133 * self.images["24mu"].frames.primary + 0.681 * self.images["70mu"].frames.primary + 1.125 * self.images["160mu"].frames.primary

        # Return the TIR map (in solar units)
        return tir_data

    # -----------------------------------------------------------------

    def get_fuv_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the FUV emission ...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        flux_fuv = np.sum(self.images["FUV"].frames.primary)

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.2 * flux_fuv / np.sum(self.disk)

        # Subtract the disk contribution to the FUV image
        new_fuv = self.images["FUV"].frames.primary - factor * self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_fuv[new_fuv < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        #new_fuv[self.fuv < self.config.non_ionizing_stars.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the new FUV frame
        return new_fuv

    # -----------------------------------------------------------------

    def get_mips_young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Subtracting the old stellar contribution from the 24 micron emission ...")

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        ## MIPS HAS BEEN CONVERTED TO LSUN (ABOVE)

        flux_mips = np.sum(self.images["24mu"].frames.primary)

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.48 * flux_mips / np.sum(self.disk)

        # Subtract the disk contribution to the 24 micron image
        new_mips = self.images["24mu"].frames.primary - factor * self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_mips[new_mips < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        #new_mips[self.mips < self.config.ionizing_stars.mips_young_stars.mips_snr_level*self.mips_errors] = 0.0

        # Return the new 24 micron frame
        return new_mips

    # -----------------------------------------------------------------

    def create_afuv_buat(self):

        """
        This function ...
        :return:
        """

        #a_fuv_buat = (-0.0333*x3) + (0.3522*x2) + (1.1960*x) + 0.4967

        pass

    # -----------------------------------------------------------------

    def create_afuv_cortese(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour: "FUV-H", "FUV-i", "FUV-r", "FUV-g" or "FUV-B"
        :return:
        """

        # Inform the user
        log.info("Creating the A(FUV) map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

        # CALCULATE FUV AND TIR MAP IN W/M2 UNIT

        # Convert the FUV map from MJy/sr to W/m2
        factor = - 20.0 + np.log10(3.e8) - np.log10(0.153e-6) + (2*np.log10(2.85/206264.806247))
        fuv_converted = self.images["FUV"] * 10.0**factor
        fuv_converted.unit = Unit("W/m2")

        # Save
        fuv_converted_path = filesystem.join(self.maps_intermediate_path, "FUV Wpm2.fits")
        fuv_converted.save(fuv_converted_path)

        # Get the TIR map in solar units
        tir_map = self.get_tir_map()

        # Convert the TIR frame from solar units to W/m2
        exponent = np.log10(3.846e26) - np.log10(4*np.pi) - (2.0*np.log10(self.distance_mpc*3.08567758e22))
        tir_map *= 10.0**exponent
        tir_map.unit = Unit("W/m2")

        # Save
        tir_path = filesystem.join(self.maps_intermediate_path, "TIR.fits")
        tir_map.save(tir_path)

        # CALCULATE TIR TO FUV RATIO

        # The ratio of TIR and FUV
        tir_to_fuv = np.log10(tir_map / fuv_converted.frames.primary)

        # Save TIR to FUV ratio map
        tir_to_fuv_path = filesystem.join(self.maps_intermediate_path, "TIRtoFUV.fits")
        tir_to_fuv.save(tir_to_fuv_path)

        # Get the sSFR map
        if ssfr_colour == "FUV-H": ssfr = self.get_fuv_h()
        elif ssfr_colour == "FUV-i": ssfr = self.get_fuv_i()
        elif ssfr_colour == "FUV-r": ssfr = self.get_fuv_r()
        elif ssfr_colour == "FUV-g": ssfr = self.get_fuv_g()
        elif ssfr_colour == "FUV-B": ssfr = self.get_fuv_b()
        else: raise ValueError("Invalid sSFR colour")

        # Calculate powers of tir_to_fuv
        tir_to_fuv2 = np.power(tir_to_fuv, 2.0)
        tir_to_fuv3 = np.power(tir_to_fuv, 3.0)
        tir_to_fuv4 = np.power(tir_to_fuv, 4.0)

        # Create an empty image
        a_fuv_cortese = Frame.zeros_like(tir_to_fuv)

        limits = []
        a1_list = []
        a2_list = []
        a3_list = []
        a4_list = []
        a5_list = []

        # Loop over all entries in the Cortese et. al
        for i in range(len(self.cortese)):

            upper = self.cortese[ssfr_colour][i]
            if i == len(self.cortese) - 1: lower = None
            else: lower = self.cortese[ssfr_colour][i+1]

            limits.append((lower, upper))

            a1 = self.cortese["a1"][i]
            a2 = self.cortese["a2"][i]
            a3 = self.cortese["a3"][i]
            a4 = self.cortese["a4"][i]
            a5 = self.cortese["a5"][i]

            a1_list.append(a1)
            a2_list.append(a2)
            a3_list.append(a3)
            a4_list.append(a4)
            a5_list.append(a5)

        # Create the FUV attenuation map
        for i in range(len(limits)):

            if limits[i][0] is None: where = ssfr < limits[i][1]
            elif limits[i][1] is None: where = ssfr > limits[i][0]
            else: where = (ssfr >= limits[i][0]) * (ssfr < limits[i][1])

            # Set the appropriate pixels
            a_fuv_cortese[where] = a1_list[i] + a2_list[i]*tir_to_fuv[where] + a3_list[i]*tir_to_fuv2[where] + a4_list[i]*tir_to_fuv3[where] + a5_list[i]*tir_to_fuv4[where]

        # Set attenuation to zero where tir_to_fuv is NaN
        a_fuv_cortese[np.isnan(tir_to_fuv)] = 0.0

        # Set attenuation to zero where sSFR is smaller than zero
        a_fuv_cortese[ssfr < 0.0] = 0.0

        # Set attenuation to zero where sSFR is greater than 10.5
        a_fuv_cortese[ssfr >= 10.5] = 0.0

        # Return the A(FUV) map
        return a_fuv_cortese

    # -----------------------------------------------------------------

    def get_fuv_h(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the FUV-H colour map ...")

        # Calculate the colour map
        fuv_h = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["H"].frames.primary)

        # Replace NaNs by zeros
        fuv_h.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        #fuv_h[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_h[fuv_h < 0.0] = 0.0

        # Mask low sigal-to-noise pixels in the fuv map, if requested
        #if self.config.ssfr.mask_low_fuv_snr: fuv_h[self.fuv < self.config.ssfr.fuv_snr_level*self.fuv_errors] = 0.0

        # Return the FUV-H colour map
        return fuv_h

    # -----------------------------------------------------------------

    def get_fuv_i(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the FUV-i colour map ...")

        # Calculate the colour map
        fuv_i = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["i"].frames.primary)

        # Replace NaNs by zeros
        fuv_i.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        #fuv_i[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_i[fuv_i < 0.0] = 0.0

        # Return the FUV-i colour map
        return fuv_i

    # -----------------------------------------------------------------

    def get_fuv_r(self):

        """
        This function ...
        :return:
        """

        # Calculate the colour map
        fuv_r = -2.5 * np.log10(self.images["FUV"].frames.primary / self.images["r"].frames.primary)

        # Replace NaNs by zeros
        fuv_r.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        #fuv_r[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_r[fuv_r < 0.0] = 0.0

        # Return the FUV-r colour map
        return fuv_r

    # -----------------------------------------------------------------

    def get_fuv_g(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def get_fuv_b(self):

        """
        This function ...
        :return:
        """

        pass

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
        dust_path = filesystem.join(self.maps_path, "dust.fits")

        # Save the dust map
        self.dust.save(dust_path)

        # Determine the path to the old stars map
        old_stars_path = filesystem.join(self.maps_path, "old_stars.fits")

        # Save the old stars map
        self.old_stars.save(old_stars_path)

        # Determine the path to the young stars map
        young_stars_path = filesystem.join(self.maps_path, "young_stars.fits")

        # Save the young stars map
        self.young_stars.save(young_stars_path)

        # Determine the path to the ionizing stars map
        ionizing_stars_path = filesystem.join(self.maps_path, "ionizing_stars.fits")

        # Save the ionizing stars map
        self.ionizing_stars.save(ionizing_stars_path)

# -----------------------------------------------------------------
