#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import Python 3 functionality
from __future__ import (absolute_import, division, print_function)

# Import standard modules
import os.path
import inspect
import numpy as np
from skimage import morphology

# Import astronomical modules
import astropy.units as u
from astropy import log
import astropy.logger

# Import Astromagic modules
from astromagic import Image
from astromagic.tools import configuration
from astromagic.core.frames import Frame
from astromagic.core.masks import Mask

# *****************************************************************

class MapMaker(object):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Determine the path to the default configuration file
        directory = os.path.dirname(os.path.dirname(inspect.getfile(inspect.currentframe())))
        default_config = os.path.join(directory, "config", "mapmaker.cfg")

        # Open the default configuration if no configuration file is specified, otherwise adjust the default
        # settings according to the user defined configuration file
        if config is None: self.config = configuration.open(default_config)
        else: self.config = configuration.open(config, default_config)

        ### SET-UP LOGGING SYSTEM

        # Set the log level
        log.setLevel(self.config.logging.level)

        # Set log file path
        if self.config.logging.path is not None: astropy.logger.conf.log_file_path = self.config.logging.path.decode('unicode--escape')

        ### TEMPORARY

        # ...

        ### INITIALIZE ATTRIBUTES

        # Input images
        self.h = None
        self.fuv = None
        self.ha = None
        self.irac = None
        self.mips = None
        self.pacsblue = None
        self.pacsred = None

        # Error frames
        self.fuv_errors = None
        self.ha_errors = None
        self.irac_errors = None
        self.mips_errors = None

        # Bulge and disk
        self.disk = None
        self.bulge = None

        # Set the low signal-to-noise mask to None initially
        self.mask = None

    # *****************************************************************

    def run(self):

        """
        This function ...
        :param image:
        :return:
        """

        # Load the input images
        self.load_images()

        # Cut-off the low signal-to-noise pixels
        self.cutoff_low_snr()

        # If requested, save the maps with masked
        if self.config.save_cutoff_maps: self.save_cutoff_maps()

        # Make the specific star formation rate map
        self.make_ssfr_map()

        # Make the FUV attenuation map
        self.make_attenuation_map()

        # Make the old stars map
        self.make_oldstars_map()

        # Make FUV and MIPS 24 emission maps
        self.make_fuv_and_24_maps()

        # Make ionizing stars map
        self.make_ionizing_stars_map()

    # *****************************************************************

    def load_images(self):

        """
        This function ...
        :return:
        """

        ### H - BAND IMAGE

        # Check whether the H-band image is present
        if os.path.isfile(self.config.h_path):

            # Open the H-band image
            self.h = Frame.from_file(self.config.h_path)

        else: raise IOError("Could not find the H-band image")

        ### FUV IMAGE

        # Check whether the FUV image is present
        if os.path.isfile(self.config.fuv_path):

            # Open the FUV image
            self.fuv = Frame.from_file(self.config.h_path)
            self.fuv_errors = Frame.from_file(self.config.h_path, plane=self.config.errors)

        else: raise IOError("Could not find the FUV image")

        ### H ALPHA IMAGE

        # Check whether the H Alpha image is present
        if os.path.isfile(self.config.ha_path):

            # Open the H Alpha image
            self.ha = Frame.from_file(self.config.ha_path)
            self.ha_errors = Frame.from_file(self.config.ha_path, plane=self.config.errors)

        else: raise IOError("Could not find the H Alpha image")

        ### 24 MICRON IMAGE

        # Check whether the 24 micron image is present
        if os.path.isfile(self.config.mips_path):

            # Open the 24 micron image
            self.mips = Frame.from_file(self.config.mips_path)
            self.mips_errors = Frame.from_file(self.config.mips_path, plane=self.config.errors)

        else: raise IOError("Could not find the 24 micron image")

        ### 3.6 MICRON IMAGE

        # Check whether the 3.6 micron image is present
        if os.path.isfile(self.config.irac_path):

            # Open the 3.6 micron image
            self.irac = Frame.from_file(self.config.irac_path)
            self.irac_errors = Frame.from_file(self.config.irac_path, plane=self.config.errors)

        else: raise IOError("Could not find the 3.6 micron image")

        ### 70 MICRON IMAGE

        # Check whether the 70 micron image is present
        if os.path.isfile(self.config.pacsblue_path):

            # Open the 70 micron image
            self.pacsblue = Frame.from_file(self.config.pacsblue_path)

        ### 160 MICRON IMAGE

        # Check whether the 160 micron image is present
        if os.path.isfile(self.config.pacsred_path):

            # Open the 160 micron image
            self.pacsred = Frame.from_file(self.config.pacsred_path)

        else: raise IOError("Could not find the 160 micron image")

        ### DISK AND BULGE

        # Check whether the disk image is present
        if os.path.isfile(self.config.disk_path):

            # Open the disk image
            self.disk = Frame.from_file(self.config.disk_path)

        else: raise IOError("Could not find the disk image")

        # Check whether the bulge image is present
        if os.path.isfile(self.config.bulge_path):

            # Open the bulge image
            self.bulge = Frame.from_file(self.config.bulge_path)

        else: raise IOError("Could not find the bulge image")

    # *****************************************************************

    def cutoff_low_snr(self):

        """
        This function ...
        :return:
        """

        # Check whether the reference image is present
        if os.path.isfile(self.config.cutoff.reference_path):

            # Open the reference image
            reference = Image(self.config.cutoff.reference_path)

            # Check whether the errors frame is present
            assert self.config.errors in reference.frames, "An error map could not be found for the reference image"

            # Create a mask for the pixels with a signal-to-noise ratio of 5 or less
            data = reference.frames[self.config.primary] < self.config.cutoff.level*reference.frames[self.config.errors]
            self.mask = Mask(data)

            # If requested, perform binary opening to remove small patches covering inner parts of the galaxy
            if self.config.cutoff.opening:

                disk_structure = morphology.disk(4)

                Frame(self.mask.astype(float)).save(self.config.saving.cutoff_mask_before_opening_path)
                self.mask = self.mask.opening(structure=disk_structure, iterations=1)

            # Save the mask as a FITS file
            Frame(self.mask.astype(float)).save(self.config.saving.cutoff_mask_path)

        # If not, raise an error
        else: raise IOError("The prepared reference image could not be found")

        # Cut-off the input images at the same contour level
        self.h[self.mask] = 0.0
        self.fuv[self.mask] = 0.0
        self.ha[self.mask] = 0.0
        self.irac[self.mask] = 0.0
        self.mips[self.mask] = 0.0
        self.pacsblue[self.mask] = 0.0
        self.pacsred[self.mask] = 0.0

        # Cut-off the bulge and disk images at the same contour level
        self.disk[self.mask] = 0.0
        self.bulge[self.mask] = 0.0

    # *****************************************************************

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

    # *****************************************************************

    def make_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the sSFR map")

        # Young non-ionizing stars (specific star formation rate) = GALEXFUV - H
        fuv_h = -2.5*(np.log10(self.fuv) - np.log10(self.h))

        # Mask nans in the sSFR map
        fuv_h.replace_nans(0.0)

        # Mask pixels outside of the low signal-to-noise contour
        fuv_h[self.mask] = 0.0

        # Set negative pixels to zero
        fuv_h[fuv_h < 0.0] = 0.0

        # Mask low sigal-to-noise pixels in the fuv map, if requested
        if self.config.ssfr.mask_low_fuv_snr: fuv_h[self.fuv < self.config.ssfr.level*self.fuv_errors] = 0.0

        # Save the resulting sSFR map
        fuv_h.save(self.config.ssfr.output_path)

        # Cache
        self.ssfr = fuv_h

    # *****************************************************************

    @property
    def tir(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the TIR map")

        D_M81 = 3.6

        factor24 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/24e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor70 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/70e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor160 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/160e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        # Convert from MJy/sr to L_sun
        mips = self.mips * 10.0**factor24
        pacsblue = self.pacsblue * 10.0**factor70
        pacsred = self.pacsred * 10.0**factor160  ## ERROR ?

        # Galametz+2013 formula for Lsun units
        tir_data = (2.133*mips) + (0.681*pacsblue) + (1.125*pacsred)

        # Convert TIR from Lsun to W/m2
        factor = np.log10(3.846e26) - np.log10(4*np.pi) - (2.0*np.log10(D_M81*3.08567758e22))
        tir_data *= 10.0**factor

        # Return the TIR map
        return tir_data

    # *****************************************************************

    def make_attenuation_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating dust attenuation map")

        # Dust = FUV attenuation = ratio of TIR and FUV luminosity

        # Convert the FUV map
        factor = - 20 + np.log10(3e8) - np.log10(0.153e-6) + (2*np.log10(2.85/206264.806247))
        fuv_converted = self.fuv * 10**factor

        # The ratio of TIR and FUV
        tir_to_fuv = np.log10(self.tir/fuv_converted)

        # Calculate powers of tir_to_fuv
        tir_to_fuv2 = np.power(tir_to_fuv, 2.0)
        tir_to_fuv3 = np.power(tir_to_fuv, 3.0)
        tir_to_fuv4 = np.power(tir_to_fuv, 4.0)

        #a_fuv_buat = (-0.0333*x3) + (0.3522*x2) + (1.1960*x) + 0.4967

        # Create an empty image
        a_fuv_cortese = Frame.zeros_like(tir_to_fuv)

        limits = []
        a = []
        b = []
        c = []
        d = []
        e = []

        limits.append((None, 4.0))
        a.append(0.50994)
        b.append(0.88311)
        c.append(0.53315)
        d.append(0.04004)
        e.append(0.04883)

        limits.append((4.0, 4.2))
        a.append(0.49867)
        b.append(0.86377)
        c.append(0.51952)
        d.append(0.04038)
        e.append(0.04624)

        limits.append((4.2, 4.6))
        a.append(0.49167)
        b.append(0.85201)
        c.append(0.51152)
        d.append(0.04060)
        e.append(0.04475)

        limits.append((4.6, 5.0))
        a.append(0.48223)
        b.append(0.83642)
        c.append(0.50127)
        d.append(0.04092)
        e.append(0.04288)

        limits.append((5.0, 5.4))
        a.append(0.46909)
        b.append(0.81520)
        c.append(0.48787)
        d.append(0.04138)
        e.append(0.04050)

        limits.append((5.4, 5.8))
        a.append(0.45013)
        b.append(0.78536)
        c.append(0.47009)
        d.append(0.04210)
        e.append(0.03745)

        limits.append((5.8, 6.3))
        a.append(0.42168)
        b.append(0.74191)
        c.append(0.44624)
        d.append(0.04332)
        e.append(0.03362)

        limits.append((6.3, 6.6))
        a.append(0.40210)
        b.append(0.71272)
        c.append(0.43139)
        d.append(0.04426)
        e.append(0.03140)

        limits.append((6.6, 6.9))
        a.append(0.37760)
        b.append(0.67674)
        c.append(0.41420)
        d.append(0.04555)
        e.append(0.02900)

        limits.append((6.9, 7.2))
        a.append(0.34695)
        b.append(0.63224)
        c.append(0.39438)
        d.append(0.04739)
        e.append(0.02650)

        limits.append((7.2, 7.5))
        a.append(0.30899)
        b.append(0.57732)
        c.append(0.37157)
        d.append(0.05000)
        e.append(0.02399)

        limits.append((7.5, 7.8))
        a.append(0.26302)
        b.append(0.51013)
        c.append(0.34522)
        d.append(0.05377)
        e.append(0.02164)

        limits.append((7.8, 8.1))
        a.append(0.20982)
        b.append(0.42980)
        c.append(0.31431)
        d.append(0.05909)
        e.append(0.01957)

        limits.append((8.1, 8.4))
        a.append(0.15293)
        b.append(0.33799)
        c.append(0.27713)
        d.append(0.06638)
        e.append(0.01792)

        limits.append((8.4, 8.8))
        a.append(0.09944)
        b.append(0.24160)
        c.append(0.23161)
        d.append(0.07580)
        e.append(0.01671)

        limits.append((8.8, 9.2))
        a.append(0.05822)
        b.append(0.15524)
        c.append(0.17801)
        d.append(0.08664)
        e.append(0.01593)

        limits.append((9.2, 9.6))
        a.append(0.03404)
        b.append(0.09645)
        c.append(0.12452)
        d.append(0.09679)
        e.append(0.01548)

        limits.append((9.6, 10.0))
        a.append(0.02355)
        b.append(0.06934)
        c.append(0.08725)
        d.append(0.10339)
        e.append(0.01526)

        limits.append((10.0, 10.5))
        a.append(0.02025)
        b.append(0.06107)
        c.append(0.07212)
        d.append(0.10588)
        e.append(0.01517)

        # Create the FUV attenuation map
        for i in range(len(limits)):

            if limits[i][0] is None: where = self.ssfr < limits[i][1]
            elif limits[i][1] is None: where = self.ssfr > limits[i][0]
            else: where = (self.ssfr >= limits[i][0]) * (self.ssfr < limits[i][1])

            # Set the appropriate pixels
            a_fuv_cortese[where] = a[i] + b[i]*tir_to_fuv[where] + c[i]*tir_to_fuv2[where] + d[i]*tir_to_fuv3[where] - e[i]*tir_to_fuv4[where]

        # Set attenuation to zero where sSFR is smaller than zero
        a_fuv_cortese[self.ssfr < 0.0] = 0.0

        # Set attenuation to zero where sSFR is greater than 10.5
        a_fuv_cortese[self.ssfr >= 10.5] = 0.0

        # Set attenuation to zero where FUV is smaller than zero
        a_fuv_cortese[self.fuv <= 0.0] = 0.0

        # Mask pixels outside of the low signal-to-noise contour
        a_fuv_cortese[self.mask] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        a_fuv_cortese[a_fuv_cortese < 0.0] = 0.0

        # Save FUV attenuation map as a FITS file
        a_fuv_cortese.save(self.config.attenuation.output_path)

    # *****************************************************************

    def make_oldstars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating old stars map")

        # Old stars = IRAC3.6 - bulge
        # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

        # Calculate factor
        flux3_6 = 60541.038
        factor = 0.54 * flux3_6 / self.bulge.sum()

        # Create the old stars map
        oldstars = self.irac - factor*self.bulge

        # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
        oldstars[self.irac < 10.0*self.irac_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        oldstars[oldstars < 0.0] = 0.0

        # Save the old stars map as a FITS file
        oldstars.save(self.config.oldstars.output_path)

    # *****************************************************************

    @property
    def fuv_young_stars(self):

        """
        This function ...
        :return:
        """

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        # TODO: calculate this flux value
        flux_fuv = 855.503

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.2*flux_fuv/self.disk.sum()

        # Subtract the disk contribution to the FUV image
        new_fuv = self.fuv - factor*self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_fuv[new_fuv < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_fuv[self.fuv < 10.0*self.fuv_errors] = 0.0

        # Return the new FUV frame
        return new_fuv

    # *****************************************************************

    @property
    def mips_young_stars(self):

        """
        This function ...
        :return:
        """

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        # TODO: calculate this flux value
        flux_mips = 27790.448

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factor = 0.48*flux_mips/self.disk.sum()

        # Subtract the disk contribution to the 24 micron image
        new_mips = self.mips - factor*self.disk

        # Make sure all pixels of the disk-subtracted maps are larger than or equal to zero
        new_mips[new_mips < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_mips[self.mips < 10.0*self.mips_errors] = 0.0

        # Return the new 24 micron frame
        return new_mips

    # *****************************************************************

    def make_ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ionizing young stars map")

        #Young ionizing stars = Ha + 0.031 x MIPS24

        D_M81 = 3.6

        # Convert to Lsun
        factor_ha = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/0.657894736e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        # Multiply
        ha = self.ha * 10**factor_ha

        #new_mips_path = os.path.join(self.in_path, "mips.fits")
        #new_mips_image = iu.open(new_mips_path)

        ionizing = ha + 0.031*self.mips

        ionizing_ratio = ha / (0.031*self.mips)

        low_snr = () + ()

        # Set pixels to zero with low signal-to-noise in the H Alpha image
        ionizing[ha < 10.0*self.ha_errors] = 0.0
        ionizing_ratio[ha < 10.0*self.ha_errors] = 0.0

        # Set pixels to zero with low signal-to-noise in the 24 micron image
        ionizing[self.mips < 10.0*self.mips_errors] = 0.0
        ionizing_ratio[self.mips < 10.0*self.mips_errors] = 0.0

        # Make sure all pixel values are larger than or equal to zero
        ionizing[ionizing < 0.0] = 0.0
        ionizing_ratio[ionizing < 0.0] = 0.0

        # Create a new image
        image = Image()
        image.name = "ionizingstars"

        # Add the ionizing stars map and the ratio map as image frames
        image.add_frame(ionizing, "primary")
        image.add_frame(ionizing_ratio, "ratio")

        # Select both frames
        image.frames.primary.select()
        image.frames.ratio.select()

        # Save the ionizing stars image as a FITS file
        image.save(self.config.ionizing.output_path)

# *****************************************************************
