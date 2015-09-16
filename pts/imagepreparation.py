#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.imagepreparation Preparing astronomical images as input for
#  SKIRT radiative transfer simulations
#
# An instance of the ImagePreparation class in this module is responsible for taking reduced astronomical image data
# of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
# dust, star formation and old stars.

# *****************************************************************

# Import standard modules
import os.path
import numpy as np

# Import the relevant PTS modules
import astromagic.utilities as iu

# *****************************************************************

# Define the full-width-half-maxima of the different images (in pixels)
fwhmax = {"2MASSH":   None,
          "GALEXFUV": 3.0,
          "Ha":       None,
          "IRAC":     2.5333738673,
          "MIPS24":   4.28666,
          "PACS70":   4.05,
          "PACS160":  3.9228070175}

# Define whether there are stars present in the different images
stars = {"2MASSH": True,
         "GALEXFUV": True,
         "Ha": True,
         "IRAC": True,
         "MIPS24": False,
         "PACS70": False,
         "PACS160": False}

# Define whether saturated stars should be recognized and removed in the image
remove_saturation = {"2MASSH": True,
                     "GALEXFUV": False,
                     "Ha": False,
                     "IRAC": True,
                     "MIPS24": False,
                     "PACS70": False,
                     "PACS160": False}

# Define whether the edges should be masked or not
edges = {"2MASSH":   True,
         "GALEXFUV": True,
         "Ha":       True,
         "IRAC":     False,
         "MIPS24":   True,
         "PACS70":   True,
         "PACS160":  True}

# Define whether we have extra regions to be masked for the different filters
extra = {"2MASSH":   True,
         "GALEXFUV": True,
         "Ha":       True,
         "IRAC":     False,
         "MIPS24":   False,
         "PACS70":   False,
         "PACS160":  False}

# Define which images should still be sky-subtracted
sky_subtracted = {"2MASSH":   False,
                 "GALEXFUV": False,
                 "Ha":       False,
                 "IRAC":     True,
                 "MIPS24":   True,
                 "PACS70":   False,
                 "PACS160":  False}

# Define the kernel files for the different filters
kernels = {"2MASSH":    "Kernel_HiRes_Gauss_03.0_to_PACS_160.fits",
           "GALEXFUV":  "Kernel_HiRes_GALEX_FUV_to_PACS_160.fits",
           "Ha":        "Kernel_HiRes_Gauss_03.0_to_PACS_160.fits",
           "IRAC":      "Kernel_HiRes_IRAC_3.6_to_PACS_160.fits",
           "MIPS24":    "Kernel_HiRes_MIPS_24_to_PACS_160.fits",
           "PACS70":    "Kernel_HiRes_PACS_70_to_PACS_160.fits",
           "PACS160":   None}

# Define the galactic attenuations for the different filters (zero for no attenuation)
attenuations = {"2MASSH":   0.036,
                "GALEXFUV": 0.5606,
                "Ha":       0.174,
                "IRAC":     0.0,
                "MIPS24":   0.0,
                "PACS70":   0.0,
                "PACS160":  0.0}

# *****************************************************************

class ImagePreparation(object):

    """
    An instance of the ImagePreparation class in this module is responsible for taking reduced astronomical image data
    of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
    dust, star formation and old stars.
    OLD STARS: IRAC, YOUNG NI STARS: FUV, YOUNG I STARS: Ha + 24micron, DUST: H + PACS70 + PACS160
    """

    # *****************************************************************

    def __init__(self, directory, filter_name=None, plot=False, save=False):

        """
        This constructor ...
        :param directory:
        :param filter_name:
        :param plot_intermediate:
        :param save_intermediate:
        :return:
        """

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = os.path.basename(directory)

        # Get the full path to the 'data', 'prep' and 'in' directories
        self.data_path = os.path.join(directory, "data")
        self.prep_path = os.path.join(directory, "prep")
        self.in_path = os.path.join(directory, "in")

        # Create the preparation and input directories if they were not yet present
        try: os.mkdir(self.prep_path)
        except OSError: pass
        try: os.mkdir(self.in_path)
        except OSError: pass

        # Set the 'plot_intermediate' and 'save_intermediate' flags
        self.plot = plot
        self.save = save

        # Get a list of files in the data directory
        files = [f for f in os.listdir(self.data_path) if os.path.isfile(os.path.join(self.data_path,f))]

        # Create a dictionary holding the path of each valid FITS file with a key that represents the filter
        self.image_paths = dict()

        # Loop over all files in the data directory
        for filename in files:

            # Ignore non-FITS files or hidden files
            if not filename.endswith(".fits") or filename.startswith("."): continue

            # Get the name of the file without the extension
            base_filename = os.path.splitext(filename)[0]

            # If a filtername was specified, only add the file that corresponds to this filter
            if filter_name is not None:

                if filter_name.lower() == base_filename.lower():

                    self.image_paths[filter_name] = os.path.join(self.data_path, filename)
                    break

            # If no filtername was specified, add each FITS file found in the data directory to the dictionary
            else: self.image_paths[base_filename] = os.path.join(self.data_path, filename)

            # If intermediate results should be saved, create a seperate directory for each filter
            if self.save:

                # Create the directory if it was not yet present
                try: os.mkdir(os.path.join(self.prep_path, base_filename))
                except OSError: pass

    # *****************************************************************

    def run(self):

        """
        This function runs the image preparation procedure
        :return:
        """

        # 1. Get galactic extinction coefficients
        #attenuations = iu.get_attenuations(self.galaxy_name, self.image_paths.keys())

        # 2. Prepare
        self.prepare()

        # 3. Fit bulge and disk
        self.bulge_and_disk()

        # 3. Make maps
        self.make_maps()

    # *****************************************************************

    def prepare(self):

        """
        This function prepares the images
        """

        # Loop over all filters for which we have an image
        for filter_name, path in self.image_paths.items():

            filter_prep_path = os.path.join(self.prep_path, filter_name)
            output_path = filter_prep_path if self.save else None

            # Open the image
            image = iu.open(path)

            # Set the fwhm of the image, if it is not None
            if fwhmax[filter_name] is not None: image.set_fwhm(fwhmax[filter_name])

            # Mask NaNs, edges and extra user-defined regions
            extra_reg = os.path.join(self.data_path, 'extra', image.name + '.reg') if extra[filter_name] else None
            iu.mask(image, edges=edges[filter_name], extra=extra_reg)

            # Interpolate over the stars indicated by the user (if the FWHM is None; the PSF will be fitted)
            if stars[filter_name]: iu.remove_stars(image, self.galaxy_name, model_stars=False, remove_saturation=remove_saturation[filter_name], output_path=output_path)

            # Subtract the sky
            if not sky_subtracted[filter_name]: iu.subtract_sky(image, self.galaxy_name, plot=self.plot, output_path=output_path, downsample_factor=int(round(2*4*image.fwhm)))

            # Convert the units to MJy / sr
            iu.convert_units(image, filter_name, attenuations)

            # Convolve this image to the PACS 160 micron resolution
            if filter_name != "PACS160": iu.convolve(image, kernels[filter_name])

            # Rebin this image to the PACS 160 micron pixel grid
            if filter_name != 'PACS160': iu.rebin(image, self.data_path, 'PACS160.fits')

            # Set the uncertainties
            if filter_name != "2MASSH": iu.set_uncertainty(image, self.prep_path, "noise.reg")

            # If requested, save the result
            if self.save: iu.save(image, filter_prep_path, 'convolved_rebinned.fits')

            # Crop the interesting part of the image
            if image.frames.errors is not None: image.frames.errors.select()
            image.crop(350, 725, 300, 825)

            # If requested, save the result
            iu.save(image, filter_prep_path, 'final.fits')

    # *****************************************************************

    def bulge_and_disk(self):

        # TODO: do the bulge/disk fitting here

        # Set path of bulge and disk images
        bulge_path = os.path.join(self.prep_path, "bulge", "M81_bulge_i59_total.fits")
        disk_path = os.path.join(self.prep_path, "disk", "M81_disk_i59_total.fits")

        # Create list
        paths = {"bulge": bulge_path, "disk": disk_path}

        # For bulge and disk ...
        for name, path in paths.items():

            # Open the image
            image = iu.open(path)

            # Set the header of the image
            image.header["EQUINOX"] = 2000.0
            image.header["NAXIS"] = 2
            image.header["NAXIS1"] = 1000
            image.header["NAXIS2"] = 1000
            image.header["CRPIX1"] = 500.5
            image.header["CRPIX2"] = 500.5
            image.header["CRVAL1"] = 148.8883333
            image.header["CRVAL2"] = +69.06527778
            image.header["CD1_1"] = -4.77942772e-4
            image.header["CD1_2"] = 0.0
            image.header["CD2_1"] = 0.0
            image.header["CD2_2"] = 4.77942772e-4
            image.header["CROTA1"] = 0.0
            image.header["CROTA2"] = 0.0
            image.header["CTYPE1"] = 'RA---TAN'
            image.header["CTYPE1"] = 'DEC--TAN'

            # Convolve the bulge image to the PACS 160 resolution
            iu.convolve(image, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")

            # Rebin the convolved bulge image to the frame of the prepared images (does not matter which one)
            IRAC_dir = os.path.join(self.prep_path, "IRAC")
            iu.rebin(image, IRAC_dir, "convolved_rebinned.fits")
            image.crop(350, 725, 300, 825)

            # Save the convolved, rebinned and cropped bulge or disk image
            iu.save(image, os.path.join(self.prep_path, 'name'), 'final.fits')

    # *****************************************************************

    def make_maps(self):

        """
        This function makes the maps of dust and stars ...
        """

        # Make the specific star formation rate map
        self.make_ssfr_map()

        # Make the total infrared map
        self.make_tir_map()

        # Make the FUV attenuation map
        self.make_attenuation_map()

        # Make the old stars map
        self.make_oldstars_map()

        # Make FUV and MIPS 24 emission maps
        self.make_fuv_and_24_maps()

        # Make ionizing stars map
        self.make_ionizing_stars_map()

    # *****************************************************************

    def make_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Young non-ionizing stars = GALEXFUV - H (sSFR)

        fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        h_path = os.path.join(self.prep_path, "2MASSH", "final.fits")

        fuv_image = iu.open(fuv_path)
        fuv_image.import_datacube(h_path, "H")

        fuv_h_data = -2.5*(np.log10(fuv_image.frames.primary.data) - np.log10(fuv_image.frames.H.data))

        fuv_image._add_frame(fuv_h_data, fuv_image.frames.primary.coordinates, "ssfr")

        fuv_image.deselect_all()

        fuv_image.frames.ssfr.select()

        fuv_image.mask_nans()

        fuv_image.masks.nans.select()

        fuv_image.apply_masks(0.0)

        iu.reset_selection(fuv_image)

        #maak alles 0 waar PACS160 SN<5
        #FUV_H_SN5 = fltarr(376,526)
        #for i=0,375 do begin
           #for j=0,525 do begin
              #if ((imaH_rebin[i,j] GT 0) && (imaFUV_rebin[i,j] GE (10*imaFUVerr_rebin[i,j]))) then begin
                 #FUV_H_SN5[i,j] = FUV_H[i,j]
              #endif else begin
                 #FUV_H_SN5[i,j] = 0
              #endelse
           #endfor
        #endfor
        #fits_write,'M81_FUV_H_StoN_new.fits',FUV_H_SN5,himaFUV_rebin

    # *****************************************************************

    def make_tir_map(self):

        """
        This function ...
        :return:
        """

        D_M81 = 3.6

        factor24 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/24e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor70 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/70e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)
        factor160 = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/160e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        mips24_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        pacs70_path = os.path.join(self.prep_path, "PACS70", "final.fits")
        pacs160_path = os.path.join(self.prep_path, "PACS160", "final.fits")

        # Open the images
        mips24_image = iu.open(mips24_path)
        pacs70_image = iu.open(pacs70_path)
        pacs160_image = iu.open(pacs160_path)

        # Convert from MJy/sr to L_sun
        mips24_image.multiply(10.0**factor24)
        pacs70_image.multiply(10.0**factor70)
        pacs160_image.frames.errors.select()
        pacs160_image.multiply(10.0**factor160)

        # Galametz+2013 formula for Lsun units
        tir_data = (2.133*mips24_image.frames.primary.data) + \
                   (0.681*pacs70_image.frames.primary.data) + \
                   (1.125*pacs160_image.frames.primary.data)

        # Convert TIR from Lsun to W/m2
        factor =  np.log10(3.846e26) - np.log10(4*np.pi) - (2*np.log10(D_M81*3.08567758e22))
        tir_data = tir_data*(10.0**factor)

    # *****************************************************************

    def make_attenuation_map(self):

        # Dust = FUV attenuation = ratio of TIR and FUV luminosity
        # where TIR = 24, 70 and 160 micron

        pass

    # *****************************************************************

    def make_oldstars_map(self):

        # Old stars = IRAC3.6 - bulge
        # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

        irac_path = os.path.join(self.prep_path, "IRAC", "final.fits")
        irac_image = iu.open(irac_path)

        bulge_path = os.path.join(self.prep_path, "bulge", "M81_bulge_i59_total.fits")
        irac_image.import_datacube(bulge_path, "bulge")

        irac_image.frames.bulge.select()

        totalbulge = irac_image.frames.bulge.sum
        flux3_6 = 60541.038
        factor = 0.54 * flux3_6 / totalbulge

        # Multiply the bulge frame with the factor
        #irac_image.multiply(factor)

        # Do the subtraction
        #irac_image.subtract()

        oldstars_data = irac_image.frames.primary.data - factor*irac_image.frames.bulge.data

        for x in range(oldstars_data.shape[1]):
            for y in range(oldstars_data.shape[0]):

                if irac_image.frames.primary.data[y,x] < 10.0*irac_image.frames.errors.data[y,x] or oldstars_data[y,x] < 0.0:

                    oldstars_data[y,x] = 0.0

        # Add old stars frame
        irac_image._add_frame(oldstars_data, irac_image.frames.primary.coordinates, "oldstars")

        irac_image.deselect_all()

        irac_image.frames.oldstars.select()

        # Save old stars map as FITS file
        iu.save(irac_image, self.in_path, "old_stars.fits")

    # *****************************************************************

    def make_fuv_and_24_maps(self):

        """
        This function ...
        :return:
        """

        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        disk_path = os.path.join(self.prep_path, "disk", "M81_disk_i59_total.fits")

        disk_image = iu.open(disk_path)

        totaldisk = disk_image.frames.primary.sum

        fluxFUV = 855.503
        flux24 = 27790.448

        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014

        factorFUV = 0.2*fluxFUV/totaldisk
        factor24 = 0.48*flux24/totaldisk

        fuv_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        fuv_image = iu.open(fuv_path)

        mips_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        mips_image = iu.open(mips_path)

        new_fuv_data = fuv_image.frames.primary.data - factorFUV*disk_image.frames.primary.data
        new_mips_data = mips_image.frames.primary.data - factor24*disk_image.frames.primary.data

        # Set zero where negative
        new_fuv_data[new_fuv_data < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_fuv_data[fuv_image.frames.primary.data < 10.0*fuv_image.frames.errors.data] = 0.0

        # Set zero where negative
        new_mips_data[new_mips_data < 0.0] = 0.0

        # Set zero where low signal-to-noise ratio
        new_mips_data[mips_image.frames.primary.data < 10.0*mips_image.frames.errors.data] = 0.0

        new_fuv_image = iu.new("fuv")
        new_mips_image = iu.new("mips")

        new_fuv_image._add_frame(new_fuv_data, None, "primary")
        new_mips_image._add_frame(new_mips_data, None, "primary")

        new_fuv_image.frames.primary.select()
        new_mips_image.frames.primary.select()

        # Save the new images
        iu.save(new_fuv_image, self.in_path, "fuv.fits")
        iu.save(new_mips_image, self.in_path, "mips.fits")

    # *****************************************************************

    def make_ionizing_stars_map(self):

        #Young ionizing stars = Ha + 0.031 x MIPS24

        D_M81 = 3.6

        # Convert to Lsun
        factorHa = (2*np.log10(2.85/206264.806247)) - 20 + np.log10(3e8/0.657894736e-6) + np.log10(4*np.pi) + (2*np.log10(D_M81*3.08567758e22)) - np.log10(3.846e26)

        ha_path = os.path.join(self.prep_path, "Ha", "final.fits")
        ha_image = iu.open(ha_path)

        ha_image.multiply(10**factorHa)

        new_mips_path = os.path.join(self.in_path, "mips.fits")
        new_mips_image = iu.open(new_mips_path)

        ionizing_stars_data = ha_image.frames.primary.data + 0.031*new_mips_image.frames.primary.data
        ionizing_stars_ratio = ha_image.frames.primary.data / (0.031*new_mips_image.frames.primary.data)

        low_snr = ha_image.frames.primary.data < 10.0*ha_image.frames.errors.data or new_mips_image.frames.primary.data < 10.0*new_mips_image.frames.errors.data

        ionizing_stars_data[low_snr] = 0.0
        ionizing_stars_ratio[low_snr] = 0.0

        ionizing_stars_data[ionizing_stars_data < 0.0] = 0.0
        ionizing_stars_ratio[ionizing_stars_data < 0.0] = 0.0

        ionizing_stars_image = iu.new("ionizingstars")

        ionizing_stars_image._add_frame(ionizing_stars_data, None, "primary")
        ionizing_stars_image._add_frame(ionizing_stars_ratio, None, "ratio")

        ionizing_stars_image.frames.primary.select()
        ionizing_stars_image.frames.ratio.select()

        # Save the new image
        iu.save(ionizing_stars_image, self.in_path, "ionizingstars.fits")

# *****************************************************************



