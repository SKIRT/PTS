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
            if filter != "PACS160": iu.convolve(image, kernels[filter_name])

            # Rebin this image to the PACS 160 micron pixel grid
            if filter != 'PACS160': iu.rebin(image, self.data_path, 'PACS160.fits')

            # NIET voor 2MASSH
            # Determine background noise in the convolved image
            # Determine mean pixel value and sigma in different apertures of the background
            # => list means = [] and sigmas = []
            # a = robust_sigma(means) = standard deviation of mean background derived for different background regions
            # b = median(sigmas) = mean of the standard deviation of pixel-by-pixel variations in different background regions
            #
            # background_uncertainty = sqrt(a^2 + b^2)
            #
            # Construct error maps
            # ima3_6err_rebin[i,j] = (ima3_6_rebin[i,j]/ima3_6_rebin[i,j])*backunc_3_6
            # imaFUVerr_rebin[i,j] = (imaFUV_rebin[i,j]/imaFUV_rebin[i,j])*backunc_FUV
            # ima24err_rebin[i,j] = sqrt((ima24err_rebin[i,j]*ima24err_rebin[i,j])+((ima24_rebin[i,j]/ima24_rebin[i,j])*backunc_24*backunc_24))
            # ima70err_rebin[i,j] = sqrt((ima70err_rebin[i,j]*ima70err_rebin[i,j])+((ima70_rebin[i,j]/ima70_rebin[i,j])*backunc_70*backunc_70))
            # ima160err[i,j] = sqrt((ima160err[i,j]*ima160err[i,j])+((ima160[i,j]/ima160[i,j])*backunc_160*backunc_160))
            # NIET: imaHerr_rebin[i,j] = (imaH_rebin[i,j]/imaH_rebin[i,j])*backunc_H
            # imaHaerr_rebin[i,j] = (imaHa_rebin[i,j]/imaHa_rebin[i,j])*backunc_Ha

            # TODO: add calibration uncertainty!

            # hextract, ima160, hima160, ima160_rebin, hima160_rebin, 350, 725, 300, 825
            # hextract, ima160err, hima160, ima160err_rebin, hima160_rebin, 350, 725, 300, 825

            # If requested, save the result
            if self.save: iu.save(image, filter_prep_path, 'convolved_rebinned.fits')

            # Crop the interesting part of the image
            image.crop(350, 725, 300, 825)

            # If requested, save the result
            if self.save: iu.save(image, filter_prep_path, 'final.fits')

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



        ## a. Old stars = IRAC3.6 - bulge
        #     From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

        IRAC_path = os.path.join(self.prep_path, "IRAC", "final.fits")
        IRAC_image = iu.open(IRAC_path)

        bulge_path = os.path.join(self.prep_path, "bulge", "M81_bulge_i59_total.fits")
        IRAC_image.import_datacube(bulge_path, "bulge")

        IRAC_image.frames.bulge.select()

        totalbulge = IRAC_image.frames.bulge.sum
        flux3_6 = 60541.038
        factor = 0.54 * flux3_6 / totalbulge

        # Multiply the bulge frame with the factor
        IRAC_image.multiply(factor)

        # Do the subtraction
        IRAC_image.subtract()

        # TODO:

        # for pixels where IRAC3.6 >= 10 * IRAC3.6_ERROR && oldstars > 0: OK
        # other pixels: oldstars = 0.0

        # Select the primary frame (= now the old stars map)
        iu.reset_selection(IRAC_image)

        # Save old stars map as FITS file
        iu.save(IRAC_image, self.in_path, "old_stars.fits")




        ## b. Young non-ionizing stars = GALEXFUV - H (sSFR)

        FUV_path = os.path.join(self.prep_path, "GALEXFUV", "final.fits")
        H_path = os.path.join(self.prep_path, "2MASSH", "final.fits")

        FUV_image = iu.open(FUV_path)
        FUV_image.import_datacube(H_path, "2MASSH")

        #FUV_H = -2.5*(alog10(imaFUV_rebin) - alog10(imaH_rebin))

        # TODO: put nans = 0
        # TODO:

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





        ## c. Young ionizing stars = Ha + 0.031 x MIPS24


        ## Compute TIR map


        ## Compute FUV-attenuation map



        ## d. Dust = FUV attenuation = ratio of TIR and FUV luminosity
        ##    where TIR = 24, 70 and 160 micron

        MIPS24_path = os.path.join(self.prep_path, "MIPS24", "final.fits")
        PACS70_path = os.path.join(self.prep_path, "PACS70", "final.fits")
        PACS160_path = os.path.join(self.prep_path, "PACS160", "final.fits")

        MIPS24_image = iu.open(MIPS24_path)
        MIPS24_image.import_datacube(PACS70_path, "PACS70")
        MIPS24_image.import_datacube(PACS160_path, "PACS160")

        #zet alles om van MJy/sr naar Lsun
        #factor24 = (2*alog10(2.85/206264.806247)) - 20 + alog10(3e8/24e-6) + alog10(4*!pi) + (2*alog10(D_M81*3.08567758e22)) - alog10(3.846e26)
        #factor70 = (2*alog10(2.85/206264.806247)) - 20 + alog10(3e8/70e-6) + alog10(4*!pi) + (2*alog10(D_M81*3.08567758e22)) - alog10(3.846e26)
        #factor160 = (2*alog10(2.85/206264.806247)) - 20 + alog10(3e8/160e-6) + alog10(4*!pi) + (2*alog10(D_M81*3.08567758e22)) - alog10(3.846e26)
        #ima24_rebin = ima24_rebin*(10^factor24)
        #ima70_rebin = ima70_rebin*(10^factor70)
        #ima160_rebin = ima160_rebin*(10^factor160)
        #ima160err_rebin = ima160err_rebin*(10^factor160)

        #Galametz+2013 formula for Lsun units
        #LTIR = fltarr(376,526)
        #LTIR = (2.133*ima24_rebin) + (0.681*ima70_rebin) + (1.125*ima160_rebin)




        ## Subtract old stellar contribution from FUV and MIPS 24 emission

        #     From the FUV and 24 micron maps we must subtract the diffuse radiation (old stellar contribution),
        #     for this we typically use an exponential disk
        #     (scale length detemermined by GALFIT)

        #totaldisk = total(ima_disk_rebin)
        #compute total stellar emission in IRAC 3.6 image
        #fluxFUV = 855.503
        #flux24 = 27790.448
        #typisch 20% en 35% respectievelijk
        #48% voor MIPS 24 komt van Lu et al. 2014
        #factorFUV = 0.2*fluxFUV/totaldisk
        #factor24 = 0.48*flux24/totaldisk
        #M81_FUV_stars = imaFUV_rebin - (factorFUV*ima_disk_rebin) ;no subtraction!!!
        #M81_24_stars = ima24_rebin - (factor24*ima_disk_rebin)
        #fits_write,'M51_stellarmass_10000_0Myr_minusbulge.fits',Mstar_10000_0Myr_minusbulge,hima_SDSS

        # TODO:

        #for i=0,375 do begin
           #for j=0,525 do begin
              #if ((imaFUV_rebin[i,j] GE (10*imaFUVerr_rebin[i,j])) && (M81_FUV_stars[i,j] GE 0)) then begin
                 #M81_FUV_stars[i,j] =  M81_FUV_stars[i,j]
              #endif else begin
                 #M81_FUV_stars[i,j] = 0
              #endelse
           #endfor
        #endfor

        #for i=0,375 do begin
           #for j=0,525 do begin
              #if ((ima24_rebin[i,j] GE (10*ima24err_rebin[i,j])) && (M81_24_stars[i,j] GE 0)) then begin
                 #M81_24_stars[i,j] =  M81_24_stars[i,j]
              #endif else begin
                 #M81_24_stars[i,j] = 0
              #endelse
           #endfor
        #endfor

        #fits_write,'M81_FUV_inputSKIRT.fits', M81_FUV_stars, hima3_6_rebin
        #fits_write,'M81_24_inputSKIRT.fits', M81_24_stars, hima3_6_rebin




        # Combineer Ha+24micron tracers


        #omzetten nr Lsun
        #factorHa = (2*alog10(2.85/206264.806247)) - 20 + alog10(3e8/0.657894736e-6) + alog10(4*!pi) + (2*alog10(D_M81*3.08567758e22)) - alog10(3.846e26)
        #imaHa_rebin = imaHa_rebin*(10^factorHa)
        #M81_ionizingstars = imaHa_rebin + (0.031*M81_24_stars)
        #M81_ionizingstars_StoN = fltarr(376,526)
        #M81_ionizingstars_ratio = fltarr(376,526)

        #for i=0,375 do begin
           #for j=0,525 do begin
              #if ((imaHa_rebin[i,j] GE (10*imaHaerr_rebin[i,j])) && (ima24_rebin[i,j] GE (10*ima24err_rebin[i,j]))) then begin
                 #M81_ionizingstars_StoN[i,j] = M81_ionizingstars[i,j]
                 #M81_ionizingstars_ratio[i,j] = imaHa_rebin[i,j] / (0.031*M81_24_stars[i,j])
              #endif else begin
                 #M81_ionizingstars_StoN[i,j] = 0
                 #M81_ionizingstars_ratio[i,j] = 0
              #endelse
           #endfor
        #endfor

        #controleer nog eens dat alles groter dan 0 is
        #for i=0,375 do begin
           #for j=0,525 do begin
              #if (M81_ionizingstars_StoN[i,j] LT 0) then begin
                 #M81_ionizingstars_StoN[i,j] = 0
                 #M81_ionizingstars_ratio[i,j] = 0
              #endif
           #endfor
        #endfor

        #fits_write,'M81_ionizingstars_inputSKIRT.fits', M81_ionizingstars_StoN, hima160_rebin
        #fits_write,'M81_ionizingstars_ratio.fits', M81_ionizingstars_ratio, hima160_rebin

# *****************************************************************
