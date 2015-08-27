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

# Import astronomical modules
from astropy import log

# Import the relevant PTS modules
import imageutilities as iu

# *****************************************************************

# Define the full-width-half-maxima of the different images (in pixels)
fwhmax = {"2MASSH":   None,
          "GALEXFUV": 3.0,
          "Ha":       None,
          "IRAC":     2.5333738673,
          "MIPS24":   4.28666,
          "PACS70":   4.05,
          "PACS160":  3.9228070175}

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

        # 2. Prepare the images
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
            iu.remove_stars(image, determine_fwhm=True, output_path=output_path)

            # Subtract the sky
            if not sky_subtracted[filter_name]: iu.subtract_sky(image, self.galaxy_name, plot=self.plot, output_path=output_path, downsample_factor=int(round(2*4*image.fwhm)))

            # Convert the units to MJy / sr
            iu.convert_units(image, filter_name, attenuations)

            # Convolve this image to the PACS 160 micron resolution
            #if filter != "PACS160": iu.convolve(image, kernels[filter_name])

            # Rebin this image to the PACS 160 micron pixel grid
            if filter != 'PACS160': iu.rebin(image, self.data_path, 'PACS160.fits')

            # If requested, save the result
            if self.save: iu.save(image, filter_prep_path, 'final.fits')

        # 2. Make maps of stars and dust


# *****************************************************************
