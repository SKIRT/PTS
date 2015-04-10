#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.datapreparation Preparing astronomical images for input in SKIRT
#
# An instance of the DataPreparation class in this module is responsible for taking reduced astronomical image data
# of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
# dust, star formation and old stars.

# -----------------------------------------------------------------

# Import standard modules
import os.path
import numpy as np

# Import astronomical modules
from astropy import units as u

# Import relevant PTS modules
from pts.image import Image
from pts.log import Log

# -----------------------------------------------------------------

# Define the full-width-half-maxima of the different images (in pixels)
fwhmax = {"2MASSH":   None, # FIT THIS TO THE STARS!
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
notsubtracted = {"2MASSH":   True,
                 "GALEXFUV": True,
                 "Ha":       True,
                 "IRAC":     False,
                 "MIPS24":   False,
                 "PACS70":   True,
                 "PACS160":  True}

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

# -----------------------------------------------------------------
#  DataPreparation class
# -----------------------------------------------------------------

## An instance of the DataPreparation class in this module is responsible for taking reduced astronomical image data
#  of a certain galaxy in different photometric filters and creating maps that represent the 2D distribution of
#  dust, star formation and old stars.
## OLD STARS: IRAC, YOUNG NI STARS: FUV, YOUNG I STARS: Ha + 24micron, DUST: H + PACS70 + PACS160
#
class DataPreparation(object):

    ## The constructor
    def __init__(self, basepath, filtername="", plot=False, save=True):

        # Get the name of the galaxy (the name of the base directory)
        self._galaxyname = os.path.basename(basepath)

        # Get the full path to the 'data', 'prep' and 'in' directories
        self._datapath = os.path.join(basepath, "data")
        self._preppath = os.path.join(basepath, "prep")
        self._inpath = os.path.join(basepath, "in")

        # Create the logging mechanism
        self._log = Log(basepath, self._galaxyname)

        # Set the 'plot' and 'save' booleans
        self._plot = plot
        self._save = save

        # Get a list of files in the data directory
        files = [f for f in os.listdir(self._datapath) if os.path.isfile(os.path.join(self._datapath,f)) ]

        # Create a dictionary holding the path of each valid FITS file with a key that represents the filter
        self._filters = dict()

        # Loop over all files in the data directory
        for filename in files:

            # Ignore non-FITS files
            if not filename.endswith(".fits"): continue

            # Add an entry to the filters dictionary
            self._filters[os.path.splitext(filename)[0]] = os.path.join(self._datapath, filename)

    ## This function ...
    def run(self):

        # For each filter
        for filter, filepath in self._filters.items():

            # Import the image
            image = self.importimage(filter)

            # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
            if not image.frames.errors:

                try:
                    path = os.path.join(self._datapath, "error", filter + ".fits")
                    image.importframe(path, "errors")
                except IOError:
                    self._log.warning("No error data found for " + filter)

            # Set the fwhm of the image, if it is not None
            if fwhmax[filter]: image.setfwhm(fwhmax[filter])

            # Mask NaNs, edges and extra user-defined regions
            self.mask(image, edges=edges[filter], extra=extra[filter])

            # Interpolate over the stars indicated by the user (if the fwhmax is None; the PSF will be fitted)
            self.interpolatestars(image, fitpsf=(not fwhmax[filter]))

            # Convert the image into units of MJy / sr
            #toMJysr(image)

            # Subtract the sky
            if filter != "GALEXFUV":

                if notsubtracted[filter]: self.subtractsky(image, image.fwhm)

            else:

                subtractskyFUV(image)

            # Calculate the conversion factor
            if filter == "GALEXFUV":

                totalfactor = conversionfactorFUV(image)

            if filter == "PACS70":

                totalfactor = conversionfactorP70(image)

            if filter == "2MASSH":

                # Convert to magnitudes
                m0 = 20.6473999

                # Convert back to fluxes
                F0 = 1024.0

                totalfactor = conversionfactorH(image)

                totalfactor = totalfactor * F0 * np.power(10.0, -m0/2.5)

            if filter == "IRAC":

                # No corrections for IRAC
                totalfactor = 1.0

            if filter == "MIPS24":

                # No corrections for MIPS 24
                totalfactor = 1.0

            if filter == "PACS160":

                totalfactor = conversionfactorP160(image)

            if filter == "Ha":

                totalfactor = conversionfactorHa(image)

            # Scale the image by the total conversion factor
            self.scale(image, totalfactor)

            # Convolve this image to the PACS 160 micron resolution
            if filter != "PACS160": self.convolve(image, kernels[filter])

            # Rebin this image to the PACS 160 micron pixel grid
            if filter != "PACS160": self.rebin(image)

    # -----------------------------------------------------------------

    ## This function ...
    def importimage(self, filter):

        # Show which image we are importing
        self._log.info("Importing image for " + filter + " band")

        # Determine the path to this FITS file
        filepath = os.path.join(self._datapath, filter + ".fits")

        # Create an image object from this FITS file and return it
        return Image(filepath)

    # -----------------------------------------------------------------

    ## This function masks the nans, edges
    def mask(self, image, edges=True, extra=False):

        # Create a mask for the pixels in the primary image that have a NaN value
        image.masknans()

        # Select the nans mask for later
        image.masks.nans.select()

        if edges:

            # Assuming the nans are in the outer parts of the image, make an expanded mask that also covers the edges
            image.expandmasks("edges")

            # Deselect the nans mask
            image.masks.nans.deselect()

            # Select the edges mask for later
            image.masks.edges.select()

        if extra:

            # Import 'extra' region
            filepath = "data/extra/" + image.name + ".reg"
            image.importregion(filepath, "extra")

            # Select the extra region
            image.regions.extra.select()

            # Create a mask 'extra' from the 'extra' region
            image.createmask()

            # Deselect the extra region
            image.regions.extra.deselect()

            # Select the extra mask for later
            image.masks.extra.select()

        # Combine the edge and extra mask (the currently active masks)
        image.combinemasks("total")

        # Select the total mask
        image.masks.deselectall()
        image.masks.total.select()

        # Apply the total mask to the primary image
        image.applymasks()

        # If requested, save the masked primary image
        if self._save:

            path = os.path.join(self._preppath, image.name, "masked.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

    # -----------------------------------------------------------------

    ## This function interpolates over the stars
    def interpolatestars(self, image, fitpsf=False):

        # Create a region object for the stars create a mask from it
        filepath = "data/stars/" + image.name + ".reg"
        image.importregion(filepath, "stars")

        # Select the stars region
        image.regions.stars.select()

        if fitpsf:

            # Estimate the FWHM of the PSF
            fwhm_x, fwhm_y = image.estimatepsf_fitskirt()

            # Circular approximation
            fwhm = 0.5 * (fwhm_x + fwhm_y)

            # Set the fwhm
            image.setfwhm(fwhm)

        # Create a mask from the stars region
        image.createmask()

        # Select the stars mask
        image.masks.stars.select()

        # Interpolate the image within the stars mask
        image.interpolate()

        # If requested, save the interpolated image
        if self._save:

            path = os.path.join(self._preppath, image.name, "interpolated.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

    # -----------------------------------------------------------------

    ## This function subtracts the sky of the specified image
    def subtractsky(self, image, fwhm):

        # Find the galaxy in the primary image (creates a region called 'galaxy')
        image.findgalaxy()

        # Select the galaxy region
        image.regions.galaxy.select()

        # Create a mask covering the galaxy ==> will also get the name "galaxy"
        image.createmask()

        # Select the total mask and the stars mask
        image.masks.total.select()
        image.masks.stars.select()

        # Find a map that represents the sky, ignoring pixels under either of these 4 masks:
        #  - 2 masks that are currently selected: the total mask and the stars mask
        #  - the galaxy mask
        #  - a mask composed of all pixels outside an annulus around the galaxy
        # and using sigma-clipping to ignore other pixels with extreme values.
        # Estimate the mean, median and standard deviation for the sky
        mean, median, error = image.findsky()

        # Show the mean, median and standard deviation
        self._log.info("Mean sky = " + str(mean))
        self._log.info("Median sky = " + str(median))
        self._log.info("Standard deviation = " + str(error))

        if self._save:

            # Deselect the primary frame and select the sky frame
            image.frames.deselectall()
            image.frames.sky.select()

            # Save the (sigma-clipped) sky as a FITS file
            path = os.path.join(self._preppath, image.name, "sky.fits")
            image.save(path)

        # Select the primary image
        image.frames.deselectall()
        image.frames.primary.select()

        # Fit the sky with a polynomial, using the FHWM
        image.fitsky(fwhm)

        # Deselect all frames and select the fittedsky frame
        image.frames.deselectall()
        image.frames.fittedsky.select()

        if self._save:

            # Save the 2D fitted sky
            path = os.path.join(self._preppath, image.name, "fittedsky.fits")
            image.save(path)

        # Select the total mask
        image.masks.total.select()

        # Subtract the fitted sky (this frame is selected) from the primary image (this frame is deselected)
        image.subtract()

        # Select the primary image frame
        image.frames.deselectall()
        image.frames.primary.select()

        if self._save:

            # Save the sky-subtracted primary image
            path = os.path.join(self._preppath, image.name, "subtracted.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

    # -----------------------------------------------------------------

    ## This function scales the image by a certain factor
    def scale(self, image, factor):

        # Multiply the primary frame by the conversion factor
        image.multiply(factor)

        # If requested, save the scaled image
        if self._save:

            path = os.path.join(self._preppath, image.name, "scaled.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

    # -----------------------------------------------------------------

    ## This function convolves the image with the PSF of the PACS 160 image
    def convolve(self, image, kernel):

        # Convolve to the PACS 160 resolution
        image.convolve(kernel)

        # If requested, save the convolved image
        if self._save:

            path = os.path.join(self._preppath, image.name, "convolved.fits")
            image.save(path)

        # If requested, save the new kernel
        if self._save:

            # Select the kernel frame
            image.frames.deselectall()
            image.frames.kernel.select()

            path = os.path.join(self._preppath, image.name, "newkernel.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

    # -----------------------------------------------------------------

    ## This function rebins the image to the resolution of the Pacs 160 micron image
    def rebin(self, image):

        # Do the rebinning
        pacs160path = os.path.join(self._datapath, "PACS160.fits")
        image.rebin(pacs160path)

        # If requested, save the rebinned image
        if self._save:

            path = os.path.join(self._preppath, image.name, "rebinned.fits")
            image.save(path)

        # Deselect all masks and regions
        image.regions.deselectall()
        image.masks.deselectall()

# -----------------------------------------------------------------

## This function ..
def subtractskyFUV(image):

    log = Log()

    # Automatically finding the galaxy in the GALEX.FUV image fails!
    filepath = "data/galaxy/GALEXFUV.reg"
    image.importregion(filepath, "galaxy")

    # Select the galaxy region
    image.regions.deselectall()
    image.regions.galaxy.select()

    # Create a mask covering the galaxy
    image.createmask()

    ####### SPECIAL THINGS TO DEFINE THE ORIENTATION OF THIS IMAGE WITHOUT HAVING TO USE THE FIND_GALAXY FUNCTION

    # Define the orientation of the galaxy
    class Orientation(object):

        def __init__(self, xpeak, ypeak, major, eps, theta):

            self.xpeak = xpeak
            self.ypeak = ypeak
            self.majoraxis = major
            self.eps = eps
            self.theta = theta

    ypeak = image.regions["galaxy"]._region[0].coord_list[0]
    xpeak = image.regions["galaxy"]._region[0].coord_list[1]
    width = image.regions["galaxy"]._region[0].coord_list[2]
    height = image.regions["galaxy"]._region[0].coord_list[3]
    theta = image.regions["galaxy"]._region[0].coord_list[4]

    major = width / (3.0 * 1.5)
    eps = 1.0 - height / width

    orientation = Orientation(xpeak, ypeak, major, eps, theta)
    image.setorientation(orientation)

    ####### END OF SPECIAL THINGS

    # Select the total mask (edges and extra)
    image.masks.total.select()

    # Find a map that represents the sky, ignoring pixels under either of these 3 masks:
    #  - the mask that is currently selected
    #  - the galaxy mask
    #  - a mask composed of all pixels outside an annulus around the galaxy
    # and using sigma-clipping to ignore other pixels with extreme values.
    # Estimate the mean, median and standard deviation for the sky
    mean, median, error = image.findsky()

    # Show the mean, median and standard deviation
    log.info("Mean sky = " + str(mean))
    log.info("Median sky = " + str(median))
    log.info("Standard deviation = " + str(error))

    # Deselect the primary frame and select the primary_masked_sky frame
    image.frames.deselectall()
    image.frames.sky.select()

    preppathFUV = "prep/GALEXFUV"

    # Save the (sigma-clipped) sky as a FITS file
    path = os.path.join(preppathFUV, "sky.fits")
    image.save(path)

    # Select the primary image
    image.frames.deselectall()
    image.frames.primary.select()

    # Fit the sky with a polynomial, using the FHWM
    image.fitsky(fwhmax["GALEXFUV"])

    # Deselect all frames and select the fittedsky_2D frame
    image.frames.deselectall()
    image.frames.fittedsky.select()

    # Save the 2D fitted sky
    path = os.path.join(preppathFUV, "fittedsky.fits")
    image.save(path)

    # Select the total mask
    image.masks.total.select()

    # Subtract the fitted sky (this frame is selected) from the primary image (this frame is deselected)
    image.subtract()

    # Select the primary image frame
    image.frames.deselectall()
    image.frames.primary.select()

    # Save the sky-subtracted primary image
    path = os.path.join(preppathFUV, "subtracted.fits")
    image.save(path)

# Calculate the conversion factor for the FUV image
def conversionfactorFUV(image):

    log = Log()

    # Get the pixel scale for this filter
    pixelscale = image.pixelscale

    wavelength = image.filter.meanwavelength()
    wavelength = (wavelength * u.micron).to(u.AA)

    # Speed of light in Angstrom per seconds
    c = (299792458 * u.m / u.s).to(u.AA / u.s)
    spectralfactor = wavelength.value**2 / c.value
    #spectralfactor = wavelength**2 / c.value
    pixelfactor = (206264.806247 / pixelscale)**2
    factor = 1.4e-15 * spectralfactor * 1e17 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    # Get the attenuation for this filter
    A = attenuations["GALEXFUV"]

    factor2 = 10**(0.4*A)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

def conversionfactorP70(image):

    log = Log()

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    return factor

def conversionfactorH(image):

    log = Log()

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 10e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    A = attenuations["2MASSH"]

    factor2 = 10**(0.4*A)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

def conversionfactorP160(image):

    log = Log()

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    return factor

def conversionfactorHa(image):

    log = Log()

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    # The wavelength (in meter)
    wavelength = 0.657894736 * 1e-6

    # The speed of light (in m / s)
    c = 299792458

    # The frequency (in Hz)
    frequency = c / wavelength

    factor = 1e23 * 1e-6 * pixelfactor / frequency

    log.info("Converting units with factor = " + str(factor))

    A = attenuations["Ha"]

    factor2 = 10**(0.4*A)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor