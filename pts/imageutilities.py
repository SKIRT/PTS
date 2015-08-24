#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import standard modules
import os.path
import logging
import numpy as np

# Import the Image class
from image import Image

# Import astronomical modules
from astropy import units as u

# *****************************************************************

def open(path):

    """
    This function ...
    :param filter_name:
    :return:
    """

    # Inform the user
    logging.info("Importing image " + path)

    # Create the image
    image = Image(path)

    # If no error map was found in the FITS file, try to find a seperate FITS file containing error data
    if image.frames.errors is None:

        try:

            error_path = os.path.join(os.path.dirname(path), "error", os.path.basename(path))
            image.import_datacube(error_path, "errors")

        except IOError: logging.warning("No error data found for " + path)

    # Return the image
    return image

# *****************************************************************

def save(image, directory, name):

    """
    This function ...
    :param image:
    :param directory:
    :param name:
    :return:
    """

    path = os.path.join(directory, name)
    image.export_datacube(path)

# *****************************************************************

def export_region(image, directory, name):

    """
    This function ...
    :param image:
    :param directory:
    :param name:
    :return:
    """

    path = os.path.join(directory, name)
    image.export_region(path)

# *****************************************************************

def mask(image, edges=True, extra=None, plot=False):

    """
    This function masks the 'NaN's, edges ...
    :param image:
    :param edges:
    :param extra:
    :return:
    """

    # Deselect all regions, masks and frames except for the primary frame
    reset_selection(image)

    # Create a mask for the pixels in the primary image that have a NaN value
    image.mask_nans()

    # Select the nans mask for later
    image.masks.nans.select()

    if edges:

        # Assuming the nans are in the outer parts of the image, make an expanded mask that also covers the edges
        image.expand_masks("edges")

        # Deselect the nans mask
        image.masks.nans.deselect()

        # Select the edges mask for later
        image.masks.edges.select()

    if extra is not None:

        # Import 'extra' region
        image.import_region(extra, "extra")

        # Select the extra region
        image.regions.extra.select()

        # Create a mask 'extra' from the 'extra' region
        image.create_mask()

        # Deselect the extra region
        image.regions.extra.deselect()

        # Select the extra mask for later
        image.masks.extra.select()

    # Combine the edge and extra mask (the currently active masks)
    image.combine_masks("total")

    # Select the total mask
    image.masks.deselect_all()
    image.masks.total.select()

    # Apply the total mask to the primary image
    image.apply_masks()

    # If requested, plot the masked primary image
    if plot: image.plot()

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def remove_stars(image, determine_fwhm=False, region_file=None, plot=False, output_path=None):

    """
    This function removes the stars ...
    :param image:
    :param fit_psf:
    :param manual:
    :param plot:
    :param output_path:
    :return:
    """

    # Deselect all regions, masks and frames except for the primary frame
    reset_selection(image)

    # Check if a region file is specified. If not, fetch the star positions automatically from the web
    if region_file is not None: image.import_region(region_file, "stars")
    else: image.find_stars()

    # Select the stars region
    image.regions.stars.select()

    # If requested, save the stars and ufos regions
    if output_path is not None: export_region(image, output_path, "stars.reg")
    image.regions.stars.deselect()
    image.regions.ufos.select()
    if output_path is not None: export_region(image, output_path, "ufos.reg")
    image.regions.deselect_all()

    # Model the stars
    image.regions.stars.select()
    image.model_stars(plot=False, upsample_factor=2.0)

    # If requested, save the modeled and unmodeled stars regions
    image.regions.deselect_all()
    image.regions.modeled_stars.select()
    if output_path is not None: export_region(image, output_path, "modeled_stars.reg")
    image.regions.modeled_stars.deselect()
    image.regions.unmodeled_stars.select()
    if output_path is not None: export_region(image, output_path, "unmodeled_stars.reg")
    image.regions.deselect_all()

    # Select the stars frame
    image.frames.primary.deselect()
    image.frames.stars.select()

    # If requested, save the stars frame
    if output_path is not None: save(image, output_path, "stars.fits")

    # Subtract the stars frame from the primary frame
    image.subtract()

    # Deselect all regions, masks and frames except for the primary frame
    reset_selection(image)

    # If requested, save the star-subtracted image
    if output_path is not None: save(image, output_path, "subtracted_stars.fits")

    # From the region of unmodeled stars, create a new region that represents the 6-sigma contours of these objects
    image.regions.unmodeled_stars.select()
    image.expand_regions(factor=6.0)
    image.regions.deselect_all()
    image.regions.unmodeled_stars_expanded.select()
    image.rename_region("unmodeled_inner")

    # Create a mask from the 6-sigma contours of unmodeled stars and select it
    image.create_mask()
    image.masks.unmodeled.select()

    # Select the unmodeled stars region again, now creating 9-sigma contours
    image.regions.deselect_all()
    image.regions.unmodeled_stars.select()
    image.expand_regions(factor=9.0)
    image.regions.deselect_all()
    image.regions.unmodeled_stars_expanded.select()
    image.rename_region("unmodeled_outer")

    # Deselect all regions, masks and frames except for the primary frame
    reset_selection(image)

    # Select the unmodeled_outer region and the unmodeled_inner mask and interpolate in the primary frame
    image.regions.unmodeled_outer.select()
    image.masks.unmodeled_inner.select()
    image.interpolate_in_regions()

    # If requested, save the star-subtracted image
    if output_path is not None: save(image, output_path, "removed_stars.fits")

    # If requested, obtain and set the FWHM of the PSF
    if determine_fwhm and image.fwhm is None:

        # Select the modeled stars region
        image.deselect_all()
        image.regions.modeled_stars.select()

        # Calculate the fwhm
        sigma = image.mean_radius()
        fwhm = 2.355 * sigma

        # Set the fwhm
        image.set_fwhm(fwhm)

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

    # If requested, plot the star-removed image
    if plot: image.plot()

# *****************************************************************

def subtract_sky(image, plot=False, output_path=None):

    """
    This function subtracts the sky in the image
    :param image:
    :param fwhm:
    :return:
    """

    # Find the galaxy in the primary image (creates a region called 'galaxy')
    image.find_galaxy(plot=plot)

    # Select the galaxy region
    image.regions.galaxy.select()

    # Create a mask covering the galaxy
    image.create_mask()

    # Find a map that represents the sky, ignoring pixels covered by
    image.masks.total.select()
    image.masks.leftovers_inner.select()
    image.masks.ufos.select()
    image.masks.modeled_stars.select()

    # Find a map that represents the sky, ignoring pixels under either of these 4 masks:
    #  - 2 masks that are currently selected: the total mask and the stars mask
    #  - the galaxy mask
    #  - a mask composed of all pixels outside an annulus around the galaxy
    # and using sigma-clipping to ignore other pixels with extreme values.
    image.find_sky()

    # If requested, save the (sigma-clipped) sky frame as a FITS file
    if output_path is not None: save(image, output_path, "sky.fits")

    # Select the primary image
    image.frames.deselect_all()
    image.frames.primary.select()

    # Fit the sky with a 2D polynomial function
    image.fit_polynomial()

    # Deselect all frames and select the fittedsky frame
    image.frames.deselect_all()
    image.frames.primary_polynomial.select()

    # Select the total mask
    image.masks.total.select()

    # Subtract the fitted sky (this frame is selected) from the primary image (this frame is deselected)
    image.subtract()

    # Select the primary image frame
    image.frames.deselect_all()
    image.frames.primary.select()

    # Deselect all masks and regions
    image.regions.deselect_all()
    image.masks.deselect_all()

    # If requested, plot the sky-subtracted primary image
    if plot: image.plot()

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def convert_units(image, filter_name, attenuations):

    """
    This function ...
    :param image:
    :param filter_name:
    :return:
    """

    # Get the attenuation
    attenuation = attenuations[filter_name]

    # Set default conversion factor
    factor = 1.0

    # Calculate the conversion factor
    if filter_name == "GALEXFUV": factor = conversionfactorFUV(image, attenuation)
    if filter_name == "PACS70": factor = conversionfactorP70(image, attenuation)
    if filter_name == "2MASSH":

        # Convert to magnitudes
        m0 = 20.6473999

        # Convert back to fluxes
        F0 = 1024.0

        # Calculate the conversion factor
        factor = conversionfactorH(image, attenuation)
        factor *= F0 * np.power(10.0, -m0/2.5)

    if filter_name == "PACS160": factor = conversionfactorP160(image, attenuation)
    if filter_name == "Ha": factor = conversionfactorHa(image, attenuation)

    # Multiply the primary frame by the conversion factor
    image.multiply(factor)

    # Deselect all masks and regions
    image.regions.deselect_all()
    image.masks.deselect_all()

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def convolve(image, kernel, plot=False):

    """
    This function convolves the image ...
    :param image:
    :param kernel:
    :param plot:
    :param save:
    :return:
    """

    # Convolve to the PACS 160 resolution
    image.convolve(kernel)

    # If requested, save the new kernel
    #if save:

        # Select the kernel frame
        #image.frames.deselect_all()
        #image.frames.kernel.select()

        #path = os.path.join(self.preppath, image.name, "newkernel.fits")
        #image.export_datacube(path)

        # Select the primary frame again
        #image.frames.deselect_all()
        #image.frames.primary.select()

    # Deselect all masks and regions
    image.regions.deselect_all()
    image.masks.deselect_all()

    # If requested, plot the convolved primary image
    if plot: image.plot()

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def rebin(image, directory, name, plot=False):

    """
    This function rebins the image to the resolution of the PACS 160 micron image
    """

    # Do the rebinning
    ref_image_path = os.path.join(directory, name)
    image.rebin(ref_image_path)

    # Deselect all masks and regions
    image.regions.deselect_all()
    image.masks.deselect_all()

    # If requested, plot the rebinned primary image
    if plot: image.plot()

# *****************************************************************

def reset_selection(image):

    """
    This function ...
    :param image:
    :return:
    """

    image.deselect_all()
    image.frames.primary.select()

# *****************************************************************

def subtract_sky_FUV(image):

    """
    This function ...
    :param image:
    :return:
    """

    # Automatically finding the galaxy in the GALEX.FUV image fails!
    filepath = "data/galaxy/GALEXFUV.reg"
    image.import_region(filepath, "galaxy")

    # Select the galaxy region
    image.regions.deselect_all()
    image.regions.galaxy.select()

    # Create a mask covering the galaxy
    image.create_mask()

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
    image.set_orientation(orientation)

    ####### END OF SPECIAL THINGS

    # Select the total mask (edges and extra)
    image.masks.total.select()

    # Find a map that represents the sky, ignoring pixels under either of these 3 masks:
    #  - the mask that is currently selected
    #  - the galaxy mask
    #  - a mask composed of all pixels outside an annulus around the galaxy
    # and using sigma-clipping to ignore other pixels with extreme values.
    image.find_sky()

    # Deselect the primary frame and select the primary_masked_sky frame
    image.frames.deselect_all()
    image.frames.sky.select()

    preppathFUV = "prep/GALEXFUV"

    # Save the (sigma-clipped) sky as a FITS file
    path = os.path.join(preppathFUV, "sky.fits")
    image.export_datacube(path)

    # Select the primary image
    image.frames.deselect_all()
    image.frames.primary.select()

    # Fit the sky with a polynomial, using the FHWM
    #image.fit_sky(fwhmax["GALEXFUV"])

    image.fit_polynomial()

    # Deselect all frames and select the fittedsky_2D frame
    image.frames.deselect_all()
    image.frames.primary_polynomial.select()

    # Save the 2D fitted sky
    path = os.path.join(preppathFUV, "fittedsky.fits")
    image.export_datacube(path)

    # Select the total mask
    image.masks.total.select()

    # Subtract the fitted sky (this frame is selected) from the primary image (this frame is deselected)
    image.subtract()

    # Select the primary image frame
    image.frames.deselect_all()
    image.frames.primary.select()

    # Save the sky-subtracted primary image
    path = os.path.join(preppathFUV, "subtracted.fits")
    image.export_datacube(path)

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def conversionfactorFUV(image, attenuation):

    """
    This function calculates the conversion factor for the FUV image
    :param image:
    :return:
    """

    # Get the pixel scale for this filter
    pixelscale = image.pixelscale

    wavelength = image.filter.centerwavelength()
    wavelength = (wavelength * u.micron).to(u.AA)

    # Speed of light in Angstrom per seconds
    c = (299792458 * u.m / u.s).to(u.AA / u.s)
    spectralfactor = wavelength.value**2 / c.value
    pixelfactor = (206264.806247 / pixelscale)**2
    factor = 1.4e-15 * spectralfactor * 1e17 * pixelfactor

    logging.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    logging.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

# *****************************************************************

def conversionfactorP70(image, attenuation):

    """
    This function calculates the conversion factor for the PACS 70 image
    :param image:
    :return:
    """

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    logging.info("Converting units with factor = " + str(factor))

    return factor

# *****************************************************************

def conversionfactorH(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 10e-6 * pixelfactor

    logging.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    logging.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

# *****************************************************************

def conversionfactorP160(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    logging.info("Converting units with factor = " + str(factor))

    return factor

# *****************************************************************

def conversionfactorHa(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    # The wavelength (in meter)
    wavelength = 0.657894736 * 1e-6

    # The speed of light (in m / s)
    c = 299792458

    # The frequency (in Hz)
    frequency = c / wavelength

    factor = 1e23 * 1e-6 * pixelfactor / frequency

    logging.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    logging.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

# *****************************************************************

