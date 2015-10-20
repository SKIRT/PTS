#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import os.path
import numpy as np

# Import the Image class
from . import Image
from .core.frames import Frame

# Import astronomical modules
from astropy import units as u
from astroquery.irsa_dust import IrsaDust
from astropy import log



def print_status(image):

    # Print the status of the image frames
    log.info("Frames:")
    frames_state = image.frames.get_state()
    for frame_name, selected in frames_state.items():

        log.info("  " + frame_name + ": " + str(selected))

    # Print the status of the image regions
    log.info("Regions:")
    regions_state = image.regions.get_state()
    for region_name, selected in regions_state.items():

        log.info("  " + region_name + ": " + str(selected))

    # Print the status of the image masks
    log.info("Masks:")
    masks_state = image.masks.get_state()
    for mask_name, selected in masks_state.items():

        log.info("  " + mask_name + ": " + str(selected))

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

def remove_stars(image, galaxy_name, region_file=None, model_stars=False, remove_saturation=False, plot=False, output_path=None):

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
    elif image.fwhm is not None:

        sigma = image.fwhm / 2.355 # in arcseconds
        image.fetch_stars(sigma/image.pixelscale, galaxy_name=galaxy_name)

    else:

        image.find_stars(galaxy_name, failed_stars_method="mean")

    if remove_saturation:

        image.regions.stars.select()
        image.split_region(criterium="flux", method="percentage", percentage=0.08)

        image.regions.stars.deselect()
        image.regions.bright.select()

        image.rename_region("bright_stars")

    # Select the stars region
    image.regions.deselect_all()
    image.regions.stars.select()

    # If requested, save the stars region
    if output_path is not None: export_region(image, output_path, "stars.reg")
    image.regions.stars.deselect()

    # If requested and present, save the ufos region
    #if image.regions.ufos is not None:

        #image.regions.ufos.select()
        #if output_path is not None: export_region(image, output_path, "ufos.reg")
        #image.regions.deselect_all()

    # If requested, model the stars
    #if model_stars:

        #image.regions.stars.select()
        #image.model_stars(plot=False, upsample_factor=2.0)

        # If requested, save the modeled and unmodeled stars regions
        #image.regions.deselect_all()
        #image.regions.modeled_stars.select()
        #if output_path is not None: export_region(image, output_path, "modeled_stars.reg")
        #image.regions.modeled_stars.deselect()
        #image.regions.unmodeled_stars.select()
        #if output_path is not None: export_region(image, output_path, "unmodeled_stars.reg")
        #image.regions.deselect_all()

        # Select the stars frame
        #image.frames.primary.deselect()
        #image.frames.stars.select()

        # If requested, save the stars frame
        #if output_path is not None: save(image, output_path, "stars.fits")

        # Subtract the stars frame from the primary frame
        #image.subtract()

        # Deselect all regions, masks and frames except for the primary frame
        #reset_selection(image)

        # If requested, save the star-subtracted image
        #if output_path is not None: save(image, output_path, "subtracted_stars.fits")

        # Set the name of the region within which we are going to interpolate
        #region_for_interpolation = "unmodeled_stars"
        #region_for_fwhm = "modeled_stars"

    #else:

    region_for_interpolation = "stars"
    region_for_fwhm = "stars"

    # From the region of (unmodeled) stars, create a new region that represents the 6-sigma contours of these objects
    image.regions[region_for_interpolation].select()
    image.expand_regions(factor=6.0)
    image.regions.deselect_all()
    image.regions[region_for_interpolation+"_expanded"].select()
    image.rename_region("inner")

    # Create a mask from the 6-sigma contours of unmodeled stars and select it
    image.create_mask()
    image.masks.inner.select()

    # Select the unmodeled stars region again, now creating 9-sigma contours
    image.regions.deselect_all()
    image.regions[region_for_interpolation].select()
    image.expand_regions(factor=9.0)
    image.regions.deselect_all()
    image.regions[region_for_interpolation+"_expanded"].select()
    image.rename_region("outer")

    # Deselect all regions, masks and frames except for the primary frame
    reset_selection(image)

    # Select the unmodeled_outer region and the unmodeled_inner mask and interpolate in the primary frame
    image.regions.outer.select()
    image.masks.inner.select()
    image.interpolate_in_regions()

    # Obtain and set the FWHM of the PSF for the image
    if image.fwhm is None:

        # Select the modeled stars region
        reset_selection(image)
        image.regions[region_for_fwhm].select()

        # Calculate the fwhm
        sigma = image.mean_radius()
        fwhm = 2.355 * sigma # in pixels

        # Set the fwhm
        image.set_fwhm(fwhm*image.pixelscale)

    # Remove saturated stars
    if remove_saturation and image.regions.bright_stars is not None:

        # Deselect all regions, masks and frames except for the primary frame
        reset_selection(image)

        # If requested, save the star-subtracted image
        if output_path is not None: save(image, output_path, "removed_with_saturation.fits")

        # Create a region of 20-sigma contours around the brightest stars and select it
        image.regions.bright_stars.select()
        image.expand_regions(factor=20.0)
        image.regions.deselect_all()
        image.regions.bright_stars_expanded.select()

        if output_path is not None: export_region(image, output_path, "bright_stars.reg")

        # Create a mask for all segments found within the ellipses of that region and select it
        image.create_segmentation_mask(image.fwhm, int(image.fwhm*2.0))
        image.regions.deselect_all()
        image.masks.segments.select()
        image.rename_mask("saturation")

        #print image.masks.get_state()

        # Interpolate the primary frame within the masked pixels
        image.interpolate()
        image.masks.deselect_all()

    # If requested, save the star-subtracted image
    if output_path is not None: save(image, output_path, "removed_stars.fits")

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

    # If requested, plot the star-removed image
    if plot: image.plot()

# *****************************************************************

def subtract_sky(image, galaxy_name, plot=False, output_path=None, downsample_factor=1):

    """
    This function subtracts the sky in the image
    :param image:
    :param fwhm:
    :return:
    """

    # Find the galaxy in the primary image (creates a region called 'galaxy') and select the resulting region
    image.find_galaxy(galaxy_name, plot=plot)
    image.regions.galaxy.select()

    # If an output path is provided, save the galaxy region
    #if output_path is not None: export_region(image, output_path, "galaxy.reg")

    # Expand the galaxy region by a factor of 5.0 and create a mask from this new region
    image.expand_regions(factor=5.0)
    image.regions.deselect_all()
    image.regions.galaxy_expanded.select()
    image.create_mask()

    # If an output path is provided, save the expanded galaxy region
    if output_path is not None: export_region(image, output_path, "galaxy_expanded.reg")

    # If present, create a mask from the modeled stars region and select this mask
    #if image.regions.modeled_stars is not None:

        #image.regions.deselect_all()
        #image.regions.modeled_stars.select()
        #image.expand_regions(factor=6.0)
        #image.regions.deselect_all()
        #image.regions.modeled_stars_expanded.select()
        #image.create_mask()     # Created masks/modeled_stars_expanded
        #image.regions.deselect_all()
        #image.masks.modeled_stars_expanded.select() # Select the new mask

    # Select all other masks that cover parts not suitable for fitting the sky (galaxy, stars)
    if image.masks.inner is not None: image.masks.inner.select()                  # stars
    image.masks.galaxy_expanded.select() # galaxy

    if image.masks.saturation is not None: image.masks.saturation.select() # saturation
    image.masks.total.select()

    # Make a sky map
    image.frames.deselect_all()
    image.frames.primary.select()
    image.copy_frames()
    image.frames.deselect_all()
    image.frames.primary_copy.select()
    image.rename_frame("sky")
    image.apply_masks()

    # If requested, save the sky map as a FITS file
    if output_path is not None: save(image, output_path, "sky.fits")

    # Select the primary image
    image.frames.deselect_all()
    image.regions.deselect_all()
    image.frames.primary.select()

    # Fit the sky with a polynomial function (ignoring pixels covered by any of the selected masks)
    if downsample_factor == 1:

        image.fit_polynomial(plot=plot)
        image.frames.deselect_all()
        image.frames.primary_polynomial.select()
        image.rename_frame("background")

    else: image.estimate_background(downsample_factor=downsample_factor, plot=plot)

    # Select the new frame and deselect all masks
    image.frames.deselect_all()
    image.frames.background.select()
    image.masks.deselect_all()

    # If requested, save the fitted sky map as a FITS file
    if output_path is not None: save(image, output_path, "sky_polynomial.fits")

    # Subtract the fitted sky (this frame is selected) from the primary image (this frame is deselected)
    image.subtract()

    # Select the primary image frame and deselect all masks and regions
    reset_selection(image)

    # If requested, plot the sky-subtracted primary image
    if plot: image.plot()

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

    # Select the errors frame if present
    if image.frames.errors is not None: image.frames.errors.select()

    # Multiply the primary frame (and the errors frame) by the conversion factor
    image.multiply(factor)

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

    # Select the errors frame if present
    if image.frames.errors is not None: image.frames.errors.select()

    # Convolve to the PACS 160 resolution
    image.convolve_fits(kernel)

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

    # Select the errors frame if present
    if image.frames.errors is not None: image.frames.errors.select()

    # Do the rebinning
    ref_image_path = os.path.join(directory, name)
    image.rebin(ref_image_path)

    # Deselect all masks and regions
    image.regions.deselect_all()
    image.masks.deselect_all()

    # If requested, plot the rebinned primary image
    if plot: image.plot()

    # Deselect all regions, masks and frames (except the primary frame)
    reset_selection(image)

# *****************************************************************

def set_uncertainty(image, directory, name):

    """
    This function ...
    :param image:
    :param directory:
    :param name:
    :return:
    """

    # NIET voor 2MASSH
    # Determine background noise in the convolved image
    # Determine mean pixel value and sigma in different apertures of the background
    # => list means = [] and sigmas = []
    # a = robust_sigma(means) = standard deviation of mean background derived for different background regions
    # b = median(sigmas) = mean of the standard deviation of pixel-by-pixel variations in different background regions
    #
    # background_uncertainty = sqrt(a^2 + b^2)

    # TODO: add calibration uncertainty!

    import_region(image, directory, name, "noise")

    image.regions.noise.select()

    uncertainty = image.uncertainty_from_regions()

    if image.frames.errors is None:

        image.add_frame(Frame(np.full(image.frames.primary.shape, uncertainty), image.frames.primary.wcs), "errors")

    else:

        image.frames.errors = np.sqrt(np.power(image.frames.errors, 2) + uncertainty**2)

    # Deselect ...
    reset_selection(image)

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

    log.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

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

    log.info("Converting units with factor = " + str(factor))

    return factor

# *****************************************************************

def conversionfactorH(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 10e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

# *****************************************************************

def conversionfactorP160(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

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

    log.info("Converting units with factor = " + str(factor))

    factor2 = 10**(0.4*attenuation)

    log.info("Correcting for galactic extinction with factor = " + str(factor2))

    # Calculate the total conversion factor
    totalfactor = factor * factor2

    return totalfactor

# *****************************************************************

def get_attenuations(galaxy_name, filter_names):

    """
    This function ...
    :param galaxy_name:
    :param filter_names:
    :return:
    """

    # Get the table of extinctions in various bands
    table = IrsaDust.get_extinction_table(galaxy_name)

    # For each filter
    for filter_name in filter_names:

        pass

# *****************************************************************
