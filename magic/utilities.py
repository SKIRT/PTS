#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       Astromagic -- the image editor for Astronomers        **
# *****************************************************************

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy import units as u
from astroquery.irsa_dust import IrsaDust
from astropy import log

# -----------------------------------------------------------------

def print_status(image):

    """
    This function ...
    :param image:
    :return:
    """

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

def reset_selection(image):

    """
    This function ...
    :param image:
    :return:
    """

    image.deselect_all()
    image.frames.primary.select()

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

def conversionfactorP160(image, attenuation):

    pixelscale = image.pixelscale

    pixelfactor = (206264.806247 / pixelscale)**2

    factor = 1e-6 * pixelfactor

    log.info("Converting units with factor = " + str(factor))

    return factor

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------

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

# -----------------------------------------------------------------
