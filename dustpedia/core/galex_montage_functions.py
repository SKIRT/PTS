#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.galex_montage_functions Functions for the mosaicing of GALEX tiles.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import math
import multiprocessing as mp
import numpy as np
import scipy.spatial
import scipy.ndimage
import lmfit
import shutil

# Import astronomical modules
from astropy.io import fits
import astropy.io.votable
import astropy.convolution
import montage_wrapper as montage
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame, sum_frames, sum_frames_quadratically
from ...magic.core.image import Image
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.misc import chrisfuncs

# -----------------------------------------------------------------

# Flux zero point for converting AB magnitudes to Jansky
ab_mag_zero_point = 3631. * Unit("Jy")

# Zero points for conversion from GALEX count/s to AB magnitude system
galex_fuv_zero_point = 18.82
galex_nuv_zero_point = 20.08

# -----------------------------------------------------------------

def level_chi_squared(level_params, image):

    """
    Fit flat plane to the image to find level
    :param level_params:
    :param image:
    :return:
    """

    level = level_params['level'].value
    chi = image - level
    chisq = chi**2.0

    return chisq

# -----------------------------------------------------------------

def level_galex_maps(fitsfile_dir, convfile_dir, target_suffix):

    """
    original name: GALEX_Zero
    Set a set of maps to the same level
    :param fitsfile_dir:
    :param convfile_dir:
    :param target_suffix:
    :return:
    """

    # Inform the user
    log.info("Inside the 'level_galex_maps' function with:")
    log.info(" - fitsfile_dir = " + fitsfile_dir)
    log.info(" - convfile_dir = " + convfile_dir)
    log.info(" - target_suffix = " + target_suffix)

    # Make list of files in target directory that have target suffix
    allfile_list = os.listdir(fitsfile_dir)
    fitsfile_list = []
    for allfile in allfile_list:
        if target_suffix in allfile:
            fitsfile_list.append(allfile)

    # Loop over each file
    for i in range(0, len(fitsfile_list)):

        # Inform the user
        log.info('Matching background of map ' + fitsfile_list[i] + " ...")

        # Read in corresponding map from directory containing convolved images
        fitsdata_conv = fits.open(convfile_dir + '/' + fitsfile_list[i])
        image_conv = fitsdata_conv[0].data
        fitsdata_conv.close()

        # Fit to level of image; save if first image, otherwise calculate appropriate offset
        level_params = lmfit.Parameters()

        #level_params.add('level', value=np.nanmedian(image_conv), vary=True) ## FOR NUMPY VERSION 1.9.0 AND ABOVE
        image_conv_nonans = image_conv[np.logical_not(np.isnan(image_conv))]
        level_params.add('level', value=np.median(image_conv_nonans), vary=True) # BELOW NUMPY VERSION 1.9.0
        #

        image_conv_clipped = chrisfuncs.SigmaClip(image_conv, tolerance=0.005, median=False, sigma_thresh=3.0)[2]
        level_result = lmfit.minimize(level_chi_squared, level_params, args=(image_conv_clipped.flatten(),))
        level = level_result.params['level'].value
        if i==0:
            level_ref = level
            continue
        average_offset = level_ref - level
        #print 'Applying offset of '+str(average_offset)+' to '+fitsfile_list[i]

        """
        # Save floor and peak values
        floor_value = np.nanmin(image_conv)
        peak_value = chrisfuncs.SigmaClip( image_conv, tolerance=0.00025, median=False, sigma_thresh=3.0)[1]
        floor_value_list.append(floor_value)
        peak_value_list.append(peak_value)
        if i==0:
            floor_value_ref = floor_value
            peak_value_ref = peak_value
            continue

        # Calculate offsets
        floor_offset = floor_value_ref - floor_value
        peak_offset = peak_value_ref - peak_value
        average_offset = peak_offset#np.mean([ floor_offset, peak_offset ])
        """

        # Read in unconvolved file, and apply offset
        fitsdata_in = fits.open(fitsfile_dir+'/'+fitsfile_list[i])
        image_in = fitsdata_in[0].data
        header_in = fitsdata_in[0].header
        fitsdata_in.close()
        image_out = image_in + average_offset
        #print 'Map mean of '+fitsfile_list[i]+' changed from '+str(np.nanmean(image_in))+' to '+str(np.nanmean(image_out))

        # Save corrected file
        image_out_hdu = fits.PrimaryHDU(data=image_out, header=header_in)
        image_out_hdulist = fits.HDUList([image_out_hdu])
        image_out_hdulist.writeto(fitsfile_dir+'/'+fitsfile_list[i], clobber=True)

# -----------------------------------------------------------------

def mosaic_galex(galaxy_name, ra, dec, width, band_dict, working_path, temp_path, meta_path, output_path, nprocesses=6):

    """
    Function to SWarp together GALEX tiles of a given source
    :param galaxy_name:
    :param ra:
    :param dec:
    :param width:
    :param band_dict:
    :param working_path:
    :param temp_path:
    :param meta_path:
    :param output_path:
    :param nprocesses: number of parallel processes
    :return:
    """

    # Inform the user
    log.info("Starting the 'mosaic_galex' function for: ")
    log.info(" - galaxy name = " + galaxy_name)
    log.info(" - ra = " + str(ra))
    log.info(" - dec = " + str(dec))
    log.info(" - width = " + str(width))
    log.info(" - band_dict = " + str(band_dict))
    log.info(" - working_path = " + working_path)
    log.info(" - temp_path = " + temp_path)
    log.info(" - meta_path = " + meta_path)
    log.info(" - output_path = " + output_path)

    ra_deg = ra.to("deg").value
    dec_deg = dec.to("deg").value
    width_deg = width.to("deg").value

    # Declare directories
    id_string = galaxy_name + '_GALEX_' + band_dict['band_long']

    # Temporary directory
    temp_path_band = fs.join(temp_path, band_dict["band_long"])

    # Raw directory in temporary directory
    temp_raw_path = fs.join(temp_path_band, "raw")
    fs.create_directory(temp_raw_path)

    # Diffs_temp directory in temporary directory
    temp_diffs_path = fs.join(temp_path_band, "diffs")
    fs.create_directory(temp_diffs_path)

    # Backsub_temp directory in temporary directory
    temp_backsub_path = fs.join(temp_path_band, "backsub")
    fs.create_directory(temp_backsub_path)

    # SWarp temp directory in temporary directory
    temp_swarp_path = fs.join(temp_path_band, "swarp")
    fs.create_directory(temp_swarp_path)

    # Reproject temp directory in temporary directory
    temp_reproject_path = fs.join(temp_path_band, "reproject")
    fs.create_directory(temp_reproject_path)

    # Convolve temp directory in temporary directory
    temp_convolve_path = fs.join(temp_path_band, "convolve")
    fs.create_directory(temp_convolve_path)

    ### CHANGE WORKING DIRECTORY TO RAW
    os.chdir(temp_raw_path)
    ###

    # Path to the overlap table
    overlap_path = fs.join(temp_path_band, "overlap_table.dat")


    # FILTER GALEX tiles for this band
    raw_files = filter_galex_tiles(galaxy_name, meta_path, overlap_path, ra, dec, width_deg, temp_raw_path, band_dict)


    # CLEAN GALEX TILES
    clean_galex_tiles(raw_files, id_string, nprocesses, working_path, temp_path_band, temp_reproject_path, band_dict)


    # Create Montage FITS header
    location_string = str(ra_deg) + ' ' + str(dec_deg)
    pix_size = 3.2

    # Make header
    header_path = fs.join(temp_path_band, id_string + ".hdr")
    montage.commands.mHdr(location_string, width_deg, header_path, pix_size=pix_size)

    # Count image files, and move to reprojection directory
    mosaic_count = 0
    for filename in os.listdir(temp_raw_path):
        if '.fits' in filename:
            mosaic_count += 1
            new_path = fs.join(temp_reproject_path, filename)  ## NEW
            if fs.is_file(new_path): fs.remove_file(new_path)  ## NEW
            shutil.move(filename, temp_reproject_path)

    # If more than one image file, commence background-matching
    if mosaic_count > 1:

        # Inform the user
        log.info("Matching background of " + id_string + " maps ...")

        # ...
        level_galex_maps(temp_reproject_path, temp_convolve_path, 'int.fits')

    # The path to the metadata table
    metatable_path = meta_path

    # Create a dictionary for the exposure times for each image
    exposure_times = dict()

    # Determine how the files are named
    filename_ends = "-" + band_dict['band_short'] + "-int"   # -fd-int for FUV, -nd-int for NUV
    filename_ends_no_int = "-" + band_dict['band_short'] # -fd for FUV, -nd for NUV

    # Get the exposure time for each image
    for path, name in fs.files_in_path(temp_reproject_path, extension="fits", contains=filename_ends, returns=["path", "name"]):

        # Open the header
        header = fits.getheader(path)

        # Search for the exposure time
        exp_time = get_total_exposure_time(header)

        # Set the exposure time
        image_name = name.split(filename_ends)[0]
        exposure_times[image_name] = exp_time

    # Reproject image and weight prior to coaddition
    montage.commands.mImgtbl(temp_reproject_path, metatable_path, corners=True)
    proj_stats_path = fs.join(temp_path_band, id_string + "_Proj_Stats.txt")
    montage.commands.mProjExec(metatable_path, header_path, temp_swarp_path, proj_stats_path, raw_dir=temp_reproject_path, debug=False, exact=True, whole=False)
    # WHOLE IS IMPORTANT HERE
    # WE ACTUALLY DON'T WANT TO REPROJECT TO THE EXACT PIXELGRID DEFINED BY THE HEADER HERE,
    # BUT RATHER CUT OFF THE ORIGINAL MAPS WHERE THEY ARE OUTSIDE OF THE FIELD OF VIEW OF THE HEADER DEFINED AREA
    # SO NO ADDING OF NEW PIXELS IS DONE TO HAVE THE EXACT PIXELGRID DEFINED BY THE HEADER
    # THIS IS PRESUMABLY DONE JUST TO MAKE THE SWARPING MORE EFFICIENT ??

    # Rename reprojected files for SWarp
    for listfile in os.listdir(temp_swarp_path):
        if '_area.fits' in listfile:
            os.remove(fs.join(temp_swarp_path, listfile))
        elif 'hdu0_' in listfile:
            os.rename(fs.join(temp_swarp_path, listfile), fs.join(temp_swarp_path, listfile.replace('hdu0_','')))

    # MOSAIC WITH SWARP
    swarp_result_path = mosaic_with_swarp(id_string, width_deg, pix_size, temp_swarp_path, ra_deg, dec_deg)



    # CONVERT BACK TO JUST COUNTS/S
    #convert_mosaic_and_error_map_to_ct_per_s(id_string, temp_mosaic_path, output_path)

    #######################

    #### LOAD SWARP RESULT, PRESUMABLY IN COUNTS/S

    # Load the resulting frame
    out_image = Frame.from_file(swarp_result_path)
    out_image.unit = "count/s"

    out_image[out_image == 0] = np.NaN
    out_image[out_image < -1E3] = np.NaN
    out_image[out_image <= 1E-8] = 0

    # Write the swarped image
    swarp_output_path = fs.join(output_path, id_string + "_swarp.fits")
    out_image.saveto(swarp_output_path)

    ### WRITE THE MOSAIC, ERROR MAP AND RELATIVE ERROR MAP IN COUNTS/S

    # SOME THINGS
    mosaic[mosaic.data == 0] = np.NaN
    mosaic[mosaic.data < -1E3] = np.NaN
    mosaic[mosaic.data <= 1E-8] = 0

    # CALCULATE RELATIVE POISSON ERROR MAP
    relerrors = errors / mosaic
    relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
    relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

    # Save mosaic as FITS file
    mosaic_output_path = fs.join(output_path, id_string + ".fits")
    mosaic.saveto(mosaic_output_path)

    # Save error map as FITS file
    errors_output_path = fs.join(output_path, id_string + "_errors.fits")
    errors.saveto(errors_output_path)

    # Save relative error map as FITS file
    relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
    relerrors.saveto(relerrors_output_path)
    
# -----------------------------------------------------------------



# -----------------------------------------------------------------



# -----------------------------------------------------------------

def old_end_of_mosaic_galex_function():

    """
    This function ...
    :return:
    """

    # Poisson temp directory in temporary directory
    temp_poisson_path = fs.join(temp_path_band, "poisson2")
    fs.create_directory(temp_poisson_path)

    # Some directories
    temp_poisson_count_path = fs.create_directory_in(temp_poisson_path, "count")
    temp_poisson_countsr_path = fs.create_directory_in(temp_poisson_path, "countsr")
    temp_poisson_rebin_path = fs.create_directory_in(temp_poisson_path, "rebin-count")
    temp_poisson_footprint_path = fs.create_directory_in(temp_poisson_path, "footprint")
    temp_poisson_weights_path = fs.create_directory_in(temp_poisson_path, "weights")
    temp_poisson_result_path = fs.create_directory_in(temp_poisson_path, "result")

    # CONVERT TO JANSKY / PIX

    # FROM COUNT / S TO AB MAG:
    # mag_AB = ZP - (2.5 * log10(CpS))
    # FROM AB MAG TO FLUX (JANSKY):
    # mag_AB = -2.5 log (Fv / 3631 Jy) => Fv[Jy] = ...

    # Calculate the conversion factor
    conversion_factor = 1.0
    conversion_factor *= ab_mag_zero_point.to("Jy").value

    if band_dict['band_long'] == "FUV":
        conversion_factor *= 10. ** (galex_fuv_zero_point / 2.5)
    elif band_dict['band_long'] == "NUV":
        conversion_factor *= 10. ** (galex_nuv_zero_point / 2.5)
    else: raise ValueError("Invalid band name: " + band_dict['band_long'])

    # DO THE CONVERSION

    # Convert and set the new unit
    out_image *= conversion_factor
    out_image.unit = "Jy/pix"

    #######################

    #### NO, BECAUSE SWARP DECIDES ITS OWN COORDINATE SYSTEM BASED ON WHAT IT HAS AS INPUT IMAGES, SO
    #### THE SWARP MOSAIC MAY NOT CORRESPOND EXACTLY TO THE TARGET HEADER OR WCS THAT WE CREATED TO
    # Load header
    #rebin_header = Header.fromtextfile(header_path)
    # To coordinate system
    #rebin_wcs = CoordinateSystem(rebin_header)
    rebin_wcs = out_image.wcs

    ## CALCULATION OF POISSON

    ## REBINNING AND CONVERSION TO COUNT

    #print(fs.files_in_path(counts_path_band))

    # Open the -int images in the temp_swarp_path that are used to make the mosaic, convert them to counts
    nswarp_images = 0
    for filename in os.listdir(temp_swarp_path):
    #for filename in os.listdir(counts_path_band):

        if not filename.endswith("-int.fits"): continue

        # Determine the path to the cleaned rebinned image inside temp_swarp_path
        # Setting the frame NAN where the corresponding cleaned image is also NAN
        cleaned_path = fs.join(temp_swarp_path, filename)

        # Determine the path to the corresponding weight map
        weight_path = cleaned_path.replace("-int", "-int.wgt")

        #print(filename_ends, filename, exposure_times.keys())

        # Get the image name
        image_name = filename.split(filename_ends)[0]

        # Increment the counter
        nswarp_images += 1

        # Determine filepath
        #filepath = fs.join(temp_swarp_path, filename)

        filepath = fs.join(counts_path_band, image_name + "-" + band_dict["band_short"] + "-cnt.fits")

        # Debugging
        log.debug("Loading the " + image_name + " frame ...")

        # Load the frame
        frame = Frame.from_file(filepath)

        # Debugging
        log.debug("Converting unit to count / sr ...")

        # Get the exposure time for this image
        exposure_time = exposure_times[image_name]

        # Debugging
        log.debug("The exposure time for this image is " + str(exposure_time))

        # Convert the frame FROM COUNT/S to COUNT
        #frame *= exposure_times[image_name]
        frame.unit = "count" # set the unit to count

        # Save the frame to the count path
        frame.saveto(fs.join(temp_poisson_count_path, image_name + ".fits"))

        # CONVERT THE FRAME FROM COUNT TO COUNT/SR
        frame /= frame.pixelarea.to("sr").value
        frame.unit = "count/sr"

        # Save the frame to the countsr path
        frame.saveto(fs.join(temp_poisson_countsr_path, image_name + ".fits"))

        # REBIN THE FRAME TO THE COMMON PIXEL GRID
        footprint = frame.rebin(rebin_wcs)

        # CONVERT THE FRAME FROM COUNT/SR TO COUNT
        frame *= frame.pixelarea.to("sr").value
        frame.unit = "count"

        # SET NANS IN THE REBINNED FRAME WHERE THE CLEANED IMAGE IN THE TEMP SWARP PATH IS ALSO NAN
        ## FIRST REBIN THE CLEANED IMAGE ALSO TO THE FINAL WCS AS DETERMINED BY SWARP
        #cleaned_frame = Frame.from_file(cleaned_path)
        weight_map = Frame.from_file(weight_path)
        weight_map.rebin(rebin_wcs)
        frame[weight_map.nans()] = np.NaN

        # Save the rebinned frame
        frame.saveto(fs.join(temp_poisson_rebin_path, image_name + ".fits"))

        # SET NANS IN THE FOOTPRINT WHERE THE WEIGHT MAP IN THE TEMP SWARP PATH IS ALSO NAN
        footprint[weight_map.nans()] = np.NaN

        # Save the (rebinned) footprint
        footprint.saveto(fs.join(temp_poisson_footprint_path, image_name + ".fits"))

    # Initialize a list to contain the frames to be summed
    a_frames = [] # the rebinned maps in counts
    ab_frames = []
    weight_frames = []

    # Loop over the files in the temp poisson rebin directory
    for path, name in fs.files_in_path(temp_poisson_rebin_path, extension="fits", returns=["path", "name"]):

        # Open the rebinned frame IN COUNTS
        a = Frame.from_file(path)

        # Set NaNs to zero
        a.replace_nans(0.0) # if we don't do this the entire combined poisson error frame is NaN

        # Get footprint
        footprint_path = fs.join(temp_poisson_footprint_path, name + ".fits")
        b = Frame.from_file(footprint_path)

        # Set NaNs to zero
        b.replace_nans(0.0)

        # Add product of primary and footprint and footprint to the appropriate list
        ab = a * b
        a_frames.append(a)
        ab_frames.append(ab)

        # Calculate weight frame
        weight_frame = b * math.sqrt(exposure_times[name])
        weight_frames.append(weight_frame)

        # Save weight frame
        weight_frame.saveto(fs.join(temp_poisson_weights_path, name + ".fits"))

    # Take the sums
    ab_sum = sum_frames(*ab_frames)

    # AB SUM SHOULD ACTUALLY BE THE SAME AS A SUM
    a_sum = sum_frames(*a_frames) # SUM ALL THE COUNT MAPS
    # Save the total count map (a_sum)
    a_sum.saveto(fs.join(temp_poisson_result_path, "total_counts.fits"))

    # Calculate the relative poisson errors
    rel_poisson_frame = ab_sum ** (-0.5)

    # Calculate the total weight map
    total_weight_map = sum_frames(*weight_frames)

    # Save rel poisson frame and total weight map
    rel_poisson_frame.saveto(fs.join(temp_poisson_result_path, "rel_poisson.fits"))
    total_weight_map.saveto(fs.join(temp_poisson_result_path, "weights.fits"))

    # Write the Poisson error frame also to the output directory
    #poisson_path = fs.join(output_path, id_string + "_relpoisson.fits")
    #rel_poisson_frame.saveto(poisson_path)

    ################ WRITE RESULT

    # Determine output image path
    out_image_path = fs.join(output_path, id_string + ".fits")

    image = Image()
    image.add_frame(out_image, "primary") # the mosaic
    image.add_frame(rel_poisson_frame, "rel_poisson") # has no unit, but Image will be saved with unit. Problem?

    # Save the image
    image.saveto(out_image_path)

    ################

    # Clean up
    log.success("Completed creating the mosaic and poisson noise map for " + id_string)
    #gc.collect()
    #shutil.rmtree(temp_dir)

# -----------------------------------------------------------------

def make_noise_maps_in_cps(band_dict, image_names_for_mosaic, counts_path_band, temp_noise_path, exposure_times):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param counts_path_band:
    :param temp_noise_path:
    :param exposure_times:
    :return:
    """

    # Inform the user
    log.info("Creating maps of the poisson noise for each GALEX tile in counts per second ...")

    # Loop over the file names
    for image_name in image_names_for_mosaic:

        # Get the path to the file in counts
        counts_filepath = fs.join(counts_path_band, image_name + "-" + band_dict["band_short"] + "-cnt.fits")

        # Load the counts map
        counts_frame = Frame.from_file(counts_filepath)
        counts_frame.unit = "ct"

        # Calculate the poisson frame
        poisson = Frame(np.sqrt(counts_frame.data))   # calculate the poisson error (in counts) in every pixel
        poisson.wcs = counts_frame.wcs
        poisson.unit = "ct"

        # Get the exposure time in seconds
        exposure_time = exposure_times[image_name]

        # Convert the poisson frame to counts/second
        poisson /= exposure_time
        poisson.unit = "ct/s"

        # Determine the path to the poisson noise map in counts per second
        new_path = fs.join(temp_noise_path, image_name + ".fits")

        # Save the poisson noise map in counts per second
        poisson.saveto(new_path)

# -----------------------------------------------------------------

def convert_frames_and_error_maps_to_per_solid_angle(band_dict, image_names_for_mosaic, temp_reproject_path,
                                                     temp_noise_path, temp_converted_path):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param temp_reproject_path:
    :param temp_noise_path:
    :param temp_converted_path:
    :return:
    """

    # Inform the user
    log.info("Converting the units to luminosity per solid angle ...")

    # Determine how the files are named
    filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
    filename_ends_no_int = "-" + band_dict['band_short']    # -fd for FUV, -nd for NUV

    # DO THE UNIT CONVERSION
    for image_name in image_names_for_mosaic:

        # Debugging
        log.debug("Converting " + image_name + " frame and error map to count / s / sr...")

        # Determine the full path to the tile in temp_reproject_path
        filepath = fs.join(temp_reproject_path, image_name + filename_ends + ".fits")

        # Load the frame
        frame = Frame.from_file(filepath)
        frame.unit = "count/s"

        # NUMBER OF SR PER PIXEL
        pixelsr = frame.pixelarea.to("sr").value

        # CONVERT FRAME
        frame /= pixelsr  # IN COUNTS PER SECOND PER SR NOW
        frame.unit = "count/(s*sr)" # set new unit

        # Determine the full path to the poisson noise map for this tile
        poisson_filepath = fs.join(temp_noise_path, image_name + ".fits")

        # Load the poisson noise frame
        poisson = Frame.from_file(poisson_filepath)
        # the unit should be set already to count/s

        # normally, this shouldn't be necessary
        poisson_pixelsr = poisson.pixelarea.to("sr").value
        assert np.isclose(pixelsr, poisson_pixelsr)

        # CONVERT ERROR MAP TO COUNTS PER SECOND PER SR
        poisson /= pixelsr  # IN COUNTS PER SECOND PER SR NOW
        poisson.unit = "count/(s*sr)" # set new unit

        # SAVE THE TILE IN COUNTS/S/SR
        frame_path = fs.join(temp_converted_path, image_name + ".fits")
        frame.saveto(frame_path)

        # SAVE THE ERROR MAP OF THE TILE IN COUNTS/S/SR
        errors_path = fs.join(temp_converted_path, image_name + "_error.fits")
        poisson.saveto(errors_path)

# -----------------------------------------------------------------

def rebin_frames_and_error_maps(temp_converted_path, temp_rebinned_path, header_path, metatable_path, proj_stats_path):

    """
    This function ...
    :param temp_converted_path:
    :param temp_rebinned_path:
    :param header_path:
    :param metatable_path:
    :param proj_stats_path:
    :return:
    """

    # Inform the user
    log.info("Rebinning the frames and poisson noise maps in counts / second to the target coordinate system ...")

    # Rebin the batch of images (tiles and noise maps) in counts/s/sr
    #rebin_batch_with_montage(temp_converted_path, temp_rebinned_path, metatable_path, header_path, proj_stats_path)

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    rebin_wcs = CoordinateSystem(rebin_header)

    # Loop over all FITS files
    for path, name in fs.files_in_path(temp_converted_path, extension="fits", returns=["path", "name"]):

        # Load the image
        frame = Frame.from_file(path)

        # Rebin
        frame.rebin(rebin_wcs)

        # Determine the new path
        new_path = fs.join(temp_rebinned_path, name + ".fits")

        # Save
        frame.saveto(new_path)

# -----------------------------------------------------------------

def rebin_weight_maps(band_dict, image_names_for_mosaic, temp_reproject_path, temp_rebinned_path, header_path):

    """
    This function ...
    :param band_dict:
    :param image_names_for_mosaic:
    :param temp_reproject_path:
    :param temp_rebinned_path:
    :param header_path:
    :return:
    """

    # Inform the user
    log.info("Rebinning weight maps ...")

    # Determine how the files are named
    filename_ends = "-" + band_dict['band_short'] + "-int"  # -fd-int for FUV, -nd-int for NUV
    filename_ends_no_int = "-" + band_dict['band_short']  # -fd for FUV, -nd for NUV

    # GET REBIN WCS
    rebin_header = Header.fromtextfile(header_path)
    rebin_wcs = CoordinateSystem(rebin_header)

    # Loop over the image_names
    for image_name in image_names_for_mosaic:

        # Determine the path to the weight map
        weight_path = fs.join(temp_reproject_path, image_name + filename_ends + ".wgt.fits")

        # Load the weight map
        weights = Frame.from_file(weight_path)

        # Rebin
        weights.rebin(rebin_wcs)

        # Determine the new path
        new_path = fs.join(temp_rebinned_path, image_name + "_weight.fits")

        # Save
        weights.saveto(new_path)

# -----------------------------------------------------------------

def combine_frames_and_error_maps(image_names_for_mosaic, temp_rebinned_path, temp_mosaic_path, wcs):

    """
    This function ...
    :param image_names_for_mosaic:
    :param temp_rebinned_path:
    :param temp_mosaic_path:
    :param wcs:
    :return:
    """

    primary_frames = []
    error_frames = []
    weight_frames = []

    for image_name in image_names_for_mosaic:

        # Determine paths
        filepath = fs.join(temp_rebinned_path, image_name + ".fits")
        errorpath = fs.join(temp_rebinned_path, image_name + "_error.fits")
        weightpath = fs.join(temp_rebinned_path, image_name + "_weight.fits")

        # Load the frames
        frame = Frame.from_file(filepath)
        errors = Frame.from_file(errorpath)
        weights = Frame.from_file(weightpath)

        # Calculate the weighted frame
        frame_weighted = frame * weights
        errors_weighted = errors / weights

        # Convert the units from counts/s/sr to counts/s
        pixelsr = frame_weighted.pixelarea.to("sr").value
        frame_weighted *= pixelsr
        errors_weighted *= pixelsr

        # Create mask where the weights are nans
        mask = weights.nans()

        # Set zero
        frame_weighted[mask] = 0.0
        errors_weighted[mask] = 0.0
        weights[mask] = 0.0

        # Add to the list
        primary_frames.append(frame_weighted)
        error_frames.append(errors_weighted)
        weight_frames.append(weights)

    #  Calculate denominator for weighted average (mosaicing)
    normalization = sum_frames(*weight_frames)

    # CALCULATE THE MOSAIC FRAME IN COUNTS/S
    mosaic_frame = sum_frames(*primary_frames) / normalization
    mosaic_frame.wcs = wcs

    # CALCULATE THE MOSAIC ERROR MAP IN COUNTS/S
    mosaic_errormap = sum_frames_quadratically(*error_frames) / normalization
    mosaic_errormap.wcs = wcs

    ## DONE

    # SAVE THE MOSAIC IN COUNTS/S
    mosaic_path = fs.join(temp_mosaic_path, "mosaic.fits")
    mosaic_frame.saveto(mosaic_path)

    # SAVE THE MOSAIC ERROR MAP IN COUNTS/S
    mosaic_error_path = fs.join(temp_mosaic_path, "mosaic_errors.fits")
    mosaic_errormap.saveto(mosaic_error_path)

    # Return the mosaic and the mosaic error map in counts/s
    return mosaic_frame, mosaic_errormap

# -----------------------------------------------------------------

def convert_mosaic_and_error_map_to_ct_per_s(id_string, temp_mosaic_path, output_path):

    """
    This function ...
    :param id_string:
    :param temp_mosaic_path:
    :param output_path:
    :return:
    """

    # Inform the user
    log.info("Converting the mosaic and error map back to counts per seconds ...")

    mosaic_path = fs.join(temp_mosaic_path, "mosaic.fits")
    mosaic_error_path = fs.join(temp_mosaic_path, "mosaic_errors.fits")

    mosaic = Frame.from_file(mosaic_path)
    errors = Frame.from_file(mosaic_error_path)

    # SOME THINGS
    mosaic[mosaic.data == 0] = np.NaN
    mosaic[mosaic.data < -1E3] = np.NaN
    mosaic[mosaic.data <= 1E-8] = 0

    # NUMBER OF SR PER PIXEL
    pixelsr = mosaic.pixelarea.to("sr").value

    # CONVERT TO COUNTS/S
    mosaic *= pixelsr
    mosaic.unit = "count/s"

    # CONVERT TO COUNTS/S
    errors *= pixelsr
    errors.unit = "count/s"


    # CALCULATE RELATIVE POISSON ERROR MAP
    relerrors = errors / mosaic
    relerrors[relerrors < 0.] = 0.0  # set negative values for relative error map to zero
    relerrors.replace_nans(0.0)  # set NaN values (because mosaic was zero) to zero

    ### SAVE

    # Save mosaic as FITS file
    mosaic_output_path = fs.join(output_path, id_string + ".fits")
    mosaic.saveto(mosaic_output_path)

    # Save error map as FITS file
    errors_output_path = fs.join(output_path, id_string + "_errors.fits")
    errors.saveto(errors_output_path)

    # Save relative error map as FITS file
    relerrors_output_path = fs.join(output_path, id_string + "_relerrors.fits")
    relerrors.saveto(relerrors_output_path)

    ###

# -----------------------------------------------------------------

def mosaic_with_swarp(id_string, width_deg, pix_size, temp_swarp_path, ra_deg, dec_deg):

    """
    This function ...
    :return:
    """

    # Use SWarp to co-add images weighted by their error maps
    log.info("Co-adding " + id_string + " maps ...")

    # Determine image width in pixels
    image_width_pixels = str(int((float(width_deg) * 3600.) / pix_size))

    # Change directories
    os.chdir(temp_swarp_path)

    # EXECUTE SWARP
    swarp_command_string = 'swarp *int.fits -IMAGEOUT_NAME ' + id_string + '_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER ' + str(ra_deg) + ',' + str(dec_deg) + ' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE ' + image_width_pixels + ',' + image_width_pixels + ' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT'
    os.system(swarp_command_string)

    # Swarp result path
    swarp_result_path = fs.join(temp_swarp_path, id_string + "_SWarp.fits")

    # Return the path to the resulting mosaic
    return swarp_result_path

# -----------------------------------------------------------------

def get_total_exposure_time(header):

    """
    This function ....
    :param header:
    :return:
    """

    # If the 'EXPTIME' keyword is defined in the header, return the value
    if "EXPTIME" in header: return float(header["EXPTIME"])

    # Otherwise, look at individual 'EXPTIME' keywords
    else:

        total = 0.0
        for key in header:

            # Add the values of all keys that state an exposure time (EXPT0001, EXPT0002, ...)
            if key.startswith("EXPT"): total += float(header[key])

        # Return the total exposure time
        return total

# -----------------------------------------------------------------

def rebin_batch_with_montage(in_path, out_path, metatable_path, header_path, proj_stats_path, remove_area=True):

    """
    This function ...
    :param in_path:
    :param out_path:
    :param metatable_path:
    :param header_path:
    :param proj_stats_path:
    :param remove_area:
    :return:
    """

    # Inform the user
    log.info("Rebinning a batch of images in " + in_path + " to " + out_path + " ...")

    # REPROJECT ALL INPUT MAPS TO TARGET HEADER
    montage.commands.mProjExec(metatable_path, header_path, out_path, proj_stats_path, raw_dir=in_path, debug=False, exact=True, whole=True)

    # Rename reprojected files
    for listfile in os.listdir(out_path):

        if '_area.fits' in listfile:
            if remove_area: os.remove(fs.join(out_path, listfile)) # remove area if requested
        elif 'hdu0_' in listfile:
            os.rename(fs.join(out_path, listfile), fs.join(out_path, listfile.replace('hdu0_', '')))

# -----------------------------------------------------------------
