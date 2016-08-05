#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.galex_montage_functions Functions for ...

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
import gc

# Import astronomical modules
from astropy.io import fits
import astropy.io.votable
import astropy.convolution
import montage_wrapper as montage
from astropy.wcs import WCS
from astropy.io.fits import Header
from astropy.units import Unit

# Import Chris' package
import ChrisFuncs

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...magic.core.frame import Frame, sum_frames
from ...magic.core.image import Image
from ...magic.basics.coordinatesystem import CoordinateSystem

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

        image_conv_clipped = ChrisFuncs.SigmaClip(image_conv, tolerance=0.005, median=False, sigma_thresh=3.0)[2]
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
        peak_value = ChrisFuncs.SigmaClip( image_conv, tolerance=0.00025, median=False, sigma_thresh=3.0)[1]
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

def clean_galex_tile(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict):

    """
    Function to clean GALEX tiles and create exposure maps
    :return:
    """

    # Inform the user ...
    log.info("Cleaning map " + raw_file + " ...")
    log.info(" - raw_file = " + raw_file)
    log.info(" - working_path = " + working_path)
    log.info(" - temp_path_band = " + temp_path_band)
    log.info(" - temp_reproject_path = " + temp_reproject_path)
    log.info(" - band_dict = " + str(band_dict))

    # Response and background paths for this band
    response_path = fs.join(working_path, "response", band_dict['band_long'])
    background_path = fs.join(working_path, "background", band_dict['band_long'])

    temp_raw_path = fs.join(temp_path_band, "raw")

    # Read in image
    in_fitsdata = fits.open(fs.join(temp_raw_path, raw_file))
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()

    # Load and align response map
    rr_path = fs.join(response_path, raw_file.replace('-int.fits','-rr.fits'))
    rr_fitsdata = fits.open(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image <= 1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = fs.join(background_path, raw_file.replace('-int.fits','-skybg.fits'))
    bg_fitsdata = fits.open(bg_path)
    bg_image = bg_fitsdata[0].data
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = scipy.ndimage.interpolation.zoom(bg_image, bg_zoom, order=0)

    # Clean image using sky background map
    out_image[ np.where( bg_image <= 1E-10 ) ] = np.NaN

    """
    # Load and align flag map
    flag_path = root_dir+'Flags/'+band_dict['band_long']+'/'+raw_file.replace('-int.fits','-flags.fits.gz')
    flag_fitsdata = fits.open(flag_path)
    flag_image = flag_fitsdata[0].data
    flag_zoom = np.float(out_image.shape[0]) / np.float(flag_image.shape[0])
    flag_image = scipy.ndimage.interpolation.zoom(flag_image, flag_zoom, order=0)

    # Nullify pixels where the bitwise flags indicate dichoric reflections, window refections, variable-pixel masks or hot-pixel masks
    out_image[ np.where( flag_image.astype(int) & (1<<1) > 0 ) ] = np.NaN # Dichoric reflection
    out_image[ np.where( flag_image.astype(int) & (1<<2) > 0 ) ] = np.NaN # Window reflection
    out_image[ np.where( flag_image.astype(int) & (1<<7) > 0 ) ] = np.NaN # Variable-pixel mask
    out_image[ np.where( flag_image.astype(int) & (1<<8) > 0 ) ] = np.NaN # Hot-pixel mask
    """

    # Set all remaining, and hence "true", zero pixels to be ever-so-slighly non-zero
    out_image += 1E-8

    # Find centre of coverage area
    cov_i = ((np.where( np.isnan(out_image)==False ))[0])
    cov_j = ((np.where( np.isnan(out_image)==False ))[1])
    cov_ellipse = ChrisFuncs.EllipseFit(cov_i, cov_j)
    cov_centre = cov_ellipse[0]
    cov_centre_i, cov_centre_j = cov_centre[0], cov_centre[1]

    # Set all pixels more than 35 arcmin (1400 pizels) from centre to be NaN, as these are typically low-quality
    cov_trim_mask = ChrisFuncs.EllipseMask(out_image, 1400, 1.0, 0.0, cov_centre_i, cov_centre_j)
    out_image[ np.where(cov_trim_mask==0) ] = np.NaN

    # Save cleaned image
    out_hdu = fits.PrimaryHDU(data=out_image, header=in_header)
    out_hdulist = fits.HDUList([out_hdu])
    out_hdulist.writeto(fs.join(temp_raw_path, raw_file), clobber=True)

    # Create convolved version of map, for later use in background-matching
    """
    if np.isnan(out_image).sum()==0:
        conv_image = scipy.ndimage.filters.gaussian_filter(out_image, 20)
    else:
    """

    temp_convolve_path = fs.join(temp_path_band, "convolve")

    kernel = astropy.convolution.kernels.Tophat2DKernel(10)
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True) #, interpolate_nan=True, normalize_kernel=True)

    # Write
    temp_convolve_image_path = fs.join(temp_convolve_path, raw_file)                  ## NEW
    if fs.is_file(temp_convolve_image_path): fs.remove_file(temp_convolve_image_path) ## NEW
    fits.writeto(temp_convolve_image_path, conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5 # SQUARE ROOT OF THE EXPOSURE TIME
    exp_hdu = fits.PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = fits.HDUList([exp_hdu])

    # Write
    temp_reproject_image_path = fs.join(temp_reproject_path, raw_file.replace('.fits','.wgt.fits'))  ## NEW
    if fs.is_file(temp_reproject_image_path): fs.remove_file(temp_reproject_image_path)              ## NEW
    exp_hdulist.writeto(temp_reproject_image_path)

# -----------------------------------------------------------------

def mosaic_galex(name, ra, dec, width, band_dict, working_path, temp_path, meta_path, output_path, nprocesses=6):

    """
    Function to SWarp together GALEX tiles of a given source
    :param name:
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
    log.info(" - name = " + name)
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
    id_string = name + '_GALEX_' + band_dict['band_long']

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

    # Poisson temp directory in temporary directory
    temp_poisson_path = fs.join(temp_path_band, "poisson")
    fs.create_directory(temp_poisson_path)

    # Create storage directories for Montage and SWarp (deleting any prior), and set appropriate Python working directory
    #os.mkdir(temp_dir + 'Raw')
    #os.mkdir(temp_dir + 'Diffs_Temp')
    #os.mkdir(temp_dir + 'Backsub_Temp')
    #os.mkdir(temp_dir + 'SWarp_Temp')
    #os.mkdir(temp_dir + 'Reproject_Temp')
    #os.mkdir(temp_dir + 'Convolve_Temp')
    #os.chdir(temp_dir + 'Raw')

    ### CHANGE WORKING DIRECTORY TO RAW

    os.chdir(temp_raw_path)

    ###

    # Path to the overlap table
    overlap_path = fs.join(temp_path_band, "overlap_table.dat")

    # Use Montage image metadata table to identify and retrieve which raw GALEX tiles overlap with entire region of interest (handling the case of only a single file)
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='circle', ra=ra.to("deg").value, dec=dec.to("deg").value, radius=(0.5*width_deg)*(2.0**0.5))

    # Get file paths of overlapping observations
    overlapping_file_paths = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype=str)

    if len(overlapping_file_paths.shape) == 0:
        overlapping_file_paths = [overlapping_file_paths.tolist()]
    for overlapping_file_path in overlapping_file_paths:
        shutil.copy(overlapping_file_path, temp_raw_path)

    # Uncompress .fits.gz files
    #[os.system('gunzip '+ listfile) for listfile in os.listdir(raw_in_temp_dir)]

    # Ensure that at least one of the raw GALEX tiles has actual flux coverage at location of source
    # raw_files = os.listdir(temp_raw_path) ## THERE WAS A 'deg' FILE IN THE DIRECTORY AS WELL WHICH WAS NOT A FITS FILE BUT A PLAIN TEXT HEADER FILE, AND SO THIS WASN'T WORKING
    raw_files = fs.files_in_path(temp_raw_path, extension="fits", returns="name")
    raw_files = [name + ".fits" for name in raw_files]
    coverage = False
    for raw_file in raw_files:

        # Read in map
        in_fitsdata = fits.open(fs.join(temp_raw_path, raw_file))
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()

        # Locate pixel coords
        in_wcs = WCS(in_header)
        location_pix = in_wcs.wcs_world2pix( np.array([[ np.float(ra.to("deg").value), np.float(dec.to("deg").value) ]]), 0 )[0]
        pix_i, pix_j = location_pix[1], location_pix[0]

        # Evalulate coverage at location, and proceed accordingly
        if True in [ coord <= 0 for coord in [ pix_i-10, pix_i+11, pix_j-10, pix_j+11 ] ]:
            continue
        try:
            image_slice = in_image[pix_i-10:pix_i+11, pix_j-10:pix_j+11]
        except:
            continue
        if np.where(image_slice>0)[0].shape[0]>0:
            coverage = True

    if not coverage:

        # Warning message
        log.warning('No GALEX '+ band_dict['band_long'] + ' coverage for ' + name)
        gc.collect()
        shutil.rmtree(temp_path_band)

    elif coverage:

        # Loop over raw tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
        log.info('Cleaning '+ str(len(raw_files)) + ' raw maps for ' + id_string + " ...")

        if nprocesses == 1:

            # CLEAN
            for raw_file in raw_files: clean_galex_tile(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict)

        else:

            # Create process pool
            pool = mp.Pool(processes=nprocesses)

            # EXECUTE THE LOOP IN PARALLEL
            for raw_file in raw_files: pool.apply_async(clean_galex_tile, args=(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict,))

            # CLOSE AND JOIN THE PROCESS POOL
            pool.close()
            pool.join()

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

        # Use SWarp to co-add images weighted by their error maps
        log.info("Co-adding " + id_string + " maps ...")
        image_width_pixels = str(int((float(width_deg)*3600.)/pix_size))
        os.chdir(temp_swarp_path)

        # EXECUTE SWARP
        swarp_command_string = 'swarp *int.fits -IMAGEOUT_NAME '+ id_string + '_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER ' + str(ra_deg) + ',' + str(dec_deg) + ' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE ' + image_width_pixels + ',' + image_width_pixels + ' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT'
        os.system(swarp_command_string)

        # Swarp result path
        swarp_result_path = fs.join(temp_swarp_path, id_string + "_SWarp.fits")

        # BEFORE:
        # Remove null values, and save finalised map to output directory
        #in_fitsdata = fits.open()
        #in_image = in_fitsdata[0].data
        #in_header = in_fitsdata[0].header
        #in_fitsdata.close()
        #out_image = in_image.copy()
        #out_image[ np.where( out_image == 0 ) ] = np.NaN
        #out_image[ np.where( out_image < -1E3 ) ] = np.NaN
        #out_image[ np.where( out_image <= 1E-8 ) ] = 0

        # NEW:
        #out_image = Frame(in_image, unit="count/s")
        out_image = Frame.from_file(swarp_result_path)
        out_image.unit = "count/s"

        out_image[out_image == 0] = np.NaN
        out_image[out_image < -1E3] = np.NaN
        out_image[out_image <= 1E-8] = 0

        # CONVERT TO JANSKY / PIX

        # FROM COUNT / S TO AB MAG:
        # mag_AB = ZP - (2.5 * log10(CpS))
        # FROM AB MAG TO FLUX (JANSKY):
        # mag_AB = -2.5 log (Fv / 3631 Jy) => Fv[Jy] = ...

        # Calculate the conversion factor
        conversion_factor = 1.0
        conversion_factor *= ab_mag_zero_point.to("Jy").value

        if band_dict['band_long'] == "FUV": conversion_factor *= 10.**(galex_fuv_zero_point/2.5)
        elif band_dict['band_long'] == "NUV": conversion_factor *= 10.**(galex_nuv_zero_point/2.5)
        else: raise ValueError("Invalid band name: " + band_dict['band_long'])

        # DO THE CONVERSION

        # Convert and set the new unit
        out_image *= conversion_factor
        out_image.unit = "Jy/pix"

        # BEFORE:
        #out_hdu = fits.PrimaryHDU(data=out_image, header=in_header)
        #out_hdulist = fits.HDUList([out_hdu])

        #### OUTPUT MOSAIC ####

        # Write mosaic
        #mosaic_path = fs.join(output_path, id_string + '.fits')
        #out_hdulist.writeto(mosaic_path, clobber=True)

        ## WRITE THE OUTPUT FRAME

        # Determine path and write mosaic: NOW DONE AT THE END TOGETHER IN AN IMAGE WITH THE REL POISSON FRAME
        #out_image_path = fs.join(output_path, id_string + ".fits")
        #out_image.save(out_image_path)

        #######################

        # Some directories
        temp_poisson_count_path = fs.create_directory_in(temp_poisson_path, "count")
        temp_poisson_countsr_path = fs.create_directory_in(temp_poisson_path, "countsr")
        temp_poisson_rebin_path = fs.create_directory_in(temp_poisson_path, "rebin")
        temp_poisson_footprint_path = fs.create_directory_in(temp_poisson_path, "footprint")
        temp_poisson_weights_path = fs.create_directory_in(temp_poisson_path, "weights")
        temp_poisson_result_path = fs.create_directory_in(temp_poisson_path, "result")

        # Load header
        rebin_header = Header.fromtextfile(header_path)

        # To coordinate system
        rebin_wcs = CoordinateSystem(rebin_header)

        ## CALCULATION OF POISSON

        ## REBINNING AND CONVERSION TO COUNT

        # Open the -int images in the temp_swarp_path that are used to make the mosaic, convert them to counts
        nswarp_images = 0
        for filename in os.listdir(temp_swarp_path):

            if not filename.endswith("-int.fits"): continue

            #print(filename_ends, filename, exposure_times.keys())

            # Get the image name
            image_name = filename.split(filename_ends)[0]

            # Increment the counter
            nswarp_images += 1

            # Determine filepath
            filepath = fs.join(temp_swarp_path, filename)

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
            frame *= exposure_times[image_name]
            frame.unit = "count" # set the unit to count

            # Save the frame to the count path
            frame.save(fs.join(temp_poisson_count_path, image_name + ".fits"))

            # CONVERT THE FRAME FROM COUNT TO COUNT/SR
            frame /= frame.pixelarea.to("sr").value
            frame.unit = "count/sr"

            # Save the frame to the countsr path
            frame.save(fs.join(temp_poisson_countsr_path, image_name + ".fits"))

            # REBIN THE FRAME TO THE COMMON PIXEL GRID
            footprint = frame.rebin(rebin_wcs)

            # CONVERT THE FRAME FROM COUNT/SR TO COUNT
            frame *= frame.pixelarea.to("sr").value
            frame.unit = "count"

            # Save the rebinned frame
            frame.save(fs.join(temp_poisson_rebin_path, image_name + ".fits"))

            # Save the (rebinned) footprint
            footprint.save(fs.join(temp_poisson_footprint_path, image_name + ".fits"))

        # Initialize a list to contain the frames to be summed
        ab_frames = []
        b_frames = []
        weight_frames = []

        # Loop over the files in the temp poisson rebin directory
        for path, name in fs.files_in_path(temp_poisson_rebin_path, extension="fits", returns=["path", "name"]):

            # Open the rebinned frame
            a = Frame.from_file(path)

            # Set NaNs to zero
            a.replace_nans(0.0) # if we don't do this the entire combined poisson error frame is NaN

            # Get footprint
            footprint_path = fs.join(temp_poisson_footprint_path, name + ".fits")
            b = Frame.from_file(footprint_path)

            # Add product of primary and footprint and footprint to the appropriate list
            ab = a * b
            ab_frames.append(ab)
            b_frames.append(b)

            # Calculate weight frame
            weight_frame = b * math.sqrt(exposure_times[name])
            weight_frames.append(weight_frame)

            # Save weight frame
            weight_frame.save(fs.join(temp_poisson_weights_path, name + ".fits"))

        # Take the sums
        ab_sum = sum_frames(*ab_frames)
        b_sum = sum_frames(*b_frames)

        # Calculate the relative poisson errors
        rel_poisson_frame = ab_sum ** (-0.5)

        # Calculate the total weight map
        total_weight_map = sum_frames(*weight_frames)

        # Save rel poisson frame and total weight map
        rel_poisson_frame.save(fs.join(temp_poisson_result_path, "rel_poisson.fits"))
        total_weight_map.save(fs.join(temp_poisson_result_path, "weights.fits"))

        # Write the Poisson error frame also to the output directory
        #poisson_path = fs.join(output_path, id_string + "_relpoisson.fits")
        #rel_poisson_frame.save(poisson_path)

        ################ WRITE RESULT

        # Determine output image path
        out_image_path = fs.join(output_path, id_string + ".fits")

        image = Image()
        image.add_frame(out_image, "primary") # the mosaic
        image.add_frame(rel_poisson_frame, "rel_poisson") # has no unit, but Image will be saved with unit. Problem?

        # Save the image
        image.save(out_image_path)

        ################

        # Clean up
        log.success("Completed Montaging and SWarping of " + id_string)

        #gc.collect()
        #shutil.rmtree(temp_dir)

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
