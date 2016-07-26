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
import multiprocessing as mp
import numpy as np
import scipy.spatial
import scipy.ndimage
import matplotlib.path

# Import astronomical modules
from astropy.io import fits
import astropy.io.votable
import astropy.convolution
import montage_wrapper as montage
from astropy.wcs import WCS

import lmfit
import shutil
import gc
import time

# Import Chris' package
import ChrisFuncs

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.logging import log

# -----------------------------------------------------------------

def GALEX_Level_Chisq(level_params, image):

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

def GALEX_Zero(fitsfile_dir, convfile_dir, target_suffix):

    """
    Set a set of maps to the same level
    :param fitsfile_dir:
    :param convfile_dir:
    :param target_suffix:
    :return:
    """

    # Make list of files in target directory that have target suffix
    allfile_list = os.listdir(fitsfile_dir)
    fitsfile_list = []
    for allfile in allfile_list:
        if target_suffix in allfile:
            fitsfile_list.append(allfile)

    # Loop over each file
    for i in range(0, len(fitsfile_list)):

        log.info('Matching background of map ' + fitsfile_list[i])

        # Read in corresponding map from directory containing convolved images
        fitsdata_conv = fits.open(convfile_dir+'/'+fitsfile_list[i])
        image_conv = fitsdata_conv[0].data
        fitsdata_conv.close()

        # Fit to level of image; save if first image, otherwise calculate appropriate offset
        level_params = lmfit.Parameters()
        level_params.add('level', value=np.nanmedian(image_conv), vary=True)
        image_conv_clipped = ChrisFuncs.SigmaClip(image_conv, tolerance=0.005, median=False, sigma_thresh=3.0)[2]
        level_result = lmfit.minimize(GALEX_Level_Chisq, level_params, args=(image_conv_clipped.flatten(),))
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
#                   raw_file, working_path, temp_path_band, temp_reproject_path, band_dict

    """
    Function to clean GALEX tiles and create exposure maps
    :return:
    """

    # Inform the user ...
    print('Cleaning map ' + raw_file)

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
    rr_path = fs.join(response_path, raw_file.replace('-int.fits','-rr.fits.gz'))
    rr_fitsdata = fits.open(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image <= 1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = fs.join(background_path, raw_file.replace('-int.fits','-skybg.fits.gz'))
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
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True)#, interpolate_nan=True, normalize_kernel=True)
    fits.writeto(fs.join(temp_convolve_path, raw_file), conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5
    exp_hdu = fits.PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = fits.HDUList([exp_hdu])
    exp_hdulist.writeto(fs.join(temp_reproject_path, raw_file.replace('.fits','.wgt.fits')))

# -----------------------------------------------------------------

def mosaic_galex(name, ra, dec, width, band_dict, working_path, temp_path, meta_path, output_path):

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
    :return:
    """

    ra_deg = ra.to("deg").value
    dec_deg = dec.to("deg").value

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
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='circle', ra=ra.to("deg").value, dec=dec.to("deg").value, radius=(0.5*width)*(2.0**0.5))

    # Get file paths of overlapping observations
    overlapping_file_paths = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype=str)

    if len(overlapping_file_paths.shape) == 0:
        overlapping_file_paths = [overlapping_file_paths.tolist()]
    for overlapping_file_path in overlapping_file_paths:
        shutil.copy(overlapping_file_path, temp_raw_path)

    # Uncompress .fits.gz files
    #[os.system('gunzip '+ listfile) for listfile in os.listdir(raw_in_temp_dir)]

    # Ensure that at least one of the raw GALEX tiles has actual flux coverage at location of source
    raw_files = os.listdir(temp_raw_path)
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

        print('No GALEX '+ band_dict['band_long'] + ' coverage for ' + name)
        gc.collect()
        shutil.rmtree(temp_path_band)

    elif coverage:

        # Loop over raw tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
        print('Cleaning '+ str(len(raw_files)) + ' raw maps for ' + id_string)

        # CLEAN
        for raw_file in raw_files: clean_galex_tile(raw_file, working_path, temp_path_band, temp_reproject_path, band_dict)

        # Create Montage FITS header
        location_string = str(ra_deg) + ' ' + str(dec_deg)
        pix_size = 3.2

        header_path = fs.join(temp_path_band, id_string + '_HDR')
        montage.commands.mHdr(location_string, width, header_path, pix_size=pix_size)

        # Count image files, and move to reprojection directory
        mosaic_count = 0
        for listfile in os.listdir(temp_raw_path):
            if '.fits' in listfile:
                mosaic_count += 1
        for listfile in os.listdir(temp_raw_path):
            if '.fits' in listfile:
                shutil.move(listfile, temp_reproject_path)

        # If more than one image file, commence background-matching
        if mosaic_count > 1:
            print('Matching background of '+id_string+' maps')
            GALEX_Zero(temp_reproject_path, temp_convolve_path, 'int.fits')

        #metatable_path = temp_dir + band + '_Image_Metadata_Table.dat'

        metatable_path = meta_path

        # Reproject image and weight prior to coaddition
        montage.commands.mImgtbl(temp_reproject_path, metatable_path, corners=True)
        proj_stats_path = fs.join(temp_path_band, id_string + "_Proj_Stats.txt")
        montage.commands.mProjExec(metatable_path, header_path, temp_swarp_path, proj_stats_path, raw_dir=temp_reproject_path, debug=False, exact=True, whole=False)

        # Rename reprojected files for SWarp
        for listfile in os.listdir(temp_swarp_path):
            if '_area.fits' in listfile:
                os.remove(fs.join(temp_swarp_path, listfile))
            elif 'hdu0_' in listfile:
                os.rename(fs.join(temp_swarp_path, listfile), fs.join(temp_swarp_path, listfile.replace('hdu0_','')))

        # Use SWarp to co-add images weighted by their error maps
        print('Co-adding ' + id_string + ' maps')
        image_width_pixels = str(int((float(width)*3600.)/pix_size))
        os.chdir(temp_swarp_path)

        # EXECUTE SWARP
        swarp_command_string = 'swarp *int.fits -IMAGEOUT_NAME '+ id_string + '_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER ' + str(ra_deg) + ',' + str(dec_deg) + ' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE ' + image_width_pixels + ',' + image_width_pixels + ' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT'
        os.system(swarp_command_string)

        # Remove null values, and save finalised map to output directory
        #in_fitsdata = fits.open(temp_dir+'/SWarp_Temp/'+id_string+'_SWarp.fits')
        in_fitsdata = fits.open(fs.join(temp_swarp_path, id_string + "_SWarp.fits"))
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()
        out_image = in_image.copy()
        out_image[ np.where( out_image==0 ) ] = np.NaN
        out_image[ np.where( out_image<-1E3 ) ] = np.NaN
        out_image[ np.where( out_image<=1E-8 ) ] = 0
        out_hdu = fits.PrimaryHDU(data=out_image, header=in_header)
        out_hdulist = fits.HDUList([out_hdu])

        # Output
        output_montages_path = fs.join(output_path, "montages")
        fs.create_directory(output_montages_path)

        # Write mosaic
        out_hdulist.writeto(fs.join(output_montages_path, id_string + '.fits'), clobber=True)

        # Clean up
        print('Completed Montaging and SWarping of ' + id_string)
        #gc.collect()
        #shutil.rmtree(temp_dir)

# -----------------------------------------------------------------
