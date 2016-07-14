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

def GALEX_Clean(raw_file, root_dir, temp_dir, out_dir, band_dict):

    """
    Function to clean GALEX tiles and create exposure maps
    :param raw_file:
    :param root_dir:
    :param temp_dir:
    :param out_dir:
    :param band_dict:
    :return:
    """

    print('Cleaning map ' + raw_file)

    # Read in image
    in_fitsdata = fits.open(temp_dir+'/Raw/'+raw_file)
    in_image = in_fitsdata[0].data
    in_header = in_fitsdata[0].header
    in_fitsdata.close()
    out_image = in_image.copy()

    # Load and align response map
    rr_path = root_dir + 'Response/' + band_dict['band_long'] + '/' + raw_file.replace('-int.fits','-rr.fits.gz')
    rr_fitsdata = fits.open(rr_path)
    rr_image = rr_fitsdata[0].data
    rr_zoom = np.float(out_image.shape[0]) / np.float(rr_image.shape[0])
    rr_image = scipy.ndimage.interpolation.zoom(rr_image, rr_zoom, order=0)

    # Clean image using response map
    out_image[ np.where( rr_image<=1E-10 ) ] = np.NaN

    # Load and align sky background map
    bg_path = root_dir + 'Background/' + band_dict['band_long'] + '/' + raw_file.replace('-int.fits','-skybg.fits.gz')
    bg_fitsdata = fits.open(bg_path)
    bg_image = bg_fitsdata[0].data
    bg_zoom = np.float(out_image.shape[0]) / np.float(bg_image.shape[0])
    bg_image = scipy.ndimage.interpolation.zoom(bg_image, bg_zoom, order=0)

    # Clean image using sky background map
    out_image[ np.where( bg_image<=1E-10 ) ] = np.NaN

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
    out_hdulist.writeto(temp_dir+'/Raw/'+raw_file, clobber=True)

    # Create convolved version of map, for later use in background-matching
    """
    if np.isnan(out_image).sum()==0:
        conv_image = scipy.ndimage.filters.gaussian_filter(out_image, 20)
    else:
    """
    kernel = astropy.convolution.kernels.Tophat2DKernel(10)
    conv_image = astropy.convolution.convolve_fft(out_image, kernel, interpolate_nan=False, normalize_kernel=True, ignore_edge_zeros=False, allow_huge=True)#, interpolate_nan=True, normalize_kernel=True)
    fits.writeto(temp_dir+'Convolve_Temp/'+raw_file, conv_image, in_header)

    # Load and align exposure time to create weight maps
    exp_image = out_image.copy()
    exp_image[ np.where( np.isnan(out_image)==False ) ] = (float(in_header['exptime']))**0.5
    exp_hdu = fits.PrimaryHDU(data=exp_image, header=in_header)
    exp_hdulist = fits.HDUList([exp_hdu])
    exp_hdulist.writeto(out_dir+raw_file.replace('.fits','.wgt.fits'))

# -----------------------------------------------------------------

def mosaic_galex(name, ra, dec, width, band_dict, working_path, temp_path, meta_path):

    """
    Function to SWarp together GALEX tiles of a given source
    :param name:
    :param ra:
    :param dec:
    :param width:
    :param band_dict:
    :param working_path:
    :param meta_path:
    :return:
    """

    # Declare directories
    id_string = name + '_GALEX_' + band_dict['band_long']


    temp_dir = fs.join(temp_path, "temp_" + band_dict["band_long"])

    raw_in_temp_dir = fs.join(temp_dir, "raw")
    fs.create_directory(raw_in_temp_dir)

    # Create storage directories for Montage and SWarp (deleting any prior), and set appropriate Python working directory
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
    os.mkdir(temp_dir+'Raw')
    os.mkdir(temp_dir+'Diffs_Temp')
    os.mkdir(temp_dir+'Backsub_Temp')
    os.mkdir(temp_dir+'SWarp_Temp')
    os.mkdir(temp_dir+'Reproject_Temp')
    os.mkdir(temp_dir+'Convolve_Temp')
    os.chdir(temp_dir+'Raw')

    overlap_path = fs.join(temp_dir, "overlap_table.dat")

    # Use Montage image metadata table to identify and retrieve which raw GALEX tiles overlap with entire region of interest (handling the case of only a single file)
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='circle', ra=ra, dec=dec, radius=(0.5*width)*(2.0**0.5))

    overlapping_file_paths = np.genfromtxt(overlap_path, skip_header=3, usecols=[31], dtype=('S500'))

    if len(overlapping_file_paths.shape)==0:
        overlapping_file_paths = [overlapping_file_paths.tolist()]
    for overlapping_file_path in overlapping_file_paths:
        shutil.copy(overlapping_file_path, raw_in_temp_dir)

    # Uncompress .fits.gz files
    #[os.system('gunzip '+ listfile) for listfile in os.listdir(raw_in_temp_dir)]

    # Ensure that at least one of the raw GALEX tiles has actual flux coverage at location of source
    raw_files = os.listdir(raw_in_temp_dir)
    coverage = False
    for raw_file in raw_files:

        # Read in map
        in_fitsdata = fits.open(fs.join(raw_in_temp_dir, raw_file))
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()

        # Locate pixel coords
        in_wcs = WCS(in_header)
        location_pix = in_wcs.wcs_world2pix( np.array([[ np.float(ra), np.float(dec) ]]), 0 )[0]
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
        shutil.rmtree(temp_dir)

    elif coverage:

        # Loop over raw tiles, creating exposure maps, and cleaning images to remove null pixels (also, creating convolved maps for later background fitting)
        print('Cleaning '+ str(len(raw_files)) + ' raw maps for ' + id_string)

        #pool = mp.Pool(processes=6)
        #for raw_file in raw_files:
        #    pool.apply_async(GALEX_Clean, args=(raw_file, root_dir, temp_dir, temp_dir+'Reproject_Temp/', band_dict,))
        #    #GALEX_Clean(raw_file, root_dir, temp_dir, temp_dir+'Reproject_Temp/', band_dict)
        #pool.close()
        #pool.join()



        # CLEEAN
        GALEX_Clean(raw_file, root_dir, temp_dir, temp_dir + "Reproject_Temp/", band_dict)



        # Create Montage FITS header
        location_string = str(ra) + ' ' + str(dec)
        pix_size = 3.2
        montage.commands.mHdr(location_string, width, temp_dir + id_string + '_HDR', pix_size=pix_size)



        # Count image files, and move to reprojection directory
        mosaic_count = 0
        for listfile in os.listdir(temp_dir+'/Raw'):
            if '.fits' in listfile:
                mosaic_count += 1
        for listfile in os.listdir(temp_dir+'/Raw'):
            if '.fits' in listfile:
                shutil.move(listfile, temp_dir+'Reproject_Temp')

        # If more than one image file, commence background-matching
        if mosaic_count > 1:
            print('Matching background of '+id_string+' maps')
            GALEX_Zero(temp_dir+'Reproject_Temp', temp_dir+'Convolve_Temp', 'int.fits')

        # Reproject image and weight prior to coaddition
        montage.commands.mImgtbl(temp_dir+'Reproject_Temp',  temp_dir+band+'_Image_Metadata_Table.dat', corners=True)
        montage.commands.mProjExec(temp_dir+band+'_Image_Metadata_Table.dat', temp_dir+id_string+'_HDR', temp_dir+'SWarp_Temp', temp_dir+id_string+'_Proj_Stats.txt', raw_dir=temp_dir+'Reproject_Temp', debug=False, exact=True, whole=False)

        # Rename reprojected files for SWarp
        for listfile in os.listdir(temp_dir+'SWarp_Temp'):
            if '_area.fits' in listfile:
                os.remove(temp_dir+'SWarp_Temp/'+listfile)
            elif 'hdu0_' in listfile:
                os.rename(temp_dir+'SWarp_Temp/'+listfile, temp_dir+'SWarp_Temp/'+listfile.replace('hdu0_',''))

        # Use SWarp to co-add images weighted by their error maps
        print('Co-adding '+id_string+' maps')
        image_width_pixels = str(int((float(width)*3600.)/pix_size))
        os.chdir(temp_dir+'/SWarp_Temp')
        os.system('swarp *int.fits -IMAGEOUT_NAME '+id_string+'_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER '+str(ra)+','+str(dec)+' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE '+image_width_pixels+','+image_width_pixels+' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT')

        # Remove null values, and save finalised map to output directory
        in_fitsdata = fits.open(temp_dir+'/SWarp_Temp/'+id_string+'_SWarp.fits')
        in_image = in_fitsdata[0].data
        in_header = in_fitsdata[0].header
        in_fitsdata.close()
        out_image = in_image.copy()
        out_image[ np.where( out_image==0 ) ] = np.NaN
        out_image[ np.where( out_image<-1E3 ) ] = np.NaN
        out_image[ np.where( out_image<=1E-8 ) ] = 0
        out_hdu = fits.PrimaryHDU(data=out_image, header=in_header)
        out_hdulist = fits.HDUList([out_hdu])
        out_hdulist.writeto(root_dir + 'Montages/' + id_string + '.fits', clobber=True)

        # Clean up
        print('Completed Montaging and SWarping of ' + id_string)
        gc.collect()
        shutil.rmtree(temp_dir)

# -----------------------------------------------------------------

def main():

    # Define paths
    root_dir = '/home/saruman/spx7cjc/DustPedia/GALEX/'

    # Read in source catalogue
    dustpedia_cat = np.genfromtxt(dropbox+'Work/Tables/DustPedia/DustPedia_LEDAWISE_Rehab-Herschel.csv', delimiter=',', names=True, usecols=[0,1,2,7], dtype=[('NAME','S50'), ('SOURCE_RA','f16'), ('SOURCE_DEC','f16'),('D25','f16')])
    name_list = dustpedia_cat['objname']

    # Read in list of already-Montaged sources
    alrady_processed = np.genfromtxt('/home/saruman/spx7cjc/DustPedia/GALEX/GALEX_Already_Processed_List.dat', dtype=('S50')).tolist()

    # State band information
    bands_dict = {'FUV':{'band_short':'fd','band_long':'FUV'},
                  'NUV':{'band_short':'nd','band_long':'NUV'}}

    # Identify targets not yet processed
    remaining_list = []
    for i in range(0, name_list.shape[0]):
        already_done = 0
        name = name_list[i]
        if name not in alrady_processed:
            remaining_list.append(i)

#    dustpedia_cat = dustpedia_cat[remaining_list]

    name_list = dustpedia_cat['objname']
    ra_list = dustpedia_cat['ra2000']
    dec_list = dustpedia_cat['de2000']
    d25_list = dustpedia_cat['d25']

    # If not yet done, produce Montage image table of raw tiles
    for band in bands_dict.keys():

        if not os.path.exists(root_dir+band+'_Image_Metadata_Table.dat'):

            print('Creating '+band+' image metadata table')
            montage.commands.mImgtbl(root_dir+'Raw/'+band, root_dir+band+'_Image_Metadata_Table.dat', corners=True)

    # Record time taken
    time_total = 0.0
    source_counter = 0.0
    source_total = name_list.shape[0]
    time_source_list = []

    # Loop over each source
    for i in range(0, dustpedia_cat.shape[0]):#np.where(name_list=='ESO358-059')[0]:

        name = name_list[i]
        ra = ra_list[i]
        dec = dec_list[i]
        d25 = d25_list[i]
        time_start = time.time()
        if d25<6.0:
            width = 0.5
        elif d25>=6.0:
            width = 1.0

        print('Processing source ' + name)

        # Check if source is in list of already-montaged sources
        if name in alrady_processed:
            print(name+' already processed')
            continue
        source_counter += 1.0

        # Check if any GALEX tiles have coverage of source in question; if not, continue
        bands_in_dict = {}
        for band in bands_dict.keys():
            montage.commands_extra.mCoverageCheck(root_dir+band+'_Image_Metadata_Table.dat', root_dir+'Temporary_Files/'+name+'_'+band+'_Overlap_Check.dat', mode='point', ra=ra, dec=dec)
            if sum(1 for line in open(root_dir+'Temporary_Files/'+name+'_'+band+'_Overlap_Check.dat'))<=3:
                print('No GALEX '+band+' coverage for ' + name)
            else:
                bands_in_dict[band] = bands_dict[band]
            os.remove(root_dir+'Temporary_Files/'+name+'_'+band+'_Overlap_Check.dat')
        if len(bands_in_dict)==0:
            alrady_processed_file = open('/home/saruman/spx7cjc/DustPedia/GALEX/GALEX_Already_Processed_List.dat', 'a')
            alrady_processed_file.write(name+'\n')
            alrady_processed_file.close()
            continue

        # Loop over bands, conducting SWarping function
        for band in bands_in_dict.keys():
            mosaic_galex(name, ra, dec, width, bands_dict[band], root_dir)#pool.apply_async( GALEX_Montage, args=(name, ra, dec, d25, width, bands_dict[band], root_dir+'Temporary_Files/', in_dir, out_dir,) )

        # Record that processing of souce has been compelted
        alrady_processed_file = open('/home/saruman/spx7cjc/DustPedia/GALEX/GALEX_Already_Processed_List.dat', 'a')
        alrady_processed_file.write(name+'\n')
        alrady_processed_file.close()

        # Clean memory, and return timings
        gc.collect()
        time_source = time.time() - time_start
        time_source_list.append(time_source)
        time_total += time_source
        time_source_mins = time_source / 60.0
        time_source_mean = time_total / source_counter
        time_source_median = np.median( time_source_list )
        time_remaining_est_mean = ( ( source_total - source_counter ) * time_source_mean ) / 3600.0
        time_remaining_est_median = ( ( source_total - source_counter ) * time_source_median ) / 3600.0
        print('Estimated time remaining from mean: '+str(time_remaining_est_mean)[:7]+' hours (estimate from median: '+str(time_remaining_est_median)[:7]+' hours)')

# -----------------------------------------------------------------
