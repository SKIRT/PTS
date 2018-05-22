#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.mosaic Contains mosaicing tools.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.table import Table
from astropy.io.fits import Header

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.basics.log import log
from ...core.tools import terminal
from ...core.tools import formatting as fmt

# -----------------------------------------------------------------

temp_montage_path = introspection.create_unique_temp_dir("montage")

# -----------------------------------------------------------------

def generate_meta_file(path, meta_dir_path=None):

    """
    This function ...
    :param path:
    :param meta_dir_path:
    :return:
    """

    import montage_wrapper as montage

    if meta_dir_path is None: meta_dir_path = path

    # Inform the user
    log.info("Generating meta file ...")

    # Path of the meta file
    meta_filename = "meta.dat"
    meta_path = fs.join(meta_dir_path, meta_filename)

    # Change working directory
    dirname = fs.name(path)
    cwd = fs.change_cwd(fs.directory_of(path))

    # Path of status file
    # FOR DEBUGGING
    #status_filename = "status.txt"
    #status_path = fs.join(meta_dir_path, status_filename)

    # Debugging
    log.debug("Changed working directory to '" + fs.directory_of(path) + "' ...")
    log.debug("Fields directory name: '" + dirname + "'")
    log.debug("Fields directory path: '" + path + "'")
    log.debug("Meta file path: '" + meta_path + "'")
    #log.debug("Status file path: '" + status_path + "'")

    dirname = path

    # Get the image table of which images cover a given part of the sky
    #montage.commands.mImgtbl(dirname, meta_path, corners=True, debug=log.is_debug, output_invalid=True, status_file=status_path) # debuggin
    montage.commands.mImgtbl(dirname, meta_path, corners=True)

    # Change working directory again
    fs.change_cwd(cwd)

    # Check whether the meta file contains any lines
    if not fs.contains_lines(meta_path):
        #log.debug("STATUS:")
        #for line in fs.read_lines(status_path): log.debug(line)
        raise RuntimeError("The meta table (" + meta_path + ") doesn't contain any lines")

    # Return the path to the created file
    return meta_path

# -----------------------------------------------------------------

def generate_overlap_file(path, ra, dec, meta_path, mode='box', width=None, radius=None, overlap_path=None):

    """
    This function ...
    :param path:
    :param ra:
    :param dec:
    :param meta_path:
    :param mode: 'box', 'point'
    :param width:
    :param radius:
    :return:
    """

    import montage_wrapper as montage

    # Inform the user
    log.info("Generating overlap file ...")

    # Path of the overlap file
    if overlap_path is None: overlap_path = fs.join(path, "overlap.dat")

    # Convert RA and DEC into degrees
    ra_deg = ra.to("deg").value
    dec_deg = dec.to("deg").value
    radius_deg = radius.to("deg").value if radius is not None else None
    width_deg = width.to("deg").value if width is not None else None

    # Check the coverage for our galaxy
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode=mode, ra=ra_deg, dec=dec_deg, width=width_deg, radius=radius_deg)

    # Check if there is any coverage for this galaxy and band
    if sum(1 for line in open(overlap_path)) <= 3: log.warning("No coverage")

    # Return the path to the created file
    return overlap_path

# -----------------------------------------------------------------

def generate_overlapping_file_paths(path, ra, dec, meta_path, mode="box", width=None, radius=None):

    """
    This function ...
    :param path:
    :param ra:
    :param dec:
    :param meta_path:
    :param mode:
    :param width:
    :param radius:
    :return:
    """

    # Generate the file
    overlap_path = generate_overlap_file(path, ra, dec, meta_path, mode, width, radius)

    # Return the paths of the overlapping files
    return get_overlapping_file_paths(overlap_path)

# -----------------------------------------------------------------

def get_overlapping_file_paths(overlap_path):

    """
    This function ...
    :param overlap_path: 
    :return: 
    """

    # Get file paths of overlapping observations
    overlapping_file_paths = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype=str)

    # Return the file paths
    return overlapping_file_paths

# -----------------------------------------------------------------

def get_field_table(cutout_center, cutout_width, band):

    """
    This function ...
    :return:
    """

    import montage_wrapper as montage

    # Get the coordinate range for this galaxy
    ra_deg = cutout_center.ra.to("deg").value
    dec_deg = cutout_center.dec.to("deg").value
    width_deg = cutout_width.to("deg").value

    # Determine path for the table
    path = fs.join(temp_montage_path, "fields.tbl")

    # Get the info
    montage.mArchiveList("SDSS", band, str(ra_deg) + " " + str(dec_deg), width_deg, width_deg, path)

    # Load the table
    table = Table.read(path, format="ascii")

    # Return the table
    return table

# -----------------------------------------------------------------

def make_header(ra, dec, width, pix_size, returns="header"):

    """
    This function ...
    :param ra:
    :param dec:
    :param width:
    :param pix_size:
    :param returns:
    :return:
    """

    import montage_wrapper as montage

    # ra and dec are taken from the ra2000 and de2000 of the attached DustPedia_LEDAWISE_Herschel.csv table
    # width is 0.5 degrees for galaxies with D25<6 arcmin, and 1 degree for galaxies with D25>=6 arcmin (as listed in DustPedia_LEDAWISE_Herschel.csv)
    # pix_size is 3.2 for GALEX, and 0.45 for SDSS.

    # Convert to degrees and the pixelsize in arcseconds
    ra_deg = ra.to("deg").value
    dec_deg = dec.to("deg").value
    width_deg = width.to("deg").value
    pix_size_arcsec = pix_size.to("arcsec").value

    # Determine the path to the temporary header file
    header_path = fs.join(temp_montage_path, "header.hdr")

    # Create the header
    montage.commands.mHdr(str(ra_deg) + ' ' + str(dec_deg), width_deg, header_path, pix_size=pix_size_arcsec)

    # Load the header
    if returns == "header": return Header.fromtextfile(header_path)
    elif returns == "path": return header_path
    elif returns == ["header", "path"]: return Header.fromtextfile(header_path), header_path
    elif returns == ["path", "header"]: return header_path, Header.fromtextfile(header_path)
    else: raise ValueError("Invalid option for 'returns'")

# -----------------------------------------------------------------

def filter_non_overlapping(ngc_name, band, fields_path, cutout_center, cutout_width, mode='box', path=None):

    """
    This function ...
    :param ngc_name:
    :param band:
    :param fields_path:
    :param cutout_center:
    :param cutout_width:
    :param mode:
    :param path:
    :return:
    """

    # Inform the user
    log.info("Getting the overlap of the SDSS observations for " + ngc_name + " in the " + band + " band ...")

    if path is None: path = fields_path

    # Debugging
    log.debug("Fields path: '" + fields_path + "'")

    # Generate meta file
    meta_path = generate_meta_file(fields_path, meta_dir_path=path)

    # Generate overlap file
    # path, ra, dec, meta_path, mode='box', width=None, radius=None
    overlap_path = generate_overlap_file(fields_path, cutout_center.ra, cutout_center.dec, meta_path, mode, cutout_width)

    # Get the names of the overlapping image files
    overlap_files = np.genfromtxt(overlap_path, skip_header=3, usecols=[32], dtype=str)

    ntotal = 0
    removed = []

    # Loop over the FITS files in the temp directory, remove non-overlapping
    for path in fs.files_in_path(fields_path, extension="fits"):

        ntotal += 1

        # If the path is not in the overlap_files list, remove the FITS file
        if path not in overlap_files:

            # Debugging
            log.debug("Removing the '" + fs.name(path) + "' image since it does not overlap with the target area ...")

            # Removed
            fs.remove_file(path)

            # Add to list
            removed.append(fs.name(path))

    # Debugging
    log.debug("Removed " + str(len(removed)) + " out of " + str(ntotal) + " images for the " + band + "band")

    # Show the removed images
    fmt.print_files_in_list(removed, "removed images for the " + band + "band")

    # Return the meta path and overlap path
    return meta_path, overlap_path

# -----------------------------------------------------------------

def reproject(input_path, output_path, metatable_path, header_path):

    """
    This function ...
    :return:
    """

    import montage_wrapper as montage

    proj_stats_path = fs.join(input_path, "Proj_Stats.txt")
    montage.commands.mProjExec(metatable_path, header_path, output_path, proj_stats_path, raw_dir=input_path, debug=False, exact=True, whole=False)

    # WHOLE IS IMPORTANT HERE
    # WE ACTUALLY DON'T WANT TO REPROJECT TO THE EXACT PIXELGRID DEFINED BY THE HEADER HERE,
    # BUT RATHER CUT OFF THE ORIGINAL MAPS WHERE THEY ARE OUTSIDE OF THE FIELD OF VIEW OF THE HEADER DEFINED AREA
    # SO NO ADDING OF NEW PIXELS IS DONE TO HAVE THE EXACT PIXELGRID DEFINED BY THE HEADER
    # THIS IS PRESUMABLY DONE JUST TO MAKE THE SWARPING MORE EFFICIENT ??

# -----------------------------------------------------------------

def mosaic(working_path, band, center, width_pixels):

    """
    This function ...
    :param working_path:
    :param band:
    :param center:
    :param width_pixels:
    :return:
    """

    ra_deg = center.ra.to("deg").value
    dec_deg = center.dec.to("deg").value

    # Change directories
    old_cwd = fs.change_cwd(working_path)

    # Determine command string
    swarp_command_string = 'swarp *int.fits -IMAGEOUT_NAME ' + band + '_SWarp.fits -WEIGHT_SUFFIX .wgt.fits -CENTER_TYPE MANUAL -CENTER ' + str(ra_deg) + ',' + str(dec_deg) + ' -COMBINE_TYPE WEIGHTED -COMBINE_BUFSIZE 2048 -IMAGE_SIZE ' + width_pixels + ',' + width_pixels + ' -MEM_MAX 4096 -NTHREADS 4 -RESCALE_WEIGHTS N  -RESAMPLE N -SUBTRACT_BACK N -VERBOSE_TYPE QUIET -VMEM_MAX 4095 -WEIGHT_TYPE MAP_WEIGHT'

    # Debugging
    log.debug("SWARP command: '" + swarp_command_string + "'")

    #os.system(swarp_command_string)
    terminal.execute(swarp_command_string)

    # Change back to the original working directory
    fs.change_cwd(old_cwd)

    # Swarp result path
    swarp_result_path = fs.join(working_path, band + "_SWarp.fits")

    # Return the path to the resulting mosaic
    return swarp_result_path

# -----------------------------------------------------------------
