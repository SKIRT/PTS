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

# Import astronomical modules
import montage_wrapper as montage
from astropy.table import Table

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection

# -----------------------------------------------------------------

temp_montage_path = fs.create_directory_in(introspection.pts_temp_dir, "montage")

# -----------------------------------------------------------------

def generate_meta_file(path):

    """
    This function ...
    :param path:
    :param ra:
    :param dec:
    :param width:
    :return:
    """

    # Inform the user
    log.info("Generating meta file ...")

    # Path of the meta file
    meta_path = fs.join(path, "meta.dat")

    # Get the image table of which images cover a given part of the sky
    montage.commands.mImgtbl(path, meta_path, corners=True)

    # Return the path to the created file
    return meta_path

# -----------------------------------------------------------------

def generate_overlap_file(path, ra, dec, width, meta_path):

    """
    This function ...
    :param path:
    :param ra:
    :param dec:
    :param width:
    :param meta_path:
    :return:
    """

    # Inform the user
    log.info("Generating overlap file ...")

    # Path of the overlap file
    overlap_path = fs.join(path, "overlap.dat")

    # Check the coverage for our galaxy
    montage.commands_extra.mCoverageCheck(meta_path, overlap_path, mode='box', ra=ra, dec=dec, width=width)

    # Return the path to the created file
    return overlap_path

# -----------------------------------------------------------------

def get_field_table(cutout_center, cutout_width, band):

    """
    This function ...
    :return:
    """

    # Get the coordinate range for this galaxy
    ra = cutout_center.ra.to("deg").value
    dec = cutout_center.dec.to("deg").value
    width = cutout_width.to("deg").value

    # Determine path for the table
    path = fs.join(temp_montage_path, "fields.tbl")

    # Get the info
    montage.mArchiveList("SDSS", band, str(ra) + " " + str(dec), width, width, path)

    # Load the table
    table = Table.read(path, format="ascii")

    # Return the table
    return table

# -----------------------------------------------------------------
