#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.image_fwhms

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import formatting as fmt
from pts.modeling.component.galaxy import get_data_images_path
from pts.magic.core.fits import get_header
from pts.magic.tools import headers
from pts.magic.convolution.kernels import has_variable_fwhm, get_fwhm
from pts.core.units.parsing import parse_unit as u
from pts.dustpedia.core.properties import DustPediaProperties
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

data_images_path = get_data_images_path(modeling_path)

# -----------------------------------------------------------------

fwhms = dict()

default_fwhm = 2.0 * u("arcsec")

# -----------------------------------------------------------------

# Get the DustPedia properties instance
properties = DustPediaProperties()

# Get FWHMs
official_fwhms = properties.fwhms

# -----------------------------------------------------------------

# Loop over the images
for image_path, image_name in fs.files_in_path(data_images_path, extension="fits", not_contains="poisson", returns=["path", "name"], recursive=True, recursion_level=1):

    # Load header and get filter
    header = get_header(image_path)
    fltr = headers.get_filter(image_name, header)

    # Set the FWHM if the instrument has a fixed PSF
    if has_variable_fwhm(fltr):
        if fltr in official_fwhms: fwhm = official_fwhms[fltr]
        else: fwhm = default_fwhm
    else: fwhm = get_fwhm(fltr)

    # Determine filter name
    filter_name = str(fltr)

    # Set FWHM
    fwhms[filter_name] = fwhm

# -----------------------------------------------------------------

# Sort on FWHM
sorted_filter_names = sorted(fwhms.keys(), key=lambda name: fwhms[name])

# -----------------------------------------------------------------

print("")
with fmt.print_in_columns(2) as print_row:
    for name in sorted_filter_names:
        print_row(name, fwhms[name])
print("")

# -----------------------------------------------------------------
