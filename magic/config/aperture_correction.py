#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add input section
definition.add_section("input", "the names of the input files")
definition.sections["input"].add_required("cutout", "file_path", "name of the cutout FITS file on which the photometry is performed")
definition.sections["input"].add_required("psf", "file_path", "the name of the FITS file with the PSF of the image on which the photometry is performed")

# Add required settings
definition.add_required("pix_arcsec", "real", "The width, in arscec, of the pixels in the map photometry is being performed upon (this is needed in case there is a pixel size mismatch with PSF)")
definition.add_required("semimaj_pix", "real", "Semi-major axis of photometric aperture, in pixels")
definition.add_required("axial_ratio", "real", "Axial ratio of photometric aperture")
definition.add_required("angle", "real", "Position angle of photometric aperture, in degrees")
definition.add_required("centre_i", "real", "Zero-indexed, 0th-axis coordinate (equivalent to y-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture")
definition.add_required("centre_j", "real", "Zero-indexed, 1st-axis coordinate (equivalent to x-axis one-indexed coordinates in FITS terms) of centre position of photometric aperture")
definition.add_required("semimaj_pix_annulus_outer", "real", "Semi-major axis length of the outer annulus ellipse")
definition.add_required("semimaj_pix_annulus_inner", "real", "Semi-major axis length of the inner annulus ellipse")
definition.add_required("axial_ratio_annulus", "real", "Axial ratio of annulus ellipses")
definition.add_required("annulus_angle", "real", "Position angle of annulus, in degrees")
definition.add_required("annulus_centre_i", "real", "0th-axis coordinate (y) of centre position of annulus ellipses")
definition.add_required("annulus_centre_j", "real", "1th-axis coordinate (x) of centre position of annulus ellipses")

# -----------------------------------------------------------------
