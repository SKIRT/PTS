#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.templates Contains functions to generate (or load) ski file templates.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools import introspection
from ..simulation.skifile import SkiFile
from .smile import SKIRTSmileSchema
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

# TODO: IN THE END, THESE FUNCTIONS COULD (SHOULD) USE THE SKIRTSMILESCHEMA CLASS TO GENERATE TEMPLATES,
# INSTEAD OF USING PRE-CREATED FILES

# -----------------------------------------------------------------

dat_ski_path = introspection.pts_modeling_ski_templates_path()

# -----------------------------------------------------------------

def get_pan_template():

    """
    This function ...
    :return: 
    """

    # Determine the path to the template ski file for panchromatic simulations
    pan_ski_path = fs.join(dat_ski_path, "pan.ski")

    # Load and return the ski file
    return SkiFile(pan_ski_path)

# -----------------------------------------------------------------

def get_oligo_template(wavelengths=None):

    """
    This function ...
    :param wavelengths:
    :return: 
    """

    if wavelengths is None: wavelengths = [1. * u("micron")]

    ski = get_pan_template()
    ski.to_oligochromatic(wavelengths)
    return ski

# -----------------------------------------------------------------

def get_oneparticle_template():

    """
    This function ...
    :return: 
    """

    # Determine the path
    path = fs.join(dat_ski_path, "oneparticle.ski")

    # Load and return
    return SkiFile(path)

# -----------------------------------------------------------------
