#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.tools.extinction Functions related to dust extinction.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# -----------------------------------------------------------------

# Equations:

# A = -2.5 x log ( I / I_0)

# => I = I_0 x 10^(-A/2.5)     (observed)
#      = I_0 x 10^(-0.4 x A)

# => I_0 = I x 10^(A/2.5)     (intrinsic)
#        = I x 10^(0.4 x A)

# -----------------------------------------------------------------

def attenuation(obs, intr):

    """
    This function ...
    :param obs: observed 
    :param intr: intrinsic 
    :return: 
    """

    return -2.5 * np.log10(obs / intr)

# -----------------------------------------------------------------

def observed(intr, att):

    """
    This function ...
    :param intr:
    :param att: 
    :return: 
    """

    return intr * intrinsic_to_observed_factor(att)

# -----------------------------------------------------------------

def intrinsic(obs, att):

    """
    This function ...
    :param obs: 
    :param att: 
    :return: 
    """

    return obs * observed_to_intrinsic_factor(att)

# -----------------------------------------------------------------

def intrinsic_to_observed_factor(att):

    """
    This function ...
    :param att: 
    :return: 
    """

    return 10**(-0.4 * att)

# -----------------------------------------------------------------

def observed_to_intrinsic_factor(att):

    """
    This function ...
    :param att: 
    :return: 
    """

    return 10**(0.4 * att)

# -----------------------------------------------------------------
