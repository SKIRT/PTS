#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.api Defines functions that bring the features of the map making classes out in a simple API.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .dust.attenuation import AttenuationDustMapMaker
from .attenuation.cortese import CorteseAttenuationMapMaker
from .oldstars.disk import DiskOldStellarMapMaker
from .youngstars.young import YoungStellarMapMaker
from .ionizingstars.ionizing import IonizingStellarMapMaker

# -----------------------------------------------------------------

def make_fuv_attenuation_map():

    """
    This function ...
    :return: 
    """

# -----------------------------------------------------------------

def make_dust_map():

    """
    This function ...
    :return: 
    """

# -----------------------------------------------------------------

def make_old_stellar_map():

    """
    This function ...
    :return: 
    """

# -----------------------------------------------------------------

def make_young_stellar_map():

    """
    This function ...
    :return: 
    """

    maker = YoungStellarMapMaker()

    maker.run()

# -----------------------------------------------------------------

def make_ionizing_stellar_map():

    """
    This function ...
    :return: 
    """

    maker = IonizingStellarMapMaker()

    maker.run()

# -----------------------------------------------------------------
