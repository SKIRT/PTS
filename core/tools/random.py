#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.random Defines a global random number generator and contains some functions to save and
#  load its state to/from file.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from . import serialization

# -----------------------------------------------------------------

# Based on: http://stackoverflow.com/questions/20911147/choose-random-seed-and-save-it

# -----------------------------------------------------------------

# SKIRT seed
skirt_seed = 4357

# -----------------------------------------------------------------

# The global random genereator = a singleton of this module
prng = None

# -----------------------------------------------------------------

def init_prng(seed=None):

    """
    This function ...
    :return:
    """

    global prng

    # Create a new random generator
    prng = np.random.RandomState(seed)

# -----------------------------------------------------------------

def setup_prng(seed):

    """
    This function ...
    :param seed:
    :return:
    """

    prng.seed(seed)

# -----------------------------------------------------------------

def get_state():

    """
    This function ...
    :return:
    """

    # Return the state
    return prng.get_state()

# -----------------------------------------------------------------

def set_state(state):

    """
    This function ...
    :return:
    """

    # Set the random state
    prng.set_state(state)

# -----------------------------------------------------------------

def save_state(path):

    """
    This function ...
    :return:
    """

    state = prng.get_state()
    serialization.dump(state, path)

# -----------------------------------------------------------------

def load_state(path):

    """
    This function ...
    :param path:
    :return:
    """

    state = serialization.load(path)
    prng.set_state(state)

# -----------------------------------------------------------------

# Initialize the global random generator
init_prng()

# -----------------------------------------------------------------
