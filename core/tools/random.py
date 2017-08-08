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
from ..basics.log import log

# -----------------------------------------------------------------

# Based on: http://stackoverflow.com/questions/20911147/choose-random-seed-and-save-it

# -----------------------------------------------------------------

# IMPORTANT: (from NumPy docs)

# randint (low[, high, size, dtype])	Return random integers from low (inclusive) to high (exclusive).  # EXCLUSIVE END !
# random_integers (low[, high, size])	Random integers of type np.int between low and high, inclusive. # INCLUSIVE END !

# THE LATTER IS DEPRECATED : WE ARE GOING TO USE RANDINT(MIN, MAX + 1) THEN

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

    # Inform the user
    log.info("Initializing the random number generator ...")

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

    # Inform the user
    log.info("Seeding the random number generator ...")

    # Give seed
    prng.seed(seed)

    # Return the instance
    return prng

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

    # Inform the user
    log.info("Saving the state of the random number generator to '" + path + "' ...")

    state = prng.get_state()
    serialization.dump(state, path)

# -----------------------------------------------------------------

def load_state(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Inform the user
    log.info("Loading the state of the random number generator from '" + path + "' ...")

    state = serialization.load(path)
    prng.set_state(state)

# -----------------------------------------------------------------

# Initialize the global random generator
init_prng()

# -----------------------------------------------------------------
