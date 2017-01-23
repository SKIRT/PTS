#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.parallelization Functions related to parallelization.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import psutil
from multiprocessing import Pool

# Import astronomical modules
from astropy.units import Unit

# -----------------------------------------------------------------

class ParallelTarget(object):

    """
    This function ...
    """

    def __init__(self, target, nprocesses):

        """
        This function ...
        :param target:
        :param nprocesses:
        """

        # Set the target
        self.target = target

        # Get the process pool
        self.nprocesses = nprocesses

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        # Initialize the process pool
        self.pool = Pool(processes=self.nprocesses)

    # -----------------------------------------------------------------

    def __call__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Launch
        self.pool.apply_async(self.target, args=tuple(args), kwds=kwargs)

    # -----------------------------------------------------------------

    def __exit__(self):

        """
        This function ...
        :return:
        """

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

# -----------------------------------------------------------------

def ncores():

    """
    This function ...
    :return:
    """

    return psutil.cpu_count(logical=False)

# -----------------------------------------------------------------

def nthreads_per_core():

    """
    This function ...
    :return:
    """

    return psutil.cpu_count() / ncores()

# -----------------------------------------------------------------

def virtual_memory():

    """
    This function ...
    :return:
    """

    return float(psutil.virtual_memory().total) * Unit("byte")

# -----------------------------------------------------------------
