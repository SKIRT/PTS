#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.monitoring Provides functions that can report information about the amount of free memory
#  or cores on the system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import resource
from sys import platform
import psutil
import multiprocessing
import numpy as np

# Import the relevant PTS classes and modules
from ..units.parsing import parse_unit as u

# -----------------------------------------------------------------

def memory_usage():

    """
    This function ...
    :return: 
    """

    # If we are on linux
    if platform == "linux" or platform == "linux2":

        kilobytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss # peak memory usage (bytes on OS X, kilobytes on Linux)
        gigabytes = kilobytes * 1e-6

        return gigabytes

    # If we are on Mac OS X
    elif platform == "darwin":

        kilobytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss # peak memory usage (bytes on OS X, kilobytes on Linux)
        gigabytes = kilobytes * 1e-9

        return gigabytes

    # We don't support Windows
    elif platform == "win32": raise EnvironmentError("The Windows operating system is not supported")

    # Unrecognized platform
    else: raise EnvironmentError("Unrecognized platform")

# -----------------------------------------------------------------

def free_memory():

    """
    This function ...
    :return:
    """

    # Get the currently available virtual memory (in gigabytes)
    memory = psutil.virtual_memory().available / 1e9
    return memory * u("GB")

# -----------------------------------------------------------------

def free_cpus():

    """
    This function ...
    :return:
    """

    # Get the total number of processors on this system
    total = multiprocessing.cpu_count()

    # Get the load of the different processors
    load = np.array(psutil.cpu_percent(percpu=True)) / 100.0

    # Calculate the the number of full processors (where 2 processors with loads x and y contribute as a processor with load x+y)
    full = np.sum(load)

    # Get the number of free processors
    free = total - full

    # Return the number of free processors (real number)
    return free

# -----------------------------------------------------------------
