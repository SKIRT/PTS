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
import psutil
import multiprocessing
import numpy as np

# -----------------------------------------------------------------

def free_memory():

    """
    This function ...
    :return:
    """

    # Get the currently available virtual memory (in gigabytes)
    memory = psutil.virtual_memory().available / 1e9
    
    return memory

# -----------------------------------------------------------------

def free_cpus():

    """
    This function ...
    :return:
    """

    # Get the total number of processors on this system
    total = multiprocessing.cpu_count()

    # Get the load of the different processors
    load = np.array(psutil.cpu_percent(percpu=True))/100.0

    # Calculate the the number of full processors (where 2 processors with loads x and y contribute as a processor with load x+y)
    full = np.sum(load)

    # Get the number of free processors
    free = total - full
    
    return free

# -----------------------------------------------------------------
