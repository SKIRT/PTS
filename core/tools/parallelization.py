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

# Import the relevant PTS classes and modules
from . import introspection
from . import terminal

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
        result = self.pool.apply_async(self.target, args=tuple(args), kwds=kwargs)
        output = PendingOutput(result)
        return output

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

class PendingOutput(object):

    """
    This function ...
    """

    def __init__(self, result):

        """
        This function ...
        :param result:
        """

        self.result = result
        self.output = None

    # -----------------------------------------------------------------

    def request(self):

        """
        This function ...
        :return:
        """

        self.output = self.result.get()

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        if self.output is None: self.request()
        for item in self.output: yield item

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        if self.output is None: self.request()
        return self.output[item]

# -----------------------------------------------------------------

def nnodes():

    """
    This function ...
    :return:
    """

    # I know no way of determining the number of nodes in the network yet
    return 1

# -----------------------------------------------------------------

def sockets_per_node():

    """
    This function ...
    :return:
    """

    if introspection.is_macos(): return 1
    elif introspection.is_linux():
        lines = terminal.execute("lscpu")
        for line in lines:
            if "socket(s):" in line.lower(): return int(line.split(": ")[1])
        raise RuntimeError("Could not determine the number of sockets")
    else: raise RuntimeError("Platforms other than MacOS and Linux are not supported")

# -----------------------------------------------------------------

def cores_per_socket():

    """
    This function ...
    :return:
    """

    if introspection.is_macos(): return ncores()
    elif introspection.is_linux():
        lines = terminal.execute("lscpu")
        for line in lines:
            if "core(s) per socket" in line.lower(): return int(line.split(": ")[1])
        raise RuntimeError("Could not determine the number of cores per socket")
    else: raise RuntimeError("Platforms other than MacOS and Linux are not supported")

# -----------------------------------------------------------------

def ncores():

    """
    This function ...
    :return:
    """

    return int(psutil.cpu_count(logical=False))

# -----------------------------------------------------------------

def nthreads_per_core():

    """
    This function ...
    :return:
    """

    return int(psutil.cpu_count() / ncores())

# -----------------------------------------------------------------

def virtual_memory():

    """
    This function ...
    :return:
    """

    return float(psutil.virtual_memory().total) * Unit("byte")

# -----------------------------------------------------------------
