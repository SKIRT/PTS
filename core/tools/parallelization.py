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

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from . import introspection
from . import terminal
from ..basics.log import log

# Check the multiprocessing state and import modules
try:
   from multiprocessing import cpu_count, Pool
   CPU_COUNT = cpu_count()
   MULTI_PROCESSING = True if CPU_COUNT > 1 else False
   log.debug("You have " + str(CPU_COUNT) + " CPU cores and multiprocessing support is present")
except ImportError:
    cpu_count = Pool = None
    CPU_COUNT = None
    MULTI_PROCESSING = False
    log.debug("You don't have multiprocessing support")

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

        # The number of tasks
        self.ntasks = 0

        # The number of completed processes
        self.ncompleted = 0

        # The process pool
        self.pool = None

    # -----------------------------------------------------------------

    def __enter__(self):

        """
        This function ...
        :return:
        """

        # Reset the number of completed processes
        self.ncompleted = 0

        # Reset the number of tasks
        self.ntasks = 0

        # Initialize the process pool if required
        if self.nprocesses > 1:
            self.pool = Pool(processes=self.nprocesses)
            self.nprocesses = self.pool._processes # make sure nprocesses is set

        # Return ourselves
        return self

    # -----------------------------------------------------------------

    def __call__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Increment the number of tasks
        self.ntasks += 1

        # Launch
        if self.pool is not None:

            result = self.pool.apply_async(self.target, args=tuple(args), kwds=kwargs, callback=self.complete)
            output = PendingOutput(result)

        else:

            output = self.target(*args, **kwargs)
            self.complete()
            output = ReadyOutput(output)

        # Return the output
        return output

    # -----------------------------------------------------------------

    def complete(self, result=None):

        """
        This function ...
        :param result:
        :return:
        """

        #self.results.extend(result)
        self.ncompleted += 1
        #print('Progress: {:.2f}%'.format((self.ncompleted / self.nprocesses) * 100))
        #print("Process " + str(self.ncompleted) + " of " + str(self.nprocesses) + " has finished a task")

        # Debugging
        log.debug("Task " + str(self.ncompleted) + " of " + str(self.ntasks) + " has been completed")

    # -----------------------------------------------------------------

    def __exit__(self, exc_type, exc_value, traceback):

        """
        This function ...
        :return:
        """

        # Error occured
        if exc_type is not None:
            log.error("A " + str(exc_type) + " occured")
            #log.error(str(exc_value))
            #print(traceback)

        # Close and join the process pool
        if self.pool is not None:
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

class ReadyOutput(object):

    """
    This class ...
    """

    def __init__(self, output):

        """
        This function ...
        :param output:
        """

        self.output = output

    # -----------------------------------------------------------------

    def request(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        This function ...
        :return:
        """

        for item in self.output: yield item

    # -----------------------------------------------------------------

    def __getitem__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

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

def has_hyperthreading():

    """
    This function ...
    :return:
    """

    return nthreads_per_core() > 1

# -----------------------------------------------------------------

def virtual_memory():

    """
    This function ...
    :return:
    """

    return float(psutil.virtual_memory().total) * Unit("byte")

# -----------------------------------------------------------------

#def set_(name):
#    sys.stdout = open(str(os.getpid()) + ".out", "w")
#    info('function f')
#    #print 'hello', name

# -----------------------------------------------------------------
