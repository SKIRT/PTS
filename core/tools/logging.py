#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.logging Provides functions for creating loggers and linking log files to them.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import sys
import logging

# Import astronomical modules
from astropy import logger

from . import inspection

# -----------------------------------------------------------------

class MemuseFilter(logging.Filter):

    def filter(self, record):
        """ This function overrides logging.Filter, adds memuse as a field
        """
        record.memuse = self.str_mem()
        return True

    # Following code from http://stackoverflow.com/a/938800/819110:
    _proc_status = '/proc/%d/status' % os.getpid()
    _scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
              'KB': 1024.0, 'MB': 1024.0*1024.0}

    def _VmB(self,VmKey):
        """Private.
        """
        # get pseudo file  /proc/<pid>/status
        try:
            t = open(self._proc_status)
            v = t.read()
            t.close()
        except:
            return 0.0  # non-Linux?
        # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
        i = v.index(VmKey)
        v = v[i:].split(None, 3)  # whitespace
        if len(v) < 3:
            return 0.0  # invalid format?
        # convert Vm value to bytes
        return float(v[1]) * self._scale[v[2]]

    def memory(self,since=0.0):
        """Return memory usage in bytes.
        """
        return self._VmB('VmSize:') - since

    def swapsize(self,since=0.0):
        """Return swap size in bytes.
        """
        return self._VmB('VmSwap:') - since

    def byte_to_mb(self,byte):
        """return size in MB (being lazy)
        """
        return byte/(1024*1024)

    def str_mem(self):
        """Return a string with the total memuse and swap size in MB
        """
        return "MemTotal:%.0fM,Swap:%.0fM"%(self.byte_to_mb(self.memory()),self.byte_to_mb(self.swapsize()) )

# -----------------------------------------------------------------

def new_log(name, level, style="astropy"):

    """
    This function ...
    :param level:
    :return:
    """

    # Create the logger
    log = logger.AstropyLogger(name)

    # Astropy logging style
    if style == "astropy":

        # Initialize the logger
        log._set_defaults()
        # Set the log level
        log.setLevel(level)

    # SKIRT logging style
    elif style == "skirt":

        # Create a stream handler
        sh = logging.StreamHandler(sys.stdout)
        if level == "INFO": sh.setLevel(logging.INFO)
        elif level == "DEBUG": sh.setLevel(logging.DEBUG)
        elif level == "WARNING": sh.setLevel(logging.WARNING)
        elif level == "ERROR": sh.setLevel(logging.ERROR)
        else: raise ValueError("Unknown log level")

        # Create and set the formatter
        formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s (%(origin)s)", "%d/%m/%Y %H:%M:%S")
        sh.setFormatter(formatter)

        # Add the stream handler
        log.addHandler(sh)

    # Return the logger
    return log

# -----------------------------------------------------------------

def new_memory_log():

    """
    This function ...
    :return:
    """

    #memory_log_conf_path = os.path.join(inspection.pts_root_dir, "memlogging.conf")
    #logging.basicConfig(memory_log_conf_path)

    #logging.basicConfig(format="%(asctime)-15s %(name)-5s %(levelname)-8s %(memuse)-22s %(message)s")

    log = logging.getLogger('')               # Get root logger

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)-15s %(name)-5s %(levelname)-8s %(memuse)-22s %(message)s")
    sh.setFormatter(formatter)

    log.addHandler(sh)

    f = MemuseFilter()                        # Create filter
    log.handlers[0].addFilter(f)         # The ugly part:adding filter to handler

    return log

# -----------------------------------------------------------------

def link_file_log(log, path, level):
    
    """
    This function ...
    """
    
    # Create a file handler
    fh = logger.FileHandler(path)
    
    # Create a formatter
    formatter = logging.Formatter("%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s (%(origin)s)", "%d/%m/%Y %H:%M:%S")

    # Set the formatter
    fh.setFormatter(formatter)
    
    # Set the level 
    fh.setLevel(level)
    
    # Add the handler to the log instance
    log.addHandler(fh)

# -----------------------------------------------------------------
