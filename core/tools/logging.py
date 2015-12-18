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
import sys
import logging

# Import astronomical modules
from astropy import logger

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
