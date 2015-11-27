#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import logging

# Import astronomical modules
from astropy import logger

# -----------------------------------------------------------------

def new_log(name, level):

    """
    This function ...
    :param level:
    :return:
    """

    # Create the logger
    log = logger.AstropyLogger(name)

    # Initialize the logger
    log._set_defaults()

    # Set the log level
    log.setLevel(level)

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
