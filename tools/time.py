#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package tools.time Useful functions for generating timestamps etc.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from datetime import datetime

# -----------------------------------------------------------------

def parse(line):

    """
    This function returns the datetime object corresponding to a certain time stamp (as a string)
    :param line:
    :return:
    """

    # Parse the line
    date, time, dummy = line.split(None, 2)
    try: day, month, year = date.split('/')
    except ValueError: return None
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')

    # Convert the strings to integers
    year = int(year)
    month = int(month)
    day = int(day)
    hour = int(hour)
    minute = int(minute)
    second = int(second)
    microsecond = int(microsecond)

    # Create and return a datetime object
    return datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=second, microsecond=microsecond)

# -----------------------------------------------------------------

def timestamp(time):

    """
    This function generates a timestamp (as a string) from a datetime object
    :param time:
    :return:
    """

    # Return a timestamp accurate up to the millisecond
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S.%f")[:-3]

# -----------------------------------------------------------------

def unique_name(name):

    """
    This function ...
    :param name:
    :return:
    """

    return name + datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")

# -----------------------------------------------------------------
