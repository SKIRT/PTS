#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.time Provides functions for dealing with timestamps etc.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import time as _time
from datetime import datetime

# -----------------------------------------------------------------

intervals = (
    ('weeks', 604800),  # 60 * 60 * 24 * 7
    ('days', 86400),    # 60 * 60 * 24
    ('hours', 3600),    # 60 * 60
    ('minutes', 60),
    ('seconds', 1),
    )

# -----------------------------------------------------------------

def display_time(seconds, granularity=2):

    """
    This function ...
    :param seconds:
    :param granularity:
    :return:
    """

    result = []
    for name, count in intervals:
        if name == "seconds": value = seconds / count
        else: value = seconds // count
        if value:
            seconds -= value * count
            if value == 1:
                name = name.rstrip('s')
            result.append("{} {}".format(value, name))
    return ', '.join(result[:granularity])

# -----------------------------------------------------------------

def time():

    """
    This function ...
    :return:
    """

    return _time.time()

# -----------------------------------------------------------------

def wait(seconds):

    """
    This function ...
    :param seconds:
    :return:
    """

    _time.sleep(seconds)

# -----------------------------------------------------------------

def iterate(seconds):

    """
    This function ...
    :param seconds:
    :return:
    """

    while True:
        yield now()
        wait(seconds)

# -----------------------------------------------------------------

def parse(timestamp):

    """
    This function ...
    :param string:
    :return:
    """

    # Parse the timestamp
    date, time = timestamp.split()

    # Combine the date and time stamp to a datetime object
    return parse_date_time(date, time)

# -----------------------------------------------------------------

def parse_line(line):

    """
    This function returns the datetime object corresponding to a certain time stamp (as a string)
    :param line:
    :return:
    """

    # Parse the line
    date, time, dummy = line.split(None, 2)

    # Combine the date and time stamp to a datetime object
    return parse_date_time(date, time)

# -----------------------------------------------------------------

def parse_date_time(date, time):

    """
    This function ...
    :param date:
    :param time:
    :return:
    """

    try: day, month, year = date.split('/')
    except ValueError: return None
    hour, minute, second = time.split(':')
    second, millisecond = second.split('.')

    # Convert the strings to integers
    year = int(year)
    month = int(month)
    day = int(day)
    hour = int(hour)
    minute = int(minute)
    second = int(second)
    millisecond = int(millisecond)

    microsecond = millisecond * 1000

    # Create and return a datetime object
    return datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=second, microsecond=microsecond)

# -----------------------------------------------------------------

def now():

    """
    This function ...
    :return:
    """

    return datetime.now()

# -----------------------------------------------------------------

def timestamp():

    """
    This function generates a timestamp (as a string) from a datetime object
    :return:
    """

    # Return a timestamp accurate up to the millisecond
    return datetime.now().strftime("%d/%m/%Y %H:%M:%S.%f")[:-3]

# -----------------------------------------------------------------

def filename_timestamp():

    """
    This function ...
    :return:
    """

    return datetime.now().strftime("%Y-%m-%d--%H-%M-%S-%f")

# -----------------------------------------------------------------

def day_and_time_as_string():

    """
    This function ...
    :return:
    """

    return datetime.now().strftime("%A, %d. %B %Y %I:%M%p")

# -----------------------------------------------------------------

def is_unique_name(name):

    """
    This function ...
    :param name: 
    :return: 
    """

    try:
        name, time = get_name_and_time_from_unique_name(name)
        return True
    except Exception: return False

# -----------------------------------------------------------------

def get_name_and_time_from_unique_name(name):

    """
    This function ...
    :param name: 
    :return: 
    """

    splitted = name.split("_")
    if len(splitted) > 2:
        first = "_".join(splitted[:-1])
        first = first.strip("_") # remove leading or trailing
        last = splitted[-1]
        return first, time_from_timestamp(last)
    elif len(splitted) == 2: return splitted[0], time_from_timestamp(splitted[1])
    else: raise ValueError("Not a timestamped unique name")

# -----------------------------------------------------------------

def get_name_from_unique_name(name):

    """
    This function ...
    :param name: 
    :return: 
    """

    splitted = name.split("_")
    if len(splitted) > 2: return "_".join(splitted[:-1])
    elif len(splitted) == 2: return splitted[0]
    else: raise ValueError("Not a timestamped unique name")

# -----------------------------------------------------------------

def get_time_from_unique_name(name):

    """
    This function ...
    :param name: 
    :return: 
    """

    splitted = name.split("_")
    last = splitted[-1]
    return time_from_timestamp(last)

# -----------------------------------------------------------------

def time_from_timestamp(stamp):

    """
    This function ...
    :param stamp: 
    :return: 
    """

    # Millisecond precision
    if len(stamp) == 24: stamp = stamp + "000"
    # Microsecond precision
    elif len(stamp) == 27: pass
    # Else
    else: raise ValueError("The timestamp '" + stamp + "' could not be recognized")

    # Parse
    return datetime.strptime(stamp, "%Y-%m-%d--%H-%M-%S-%f")

# -----------------------------------------------------------------

def unique_name(name=None, separator="_", precision="milli"):

    """
    This function ...
    :param name:
    :param separator:
    :param precision:
    :return:
    """

    if precision == "milli" or precision == "millisecond": ndigits = 3
    elif precision == "micro" or precision == "microsecond": ndigits = 0
    else: raise ValueError("Invalid precision")

    # Return the timestamped name
    if name is None: return strip_last_digits(filename_timestamp(), ndigits)
    else: return name + separator + strip_last_digits(filename_timestamp(), ndigits)

# -----------------------------------------------------------------

def strip_last_digits(string, ndigits):

    """
    This function ...
    :param string:
    :param ndigits:
    :return:
    """

    if ndigits == 0: return string
    else: return string[:-ndigits]

# -----------------------------------------------------------------

def pretty_time():

    """
    This function ...
    :return:
    """

    return _time.strftime("%A %B %dth %Y %H:%M:%S", _time.localtime())

# -----------------------------------------------------------------

def pretty_lookback_time(time):

    """
    Get a datetime object or a int() Epoch timestamp and return a
    pretty string like 'an hour ago', 'Yesterday', '3 months ago',
    'just now', etc
    """

    # Get time difference
    if type(time) is int: diff = now() - datetime.fromtimestamp(time)
    elif isinstance(time, datetime): diff = now() - time
    else: raise ValueError("Invalid time object")

    second_diff = diff.seconds
    day_diff = diff.days

    if day_diff < 0: return ''

    if day_diff == 0:

        if second_diff < 10:
            return "just now"
        if second_diff < 60:
            return str(second_diff) + " seconds ago"
        if second_diff < 120:
            return "a minute ago"
        if second_diff < 3600:
            return str(second_diff / 60) + " minutes ago"
        if second_diff < 7200:
            return "an hour ago"
        if second_diff < 86400:
            return str(second_diff / 3600) + " hours ago"

    if day_diff == 1: return "Yesterday"
    if day_diff < 7:
        return str(day_diff) + " days ago"
    if day_diff < 31:
        return str(day_diff / 7) + " weeks ago"
    if day_diff < 365:
        return str(day_diff / 30) + " months ago"

    return str(day_diff / 365) + " years ago"

# -----------------------------------------------------------------

def validate_format(string):

    """
    This function ...
    :param date_text:
    :return:
    """

    try: datetime.strptime(string, '%Y-%m-%d')
    except ValueError: raise ValueError("Incorrect date format")

# -----------------------------------------------------------------
