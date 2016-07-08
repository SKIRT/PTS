#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.network Network-related functions

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from subprocess import call, check_output
import urllib
from . import filesystem as fs
from .logging import log

# -----------------------------------------------------------------

def download_file(url, path):

    """
    This function ...
    :param url
    :param path:
    :return:
    """

    # Get the name of the file
    filename = fs.name(url)

    # Determine the local path to the file
    filepath = fs.join(path, filename) if fs.is_directory(path) else path

    # Debugging
    log.debug("Downloading '" + filename + "' to '" + path + "' ...")

    # Download
    urllib.urlretrieve(url, filepath)

    # Return the file path
    return filepath

# -----------------------------------------------------------------

def download_files(urls, path):

    """
    This function ...
    :param urls:
    :param path:
    :return:
    """

    paths = []

    count = len(urls)

    # Loop over the urls
    index = 0
    for url in urls:

        filename = fs.name(url)
        filepath = fs.join(path, filename)

        # Debugging
        log.debug("Downloading '" + filename + "' to '" + path + "' ... (" + str(index+1) + " of " + str(count) + ")")

        # Download
        urllib.urlretrieve(url, filepath)

        # If succesful, add the file path to the list
        paths.append(filepath)

        index += 1

    # Return paths
    return paths

# -----------------------------------------------------------------

def dns_ips():

    """
    This function ...
    :return:
    """

    # Alternative: cat /etc/resolv.conf ?

    command = "scutil --dns | grep 'nameserver\[[0-9]*\]' | sort | uniq"
    output = check_output(command, shell=True)

    entries = [entry.strip() for entry in output.split("\n") if entry.strip()]

    return [entry.split(" : ")[1] for entry in entries]

# -----------------------------------------------------------------

def dns_search_domains():

    """
    This function ...
    :return:
    """

    # Alternative: cat /etc/resolv.conf ?

    command = "scutil --dns | grep 'search domain\[[0-9]*\]' | sort | uniq"
    output = check_output(command, shell=True)

    entries = [entry.strip() for entry in output.split("\n") if entry.strip()]

    return [entry.split(" : ")[1] for entry in entries]

# -----------------------------------------------------------------
