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
import httplib
from . import filesystem as fs
from .logging import log
from . import archive

# -----------------------------------------------------------------

def exists(url):

    """
    This function ...
    :param url:
    :return:
    """

    from urlparse import urlparse

    p = urlparse(url)
    conn = httplib.HTTPConnection(p.netloc)
    conn.request('HEAD', p.path)
    resp = conn.getresponse()
    return resp.status < 400

# -----------------------------------------------------------------

def download_and_decompress_file(url, path, remove=True, overwrite=False):

    """
    This function ...
    :param url:
    :param path:
    :param remove:
    :param overwrite:
    :return:
    """

    # Check if path is a directory
    if not fs.is_directory(path): raise ValueError("Second argument must be an existing directory")

    # Download the file and decompress
    filepath = download_file(url, path, overwrite=overwrite)
    decompressed_filepath = archive.decompress_file_in_place(filepath, remove=remove)
    return decompressed_filepath

# -----------------------------------------------------------------

def download_file(url, path, overwrite=False):

    """
    This function ...
    :param url
    :param path:
    :param overwrite:
    :return:
    """

    # Get the name of the file
    filename = fs.name(url)

    # Determine the local path to the file
    filepath = fs.join(path, filename) if fs.is_directory(path) else path

    # Check filepath
    if fs.is_file(filepath):
        if overwrite: fs.remove_file(filepath)
        else: raise IOError("File is already present: " + filepath)

    # Debugging
    log.debug("Downloading '" + filename + "' to '" + path + "' ...")
    log.debug("URL: " + url)

    # Download
    urllib.urlretrieve(url, filepath)

    # Return the file path
    return filepath

# -----------------------------------------------------------------

def download_and_decompress_files(urls, path, remove=True, overwrite=False):

    """
    This function ...
    :param urls:
    :param path:
    :param remove:
    :return:
    """

    # Debugging
    log.debug("Downloading the files to '" + path + "' ...")

    # Download the files
    paths = download_files(urls, path, overwrite=overwrite)

    # Debugging
    log.debug("Decompressing the files ...")

    # Decompress the files and remove the originals
    new_paths = archive.decompress_files(paths, remove=remove)

    # Return the paths of the decompressed files
    return new_paths

# -----------------------------------------------------------------

def download_files(urls, path, overwrite=False):

    """
    This function ...
    :param urls:
    :param path:
    :param overwrite:
    :return:
    """

    paths = []

    count = len(urls)

    # Loop over the urls
    index = 0
    for url in urls:

        filename = fs.name(url)
        filepath = fs.join(path, filename)

        # Check if the file is present
        if fs.is_file(filepath):
            if overwrite: fs.remove_file(filepath)
            else: raise IOError("File is already present: " + filepath)

        # Debugging
        log.debug("Downloading '" + filename + "' to '" + path + "' ... (" + str(index+1) + " of " + str(count) + ")")
        log.debug("URL: " + url)

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
