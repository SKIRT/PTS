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
from subprocess import check_output
import urllib
import httplib
from . import filesystem as fs
from .logging import log
from . import archive
from . import progress

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

def download_and_decompress_file(url, path, remove=True, overwrite=False, progress_bar=False):

    """
    This function ...
    :param url:
    :param path:
    :param remove:
    :param overwrite:
    :param progress_bar:
    :return:
    """

    # Check if path is a directory
    if not fs.is_directory(path): raise ValueError("Second argument must be an existing directory")

    # Download the file and decompress
    filepath = download_file(url, path, overwrite=overwrite, progress_bar=progress_bar)
    decompressed_filepath = archive.decompress_file_in_place(filepath, remove=remove)
    return decompressed_filepath

# -----------------------------------------------------------------

def download_file_no_requests(url, path, overwrite=False):

    """
    This function ...
    :param url:
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

def download_file(url, path, overwrite=False, progress_bar=False, stream=False, chunk_size=1024, session=None):

    """
    This function ...
    :param url
    :param path:
    :param overwrite:
    :param progress_bar:
    :param stream
    :param chunk_size:
    :param session:
    :return:
    """

    # Import here to enable this module to be imported with a clean python install
    import requests

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

    # Show progress bar, so stream
    if progress_bar:

        # Request
        if session is None: session = requests.session()
        r = session.get(url, stream=True)

        # Open the local file
        with open(filepath, 'wb') as f:

            total_length = int(r.headers.get('content-length'))
            for chunk in progress.bar(r.iter_content(chunk_size=chunk_size), expected_size=(total_length / chunk_size) + 1):
                if chunk:
                    f.write(chunk)
                    f.flush()

    # User wants streaming
    elif stream:

        # Request
        if session is None: session = requests.session()
        r = session.get(url, stream=True)

        # Open the local file, and load the content in it
        with open(filepath, 'wb') as f:

            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
                    # f.flush() # commented by recommendation from J.F.Sebastian

    # Regular download
    elif session is not None:

        # Request
        if session is None: session = requests.session()
        r = session.get(url)

        # Open the local file, and load the content in it
        with open(filepath, 'wb') as f:
            for chunk in r.iter_content(chunk_size=chunk_size):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)

    # Regular download , no session
    else: urllib.urlretrieve(url, filepath)

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
