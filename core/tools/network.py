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
from urllib2 import urlopen, HTTPError
from . import filesystem as fs
from ..basics.log import log
from . import archive
from . import progress

# -----------------------------------------------------------------

def exists(url, allow_redirection=True, check_open=True):

    """
    This function ...
    :param url:
    :param allow_redirection:
    :param check_open:
    :return:
    """

    from urlparse import urlparse

    p = urlparse(url)
    conn = httplib.HTTPConnection(p.netloc)

    try: conn.request('HEAD', p.path)
    except: return False

    resp = conn.getresponse()
    conn.close()

    #print(resp.status)
    #print(resp.read()) # doesn't show anything

    # successful
    # OK = 200
    # CREATED = 201
    # ACCEPTED = 202
    # NON_AUTHORITATIVE_INFORMATION = 203
    # NO_CONTENT = 204
    # RESET_CONTENT = 205
    # PARTIAL_CONTENT = 206
    # MULTI_STATUS = 207
    # IM_USED = 226
    #
    # # redirection
    # MULTIPLE_CHOICES = 300
    # MOVED_PERMANENTLY = 301
    # FOUND = 302
    # SEE_OTHER = 303
    # NOT_MODIFIED = 304
    # USE_PROXY = 305
    # TEMPORARY_REDIRECT = 307
    #
    # # client error
    # BAD_REQUEST = 400
    # UNAUTHORIZED = 401
    # PAYMENT_REQUIRED = 402
    # FORBIDDEN = 403
    # NOT_FOUND = 404
    # METHOD_NOT_ALLOWED = 405
    # NOT_ACCEPTABLE = 406
    # PROXY_AUTHENTICATION_REQUIRED = 407
    # REQUEST_TIMEOUT = 408
    # CONFLICT = 409
    # GONE = 410
    # LENGTH_REQUIRED = 411
    # PRECONDITION_FAILED = 412
    # REQUEST_ENTITY_TOO_LARGE = 413
    # REQUEST_URI_TOO_LONG = 414
    # UNSUPPORTED_MEDIA_TYPE = 415
    # REQUESTED_RANGE_NOT_SATISFIABLE = 416
    # EXPECTATION_FAILED = 417
    # UNPROCESSABLE_ENTITY = 422
    # LOCKED = 423
    # FAILED_DEPENDENCY = 424
    # UPGRADE_REQUIRED = 426
    #
    # # server error
    # INTERNAL_SERVER_ERROR = 500
    # NOT_IMPLEMENTED = 501
    # BAD_GATEWAY = 502
    # SERVICE_UNAVAILABLE = 503
    # GATEWAY_TIMEOUT = 504
    # HTTP_VERSION_NOT_SUPPORTED = 505
    # INSUFFICIENT_STORAGE = 507
    # NOT_EXTENDED = 510

    # Allow redirection?
    if allow_redirection:
        if check_open: return resp.status < 400 and can_be_opened(url)
        else: return resp.status < 400
    else:
        if check_open: return resp.status < 300 and can_be_opened(url)
        else: return resp.status < 300

# -----------------------------------------------------------------

def can_be_opened(url):

    """
    This function ...
    :param url:
    :return:
    """

    try:
        data = urlopen(url)
        return True
    except HTTPError as e:
        message = str(e)
        if "404" in message: return False # Not existing
        elif "403" in message: return False # forbidden
        else: raise ValueError("Unknown error: '" + message + "'")

# -----------------------------------------------------------------

def is_forbidden(url):

    """
    Thisf unction ...
    :param url:
    :return:
    """

    try:
        data = urlopen(url)
        return False
    except HTTPError as e:
        message = str(e)
        if "403" in message: return True
        else: raise ValueError("Unknown error: '" + message + "'")

# -----------------------------------------------------------------

def is_available(url):

    """
    This function ...
    :param url:
    :return:
    """

    return exists(url) and not is_forbidden(url)

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

def download_and_decompress_directory(url, path, remove=True, overwrite=False, progress_bar=False, into_root=False):

    """
    This function ...
    :param url:
    :param path:
    :param remove:
    :param overwrite:
    :param progress_bar:
    :param into_root:
    :return:
    """

    # Check if path is a directory
    if not fs.is_directory(path): raise ValueError("Second argument must be an existing directory")

    # Download the file and decompress into directory
    filepath = download_file(url, path, overwrite=overwrite, progress_bar=progress_bar)
    decompressed_path = archive.decompress_directory_in_place(filepath, remove=remove, into_root=into_root)
    return decompressed_path

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

def download_file(url, path, new_name=None, overwrite=False, progress_bar=False, stream=False, chunk_size=1024, session=None):

    """
    This function ...
    :param url
    :param path:
    :param new_name:
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
    if new_name is not None: filename = new_name
    else: filename = fs.name(url)

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
        r = session.get(url, stream=True, timeout=(60,600)) # (connect timeout, read timeout)

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
        r = session.get(url, stream=True, timeout=(60,600)) # (connect timeout, read timeout)

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
        r = session.get(url, timeout=(60,600)) # (connect timeout, read timeout)

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

def download_and_decompress_files(urls, path, remove=True, overwrite=False, info=None):

    """
    This function ...
    :param urls:
    :param path:
    :param remove:
    :param overwrite:
    :param info:
    :return:
    """

    # Debugging
    log.debug("Downloading the files to '" + path + "' ...")

    # Download the files
    paths = download_files(urls, path, overwrite=overwrite, info=info)

    # Debugging
    log.debug("Decompressing the archived files ...")

    # Decompress the files and remove the originals
    #new_paths = archive.decompress_files(paths, remove=remove)

    new_paths = []
    compressed_paths = []
    for path in paths:

        if archive.is_archive(path): compressed_paths.append(path)
        else: new_paths.append(path)

    # Decompress
    new_paths += archive.decompress_files(compressed_paths, remove=remove)

    # Return the paths of the decompressed files
    return new_paths

# -----------------------------------------------------------------

def download_files(urls, path, overwrite=False, info=None):

    """
    This function ...
    :param urls:
    :param path:
    :param overwrite:
    :param info:
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
        if info is not None: log.debug("Downloading '" + filename + "' to '" + path + "' ... (" + str(index+1) + " of " + str(count) + " " + info + ")")
        else: log.debug("Downloading '" + filename + "' to '" + path + "' ... (" + str(index+1) + " of " + str(count) + ")")
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
