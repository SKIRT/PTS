#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.filesystem Provides useful functions for manipulating the local file system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import shutil

# Import the relevant PTS classes and modules
from . import time

# -----------------------------------------------------------------

def join(path_a, path_b):

    """
    This function ...
    :param path_a:
    :param path_b:
    :return:
    """

    return os.path.join(path_a, path_b)

# -----------------------------------------------------------------

def is_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.isfile(path)

# -----------------------------------------------------------------

def is_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.isdir(path)

# -----------------------------------------------------------------

def create_directory(path, recursive=False):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    # Check whether the directory does not exist yet
    if not os.path.isdir(path):

        # Create the directory, recursively or not
        if recursive: os.makedirs(path)
        else: os.mkdir(path)

# -----------------------------------------------------------------

def create_directories(paths, recursive=False):
    
    """
    This function ...
    :param paths:
    :param recursive:
    """
    
    # Loop over the different paths in the list
    for path in paths: create_directory(path)

# -----------------------------------------------------------------

def create_temporary_directory(prefix=None):

    """
    This function ...
    :param prefix:
    :return:
    """

    # Add a timestamp to the prefix
    name = time.unique_name(prefix) if prefix is not None else time.unique_name("", "")

    # Set the path to the temporary directory
    path = os.path.join(os.getcwd(), name)

    # Create the directory
    create_directory(path)

    # Return the directory path
    return path

# -----------------------------------------------------------------

def remove_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    shutil.rmtree(path)

# -----------------------------------------------------------------

def remove_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    os.remove(path)

# -----------------------------------------------------------------

def files_in_path(path=None, recursive=False, ignore_hidden=True, extension=None, contains=None, not_contains=None, names=False, extensions=False):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    if path is None: path = os.getcwd()

    # Initialize a list to contain the paths of the files that are found in the given directory
    file_paths = []

    # Get the list of items
    if recursive: items = [os.path.join(dp, f) for dp, dn, fn in os.walk(path) for f in fn]
    else: items = os.listdir(path)

    # Loop over all items; get files that match the specified conditions
    for item in items:

        # Determine the full path
        item_path = os.path.join(path, item)

        # Get the file name and extension
        item_name = os.path.splitext(item)[0]
        item_extension = os.path.splitext(item)[1][1:]

        # Ignore hidden files if requested
        if ignore_hidden and item.startswith("."): continue

        # Ignore files with extension different from the one that is specified
        if extension is not None and item_extension != extension: continue

        # Ignore filenames that do not contain a certain string, if specified
        if contains is not None and contains not in item_name: continue

        # Ignore filenames that do contain a certain string that it should not contain, if specified
        if not_contains is not None and not_contains in item_name: continue

        # Check if the current item is a file; if not skip it
        if not os.path.isfile(item_path): continue

        # Add the relevant info to the list
        thing = [item_path]
        if names: thing.append(item_name)
        if extensions: thing.append(item_extension)
        file_paths.append(thing if len(thing) > 1 else thing[0])

    # Return the list of file paths
    return file_paths

# -----------------------------------------------------------------

def directories_in_path(path=None, recursive=False, ignore_hidden=True, names=False):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    if path is None: path = os.getcwd()

    # Initialize a list to contain the paths of the directories that are found in the given directory
    directory_paths = []

    # Get the list of items
    if recursive: items = [os.path.join(dp, f) for dp, dn, fn in os.walk(path) for f in fn]
    else: items = os.listdir(path)

    # List all items in the specified directory
    for item in items:

        # Determine the full path
        item_path = os.path.join(path, item)

        # Ignore hidden directories if requested
        if ignore_hidden and item.startswith("."): continue

        # Check if the current item is a directory; if not skip it
        if not os.path.isdir(item_path): continue

        # Add the directory path to the list
        thing = [item_path]
        if names: thing.append(item)
        directory_paths.append(thing if len(thing) > 1 else thing[0])

    # Return the list of directory paths
    return directory_paths

# -----------------------------------------------------------------

def copy_file(file_path, directory_path, new_name=None):

    """
    This function ...
    :return:
    """

    if new_name is not None: destination = os.path.join(directory_path, new_name)
    else: destination = directory_path

    shutil.copy(file_path, destination)

# -----------------------------------------------------------------

def copy_files(file_paths, directory_path):

    """
    This function ...
    :param file_paths:
    :param directory_path:
    :return:
    """

    for file_path in file_paths: copy_file(file_path, directory_path)

# -----------------------------------------------------------------
