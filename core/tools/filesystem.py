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

def files_in_path(path, recursive=False, ignore_hidden=True, extension=None, contains=None):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    # Initialize a list to contain the paths of the files that are found in the given directory
    file_paths = []

    # List all items in the specified directory
    for item in os.listdir(path):

        # Determine the full path
        item_path = os.path.join(path, item)

        # Ignore hidden files if requested
        if ignore_hidden and item.startswith("."): continue

        # Ignore files with extension different from the one that is specified
        if extension is not None and os.path.splitext(item)[1][1:] != extension: continue

        # Ignore filenames that do not contain a certain string, if specified
        if contains is not None and contains not in item: continue

        # Check if the current item is a file
        if os.path.isfile(item_path): file_paths.append(item_path)

    # Return the list of file paths
    return file_paths

# -----------------------------------------------------------------

def directories_in_path(path, recursive=False, ignore_hidden=True):

    """
    This function ...
    :param path:
    :param recursive:
    :return:
    """

    # Initialize a list to contain the paths of the directories that are found in the given directory
    directory_paths = []

    # List all items in the specified directory
    for item in os.listdir(path):

        # Determine the full path
        item_path = os.path.join(path, item)

        # Ignore hidden directories if requested
        if ignore_hidden and item.startswith("."): continue

        # Check if the current item is a directory
        if os.path.isdir(item_path): directory_paths.append(item_path)

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
