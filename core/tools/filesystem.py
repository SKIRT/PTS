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
import platform
import subprocess

# Import the relevant PTS classes and modules
from . import time

# -----------------------------------------------------------------

def cwd():

    """
    This function returns the path to the current working directory
    :return:
    """

    return os.getcwd()

# -----------------------------------------------------------------

def change_cwd(path):

    """
    This function ...
    :param path:
    :return:
    """

    previous = cwd()
    os.chdir(path)
    return previous

# -----------------------------------------------------------------

def home():

    """
    This function returns the full path to the home directory
    :return:
    """

    return os.path.expanduser('~')

# -----------------------------------------------------------------

def to_home_directory():

    """
    This function ...
    :param path:
    :return:
    """

    change_cwd(home())

# -----------------------------------------------------------------

def touch(path):

    """
    This function replicates the behaviour of 'touch' in terminal
    :param path:
    :return:
    """

    with open(path, 'a'):
        os.utime(path, None)

# -----------------------------------------------------------------

def open_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Check if existing
    if not is_file(path): raise ValueError("The file '" + path + "' does not exist")

    if platform.system() == "Darwin": subprocess.Popen(["open", path])
    else: subprocess.Popen(["xdg-open", path])

# -----------------------------------------------------------------

def open_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Check if existing
    if not is_directory(path): raise ValueError("The file '" + path + "' does not exist")

    if platform.system() == "Darwin": subprocess.Popen(["open", path])
    else: subprocess.Popen(["xdg-open", path])

# -----------------------------------------------------------------

def desktop():

    """
    This function returns the full path to the desktop directory
    :return:
    """

    return os.path.expanduser("~/Desktop") # On Mac OS X

# -----------------------------------------------------------------

def documents():

    """
    This function returns the full path to the documents directory
    :return:
    """

    return os.path.expanduser("~/Documents") # On Mac OS X

# -----------------------------------------------------------------

def absolute_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.abspath(os.path.expanduser(path))

# -----------------------------------------------------------------

def absolute_or_in(path, in_path):

    """
    This function ...
    :param path:
    :param in_path:
    :return:
    """

    if os.path.isabs(path): return path
    else: return join(in_path, path)

# -----------------------------------------------------------------

def is_subdirectory(path, parent_path):

    """
    This function returns whether
    :param path:
    :param parent_path:
    :return:
    """

    if not is_directory(path): raise ValueError("Not a directory: " + path)

    path = os.path.realpath(path)
    parent_path = os.path.realpath(parent_path)
    return path.startswith(parent_path)

# -----------------------------------------------------------------

def contains_file(directory, filename):

    """
    This function ...
    :param directory:
    :param filename:
    :return:
    """

    return is_file(join(directory, filename))

# -----------------------------------------------------------------

def contains_files(directory, filenames):

    """
    This function ...
    :param directory:
    :param filenames:
    :return:
    """

    # Loop over the filenames
    for filename in filenames:

        filepath = join(directory, filename)
        if not is_file(filepath): return False

    return True

# -----------------------------------------------------------------

def join(*args):

    """
    This function ...
    :param args:
    :return:
    """

    return os.path.join(*args)

# -----------------------------------------------------------------

def file_or_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    if is_file(path): return "file"
    elif is_directory(path): return "directory"
    else: return None

# -----------------------------------------------------------------

def is_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.isfile(path)

# -----------------------------------------------------------------

def has_file(directory, filename):

    """
    This function ...
    :param directory:
    :param filename:
    :return:
    """

    return is_file(join(directory, filename))

# -----------------------------------------------------------------

def is_empty(directory, ignore_hidden=True):

    """
    This function ...
    :param directory:
    :param ignore_hidden:
    :return:
    """

    items = os.listdir(directory)
    if ignore_hidden: items = [item for item in items if not item.startswith(".")]
    return len(items) == 0

# -----------------------------------------------------------------

def is_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.isdir(path)

# -----------------------------------------------------------------

def is_mount_point(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.ismount(path)

# -----------------------------------------------------------------

def directory_of(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.dirname(path)

# -----------------------------------------------------------------

def directory_and_name(path):

    """
    This function ...
    :param path:
    :return:
    """

    root, name = os.path.split(path)
    root = os.path.realpath(root)

    return root, name

# -----------------------------------------------------------------

def name(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.basename(path)

# -----------------------------------------------------------------

def strip_extension(name):

    """
    This function ...
    :param name:
    :return:
    """

    return os.path.splitext(name)[0]

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

def create_directory_in(base_path, name):

    """
    This function ...
    :param base_path:
    :param name:
    :return:
    """

    directory_path = join(base_path, name)
    create_directory(directory_path)
    return directory_path

# -----------------------------------------------------------------

def create_directories(*paths, **kwargs):
    
    """
    This function ...
    :param paths:
    :param kwargs:
    """
    
    # Loop over the different paths in the list
    for path in paths: create_directory(path, kwargs.pop("recursive", False))

# -----------------------------------------------------------------

def create_directories_in(base_path, names):

    """
    This function ...
    :param base_path:
    :param names:
    :return:
    """

    paths = []
    for name in names: paths.append(create_directory_in(base_path, name))
    return paths

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

def clear_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    for file_path in files_in_path(path): remove_file(file_path)
    for directory_path in directories_in_path(path): remove_directory(directory_path)

# -----------------------------------------------------------------

def remove_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    shutil.rmtree(path)

# -----------------------------------------------------------------

def remove_directories(paths):

    """
    This function ...
    :param paths:
    :return:
    """

    for path in paths: remove_directory(path)

# -----------------------------------------------------------------

def remove_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    os.remove(path)

# -----------------------------------------------------------------

def remove_files(paths):

    """
    This function ...
    :param paths:
    :return:
    """

    for path in paths: remove_file(path)

# -----------------------------------------------------------------

def size(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.getsize(path)

# -----------------------------------------------------------------

def ls(path=None):

    """
    This function ...
    :param path:
    :return:
    """

    if path is None: path = cwd()
    return os.listdir(path)

# -----------------------------------------------------------------

def files_in_path(path=None, recursive=False, ignore_hidden=True, extension=None, contains=None, not_contains=None,
                  extensions=False, returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None, sort=None):

    """
    This function ...
    :param path:
    :param recursive:
    :param ignore_hidden:
    :param extension:
    :param contains:
    :param not_contains:
    :param extensions:
    :param returns: a string ("path", "name", or "directory") OR a list [], with elements equal to "path", "name" or "directory" (e.g. [path, name] or [name, directory])
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :param sort: a function which determines how the files should be sorted based on their filename. Hidden items (starting with .) are placed first.
    :return:
    """

    # If no path is given, use the current working directory
    if path is None: path = os.getcwd()

    # Determine absolute path
    path = absolute_path(path)

    # Initialize a list to contain the paths of the files that are found in the given directory
    file_paths = []

    # Get the list of items
    if recursive: items = [join(dp, f) for dp, dn, fn in os.walk(path) for f in fn]
    else: items = [join(path, item) for item in os.listdir(path)]

    # If the files have to be sorted on their name
    if sort is not None:
        #sort_function = lambda x: sort(strip_extension(name(x)))
        def sort_function(x):
            item_name = strip_extension(name(x))
            if item_name.startswith("."): value = 0
            else:
                try: value = sort(item_name)
                except ValueError: value = 0
            return value
        items.sort(key=sort_function)

    # Loop over all items; get files that match the specified conditions
    for item_path in items:

        # Determine the name of the item
        item = name(item_path)

        # Get the path to the containing directory
        directory = directory_of(item_path)

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

        # Ignore filenames that do not match the exact filename, if specified
        if exact_name is not None and exact_name != item_name: continue

        # If the filename matches the 'exact not name', skip it
        if exact_not_name is not None and exact_not_name == item_name: continue

        # Ignore filenames that do not start or end with the specified strings
        if startswith is not None and not item_name.startswith(startswith): continue
        if endswith is not None and not item_name.endswith(endswith): continue

        # Check if the current item is a file; if not skip it
        if not os.path.isfile(item_path): continue

        # Create the return value
        if isinstance(returns, basestring):

            if returns == "path": thing = item_path
            elif returns == "name": thing = item_name + "." + item_extension if extensions else item_name
            elif returns == "directory": thing = directory
            else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        else: # Assume 'returns' is a list

            thing = []

            for return_value in returns:

                if return_value == "path": thing.append(item_path)
                elif return_value == "name": thing.append(item_name + "." + item_extension if extensions else item_name)
                elif return_value == "directory": thing.append(directory)
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        file_paths.append(thing)

    # Return the list of file paths
    return file_paths

# -----------------------------------------------------------------

def directories_in_path(path=None, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                        returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None, sort=None):

    """
    This function ...
    :param path:
    :param recursive:
    :param ignore_hidden:
    :param contains:
    :param not_contains:
    :param returns: a string ("path", "name", or "directory") OR a list [], with elements equal to "path", "name" or "directory" (e.g. [path, name] or [name, directory])
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :param sort: a function which determines how the directories should be sorted based on their name. Hidden items (starting with .) are placed first.
    :return:
    """

    # If no path is given, use the current working directory
    if path is None: path = os.getcwd()

    # Determine absolute path
    path = absolute_path(path)

    # Initialize a list to contain the paths of the directories that are found in the given directory
    directory_paths = []

    # Get the list of items
    if recursive: items = [join(dp, d) for dp, dn, fn in os.walk(path) for d in dn]
    else: items = [join(path, item) for item in os.listdir(path)]

    # If the directories have to be sorted on their name
    if sort is not None:
        #sort_function = lambda x: sort(strip_extension(name(x)))
        def sort_function(x):
            item_name = name(x)
            if item_name.startswith("."): value = 0
            else:
                try: value = sort(item_name)
                except ValueError: value = 0
            return value
        items.sort(key=sort_function)

    # List all items in the specified directory
    for item_path in items:

        # Get the name of the directory
        item = name(item_path)

        # Get the path to the containing directory
        directory = directory_of(item_path)

        # Ignore hidden directories if requested
        if ignore_hidden and item.startswith("."): continue

        # Ignore names that do not contain a certain string, if specified
        if contains is not None and contains not in item: continue

        # Ignore names that do contain a certain string that it should not contain, if specified
        if not_contains is not None and not_contains in item: continue

        # If the directory name does not match the exact name, skip it
        if exact_name is not None and exact_name != item: continue

        # If the directory name matches the 'exact not name', skip it
        if exact_not_name is not None:

            if isinstance(exact_not_name, basestring):
                if exact_not_name == item: continue
            elif isinstance(exact_not_name, list):
                if item in exact_not_name: continue
            else: raise ValueError("Invalid option for 'exact_not_name': must be string or list")

        # Ignore directory names that do not start or end with the specified strings
        if startswith is not None and not item.startswith(startswith): continue
        if endswith is not None and not item.endswith(endswith): continue

        # Check if the current item is a directory; if not skip it
        if not os.path.isdir(item_path): continue

        # Create the return value
        if isinstance(returns, basestring):

            if returns == "path": thing = item_path
            elif returns == "name": thing = item
            elif returns == "directory": thing = directory
            else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        else: # Assume 'returns' is a list

            thing = []

            for return_value in returns:

                if return_value == "path": thing.append(item_path)
                elif return_value == "name": thing.append(item)
                elif return_value == "directory": thing.append(directory)
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        directory_paths.append(thing)

    # Return the list of directory paths
    return directory_paths

# -----------------------------------------------------------------

def rename_file(directory, original_name, new_name):

    """
    This function ...
    :param directory:
    :param original_name:
    :param new_name:
    :return:
    """

    # Check whether the original file exists
    original_path = join(directory, original_name)
    if not is_file(original_path): raise ValueError("File '" + original_path + "' does not exist")

    # Determine new file path
    new_path = join(directory, new_name)

    # Rename the file
    os.rename(original_path, new_path)

# -----------------------------------------------------------------

def copy_file(file_path, directory_path, new_name=None):

    """
    This function ...
    :param file_path:
    :param directory_path:
    :param new_name:
    :return:
    """

    if new_name is not None: destination = join(directory_path, new_name)
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
