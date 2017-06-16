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
import sys
import shutil
import platform
import subprocess
import datetime

# Import the relevant PTS classes and modules
from . import time
from . import types
from . import sequences

# -----------------------------------------------------------------

# Check if Python 2 or higher
if (sys.version_info > (3, 0)): PYTHON_2 = False # Python 3 or higher
else: PYTHON_2 = True # Python 2

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

def open_file(path, wait=False):

    """
    This function ...
    :param path:
    :param wait:
    :return:
    """

    # Check if existing
    if not is_file(path): raise ValueError("The file '" + path + "' does not exist")

    # Determine command
    if platform.system() == "Darwin": command = ["open", path]
    else: command = ["xdg-open", path]

    # Call the command
    if wait:
        subprocess.call(command, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
        time.wait(10)
    else: subprocess.Popen(command)

# -----------------------------------------------------------------

def open_directory(path, wait=False):

    """
    This function ...
    :param path:
    :param wait:
    :return:
    """

    # Check if existing
    if not is_directory(path): raise ValueError("The directory '" + path + "' does not exist")

    # Determine command
    if platform.system() == "Darwin": command = ["open", path]
    else: command = ["xdg-open", path]

    # Call the command
    if wait:
        subprocess.call(command, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
        time.wait(10)
    else: subprocess.Popen(command)

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

def create_absolute_or_in(path, in_path):

    """
    This function ...
    :param path:
    :param in_path:
    :return:
    """

    absolute_path = absolute_or_in(path, in_path)

    if is_file(absolute_path): return absolute_path
    elif is_directory(absolute_path): return absolute_path
    elif has_extension(absolute_path): return absolute_path # file
    else: # directory
        create_directory(absolute_path)
        return absolute_path

# -----------------------------------------------------------------

def absolute_or_in_cwd(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    return absolute_or_in(path, cwd())

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

def contains_directory(directory, dirname):

    """
    This function ...
    :param directory:
    :param dirname:
    :return:
    """

    return is_directory(join(directory, dirname))

# -----------------------------------------------------------------

def contains_files(directory, filenames=None):

    """
    This function ...
    :param directory:
    :param filenames:
    :return:
    """

    if filenames is not None:

        # Loop over the filenames
        for filename in filenames:

            filepath = join(directory, filename)
            if not is_file(filepath): return False

        return True

    else: return len(files_in_path(directory)) > 0

# -----------------------------------------------------------------

def contains_directories(directory, dirnames=None):

    """
    This function ...
    :param directory:
    :param dirnames:
    :return:
    """

    if dirnames is not None:

        # Loop over the dirnames
        for dirname in dirnames:

            dirpath = join(directory, dirname)
            if not is_directory(dirpath): return False

        return True

    else: return len(directories_in_path(directory)) > 0

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

def has_extension(path):

    """
    This function ...
    :param path:
    :return:
    """

    filename = name(path)
    return "." in filename

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

def is_empty(directory, ignore_hidden=True, besides=None):

    """
    This function ...
    :param directory:
    :param ignore_hidden:
    :param besides:
    :return:
    """

    items = os.listdir(directory)
    if ignore_hidden: items = [item for item in items if not item.startswith(".")]

    if besides is not None:
        return len(items) == 0 or (len(items) == 1 and items[0] == besides)
    else: return len(items) == 0

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

def strip_extension(name_or_path, double=False):

    """
    This function ...
    :param name_or_path:
    :param double:
    :return:
    """

    name_or_path = os.path.splitext(name_or_path)[0]
    if double: name_or_path = os.path.splitext(name_or_path)[0]
    return name_or_path

# -----------------------------------------------------------------

def get_extension(name_or_path, double=False):

    """
    This function ...
    :param name_or_path:
    :param double:
    :return:
    """

    without_extension = strip_extension(name_or_path, double=double)
    if name_or_path == without_extension: return ""
    else: return name_or_path.split(without_extension + ".")[1]

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

def remove_directories_and_files(paths):

    """
    This function ...
    :param paths:
    :return:
    """

    # Loop over the paths
    for path in paths:

        if is_file(path): remove_file(path)
        elif is_directory(path): remove_directory(path)
        else: raise ValueError("Not an existing file or directory: " + path)

# -----------------------------------------------------------------

def replace_file(path, new_path):

    """
    This function ...
    :param path:
    :param new_path:
    :return:
    """

    remove_file(path)
    copy_file(new_path, directory_of(path), new_name=name(path))

# -----------------------------------------------------------------

def replace_directory(path, new_path):

    """
    This function ...
    :param path:
    :param new_path:
    :return:
    """

    remove_directory(path)
    copy_directory(new_path, directory_of(path), new_name=name(path))

# -----------------------------------------------------------------

def directory_size(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..units.parsing import parse_unit as u
    total_size = 0. * u("byte")
    for dirpath, dirnames, filenames in os.walk(path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += file_size(fp)

    return total_size

# -----------------------------------------------------------------

def file_size(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..units.parsing import parse_unit as u
    return os.path.getsize(path) * u("byte")

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

def extension_of(path, filename):

    """
    This function ....
    :param path:
    :param filename:
    :return:
    """

    paths = files_in_path(path, exact_name=filename)
    if len(paths) == 0: raise IOError("No such file: '" + filename + ''" in '" + path + "'")
    elif len(paths) > 1: raise IOError("Multiple files with the name '" + filename + "' in '" + path + "'")
    else: return get_extension(paths[0])

# -----------------------------------------------------------------

def walk_level(some_dir, level=1):
    
    """
    This function ...
    :param some_dir: 
    :param level: 
    :return: 
    """
    
    some_dir = some_dir.rstrip(os.path.sep)
    
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]

# -----------------------------------------------------------------

def files_in_path(path=None, recursive=False, ignore_hidden=True, extension=None, contains=None, not_contains=None,
                  extensions=False, returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None,
                  sort=None, contains_operator="OR", recursion_level=None, unpack=False):

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
    :param contains_operator: relevant for when 'contains' is specified as a sequence (should they be all contained or only at least one?)
    :param recursion_level:
    :param unpack:
    :return:
    """

    # If no path is given, use the current working directory
    if path is None: path = os.getcwd()

    # Determine absolute path
    path = absolute_path(path)

    # Initialize a list to contain the paths of the files that are found in the given directory
    file_paths = []

    # Get the list of items
    if recursive:
        if recursion_level is None: items = [join(dp, f) for dp, dn, fn in os.walk(path) for f in fn]
        else: items = [join(dp, f) for dp, dn, fn in walk_level(path, recursion_level) for f in fn]
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

    if returns == "dict":
        if unpack: raise ValueError("Cannot do unpacking when return value is supposed to be 'dict'")
        returns = ["name", "path"]
        return_dict = True
    else: return_dict = False

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
        if extension is not None:
            if types.is_string_type(extension):
                if item_extension != extension: continue
            elif types.is_sequence(extension):
                if item_extension not in extension: continue
            else: raise ValueError("Unknown type for 'extension': " + str(extension))

        # Ignore filenames that do not contain a certain string, if specified
        if contains is not None:
            if types.is_string_type(contains):
                if contains not in item_name: continue
            elif types.is_sequence(contains):
                if contains_operator == "OR":
                    if not sequences.any_in(contains, item_name): continue
                elif contains_operator == "AND":
                    if not sequences.all_in(contains, item_name): continue
                else: raise ValueError("Invalid contains operator")
            else: raise ValueError("contains should be string or sequence")

        # Ignore filenames that do contain a certain string that it should not contain, if specified
        if not_contains is not None:
            if types.is_string_type(not_contains):
                if not_contains in item_name: continue
            elif types.is_sequence(not_contains):
                if sequences.any_in(not_contains, item_name): continue
            else: raise ValueError("contains should be string or sequence")

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
        if types.is_string_type(returns):

            if unpack: raise ValueError("Cannot unpack with only one return value")

            if returns == "path": thing = item_path
            elif returns == "name": thing = item_name + "." + item_extension if extensions else item_name
            elif returns == "directory": thing = directory
            else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        else: # Assume 'returns' is a list or tuple

            thing = []

            for return_value in returns:

                if return_value == "path": thing.append(item_path)
                elif return_value == "name": thing.append(item_name + "." + item_extension if extensions else item_name)
                elif return_value == "directory": thing.append(directory)
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        file_paths.append(thing)

    # Return the list of file paths
    if return_dict: return dict(file_paths) # file_paths is list of tuples
    elif unpack: return sequences.unpack(file_paths, default_size=len(returns))
    else: return file_paths

# -----------------------------------------------------------------

def find_file_in_path(path, recursive=False, ignore_hidden=True, extension=None, contains=None, not_contains=None,
                      exact_name=None, exact_not_name=None, startswith=None, endswith=None):

    """
    This function ...
    :param path: 
    :param recursive:
    :param ignore_hidden:
    :param extension:
    :param contains:
    :param not_contains:
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :return: 
    """

    # Get paths
    paths = files_in_path(path, recursive=recursive, ignore_hidden=ignore_hidden, extension=extension, contains=contains,
                          not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name,
                          startswith=startswith, endswith=endswith)
    if len(paths) == 1: return paths[0]
    elif len(paths) == 0: raise ValueError("Not found")
    else: raise ValueError("Multiple files found")

# -----------------------------------------------------------------

def directories_in_path(path=None, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                        returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None, sort=None,
                        recursion_level=None, unpack=False):

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
    :param recursion_level:
    :param unpack:
    :return:
    """

    # If no path is given, use the current working directory
    if path is None: path = os.getcwd()

    # Determine absolute path
    path = absolute_path(path)

    # Initialize a list to contain the paths of the directories that are found in the given directory
    directory_paths = []

    # Get the list of items
    if recursive:
        if recursion_level is None: items = [join(dp, d) for dp, dn, fn in os.walk(path) for d in dn]
        else: items = [join(dp, d) for dp, dn, fn in walk_level(path, recursion_level) for d in dn]
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

    if returns == "dict":
        if unpack: raise ValueError("Cannot unpack when return type is 'dict'")
        returns = ["name", "path"]
        return_dict = True
    else: return_dict = False

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

            if types.is_string_type(exact_not_name):
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
        if types.is_string_type(returns):

            if unpack: raise ValueError("Cannot unpack with only one return value")

            if returns == "path": thing = item_path
            elif returns == "name": thing = item
            elif returns == "directory": thing = directory
            else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        else: # Assume 'returns' is a list or tuple

            thing = []

            for return_value in returns:

                if return_value == "path": thing.append(item_path)
                elif return_value == "name": thing.append(item)
                elif return_value == "directory": thing.append(directory)
                else: raise ValueError("Invalid option for 'returns': should be (a list of) 'path', 'name' or 'directory'")

        directory_paths.append(thing)

    # Return the list of directory paths
    if return_dict: return dict(directory_paths)
    elif unpack: return sequences.unpack(directory_paths, default_size=len(returns))
    else: return directory_paths

# -----------------------------------------------------------------

def find_directory_in_path(path, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                           exact_name=None, exact_not_name=None, startswith=None, endswith=None):

    """
    This function ...
    :param path: 
    :param recursive:
    :param ignore_hidden:
    :param contains:
    :param not_contains:
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :return: 
    """

    paths = directories_in_path(path, recursive=recursive, ignore_hidden=ignore_hidden, contains=contains,
                                not_contains=not_contains, exact_name=exact_name, exact_not_name=exact_not_name,
                                startswith=startswith, endswith=endswith)
    if len(paths) == 1: return paths[0]
    elif len(paths) == 0: raise ValueError("Not found")
    else: raise ValueError("Multiple directories found")

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

def rename_directory(parent, original_name, new_name):

    """
    This function ...
    :param parent:
    :param original_name:
    :param new_name:
    :return:
    """

    # Check whether the original directory exists
    original_path = join(parent, original_name)
    if not is_directory(original_path): raise ValueError("Directory '" + original_path + "' does not exist")

    # Determine new directory path
    new_path = join(parent, new_name)

    # Rename the directory
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

    # Copy
    shutil.copy(file_path, destination)

    if is_file(destination): return destination
    elif is_directory(destination): return join(destination, name(file_path))
    else: raise ValueError("Don't understand the destination: " + destination)

# -----------------------------------------------------------------

def copy_directory(path, directory_path, new_name=None):

    """
    This function ...
    :param path:
    :param directory_path:
    :param new_name:
    :return:
    """

    # Create the directory
    dirname = new_name if new_name is not None else name(path)
    copy_path = create_directory_in(directory_path, dirname)

    # Copy contents
    copy_from_directory(path, copy_path)

# -----------------------------------------------------------------

def copy_from_directory(from_directory, to_directory, **kwargs):

    """
    This function ..
    :param from_directory: 
    :param to_directory: 
    :param kwargs:
    :return: 
    """

    copy_files(files_in_path(from_directory, **kwargs), to_directory)

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

def copy_and_decompress_files(file_paths, directory_path):

    """
    This function ...
    :param file_paths: 
    :param directory_path: 
    :return: 
    """

    # Import here to avoid circular imports
    from . import archive

    # Loop over the files
    for file_path in file_paths:

        # Either decompress in the new directory of copy (if not an archive)
        if archive.is_archive(file_path): archive.decompress_file(file_path, directory_path)
        else: copy_file(file_path, directory_path)

# -----------------------------------------------------------------

def compressed_files_in_path(directory_path, returns="path"):

    """
    This function ...
    :param directory_path: 
    :param returns:
    :return: 
    """

    from .archive import extensions

    # Return
    return files_in_path(directory_path, returns=returns, extension=extensions)

# -----------------------------------------------------------------

def move_file(file_path, directory_path, new_name=None):

    """
    This function ...
    :param file_path:
    :param directory_path:
    :param new_name:
    :return:
    """

    if new_name is not None: destination = join(directory_path, new_name)
    else: destination = directory_path

    # Move the file
    shutil.move(file_path, destination)

# -----------------------------------------------------------------

def move_files(file_paths, directory_path):

    """
    This function ...
    :param file_paths:
    :param directory_path:
    :return:
    """

    for file_path in file_paths: copy_file(file_path, directory_path)

# -----------------------------------------------------------------

def read_lines(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Resolve path
    path = absolute_path(path)

    # Open the file
    with open(path, 'r') as fh:

        # Loop over the lines, cut off the end-of-line characters
        for line in fh.readlines(): yield line[:-1]

# -----------------------------------------------------------------

def get_lines(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    return list(read_lines(path))

# -----------------------------------------------------------------

def contains_lines(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    return len(get_lines(path)) > 0

# -----------------------------------------------------------------

def get_line(path, row):

    """
    This function ...
    :param path: 
    :param row: 
    :return: 
    """

    index = 0
    for line in read_lines(path):
        if index == row: return line
        index += 1
    return None

# -----------------------------------------------------------------

def read_last_lines(path, nlines):

    """
    This function ...
    :param path:
    :param nlines:
    :return:
    """

    for line in list(read_lines(path))[-nlines:]: yield line

# -----------------------------------------------------------------

def read_lines_reversed(path, buf_size=8192):

    """
    This function is a generator that returns the lines of a file in reverse order
    FROM http://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python
    :param path:
    :param buf_size:
    :return:
    """

    # Resolve path
    path = absolute_path(path)

    # Open the file
    with open(path) as fh:

        segment = None
        offset = 0
        if PYTHON_2: fh.seek(0, os.SEEK_END)
        else: file_size = fh.seek(0, os.SEEK_END)
        total_size = remaining_size = fh.tell()

        while remaining_size > 0:
            offset = min(total_size, offset + buf_size)
            if PYTHON_2: fh.seek(-offset, os.SEEK_END)
            else: fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # the first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # if the previous chunk starts right from the beginning of line
                # do not concact the segment to the last line of new chunk
                # instead, yield the segment first
                if buffer[-1] is not '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if len(lines[index]):
                    yield lines[index]

        yield segment

# -----------------------------------------------------------------

def write_line(filepath, line):

    """
    This function ...
    :param filepath:
    :param line:
    :return:
    """

    with open(filepath, 'w') as fh: fh.write(line + "\n")

# -----------------------------------------------------------------

def write_lines(filepath, lines):

    """
    This function ...
    :param filepath:
    :param lines:
    :return:
    """

    with open(filepath, 'w') as fh:
        for line in lines: fh.write(line + "\n")

# -----------------------------------------------------------------

def write_data(filepath, title, *cols):

    """
    This function ...
    :param filepath: 
    :param title:
    :param cols: the columns, as lists
    :return: 
    """

    # Check whether the columns have the same length
    if not sequences.equal_sizes(*cols): raise ValueError("Columns must have the same length")

    # Open the file
    with open(filepath, 'w') as fh:

        print("# " + title, file=fh)

        # Loop over the rows
        for index in range(len(cols[0])):

            # Construct the line
            values = [col[index] for col in cols]
            line = " ".join(str(value) for value in values)

            # Write the row
            print(line, file=fh)

# -----------------------------------------------------------------

def write_multi_data(filepath, *data):

    """
    THis function ...
    :param filepath: 
    :param data: title_a, column, column, ..., title_b, column, column, column, ..., title_c, ...
    :return: 
    """

    # Open the file
    with open(filepath, 'w') as fh:

        new_columns = None

        # Loop over the data
        for entry in data:

            # New title
            if types.is_string_type(entry):

                # Write previous data (if present)
                if new_columns is not None:

                    # Check whether the columns have the same length
                    if not sequences.equal_sizes(*new_columns): raise ValueError("Columns for a certain title must have the same length")

                    # Loop over the rows
                    for index in range(len(new_columns[0])):

                        # Construct the line
                        values = [col[index] for col in new_columns]
                        line = " ".join(str(value) for value in values)

                        # Write the row
                        print(line, file=fh)

                # Write new title
                print("# " + entry, file=fh)

                # Clear columns
                new_columns = []

            # Add column
            elif types.is_sequence(entry): new_columns.append(entry)

            # Invalid entry
            else: raise ValueError("Invalid input: must be strings and sequences of values")

        # Write last columns

        if new_columns is None: raise ValueError("No column data was found")

        # Check whether the columns have the same length
        if not sequences.equal_sizes(*new_columns): raise ValueError("Columns for a certain title must have the same length")

        # Loop over the rows
        for index in range(len(new_columns[0])):

            # Construct the line
            values = [col[index] for col in new_columns]
            line = " ".join(str(value) for value in values)

            # Write the row
            print(line, file=fh)

# -----------------------------------------------------------------

def append_line(filepath, line):

    """
    This function ...
    :param filepath:
    :param line:
    :return:
    """

    with open(filepath, 'a') as fh: fh.write(line + "\n")

# -----------------------------------------------------------------

def append_lines(filepath, lines):

    """
    This function ...
    :param filepath:
    :param lines:
    :return:
    """

    with open(filepath, 'a') as fh:
        for line in lines: fh.write(line + "\n")

# -----------------------------------------------------------------

def creation_date(filepath):

    """
    Try to get the date that a file was created, falling back to when it was
    last modified if that isn't possible.
    See http://stackoverflow.com/a/39501288/1709587 for explanation.
    """

    if platform.system() == 'Windows': seconds = os.path.getctime(filepath)
    else:
        stat = os.stat(filepath)
        try: seconds = stat.st_birthtime
        except AttributeError:
            # We're probably on Linux. No easy way to get creation dates here,
            # so we'll settle for when its content was last modified.
            seconds = stat.st_mtime

    # Return datetime object
    return datetime.datetime.fromtimestamp(seconds)

# -----------------------------------------------------------------

def reverse_readline(filename, buf_size=8192):

    """a generator that returns the lines of a file in reverse order"""

    with open(filename) as fh:
        for line in reverse_read_line_impl(fh, buf_size=buf_size): yield line

# -----------------------------------------------------------------

def reverse_read_line_impl(fh, buf_size=8192):

    """
    Thisf unction ...
    :param fh: 
    :param buf_size:
    :return: 
    """

    segment = None
    offset = 0
    fh.seek(0, os.SEEK_END)
    file_size = remaining_size = fh.tell()

    while remaining_size > 0:

        offset = min(file_size, offset + buf_size)
        fh.seek(file_size - offset)
        buffer = fh.read(min(remaining_size, buf_size))
        remaining_size -= buf_size
        lines = buffer.split('\n')
        # the first line of the buffer is probably not a complete line so
        # we'll save it and append it to the last line of the next buffer
        # we read
        if segment is not None:
            # if the previous chunk starts right from the beginning of line
            # do not concact the segment to the last line of new chunk
            # instead, yield the segment first
            if buffer[-1] is not '\n':
                lines[-1] += segment
            else:
                yield segment
        segment = lines[0]
        for index in range(len(lines) - 1, 0, -1):
            if len(lines[index]):
                yield lines[index]
    # Don't yield None if the file was empty
    if segment is not None:
        yield segment

# -----------------------------------------------------------------

def relative_to(path, base_path):

    """
    This function ...
    :param path:
    :param base_path:
    :return:
    """

    # Both absolute
    path = absolute_path(path)
    base_path = absolute_path(base_path)

    relative = path.split(base_path)[1]
    if relative.startswith("/"): return relative[1:]
    else: return relative

# -----------------------------------------------------------------
