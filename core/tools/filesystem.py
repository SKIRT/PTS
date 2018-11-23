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
import psutil
import shutil
import platform
import subprocess
import datetime
import filecmp
import warnings
from collections import OrderedDict
from itertools import takewhile

# Import the relevant PTS classes and modules
from . import time
from . import types
from . import sequences
from . import strings

# -----------------------------------------------------------------

# Check if Python 2 or higher
if (sys.version_info > (3, 0)): PYTHON_2 = False # Python 3 or higher
else: PYTHON_2 = True # Python 2

# -----------------------------------------------------------------

home = os.path.expanduser('~')
downloads = os.path.join(home, "Downloads")
documents = os.path.join(home, "Documents")
desktop = os.path.join(home, "Desktop")

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

def to_home_directory():
    change_cwd(home)

# -----------------------------------------------------------------

def to_desktop():
    change_cwd(desktop)

# -----------------------------------------------------------------

def to_documents():
    change_cwd(documents)

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

def show_file_in_directory(path, wait=False):

    """
    This function ...
    :param path:
    :param wait:
    :return:
    """

    # Check if existing
    if not is_file(path): raise ValueError("The file '" + path + "' does not exist")

    # Check system
    if platform.system() != "Darwin": raise NotImplementedError("Not implemented for UNIX")

    # Set command
    command = ["open", "-F", "-R", "-a", "Finder", path, "-R"] # -R for revealing file in Finder, instead of opening

    # Call the command
    if wait:
        subprocess.call(command, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
        time.wait(10)
    else: subprocess.Popen(command)

# -----------------------------------------------------------------

def absolute_path(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Starts with environment variable
    if path.startswith("$"):
        if path.startswith("$HOME"): return path.replace("$HOME", home)
        elif path.startswith("$DESKTOP"): return path.replace("$DESKTOP", desktop)
        elif path.startswith("$DOCUMENTS"): return path.replace("$DOCUMENTS", documents)
        elif path.startswith("$DOWNLOADS"): return path.replace("$DOWNLOADS", downloads)
        else: raise ValueError("Unknown environment variable path: '" + path.split("/")[0] + "'")

    # Relative path (or with ~ for home directory)
    else: return os.path.abspath(os.path.expanduser(path))

# -----------------------------------------------------------------

def is_absolute(path):

    """
    This function ...
    :param path:
    :return:
    """

    return os.path.isabs(path)

# -----------------------------------------------------------------

def absolute_or_in(path, in_path):

    """
    This function ...
    :param path:
    :param in_path:
    :return:
    """

    if is_absolute(path): return path
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

def is_subdirectory(path, parent_path, check_existing=True):

    """
    This function returns whether
    :param path:
    :param parent_path:
    :param check_existing:
    :return:
    """

    if check_existing and not is_directory(path): raise ValueError("Not a directory: " + path)

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

    # Filenames are given
    if filenames is not None:

        # Loop over the filenames
        for filename in filenames:

            # Check
            filepath = join(directory, filename)
            if not is_file(filepath): return False

        # All checks passed
        return True

    # No filenames are given
    else: return len(files_in_path(directory)) > 0

# -----------------------------------------------------------------

def contains_directories(directory, **kwargs):

    """
    This function ...
    :param directory:
    :param kwargs:
    :return:
    """

    # Dirnames are given
    if kwargs.get("dirnames", None) is not None:

        dirnames = kwargs.pop("dirnames")

        # Loop over the dirnames
        for dirname in dirnames:

            dirpath = join(directory, dirname)
            if not is_directory(dirpath): return False

        return True

    # No names are given
    else: return len(directories_in_path(directory, **kwargs)) > 0

# -----------------------------------------------------------------

def contains_path(directory, path):

    """
    This function ...
    :param directory:
    :param path:
    :return:
    """

    # Set absolute paths
    directory = absolute_path(directory)
    path = absolute_path(path)

    # Check whether directory path is contained within the other path
    #return directory in path # NO: returns TRUE for contains_path("/test/path/ol", "/test/path/ola")

    # If the paths are the same: doesn't CONTAIN
    if directory == path: return False

    # Split
    directory_parts = directory.split("/")
    path_parts = path.split("/")

    #print(directory_parts)
    #print(path_parts)

    # If the path has less parts than the directory, the directory cannot possibly contain the path
    if len(path_parts) < len(directory_parts): return False

    # Check each part
    #for index, part in enumerate(path_parts):
        #if part not in directory_parts: return False
    #return True

    # Loop over the directory parts
    for index, part in enumerate(directory_parts):
        path_part = path_parts[index]
        #print(part, path_part)
        if path_part != part: return False
    return True

# -----------------------------------------------------------------

def any_contains_path(directories, path):

    """
    This function ...
    :param directories:
    :param path:
    :return:
    """

    for directory in directories:
        #print(directory, path)
        if contains_path(directory, path): return True
    return False

# -----------------------------------------------------------------

def contains_any_path(directory, paths):

    """
    This function ...
    :param directory:
    :param paths:
    :return:
    """

    for path in paths:
        if contains_path(directory, path): return True
    return False

# -----------------------------------------------------------------

def join(*args):

    """
    This function ...
    :param args:
    :return:
    """

    return os.path.join(*args)

# -----------------------------------------------------------------

def get_filepath(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    # Construct path, check and return
    path = join(*args)
    if not is_file(path): raise IOError(kwargs.pop("error_message", "File '" + path + "' does not exist"))
    return path

# -----------------------------------------------------------------

def get_dirpath(*args):

    """
    This function ...
    :param args:
    :return:
    """

    # Construct path, check and return
    path = join(*args)
    if not is_directory(path): raise IOError("Directory '" + path + "' does not exist")
    return path

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

def has_file(directory, filename, extension=None):

    """
    This function ...
    :param directory:
    :param filename:
    :param extension:
    :return:
    """

    if extension is not None: filename += "." + extension
    return is_file(join(directory, filename))

# -----------------------------------------------------------------

def is_empty(directory, ignore_hidden=True, besides=None, recursive=False):

    """
    This function ...
    :param directory:
    :param ignore_hidden:
    :param besides:
    :param recursive:
    :return:
    """

    # Recursive: look for only files within the directory hierarchy
    if recursive:

        filepaths = files_in_path(directory, ignore_hidden=ignore_hidden, exact_not_name=besides, recursive=True)
        return len(filepaths) == 0

    # Not recursive: count directories and files within the specified directory
    else:

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

def create_directory_in(base_path, name, clear=False, recursive=False):

    """
    This function ...
    :param base_path:
    :param name:
    :param clear:
    :param recursive:
    :return:
    """

    directory_path = join(base_path, name)
    if is_directory(directory_path):
        if clear: clear_directory(directory_path)
    else: create_directory(directory_path, recursive=recursive)
    return directory_path

# -----------------------------------------------------------------

def create_directory_in_cwd(name, clear=False):

    """
    This function ...
    :param name:
    :param clear:
    :param recursive:
    :return:
    """
    
    return create_directory_in(cwd(), name, clear=clear)

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

def clear_directory(path, recursive=False):

    """
    This function ...
    :param path:
    :return:
    """

    if recursive:

        for file_path in files_in_path(path): remove_file(file_path)
        for directory_path in directories_in_path(path): clear_directory(directory_path, recursive=True)

    else:

        for file_path in files_in_path(path): remove_file(file_path)
        for directory_path in directories_in_path(path): remove_directory(directory_path)

# -----------------------------------------------------------------

def remove_directory_if_present(path):

    """
    This function ...
    :param path:
    :return:
    """

    if is_directory(path):
        remove_directory(path)
        return True
    else: return False

# -----------------------------------------------------------------

def remove_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Test
    if not is_directory(path): raise ValueError("Is not a directory: '" + path + "'")

    # Remove
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

def remove_directories_but_keep(paths, keep_path):

    """
    This function ...
    :param paths:
    :param keep_path:
    :return:
    """

    for path in paths:

        # Don't remove if it contains the keep_path (looks recursively)
        if contains_path(path, keep_path): continue

        # Remove the
        remove_directory(path)

# -----------------------------------------------------------------

def remove_file_if_present(path):

    """
    This function ...
    :param path:
    :return:
    """

    if is_file(path):
        remove_file(path)
        return True
    else: return False

# -----------------------------------------------------------------

def remove_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    if not is_file(path): raise ValueError("Not a file: '" + path + "'")
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

def remove_files_in_path(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    remove_files(files_in_path(*args, **kwargs))

# -----------------------------------------------------------------

def remove_directories_in_path(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    remove_directories(directories_in_path(*args, **kwargs))

# -----------------------------------------------------------------

def remove_directory_or_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    if is_file(path): remove_file(path)
    elif is_directory(path): remove_directory(path)
    else: raise ValueError("Not an existing file or directory: " + path)

# -----------------------------------------------------------------

def remove_directories_and_files(paths):

    """
    This function ...
    :param paths:
    :return:
    """

    # Loop over the paths
    for path in paths: remove_directory_or_file(path)

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

    return total_size.to("GB")

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

def nfiles_in_path(*args, **kwargs):

    """
    Thisf unction ...
    :param args:
    :param kwargs:
    :return:
    """

    return len(files_in_path(*args, **kwargs))

# -----------------------------------------------------------------

def has_files_in_path(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    return nfiles_in_path(*args, **kwargs) > 0

# -----------------------------------------------------------------

def ndirectories_in_path(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    return len(directories_in_path(*args, **kwargs))

# -----------------------------------------------------------------

def has_directories_in_path(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    return ndirectories_in_path(*args, **kwargs) > 0

# -----------------------------------------------------------------

def nitems_in_path(*args, **kwargs):

    """
    This function ...
    :param args: 
    :param kwargs: 
    :return: 
    """

    files = files_in_path(*args, **kwargs)
    directories = directories_in_path(*args, **kwargs)
    return len(files) + len(directories)

# -----------------------------------------------------------------

def files_in_cwd(**kwargs):

    """
    This function ...
    :param kwargs:
    :return:
    """

    return files_in_path(cwd(), **kwargs)

# -----------------------------------------------------------------

def files_in_path(path=None, recursive=False, ignore_hidden=True, extension=None, not_extension=None, contains=None,
                  not_contains=None, extensions=False, returns="path", exact_name=None, exact_not_name=None,
                  startswith=None, endswith=None, sort=None, contains_operator="OR", recursion_level=None, unpack=False,
                  convert=None, convert_split_index=None, convert_split_pattern=" ", directory_not_contains=None,
                  directory_exact_not_name=None):

    """
    This function ...
    :param path:
    :param recursive:
    :param ignore_hidden:
    :param extension:
    :param not_extension:
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
    :param convert:
    :param convert_split_index:
    :param convert_split_pattern:
    :param directory_not_contains:
    :param directory_exact_not_name:
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
                try:
                    if convert_split_index is not None:
                        try: item_name = item_name.split(convert_split_pattern)[convert_split_index]
                        except IndexError: pass
                    value = sort(item_name)
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
        directory_name = name(directory)

        # Get the file name and extension
        item_name = os.path.splitext(item)[0]
        item_extension = os.path.splitext(item)[1][1:]

        # Ignore hidden files if requested
        if ignore_hidden and item.startswith("."): continue

        # Ignore in directory?
        if directory_exact_not_name is not None:
            if types.is_string_type(directory_exact_not_name):
                if directory_name == directory_exact_not_name: continue
            elif types.is_string_sequence(directory_exact_not_name):
                if directory_name in directory_exact_not_name: continue
            else: raise ValueError("Unkown type for 'directory_exact_not_name': " + str(directory_exact_not_name))

        # Directory name cannot contain string(s)
        if directory_not_contains is not None:
            if types.is_string_type(directory_not_contains):
                if directory_name == directory_not_contains: continue
            elif types.is_string_sequence(directory_not_contains):
                if sequences.any_in(directory_not_contains, directory_name): continue
            else: raise ValueError("directory_not_contains should be string or sequence")

        # Ignore files with extension different from the one that is specified
        if extension is not None:
            if types.is_string_type(extension):
                if item_extension != extension: continue
            elif types.is_string_sequence(extension):
                if item_extension not in extension: continue
            else: raise ValueError("Unknown type for 'extension': " + str(extension))

        # Ignore files with extensions that are in not_extension
        if not_extension is not None:
            if types.is_string_type(not_extension):
                if item_extension == not_extension: continue
            elif types.is_string_sequence(not_extension):
                if item_extension in not_extension: continue
            else: raise ValueError("Unknown type for 'not_extension': " + str(not_extension))

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
        if exact_name is not None:

            # One name
            if types.is_string_type(exact_name):
                if exact_name != item_name: continue

            # Sequence of names
            elif types.is_sequence(exact_name):
                if item_name not in exact_name: continue

            # Invalid
            else: raise ValueError("Invalid value for 'exact_name'")

        # If the filename matches the 'exact not name', skip it
        if exact_not_name is not None:

            # One name
            if types.is_string_type(exact_not_name):
                if exact_not_name == item_name: continue

            # Sequence
            elif types.is_sequence(exact_not_name):
                if item_name in exact_not_name: continue

            # Invalid
            else: raise ValueError("Invalid value for 'exact_not_name'")

        # Ignore filenames that do not start or end with the specified strings
        if startswith is not None:

            if types.is_string_type(startswith):
                if not item_name.startswith(startswith): continue
            elif types.is_sequence(startswith):
                for sw in startswith:
                    if item_name.startswith(sw): break
                else: # break not encountered
                    continue
            else: raise ValueError("Invalid value for 'startswith'")

        # Ignore filenames that do not end wwith the specified string
        if endswith is not None:

            if types.is_string_type(endswith):
                if not item_name.endswith(endswith): continue
            elif types.is_sequence(endswith):
                for ew in endswith:
                    if item_name.endswith(ew): break
                else: # break not encountered
                    continue

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

        # CONVERT?
        if convert is not None:
            if convert_split_index is not None: thing = thing.split(convert_split_pattern)[convert_split_index]
            thing = convert(thing)

        # Add to the list
        file_paths.append(thing)

    # Return the list of file paths
    if return_dict:
        if sort is not None: return OrderedDict(file_paths)
        else: return dict(file_paths) # file_paths is list of tuples
    elif unpack: return sequences.unpack(file_paths, default_size=len(returns))
    else: return file_paths

# -----------------------------------------------------------------

def find_file_in_path(path, **kwargs):

    """
    This function ...
    :param path: 
    :param kwargs:
    :return: 
    """

    return_none = kwargs.pop("return_none", False)

    # Get paths
    paths = files_in_path(path, **kwargs)
    if len(paths) == 1: return paths[0]
    elif len(paths) == 0:
        if return_none: return None
        else: raise ValueError("Not found")
    else: raise ValueError("Multiple files found")

# -----------------------------------------------------------------

def directories_in_cwd(**kwargs):

    """
    This function ...
    :param kwargs:
    :return:
    """

    return directories_in_path(cwd(), **kwargs)

# -----------------------------------------------------------------

def directories_in_path(path=None, recursive=False, ignore_hidden=True, contains=None, not_contains=None,
                        returns="path", exact_name=None, exact_not_name=None, startswith=None, endswith=None, sort=None,
                        recursion_level=None, unpack=False, convert=None):

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
    :param convert:
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
        if exact_name is not None:

            if types.is_string_type(exact_name):
                if exact_name != item: continue
            elif types.is_sequence(exact_name):
                if item not in exact_name: continue
            else: raise ValueError("Invalid option for 'exact_name': must be string or list")

        # If the directory name matches the 'exact not name', skip it
        if exact_not_name is not None:

            if types.is_string_type(exact_not_name):
                if exact_not_name == item: continue
            elif types.is_sequence(exact_not_name):
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

        # Convert
        if convert is not None: thing = convert(thing)

        # Add
        directory_paths.append(thing)

    # Return the list of directory paths
    if return_dict:
        if sort is not None: return OrderedDict(directory_paths)
        else: return dict(directory_paths)
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

    # Return the new file path
    return new_path

# -----------------------------------------------------------------

def rename_file_path(file_path, new_name):

    """
    This function ...
    :param file_path:
    :param new_name:
    :return:
    """

    directory = directory_of(file_path)
    original_name = name(file_path)
    return rename_file(directory, original_name, new_name)

# -----------------------------------------------------------------

def add_prefix(file_path, prefix):

    """
    This function ...
    :param file_path:
    :param prefix:
    :return:
    """

    old_name = strip_extension(name(file_path))
    extension = get_extension(file_path)
    new_name = prefix + old_name + "." + extension
    rename_file_path(file_path, new_name)

# -----------------------------------------------------------------

def added_prefix(file_path, prefix):

    """
    This function ...
    :param file_path:
    :param prefix:
    :return:
    """

    old_name = strip_extension(name(file_path))
    extension = get_extension(file_path)
    new_name = prefix + old_name + "." + extension

    directory = directory_of(file_path)
    return join(directory, new_name)

# -----------------------------------------------------------------

def add_suffix(file_path, suffix):

    """
    This function ...
    :param file_path:
    :param suffix:
    :return:
    """

    old_name = strip_extension(name(file_path))
    extension = get_extension(file_path)
    new_name = old_name + suffix + "." + extension
    rename_file_path(file_path, new_name)

# -----------------------------------------------------------------

def added_suffix(file_path, suffix):

    """
    This function ...
    :param file_path:
    :param suffix:
    :return:
    """

    old_name = strip_extension(name(file_path))
    extension = get_extension(file_path)
    new_name = old_name + suffix + "." + extension

    directory = directory_of(file_path)
    return join(directory, new_name)

# -----------------------------------------------------------------

def replace_file_path(file_path, replace, replace_with):

    """
    This function ...
    :param file_path:
    :param replace:
    :param replace_with:
    :return:
    """

    directory = directory_of(file_path)
    original_name = name(file_path)
    new_name = original_name.replace(replace, replace_with)
    return rename_file(directory, original_name, new_name)

# -----------------------------------------------------------------

def replace_file_paths_in_path(path, replace, replace_with, **kwargs):

    """
    This function ...
    :param path:
    :param replace:
    :param replace_with:
    :param kwargs:
    :return:
    """

    # Loop over the files
    for filepath in files_in_path(path, **kwargs):

        # Replace in the filepath
        replace_file_path(filepath, replace, replace_with)

# -----------------------------------------------------------------

def replace_file_paths_in_cwd(replace, replace_with, **kwargs):

    """
    This function ...
    :param replace:
    :param replace_with:
    :param kwargs:
    :return:
    """

    # Loop over the files
    for filepath in files_in_cwd(**kwargs):

        # Replace in the filepath
        replace_file_path(filepath, replace, replace_with)

# -----------------------------------------------------------------

def remove_extension(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    filename = name(filepath)
    directory = directory_of(filepath)
    new_filename = strip_extension(filename)

    # Rename, return the new file path
    return rename_file(directory, filename, new_filename)

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

def copy_file(file_path, directory_path, new_name=None, remove=False):

    """
    This function ...
    :param file_path:
    :param directory_path:
    :param new_name:
    :param remove:
    :return:
    """

    if new_name is not None: destination = join(directory_path, new_name)
    else: destination = directory_path

    # Remove previous
    if is_file(destination) and remove: remove_file(destination)

    # Copy
    shutil.copy(file_path, destination)

    if is_file(destination): return destination
    elif is_directory(destination): return join(destination, name(file_path))
    else: raise ValueError("Don't understand the destination: " + destination)

# -----------------------------------------------------------------

def copy_directory(path, directory_path, new_name=None, replace_files=False, replace_directories=False):

    """
    This function ...
    :param path:
    :param directory_path:
    :param new_name:
    :param replace_files:
    :param replace_directories:
    :return:
    """

    # Create the directory
    dirname = new_name if new_name is not None else name(path)
    copy_path = join(directory_path, dirname)
    if is_directory(copy_path):
        if replace_directories: clear_directory(copy_path)
    else: create_directory(copy_path)

    # Copy contents
    copy_from_directory(path, copy_path, replace_files=replace_files, replace_directories=replace_directories)

    # Return the directory path
    return copy_path

# -----------------------------------------------------------------

def copy_contents(path, directory_path):

    """
    Thisf unction ...
    :param path:
    :param directory_path:
    :return:
    """

    # Create the directory
    if not is_directory(directory_path): create_directory(directory_path)

    # Copy contents
    copy_from_directory(path, directory_path)

# -----------------------------------------------------------------

def clear_and_copy_from_directory(from_directory, to_directory, **kwargs):

    """
    This function ...
    :param from_directory:
    :param to_directory:
    :param kwargs:
    :return:
    """

    clear_directory(to_directory)
    copy_from_directory(from_directory, to_directory, **kwargs)

# -----------------------------------------------------------------

def copy_from_directory(from_directory, to_directory, **kwargs):

    """
    This function ..
    :param from_directory: 
    :param to_directory: 
    :param kwargs:
    :return: 
    """

    # Replace?
    replace_files = kwargs.pop("replace_files", False)
    replace_directories = kwargs.pop("replace_directories", False)

    # Set for files
    kwargs["replace"] = replace_files

    # Copy files
    copy_files_from_directory(from_directory, to_directory, **kwargs)

    # Remove invalid parameters for directories_in_path
    if "extension" in kwargs: del kwargs["extension"]
    if "not_extension" in kwargs: del kwargs["not_extension"]
    if "extensions" in kwargs: del kwargs["extensions"]

    # Set for directories
    kwargs.pop("replace")
    kwargs["replace_files"] = replace_files
    kwargs["replace_directories"] = replace_directories

    # Copy directories
    copy_directories_from_directory(from_directory, to_directory, **kwargs)

# -----------------------------------------------------------------

def copy_files_from_directory(from_directory, to_directory, **kwargs):

    """
    This function ...
    :param from_directory:
    :param to_directory:
    :param kwargs:
    :return:
    """

    # Replace?
    replace = kwargs.pop("replace", False)

    # Copy files
    copy_files(files_in_path(from_directory, **kwargs), to_directory, replace=replace)

# -----------------------------------------------------------------

def copy_directories_from_directory(from_directory, to_directory, **kwargs):

    """
    This function ...
    :param from_directory:
    :param to_directory:
    :param kwargs:
    :return:
    """

    # Replace?
    replace_files = kwargs.pop("replace_files", False)
    replace_directories = kwargs.pop("replace_directories", False)

    # Copy directories
    copy_directories(directories_in_path(from_directory, **kwargs), to_directory, replace_files=replace_files, replace_directories=replace_directories)

# -----------------------------------------------------------------

def copy_files(file_paths, directory_path, replace=False):

    """
    This function ...
    :param file_paths:
    :param directory_path:
    :param replace:
    :return:
    """

    for file_path in file_paths: copy_file(file_path, directory_path, remove=replace)

# -----------------------------------------------------------------

def copy_directories(directory_paths, directory_path, replace_files=False, replace_directories=False):

    """
    This function ...
    :param directory_paths:
    :param directory_path:
    :param replace_files:
    :param replace_directories:
    :return:
    """

    for dirpath in directory_paths: copy_directory(dirpath, directory_path, replace_files=replace_files, replace_directories=replace_directories)

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

    for file_path in file_paths: move_file(file_path, directory_path)

# -----------------------------------------------------------------

def read_lines_filehandle(fh, newlines=False):

    """
    This function ...
    :param fh:
    :param newlines:
    :return:
    """

    # Loop over the lines, cut off the end-of-line characters
    #for line in fh.readlines(): # MUCHHHHHH SLOWER: THIS LOADS EVERYTHING IN MEMORY FIRST!!!!!!!!
    for line in fh:
        if newlines: yield line
        elif line.endswith("\n"): yield line[:-1]
        else: yield line # sometimes, the last line does not end on a newline character!

# -----------------------------------------------------------------

def get_lines_filehandle(fh, newlines=False):

    """
    This function ...
    :param fh:
    :param newlines:
    :return:
    """

    return list(read_lines_filehandle(fh, newlines=newlines))

# -----------------------------------------------------------------

def read_lines(path, newlines=False):

    """
    This function ...
    :param path:
    :param newlines:
    :return:
    """

    # Resolve path
    path = absolute_path(path)

    # Open the file
    with open(path, 'r') as fh:
        for line in read_lines_filehandle(fh, newlines=newlines): yield line

# -----------------------------------------------------------------

def read_words(path, newlines=False):

    """
    This function ...
    :param path:
    :param newlines:
    :return:
    """

    for line in read_lines(path, newlines=newlines):
        for part in line.split(" "):
            if "\n" in part:
                before, after = part.split("\n")
                before = before.strip()
                after = after.strip()
                if before: yield before
                yield "\n"
                if after: yield after
            else:
                part = part.strip()
                if part: yield part

# -----------------------------------------------------------------

def get_lines(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    return list(read_lines(path))

# -----------------------------------------------------------------

def get_nlines(path):

    """
    This function ...
    :param path:
    :return:
    """

    from subprocess import check_output
    command = "wc -l '" + path + "'"
    output = check_output(command, shell=True)
    return int(output.split()[0])

# -----------------------------------------------------------------

def get_nlines_noheader(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Return
    return get_nlines(path) - get_nheader_lines(path)

# -----------------------------------------------------------------

def get_text(path):

    """
    This function ...
    :param path:
    :return:
    """

    return "\n".join(get_lines(path))

# -----------------------------------------------------------------

def contains_lines(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    return len(get_lines(path)) > 0

# -----------------------------------------------------------------

def has_lines(path):

    """
    Thisfunction ...
    :param path:
    :return:
    """

    return contains_lines(path)

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

def tail(filepath, window=20):

    """
    Returns the last `window` lines of file `f` as a list.
    """

    fh = open(filepath)

    if window == 0: return []
    BUFSIZ = 1024
    fh.seek(0, 2)
    bytes = fh.tell()
    size = window + 1
    block = -1
    data = []

    while size > 0 and bytes > 0:

        if bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            fh.seek(block * BUFSIZ, 2)
            # read BUFFER
            data.insert(0, fh.read(BUFSIZ))
        else:
            # file too small, start from begining
            fh.seek(0,0)
            # only read what was not read
            data.insert(0, fh.read(bytes))
        linesFound = data[0].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1

    # Return
    return ''.join(data).splitlines()[-window:]

# -----------------------------------------------------------------

def get_line_values(path, nlines, cls=float, last=False, start=None):

    """
    This function ...
    :param path:
    :param nlines:
    :param cls:
    :param last:
    :param start:
    :return:
    """

    if last: return [cls(line) for line in get_last_lines(path, nlines)]
    else: return [cls(line) for line in get_first_lines(path, nlines, skip_header=True)]

# -----------------------------------------------------------------

def get_last_lines(path, nlines):

    """
    This function ...
    :param path:
    :param nlines:
    :return:
    """

    return list(read_last_lines(path, nlines))

# -----------------------------------------------------------------

def get_last_line(path):

    """
    This function ...
    :param path:
    :return:
    """

    return get_last_lines(path, 1)[0]

# -----------------------------------------------------------------

def read_first_lines(path, nlines, skip_header=False):

    """
    This function ...
    :param path:
    :param nlines:
    :param skip_header:
    :return:
    """

    count = 0
    for line in read_lines(path):
        if skip_header and line.startswith("#"): continue
        yield line
        count += 1
        if count == nlines: return

# -----------------------------------------------------------------

def get_first_lines(path, nlines, skip_header=False):

    """
    This function ...
    :param path:
    :param nlines:
    :param skip_header:
    :return:
    """

    return list(read_first_lines(path, nlines, skip_header=skip_header))

# -----------------------------------------------------------------

def get_first_line(path):

    """
    This function ...
    :param path:
    :return:
    """

    return get_first_lines(path, 1)[0]

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

def remove_first_lines(filepath, startswith):

    """
    This function ...
    :param filepath:
    :param startswith:
    :return:
    """

    lines = []

    for line in read_lines_reversed(filepath):

        # Skip empty
        if line.strip() == ".": continue
        if not line.strip(): continue

        if line.strip().startswith(startswith): break
        lines.append(line)

    # Remove
    remove_file(filepath)

    # Write
    write_lines(filepath, reversed(lines))

# -----------------------------------------------------------------

def get_header_lines(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    lines = []

    for line in read_lines(filepath):

        # We are no longer at the header
        if not line.startswith("#"): break

        line = line.split("#", 1)[1].strip()
        lines.append(line)

    # Return the lines
    return lines

# -----------------------------------------------------------------

def get_nheader_lines(filepath):
    return len(get_header_lines(filepath))

# -----------------------------------------------------------------

def get_first_header_line(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    for line in read_lines(filepath):

        # We are no longer at the header
        if line.startswith("#"):

            line = line.split("#", 1)[1].strip()
            return line

    # No header
    raise IOError("No header")

# -----------------------------------------------------------------

def get_last_header_line(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    return get_header_lines(filepath)[-1]

# -----------------------------------------------------------------

def get_header_labels(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    return get_first_header_line(filepath).split()

# -----------------------------------------------------------------

def get_column_names(filepath, lower=False, return_units=False, capitalize=False):

    """
    This function ...
    :param filepath:
    :param lower:
    :param return_units:
    :param capitalize:
    :return:
    """

    from ..units.utils import clean_unit_string

    lines = get_header_lines(filepath)
    nlines = len(lines)

    # Only one line: return splitted first line
    if nlines == 1: column_names, column_units = _get_column_names_from_line(lines[0], return_units=True)

    # Header has lines stating the column index and names
    elif strings.any_startswith(lines, "column"):

        #print("2")
        #print(lines)

        colnames = dict()
        for line in lines:
            if not line.startswith("column"): continue
            colno = int(line.split("column")[1].split(":")[0].strip().split(" ")[0])
            if ":" in line: name = line.split(":")[1].strip()
            else: name = line.split(str(colno))[1].strip()
            colnames[colno] = name

        minno = min(colnames.keys())
        maxno = max(colnames.keys())

        first_is_one = minno == 1
        if first_is_one: ncolumns = maxno
        else: ncolumns = maxno - 1

        #print(ncolumns)
        column_names = [None] * ncolumns
        column_units = [None] * ncolumns
        #print(maxno)

        # Fill in column names
        for colno in colnames:

            if first_is_one: new_colno = colno - 1
            else: new_colno = colno
            #print(new_colno, colnames[colno])

            name, unit = _get_name_and_unit(colnames[colno])
            if unit is not None: unit = clean_unit_string(unit)
            if unit == "": unit = None
            #print(name, unit)
            #names[new_colno] = colnames[colno]
            column_names[new_colno] = name #if remove_units else colnames[colno]
            column_units[new_colno] = unit

    # Return splitted last line of the header
    else:

        #print("3")
        #return lines[-1].split()
        #column_names = lines[-1].split()
        #column_units = [None] * len(column_names)

        column_names, column_units = _get_column_names_from_line(lines[-1], return_units=True)

    # Make lowercase
    if lower: column_names = [name.lower() for name in column_names]

    # Capitalize first letter?
    if capitalize: column_names = [name.capitalize() for name in column_names]

    # Return
    if return_units: return column_names, column_units
    else: return column_names

# -----------------------------------------------------------------

def get_column_units(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    names, units = get_column_names(filepath, return_units=True)
    return units

# -----------------------------------------------------------------

def get_column_index(filepath, name, return_none=False):

    """
    This function ...
    :param filepath:
    :param name:
    :param return_none:
    :return:
    """

    try: return get_column_names(filepath).index(name)
    except ValueError:
        if return_none: return None
        else: raise ValueError("No column '" + name + "' in the file")

# -----------------------------------------------------------------

def _get_column_names_from_line(line, return_units=False):

    """
    This function ...
    :param line:
    :param return_units:
    :return:
    """

    from ..units.utils import clean_unit_string

    # splitted = lines[0].split()
    splitted = strings.split_except_within_round_brackets_and_double_quotes(line, add_quotes=False)
    #print(splitted)
    # Sanitize
    names = []
    column_units = []
    # print(splitted)
    for part in splitted:
        # print(part)

        if "[" in part and "]" in part:

            if strings.is_wrapped_by_round_brackets(part): names.append(part.split("[")[0].strip() + ")")
            else: names.append(part.split("[")[0].strip())

            unit = part.split("[")[1].split("]")[0]
            unit = clean_unit_string(unit)
            if unit == "": unit = None
            # print(unit)
            column_units.append(unit)

        elif strings.is_wrapped_by_squared_brackets(part): pass  # this is just the unit of the previous column name
        else: names.append(part)

    # print(names)

    # Check
    final_names = []
    for index in range(len(names)):

        name = names[index]
        # print(name)

        if strings.is_wrapped_by_round_brackets(name):  # this part of a formula e.g. log10(F_X) of the previous column

            nname = name[1:-1]
            if "[" in nname and "]" in nname: nname = nname.split("[")[1].split("]")[0]
            final_names[-1] += "(" + nname + ")"

        else: final_names.append(name)

    # Return the names
    if return_units: return final_names, column_units
    else: return final_names

# -----------------------------------------------------------------

def _get_name_and_unit(name_and_unit, capitalize=False):

    """
    This function ...
    :param name_and_unit:
    :param capitalize:
    :return:
    """

    #if "(" in name_and_unit and ")" in name_and_unit:
    if name_and_unit.count("(") == 1 and name_and_unit.count(")") == 1:

        before = name_and_unit.split(" (")[0] # before unit
        after = name_and_unit.split(")")[1] # after unit
        name = before.capitalize() + after if capitalize else before + after
        unit = name_and_unit.split(" (")[1].split(")")[0]
        if "dimensionless" in unit: unit = None

    elif name_and_unit.count("[") == 1 and name_and_unit.count("]") == 1:

        before = name_and_unit.split(" [")[0] # before unit
        after = name_and_unit.split("]")[1] # after unit
        name = before.capitalize() + after if capitalize else before + after
        unit = name_and_unit.split(" [")[1].split("]")[0]
        if "dimensionless" in unit: unit = None

    else:
        name = name_and_unit.capitalize() if capitalize else name_and_unit
        unit = None

    if ", i.e." in name: name = name.split(", i.e.")[0]

    # Return
    return name, unit

# -----------------------------------------------------------------

def get_ncolumns(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    return len(get_column_names(filepath))

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

def add_lines(filepath, lines, create=False):

    """
    This function ...
    :param filepath:
    :param lines:
    :param create:
    :return:
    """

    # Create the file?
    if not is_file(filepath):
        if create: mode = "w"
        else: raise IOError("File '" + filepath + "' is not present")
    else: mode = "a"

    # Write
    with open(filepath, mode) as fh:
        for line in lines: fh.write(line + "\n")

# -----------------------------------------------------------------

def write_text(filepath, text):

    """
    This function ...
    :param filepath:
    :param text:
    :return:
    """

    with open(filepath, 'w') as fh: fh.write(text + "\n")

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
            # NO: Giving an error is clearer!
            #raise NotImplementedError("Getting the creation date is impossible")
            # OK-> temporarily just give the modification date

    # Return datetime object
    return datetime.datetime.fromtimestamp(seconds)

# -----------------------------------------------------------------

def first_created_path(*paths):

    """
    This function ...
    :param paths:
    :return:
    """

    return min(paths, key=creation_date)

# -----------------------------------------------------------------

def last_created_path(*paths):

    """
    This function ...
    :param paths:
    :return:
    """

    return max(paths, key=creation_date)

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

    # WARNING: at the moment, absolute_path converts empty strings into the current working directory. Is this desired?

    # Both absolute
    path = absolute_path(path)
    base_path = absolute_path(base_path)

    relative = path.split(base_path)[1]
    if relative.startswith("/"): return relative[1:]
    else: return relative

# -----------------------------------------------------------------

def in_localhost_path(path, serve_path, port=8000):

    """
    Thisf unction ...
    :param path:
    :param serve_path:
    :return:
    """

    return "http://localhost:" + str(port) + "/" + relative_to(path, serve_path)

# -----------------------------------------------------------------

def base_directory(path):

    """
    This function ...
    :param path:
    :return:
    """

    return path.split("/")[0]

# -----------------------------------------------------------------

def replace_strings(path, replacement_dict):

    """
    This function ...
    :param path:
    :param replacement_dict:
    :return:
    """

    new_lines = []

    which_system = None

    # Read the lines
    for line in read_lines(path):

        if "<dustSystem" in line: which_system = "dust"
        elif "/dustSystem" in line: which_system = None

        if "<stellarSystem" in line: which_system = "stellar"
        elif "/stellarSystem" in line: which_system = None

        # Loop over the replacements
        for from_string in replacement_dict:
            to_string = replacement_dict[from_string]

            # Determine the new line
            if from_string in line: line = line.replace(from_string, to_string)

            if 'components type=""' in line:
                if which_system == "dust": line = line.replace('components type=""', 'components type="DustComp"')
                elif which_system == "stellar": line = line.replace('components type=""', 'components type="StellarComp"')
                else: raise RuntimeError("Something went wrong")

        # Add the line
        new_lines.append(line)

    # Remove the file
    remove_file(path)

    # Write lines
    write_lines(path, new_lines)

# -----------------------------------------------------------------

def add_extension(filename, extension):

    """
    This function ...
    :param filename:
    :param extension:
    :return:
    """

    return filename + "." + extension

# -----------------------------------------------------------------

def appended_filename(filepath, append_with, sep=None):

    """
    This function ...
    :param filepath:
    :param append_with:
    :param sep:
    :return:
    """

    if sep is None: sep = ""

    filename = name(filepath)
    the_name = strip_extension(filename)
    extension = get_extension(filename)
    appended_name = the_name + sep + append_with
    return add_extension(appended_name, extension)

# -----------------------------------------------------------------

def prepended_filename(filepath, prepend_with, sep=None):

    """
    This function ...
    :param filepath:
    :param prepend_with:
    :param sep:
    :return:
    """

    if sep is None: sep = ""
    filename = name(filepath)
    return prepend_with + sep + filename

# -----------------------------------------------------------------

def prepended_and_appended_filename(filepath, prepend, append, sep=None):

    """
    This function ...
    :param filepath:
    :param prepend:
    :param append:
    :param sep:
    :return:
    """

    if sep is None: sep = ""

    filename = name(filepath)
    the_name = strip_extension(filename)
    extension = get_extension(filename)

    new_name = prepend + sep + the_name + sep + append
    return add_extension(new_name, extension)

# -----------------------------------------------------------------

def appended_filepath(filepath, append_with, sep=None):

    """
    This function ...
    :param filepath:
    :param append_with:
    :param sep:
    :return:
    """

    return join(directory_of(filepath), appended_filename(filepath, append_with, sep=sep))

# -----------------------------------------------------------------

def prepended_filepath(filepath, prepend_with, sep=None):

    """
    This function ...
    :param filepath:
    :param prepend_with:
    :param sep:
    :return:
    """

    return join(directory_of(filepath), prepended_filename(filepath, prepend_with, sep=sep))

# -----------------------------------------------------------------

def prepended_and_appended_filepath(filepath, prepend, append, sep=None):

    """
    This function ...
    :param filepath:
    :param prepend:
    :param append:
    :param sep:
    :return:
    """

    return join(directory_of(filepath), prepended_and_appended_filename(filepath, prepend, append, sep=sep))

# -----------------------------------------------------------------

def equal_files(filepath_a, filepath_b, shallow=True):

    """
    This function ...
    :param filepath_a:
    :param filepath_b:
    :param shallow:
    :return:
    """

    return filecmp.cmp(filepath_a, filepath_b, shallow=shallow)

# -----------------------------------------------------------------

def all_hidden_files(filepaths):

    """
    This function ...
    :param filepaths:
    :return:
    """

    for filepath in filepaths:
        filename = name(filepath)
        if not filename.startswith("."): return False
    return True

# -----------------------------------------------------------------

def equal_directories(directory_a, directory_b, shallow=True, create=False, report=False, ignore_hidden=True):

    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    @param dir1: First directory path
    @param dir2: Second directory path

    @return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
    """

    # Check existence
    if not is_directory(directory_a):
        if create:
            if report: print("Directory '" + directory_a + "' does not exist yet: creating ...")
            create_directory(directory_a)
        else: raise ValueError("The directory '" + directory_a + "' does not exist")
    if not is_directory(directory_b):
        if create:
            if report: print("Directory '" + directory_b + "' does not exist yet: creating ...")
            create_directory(directory_b)
        else: raise ValueError("The directory '" + directory_b + "' does not exist")

    # Do directory comparison
    dirs_cmp = filecmp.dircmp(directory_a, directory_b)
    if len(dirs_cmp.left_only) > 0:
        if ignore_hidden and all_hidden_files(dirs_cmp.left_only): pass
        else:
            if report: print("Only in '" + directory_a + "': " + str(dirs_cmp.left_only))
            return False
    if len(dirs_cmp.right_only) > 0:
        if ignore_hidden and all_hidden_files(dirs_cmp.right_only): pass
        else:
            if report: print("Only in '" + directory_b + "': " + str(dirs_cmp.right_only))
            return False
    if len(dirs_cmp.funny_files) > 0:
        if report: print("Funny files: " + str(dirs_cmp.funny_files))
        return False

    # Check files
    (_, mismatch, errors) = filecmp.cmpfiles(directory_a, directory_b, dirs_cmp.common_files, shallow=shallow)
    if len(mismatch) > 0:
        if report: print("Mismatches: " + str(mismatch))
        return False
    if len(errors) > 0:
        if report: print("Errors: " + str(errors))
        return False

    # Loop over common directories and check them
    for common_dir in dirs_cmp.common_dirs:

        # Get paths
        new_dir1 = os.path.join(directory_a, common_dir)
        new_dir2 = os.path.join(directory_b, common_dir)

        # Recursively call this function
        if not equal_directories(new_dir1, new_dir2, shallow=shallow, report=report): return False

    # All checks passed
    return True

# -----------------------------------------------------------------

def equal_number_of_files(directory_a, directory_b, recursive=True, create=False):

    """
    This function ...
    :param directory_a:
    :param directory_b:
    :param recursive:
    :param create:
    :return:
    """

    # Get the number of files in directory a
    if is_directory(directory_a): nfiles_a = nfiles_in_path(directory_a, recursive=recursive)
    else:
        if create: create_directory(directory_a)
        else: raise ValueError("The directory '" + directory_a + "' does not exist")
        nfiles_a = 0

    # Get the number of files in directory b
    if is_directory(directory_b): nfiles_b = nfiles_in_path(directory_b, recursive=recursive)
    else:
        if create: create_directory(directory_b)
        else: raise ValueError("The directory '" + directory_b + "' does not exist")
        nfiles_b = 0

    # Compare
    return nfiles_a == nfiles_b

# -----------------------------------------------------------------

def equal_number_of_directories(directory_a, directory_b, recursive=True, create=False):

    """
    This function ...
    :param directory_a:
    :param directory_b:
    :param recursive:
    :param create:
    :return:
    """

    # Get the number of directories
    if is_directory(directory_a): ndirectories_a = ndirectories_in_path(directory_a, recursive=recursive)
    else:
        if create: create_directory(directory_a)
        else: raise ValueError("The directory '" + directory_a + "' does not exist")
        ndirectories_a = 0

    # Get the number of directories
    if is_directory(directory_b): ndirectories_b = ndirectories_in_path(directory_b, recursive=recursive)
    else:
        if create: create_directory(directory_b)
        else: raise ValueError("The directory '" + directory_b + "' does not exist")
        ndirectories_b = 0

    # Compare
    return ndirectories_a == ndirectories_b

# -----------------------------------------------------------------

def equal_number_of_items(directory_a, directory_b, recursive=True, create=False):

    """
    This function ...
    :param directory_a:
    :param directory_b:
    :param recursive:
    :param create:
    :return:
    """

    # Get the number of items
    if is_directory(directory_a): nitems_a = nitems_in_path(directory_a, recursive=recursive)
    else:
        if create: create_directory(directory_a)
        else: raise ValueError("The directory '" + directory_a + "' does not exist")
        nitems_a = 0

    # Get the number of items
    if is_directory(directory_b): nitems_b = nitems_in_path(directory_b, recursive=recursive)
    else:
        if create: create_directory(directory_b)
        else: raise ValueError("The directory '" + directory_b + "' does not exist")
        nitems_b = 0

    # Compare
    return nitems_a == nitems_b

# -----------------------------------------------------------------

def update_file_in(source, directory, create=False, report=False, return_target=False):

    """
    This function ...
    :param source:
    :param directory:
    :param create:
    :param report:
    :param return_target
    :return:
    """

    filename = name(source)
    target = join(directory, filename)
    updated = update_file(source, target, create=create, report=report)

    if return_target: return updated, target
    else: return updated

# -----------------------------------------------------------------

def update_file(source, target, create=False, report=False):

    """
    This function ...
    :param source:
    :param target:
    :param create:
    :param report:
    :return:
    """

    # Source does not exist
    if not is_file(source): raise ValueError("The file '" + source + "' does not exist")

    # Target does not yet exist
    if not is_file(target):

        if create:
            if report: print("File does not yet exist: copying ...")
            copy_file(source, target)
        else: raise ValueError("The file '" + target + "' does not exist")
        return True

    # Replace if not equal
    elif not equal_files(source, target):
        if report: print("Files not equal: replacing ...")
        replace_file(target, source)
        return True

    # No update required
    return False

# -----------------------------------------------------------------

def update_directory_in(source, directory, create=False, report=False, ignore_hidden=True):

    """
    This function ...
    :param source:
    :param directory:
    :param create:
    :param report:
    :param ignore_hidden:
    :return:
    """

    dirname = name(source)
    target = join(directory, dirname)
    return update_directory(source, target, create=create, report=report, ignore_hidden=ignore_hidden)

# -----------------------------------------------------------------

def update_directory(source, target, create=False, report=False, ignore_hidden=True):

    """
    This function ...
    :param source:
    :param target:
    :param create:
    :param report:
    :param ignore_hidden:
    :return:
    """

    # Source directory does not exist
    if not is_directory(source): raise ValueError("The directory '" + source + "' does not exist")

    # if not fs.equal_number_of_files(fonts_path, mount_fonts_path, create=True):
    if not equal_directories(source, target, create=create, report=report, ignore_hidden=ignore_hidden):
        clear_directory(target)
        copy_from_directory(source, target)
        return True

    # No update required
    else: return False

# -----------------------------------------------------------------

def nallowed_open_files():

    """
    This function ...
    :return:
    """

    import resource
    return resource.getrlimit(resource.RLIMIT_NOFILE)[0]

# -----------------------------------------------------------------

def set_nallowed_open_files(number):

    """
    This function ...
    :param number:
    :return:
    """

    from .introspection import is_linux, is_macos

    if is_macos():

        import resource
        resource.setrlimit(resource.RLIMIT_NOFILE, (number, -1))

    elif is_linux():

        import warnings
        warnings.warn("On Linux, setting the number of allowed open files is currently not implemented")

    else: raise RuntimeError("Platforms other than MacOS and Linux are not supported")

# -----------------------------------------------------------------

def open_files():

    """
    This function ...
    :return:
    """

    all = []
    for proc in psutil.process_iter():
        try:
            files = proc.open_files()
            all.extend(fh.path for fh in files)
        except psutil.AccessDenied: print("Acces denied for process", proc.name())
    return all

# -----------------------------------------------------------------

def nopen_files():

    """
    This function ...
    :return: 
    """


    return len(open_files())

# -----------------------------------------------------------------

def get_backup_filepath(filepath, suffix=None, prefix=None, sep="_"):

    """
    Thisf unction ...
    :param filepath:
    :param suffix:
    :param prefix:
    :param sep:
    :return:
    """

    # Set name for backup file
    if prefix is not None and suffix is not None: backup_filepath = prepended_and_appended_filepath(filepath, prefix + sep, sep + suffix)
    elif prefix is None and suffix is not None: backup_filepath = appended_filepath(filepath, sep + suffix)
    elif suffix is None and prefix is not None: backup_filepath = prepended_filepath(filepath, prefix + sep)
    else:
        suffix = "backup"  # default
        backup_filepath = appended_filepath(filepath, sep + suffix)

    # Return
    return backup_filepath

# -----------------------------------------------------------------

def backup_file(filepath, suffix=None, prefix=None, exists="error", remove=False, sep="_", check_filepath=True,
                remove_if_exists=False):

    """
    This function ...
    :param filepath:
    :param suffix:
    :param prefix:
    :param exists:
    :param remove:
    :param sep:
    :param check_filepath:
    :param remove_if_exists:
    :return:
    """

    # Check
    if check_filepath and not is_file(filepath): raise ValueError("File '" + filepath + "' does not exist")

    # Get backup filepath
    backup_filepath = get_backup_filepath(filepath, suffix=suffix, prefix=prefix)

    # Already exists
    if is_file(backup_filepath):

        # First backup the seemingly earlier backup
        if exists == "backup":

            backup_file(backup_filepath, suffix=suffix, prefix=prefix, sep=sep, exists=exists) # backup the backup
            remove_file(backup_filepath)

        # Give error
        elif exists == "error": raise IOError("Backup file path already exists (" + backup_filepath + ")")

        # Just overwrite the previous, but give warning
        elif exists == "overwrite":
            warnings.warn("The backup '" + backup_filepath + "' already exists: will be overwritten ...")
            remove_file(backup_filepath)

        # Do nothing, trust earlier backup
        elif exists == "pass":
            warnings.warn("The backup '" + backup_filepath + "' already exists: skipping ...")
            if remove and remove_if_exists and is_file(filepath): remove_file(filepath)
            return

        # Invalid
        else: raise ValueError("Invalid value for 'exists'")

    # Actually make the backup
    copy_file(filepath, backup_filepath)

    # Remove original?
    if not is_file(backup_filepath): raise IOError("Something went wrong: backup file not created")
    if remove: remove_file(filepath)

# -----------------------------------------------------------------

def backup_directory(path, suffix="backup"):

    """
    This function ...
    :param path:
    :param suffix:
    :return:
    """

    backup_path = path + "_" + suffix
    if is_directory(backup_path): raise IOError("Backup directory already exists (" + backup_path + ")")
    copy_contents(path, backup_path)

# -----------------------------------------------------------------

def backup_files(filepaths, suffix="backup"):

    """
    This function ...
    :param filepaths:
    :param suffix:
    :return:
    """

    for filepath in filepaths: backup_file(filepath, suffix=suffix)

# -----------------------------------------------------------------

def backup_files_in_directory(path, suffix="backup", extension=None, contains=None):

    """
    This function ...
    :param path:
    :param suffix:
    :param extension:
    :param contains:
    :return:
    """

    paths = files_in_path(path, extension=extension, contains=contains)
    return backup_files(paths, suffix=suffix)

# -----------------------------------------------------------------

def read_start(path, ncharacters):

    """
    Thisfunction ...
    :param ncharacters:
    :return:
    """

    with open(path) as f: return f.read(ncharacters)

# -----------------------------------------------------------------

def get_file_hash(path, blocksize=2**20):

    """
    This function ...
    :param path:
    :param blocksize:
    :return:
    """

    import hashlib

    # Open the file
    f = open(path)

    md5 = hashlib.md5()

    # Read the file
    while True:
        data = f.read(blocksize)
        if not data: break
        md5.update(data)

    # Close the file
    f.close()

    # Return the hash
    #return md5.digest()
    return md5.hexdigest()

# -----------------------------------------------------------------

def equal_files_hash(filepath_a, filepath_b):

    """
    This function ...
    :param filepath_a:
    :param filepath_b:
    :return:
    """

    return get_file_hash(filepath_a) == get_file_hash(filepath_b)

# -----------------------------------------------------------------

def common_directory(paths, sep='/'):

    """
    This function ...
    :param paths:
    :param sep:
    :return:
    """

    bydirectorylevels = zip(*[p.split(sep) for p in paths])
    return sep.join(x[0] for x in takewhile(sequences.all_equal, bydirectorylevels))

# -----------------------------------------------------------------

def get_volume_path(volname):

    """
    This function ...
    :param volname:
    :return:
    """

    from .introspection import is_linux, is_macos, username
    if is_macos(): return join("/Volumes", volname)
    elif is_linux(): return join("/media", username(), volname)
    else: raise RuntimeError("Platforms other than MacOS and Linux are not supported")

# -----------------------------------------------------------------

def file_nbytes(path):

    """
    This function ...
    :param path:
    :return:
    """

    from .introspection import is_linux, is_macos
    #from .terminal import execute
    from subprocess import check_output

    # Mac
    if is_macos():

        command = "cat '" + path + "' | wc -c"
        #print(command)
        output = check_output(command, shell=True)
        #output = execute(command)
        #print(output)
        return int(output)

    # Linux
    elif is_linux():

        command = "du -b '" + path + "' | awk '{ print $1}' | bc"
        output = check_output(command, shell=True)
        #print(output)
        return int(output)

    # Not supported
    else: raise RuntimeError("Platforms other than MacOS and Linux are not supported")

# -----------------------------------------------------------------

def file_nbits(path):

    """
    This function ...
    :param path:
    :return:
    """

    return file_nbytes(path) * 8

# -----------------------------------------------------------------

def get_columns(filepath, method="numpy", dtype=None, indices=None):

    """
    This function ...
    :param filepath:
    :param method:
    :param dtype:
    :param indices:
    :return:
    """

    from ..basics.log import log

    # Debugging
    log.debug("Reading data (this can take a while) ...")

    # Using NumPy
    if method == "numpy":

        import numpy as np
        columns = np.loadtxt(filepath, unpack=True, ndmin=2, dtype=dtype, usecols=indices)

    # Using Pandas
    elif method == "pandas":

        import pandas as pd
        df = pd.read_csv(filepath, sep=" ", comment="#", header=None, dtype=dtype, usecols=indices)
        #print(df)
        ncolumns = len(df.columns)
        #print("ncols", ncolumns)
        if indices is not None: column_indices = indices
        else: column_indices = range(ncolumns)
        columns = [df[index].values for index in column_indices]
        #columns = [df.columns[index] for index in range(ncolumns)]

    # Invalid
    else: raise ValueError("Invalid method: must be 'numpy' or 'pandas'")

    # Return the columns
    return columns

# -----------------------------------------------------------------

def write_columns(*args, **kwargs):

    """
    This function ...
    :param args:
    :param kwargs:
    :return:
    """

    import numpy as np

    # Get filepath
    filepath = args[0]
    if not types.is_string_type(filepath): raise ValueError("First argument must be the filepath")

    # Get columns
    cols = args[1:]

    # Save
    np.savetxt(filepath, np.transpose(cols))

# -----------------------------------------------------------------

def get_2d_data(filepath, method="numpy", dtype=None, columns=None, base="column"):

    """
    Thisn function ...
    :param filepath:
    :param method:
    :param dtype:
    :param columns:
    :param base:
    :return:
    """

    from ..basics.log import log

    # Debugging
    log.debug("Reading data (this can take a while) ...")

    # Using NumPy
    if method == "numpy":

        import numpy as np
        data = np.loadtxt(filepath, ndmin=2, dtype=dtype, usecols=columns)

    # Using pandas
    elif method == "pandas":

        import pandas as pd
        df = pd.read_csv(filepath, sep=" ", comment="#", header=None, dtype=dtype, usecols=columns)
        data = df.values

    # Invalid
    else: raise ValueError("Invalid method: must be 'numpy' or 'pandas'")

    # Return the data
    if base == "column": return data.transpose()
    elif base == "row": return data
    else: raise ValueError("Invalid base: mus be 'column' or 'row'")

# -----------------------------------------------------------------

def get_column(filepath, index, dtype, method="numpy"):

    """
    This function ...
    :param filepath:
    :param index:
    :param dtype
    :param method:
    :return:
    """

    from ..basics.log import log

    # Debugging
    log.debug("Reading data (this can take a while) ...")

    # Using NumPy
    if method == "numpy":

        import numpy as np
        column = np.loadtxt(filepath, usecols=index, dtype=dtype)

    # Using Pandas
    elif method == "pandas":

        import pandas as pd
        df = pd.read_csv(filepath, sep=" ", comment="#", header=None, usecols=[index], dtype=dtype)
        column = df[index].values # HAS TO BE INDEX, NOT ZERO (PANDAS REMEMBERS THE ORIGINAL COLUMN INDEX, OTHER COLUMNS ARE INACCESSIBLE)

    # Invalid
    else: raise ValueError("Invalid method: must be 'numpy' or 'pandas'")

    # Return the column
    return column

# -----------------------------------------------------------------

def get_row(filepath, index):

    """
    This function ...
    :param filepath:
    :param index:
    :return:
    """

    from ..basics.log import log

    # Debugging
    log.debug("Reading data (this can take a while) ...")

    import pandas as pd
    df = pd.read_csv(filepath, sep=" ", comment="#", header=None, skiprows=lambda i: i != index)

    # Return the row as an array
    return df.values[0]

# -----------------------------------------------------------------
