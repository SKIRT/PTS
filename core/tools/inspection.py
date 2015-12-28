#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.inspection Contains some useful variables that store SKIRT and PTS installation directories
#  and provides functions for checking the presence and use of SKIRT and PTS dependencies.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import sys
import imp
import inspect
import subprocess
from operator import itemgetter, methodcaller
from collections import defaultdict
from contextlib import contextmanager
from distutils.spawn import find_executable
from importlib import import_module

# -----------------------------------------------------------------

# The path to the root PTS directory
pts_root_dir = inspect.getfile(inspect.currentframe()).split("/pts")[0]

# The path to the PTS package directory (PTS/pts)
pts_package_dir = os.path.join(pts_root_dir, "pts")

# The path to the PTS user directory (PTS/user)
pts_user_dir = os.path.join(pts_root_dir, "user")

# The path to the PTS do directory containing launchable scripts (PTS/pts/do)
pts_do_dir = os.path.join(pts_package_dir, "do")

# The path to the 'dat' directory for a given PTS subproject
def pts_dat_dir(subproject): return os.path.join(pts_package_dir, subproject, "dat")

# -----------------------------------------------------------------

# The path to the SKIRT executable
skirt_path = find_executable("skirt")

# The path to the root SKIRT directory
skirt_root_dir = skirt_path.split("/release")[0] if skirt_path is not None else None

# The path to the SKIRT repository
skirt_repo_dir = os.path.join(skirt_root_dir, "git") if skirt_path is not None else None

# The path to the SKIRT release directory
skirt_release_dir = os.path.join(skirt_root_dir, "release") if skirt_path is not None else None

# The path to the SKIRT run directory
skirt_run_dir = os.path.join(skirt_root_dir, "run") if skirt_path is not None else None

# -----------------------------------------------------------------

def has_mpi():

    """
    This function ...
    :return:
    """

    # Try opening the 'mpirun' executable
    try:
        devnull = open(os.devnull)
        subprocess.Popen("mpirun", stdout=devnull, stderr=devnull).communicate()
        return True
    except: return False

# -----------------------------------------------------------------

def skirt_version():

    """
    This function ...
    :return:
    """

    # Execute skirt with incorrect argument list and get its output
    process = subprocess.Popen([skirt_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0]

    # Return the relevant portion of the output
    return "SKIRT" + output.splitlines()[0].partition("SKIRT")[2]

# -----------------------------------------------------------------

def get_pip_versions():

    """
    This function ...
    :return:
    """

    packages = dict()

    # Launch the 'pip freeze' command and get the output
    output = subprocess.check_output(["pip", "freeze"])

    # Loop over the different package names and record the version number
    for entry in output.split("\n"):

        # Skip empty strings
        if not entry: continue

        # Get the package name and version
        name, version = entry.split("==")

        # Add them to the dictionary
        packages[name.lower()] = version

    # Return the dictionary
    return packages

# -----------------------------------------------------------------

@contextmanager
def ignore_site_packages_paths():

    """
    This function ...
    :return:
    """

    paths = sys.path[:]
    # remove working directory so that all
    # local imports fail
    if os.getcwd() in sys.path:
        sys.path.remove(os.getcwd())
    # remove all third-party paths
    # so that only stdlib imports will succeed
    sys.path = list(set(filter(
        None,
        filter(lambda i: all(('site-packages' not in i,
                              'python' in i or 'pypy' in i)),
               map(methodcaller('lower'), sys.path))
    )))
    yield
    sys.path = paths

# -----------------------------------------------------------------

def is_std_lib(module):

    """
    This function ...
    :param module:
    :return:
    """

    if not module:
        return False

    if module in sys.builtin_module_names:
        return True

    with ignore_site_packages_paths():
        imported_module = sys.modules.pop(module, None)
        try:
            import_module(module)
        except ImportError:
            return False
        else:
            return True
        finally:
            if imported_module:
                sys.modules[module] = imported_module

# -----------------------------------------------------------------

def get_all_dependencies():

    """
    This function ...
    :return:
    """

    # Create an empty dictionary to contain the required modules together with the places of use
    modules = defaultdict(list)

    # Recursively loop over all files inside this working directory
    for directory, subdirs, files in os.walk(pts_package_dir):

        # Loop over all files in the (sub)directory
        for filename in files:

            # If the file is not a python script, skip it
            if not filename.endswith(".py"): continue

            # Determine the full path to this file
            filepath = os.path.join(directory, filename)

            # Read the lines of the script file
            for line in open(filepath, 'r'):

                # Look for an 'import yyy' or 'from yyy import zzz' statement
                if line.startswith("import ") or (line.startswith("from ") and "import" in line):

                    # Get the name of the module
                    module = line.split()[1].split(".")[0]

                    # Add the module name to the list
                    if module: modules[module].append(filepath)

    return modules

# -----------------------------------------------------------------

def add_dependencies(dependencies, script_path, prefix=""):

    """
    This function ...
    :param path:
    :return:
    """

    # Skip files that are not python scripts
    if not script_path.endswith(".py"): raise ValueError("Not a valid script path")

    # Read the lines of the script file
    for line in open(script_path, 'r'):

        # If the current line does not contain an 'import yyy' or 'from yyy import zzz' statement, skip it
        if not (line.startswith("import ") or (line.startswith("from ") and "import" in line)): continue

        # Get the path to the modules that are being imported in the current line
        modules = get_modules(line, script_path)

        for module in modules:

            # Check if the imported module refers to a PTS module or an external package
            if module.startswith("/"): add_dependencies(dependencies, module, prefix=prefix+"  ")
            else: dependencies[module].add(script_path)

# -----------------------------------------------------------------

def is_present(package):

    """
    This function ...
    :return:
    """

    try:
        imp.find_module(package)
        return True
    except ImportError:
        return False

# -----------------------------------------------------------------

def get_modules(import_statement, script_path):

    """
    This function ...
    :param import_statement:
    :return:
    """

    splitted = import_statement.split()

    if len(splitted) <= 2: imported = []

    elif "," in splitted[3]:

        imported = [splitted[3][:-1]]

        for more in splitted[4:]:

            if "," in more: more = more[:-1]
            imported.append(more)

    else: imported = [splitted[3]]

    which = []

    # Check if this line denotes a relative import statement
    if splitted[1].startswith("."):

        after_dots = splitted[1].lstrip(".")

        number_of_dots = len(splitted[1]) - len(after_dots)

        # Determine the path to the PTS subpackage
        subpackage_dir = script_path
        for i in range(number_of_dots):
            subpackage_dir = os.path.dirname(subpackage_dir)

        subpackage_name = after_dots.split(".")[0]
        subpackage_path = os.path.join(subpackage_dir, subpackage_name)

        for name in imported: which.append(which_module(subpackage_path, name))

    # Absolute import of a pts class or module
    elif splitted[1].startswith("pts"):

        parts = splitted[1].split(".")[1:]

        subpackage_dir = pts_package_dir
        for part in parts:
            subpackage_dir = os.path.join(subpackage_dir, part)

        for name in imported: which.append(which_module(subpackage_dir, name))

    # External module
    else:

        # Get the name of the module
        module = splitted[1].split(".")[0]

        which.append(module)

    return which

# -----------------------------------------------------------------

def which_module(subpackage, name):

    """
    This function ...
    :param subpackage:
    :param name:
    :return:
    """

    if "," in name: print(name)
    if name.islower() and os.path.isfile(os.path.join(subpackage, name + ".py")): return os.path.join(subpackage, name + ".py")

    elif os.path.isfile(subpackage + ".py"):

        #print("  " + subpackage + ".py")
        return subpackage + ".py"

    elif os.path.isfile(os.path.join(subpackage, "__init__.py")):

        #print("  " + subpackage + ":")
        return os.path.join(subpackage, "__init__.py")

    else: raise ValueError("Don't know how to get further with " + subpackage + " and " + name)

# -----------------------------------------------------------------

def get_scripts():

    """
    This function ...
    :return:
    """

    # Loop over the directories within the 'do' subpackage
    scripts = []
    for item in os.listdir(pts_do_dir):

        # Get the full path to the item
        item_path = os.path.join(pts_do_dir, item)

        # Skip items that are not directories
        if not os.path.isdir(item_path): continue

        # Loop over the files in the directory
        for name in os.listdir(item_path):

            # Get the full name to the file
            file_path = os.path.join(item_path, name)

            # Skip items that are not files
            if not os.path.isfile(file_path): continue

            # Add the file path
            if file_path.endswith(".py") and "__init__" not in file_path: scripts.append((item, name))

    # Return the sorted list of script names
    return sorted(scripts, key=itemgetter(1))

# -----------------------------------------------------------------

def find_matches(scripts, name):

    """
    This function ...
    :return:
    """

    # Get a list of the script names that match the first command line argument, if there is one
    if "/" in name:

        matches = []
        dir_name = name.split("/")[0]
        script_name = name.split("/")[1]

        # Loop over all found items
        for item in scripts:
            if item[0] == dir_name and item[1].startswith(script_name): matches.append(item)

    # Return the list of matching scripts
    elif name is not None: return filter(lambda item: item[1].startswith(name), scripts)
    else: return []

# -----------------------------------------------------------------

def find_matching_script(script_name):

    """
    This function ...
    :return:
    """

    # Find matching scripts
    scripts = get_scripts()
    matches = find_matches(scripts, script_name)

    # If there is a unique match, return it
    if len(matches) == 1: return matches[0]

    # No matches
    elif len(matches) == 0:
        print("No match found. Available scripts:")

        # Sort on the 'do' subfolder name
        scripts = sorted(scripts, key=itemgetter(0))

        current_dir = None
        for script in scripts:

            if current_dir == script[0]:
                print(" "*len(current_dir) + "/" + script[1])
            else:
                print(script[0] + "/" + script[1])
                current_dir = script[0]
        return None

    # Multiple matches
    else:
        print("The command you provided is ambiguous. Possible matches:")

        # Sort on the 'do' subfolder name
        matches = sorted(matches, key=itemgetter(0))

        current_dir = None
        for script in matches:

            if current_dir == script[0]:
                print(" "*len(current_dir) + "/" + script[1])
            else:
                print(script[0] + "/" + script[1])
                current_dir = script[0]
        return None

# -----------------------------------------------------------------
