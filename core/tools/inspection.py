#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
This module ...
"""

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import imp
import inspect
import subprocess
from operator import itemgetter
from collections import defaultdict
from distutils.spawn import find_executable

# -----------------------------------------------------------------

# The path to the root PTS directory
pts_root_dir = inspect.getfile(inspect.currentframe()).split("/pts")[0]

# The path to the PTS package directory (PTS/pts)
pts_package_dir = os.path.join(pts_root_dir, "pts")

# The path to the PTS user directory (PTS/user)
pts_user_dir = os.path.join(pts_root_dir, "user")

# The path to the PTS do directory containing launchable scripts (PTS/pts/do)
pts_do_dir = os.path.join(pts_package_dir, "do")

# -----------------------------------------------------------------

# The path to the SKIRT executable
skirt_path = find_executable("skirt")

# The path to the root SKIRT directory
skirt_root_dir = skirt_path.split("/release")[0]

# The path to the SKIRT repository
skirt_repo_dir = os.path.join(skirt_root_dir, "git")

# The path to the SKIRT release directory
skirt_release_dir = os.path.join(skirt_root_dir, "release")

# The path to the SKIRT run directory
skirt_run_dir = os.path.join(skirt_root_dir, "run")

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
    process = subprocess.Popen([skirt_path, "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()[0]

    # Return the relevant portion of the output
    return "SKIRT" + output.splitlines()[0].partition("SKIRT")[2]

# -----------------------------------------------------------------

def list_dependencies(script_path):

    """
    This function ...
    :param path:
    :return:
    """

    if not script_path.endswith(".py"): raise ValueError("Not a valid script path")

    # Initialize a set to contain the required modules
    dependencies = set()

    # Read the lines of the script file
    for line in open(script_path, 'r'):

        # Look for an 'import yyy' or 'from yyy import zzz' statement
        if line.startswith("import ") or (line.startswith("from ") and "import" in line):

            splitted = line.split()

            # Check if this line denotes a relative import statement
            if splitted[1].startswith("."):

                after_dots = splitted[1].lstrip(".")

                number_of_dots = len(splitted[1]) - len(after_dots)

                # Determine the path to the PTS module
                module_dir = script_path
                for i in range(number_of_dots):
                    module_dir = os.path.dirname(module_dir)

                module_name = after_dots.split(".")[0]
                module_path = os.path.join(module_dir, module_name)

                print("module path: ", module_path)

            # Absolute import of a pts class or module
            elif splitted[1].startswith("pts"):

                parts = splitted[1].split(".")[1:]

                module_dir = pts_package_dir
                for part in parts:
                    module_dir = os.path.join(module_dir, part)

                module_name = splitted[3]

                # Check if a module with that name exists in the directory
                if os.path.isfile(os.path.join(module_dir, module_name + ".py")): print("module path: ", os.path.join(module_dir, module_name + ".py"))
                else:

                    # Open the subpackage's __init__ file
                    init_path = os.path.join(module_dir, "__init__.py")

                    # List the dependencies of the init file
                    list_dependencies(init_path)

                    #print("module path: ", module_dir, module_name)

            # External module
            else:

                # Get the name of the module
                module = line.split()[1].split(".")[0]

                # Check whether the module is present on this system
                try:
                    imp.find_module(module)
                    message = "present"
                except ImportError:
                    found = False
                    message = "not found!"

                print(module, ":", message)

            # Get the path of the script, relative to the 'PTS/git' directory
            #rel_filepath = script_path.split("PTS/pts/")[1]

            # Add the module
            #dependencies.add(module)

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
