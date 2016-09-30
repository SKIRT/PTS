#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.depends List the dependencies for a certain PTS do script.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import argparse
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.do.commandline import show_all_available, show_possible_matches
from pts.core.tools import introspection
from pts.core.tools.logging import log
from pts.core.tools import filesystem as fs

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("script", type=str, nargs="?", help="the name of the PTS do script for which to determine the dependencies")
parser.add_argument("-m", "--modules", action="store_true", help="show the PTS modules which import a given package")
parser.add_argument("-s", "--standard", action="store_true", help="show import packages from the python standard library")
parser.add_argument("-v", "--version", action="store_true", help="show the version numbers of the required packages")
parser.add_argument("-c", "--canopy", action="store_true", help="show whether the package is available through Canopy")
parser.add_argument("-p", "--pip", action="store_true", help="show whether the package is available through pip") # not working yet
parser.add_argument("-d", "--description", action="store_true", help="show a description") # not working yet

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

# If no script name is given, execute the "list_dependencies.py" script to list all dependencies of PTS and the
# PTS modules that use them
if arguments.script is None: dependencies = introspection.get_all_dependencies()

else:

    scripts = introspection.get_scripts()
    tables = introspection.get_arguments_tables()

    # Find matching 'do' commands (actual scripts or tabulated commands)
    matches = introspection.find_matches_scripts(arguments.script, scripts)
    table_matches = introspection.find_matches_tables(arguments.script, tables)

    # List the dependencies of the matching script
    dependencies = defaultdict(set)

    # No match
    if len(matches) + len(table_matches) == 0:
        show_all_available(scripts, tables)
        exit()

    # More matches
    elif len(matches) + len(table_matches) > 1:
        show_possible_matches(matches, table_matches, tables)
        exit()

    # Exactly one match from existing do script
    elif len(matches) == 1 and len(table_matches) == 0:

        # Determine the full path to the matching script
        script_path = os.path.join(introspection.pts_do_dir, matches[0][0], matches[0][1])

        introspection.add_dependencies(dependencies, script_path, set())

    # Exactly one match from tabulated command
    elif len(table_matches) == 1 and len(matches) == 0:

        # from pts.core.tools import logging
        # configuration

        # Path to class module

        table_match = table_matches[0]
        subproject = table_match[0]
        index = table_match[1]

        relative_class_module_path = tables[subproject]["Path"][index].replace(".", "/").rsplit("/", 1)[0] + ".py"

        class_module_path = fs.join(introspection.pts_subproject_dir(subproject), relative_class_module_path)

        logging_path = fs.join(introspection.pts_package_dir, "core", "tools", "logging.py")

        command_name = tables[subproject]["Command"][index]
        configuration_name = tables[subproject]["Configuration"][index]
        if configuration_name == "--": configuration_name = command_name
        configuration_module_path = fs.join(introspection.pts_root_dir, "pts/" + subproject + "/config/" + configuration_name + ".py")

        # Add dependencies
        encountered = set()
        introspection.add_dependencies(dependencies, logging_path, encountered)
        introspection.add_dependencies(dependencies, configuration_module_path, encountered)
        introspection.add_dependencies(dependencies, class_module_path, encountered)

# Get the names and versions of all installed python packages
packages = introspection.installed_python_packages()

# Get the names of the packages available through Canopy
if arguments.canopy:

    url = "https://www.enthought.com/products/canopy/package-index/"

    import requests
    from lxml import html

    page_as_string = requests.get(url).content
    tree = html.fromstring(page_as_string)

    tables = [e for e in tree.iter() if e.tag == 'table']
    table = tables[-1]

    table_rows = [e for e in table.iter() if e.tag == 'tr']
    column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

    table_rows = [e for e in table.iter() if e.tag == 'tr']
    column_headings = [e.text_content() for e in table_rows[0].iter() if e.tag == 'th']

    # return table_rows, column_headings

    #for row in table_rows: print(row.text_content())

    #objname = str(table_rows[1].text_content().split("\n")[1]).strip()

    canopy_packages = []

    for row in table_rows[1:]:

        name = row[0].text_content().strip()
        canopy_packages.append(name.lower())

else: canopy_packages = None

# Loop over the packages and report their presence
for dependency in sorted(dependencies, key=str.lower):

    # Get the list of PTS scripts for this dependency
    script_list = dependencies[dependency]

    # Skip packages from the standard library, unless the appropriate flag is enabled
    if introspection.is_std_lib(dependency) and not arguments.standard: continue

    # Check presency and version
    if dependency in packages:
        present = True
        version = packages[dependency] if arguments.version else None
    else:
        present = introspection.is_present_package(dependency)
        version = None

    # Check whether the package is available through Canopy
    in_canopy = dependency.lower() in canopy_packages if canopy_packages is not None else None
    in_canopy_string = " <<CANOPY>>" if in_canopy else ""

    # Check whether the package is available through pip
    #if arguments.pip:
    #else: pip = None

    # Check whether the current package is present
    if present:

        # Show package name, whether it's present and version number (if requested)
        if version is not None: log.success(dependency + ": present (version " + version + ")" + in_canopy_string)
        else: log.success(dependency + ": present" + in_canopy_string)

    # The package is not present
    else: log.error(dependency + ": not found" + in_canopy_string)

    # List the PTS modules that have this dependency
    if arguments.modules:
        for script in script_list: log.info("    " + script.split("PTS/pts/")[1])

# -----------------------------------------------------------------
