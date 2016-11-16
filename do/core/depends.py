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
import argparse
import requests
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.do.commandline import show_all_available, show_possible_matches
from pts.core.tools import introspection
from pts.core.tools.logging import log
from pts.core.tools import filesystem as fs
from pts.core.tools.pypi import search
from pts.core.tools import formatting as fmt

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser()
parser.add_argument("script", type=str, nargs="?", help="the name of the PTS do script for which to determine the dependencies")
parser.add_argument("-m", "--modules", action="store_true", help="show the PTS modules which import a given package")
parser.add_argument("-s", "--standard", action="store_true", help="show import packages from the python standard library")
parser.add_argument("-v", "--version", action="store_true", help="show the version numbers of the required packages")
parser.add_argument("-c", "--canopy", action="store_true", help="show whether the package is available through Canopy")
parser.add_argument("-p", "--pip", action="store_true", help="show whether the package is available through pip")
parser.add_argument("-d", "--description", action="store_true", help="show a description")
parser.add_argument("-j", "--subprojects", action="store_true", help="show the dependencies per subproject")

# Parse the command line arguments
arguments = parser.parse_args()

# -----------------------------------------------------------------

def print_dependency(dependency, data):

    script_list = data[dependency][0]
    version = data[dependency][2]
    in_canopy = data[dependency][3]
    in_pip = data[dependency][4]
    description = data[dependency][5]

    info = []
    if arguments.version: info.append(("version", version))
    if arguments.canopy: info.append(("in canopy", in_canopy))
    if arguments.pip: info.append(("in pip", in_pip))
    if arguments.description: info.append(("description", description))

    had_indent = False

    if len(info) == 0:
        print("   " + fmt.underlined + dependency + fmt.reset)
    elif len(info) == 1:
        if info[0][1] is not None: print("   " + fmt.underlined + dependency + fmt.reset + " (", info[0][1], ")")
    else:
        print("   " + fmt.underlined + dependency + fmt.reset)
        print("")
        had_indent = True
        for info_entry in info:
            if info_entry[1] is not None: print("    - " + info_entry[0] + ": ", info_entry[1])

    # List the PTS modules that have this dependency
    if arguments.modules:
        print("")
        had_indent = True
        for script in script_list: print("    * " + script.split("PTS/pts/")[1])

    if had_indent: print("")

# -----------------------------------------------------------------

dependencies_for_subproject = None

# If no script name is given, execute the "list_dependencies.py" script to list all dependencies of PTS and the
# PTS modules that use them
if arguments.script is None:

    if arguments.subprojects:

        dependencies = defaultdict(set)

        dependencies_for_subproject = dict()

        directories_for_subproject = defaultdict(list)

        # Loop over the subdirectories of the 'do' directory
        for path, name in fs.directories_in_path(introspection.pts_do_dir, returns=["path", "name"]):

            subproject = name
            directories_for_subproject[subproject].append(path)

        # Loop over the other directories in 'pts' (other than 'do' and 'doc')
        for path, name in fs.directories_in_path(introspection.pts_package_dir, returns=["path", "name"], exact_not_name=["do", "doc"]):

            subproject = name
            directories_for_subproject[subproject].append(path)

        encountered_internal_modules = set()

        for subproject in directories_for_subproject:

            #print(subproject, directories_for_subproject[subproject])

            # List the dependencies of the matching script
            dependencies_for_this = defaultdict(set)

            for dir_path in directories_for_subproject[subproject]:

                for file_path in fs.files_in_path(dir_path, extension="py", recursive=True):

                    #if subproject == "dustpedia": print(file_path)

                    introspection.add_dependencies(dependencies_for_this, file_path, encountered_internal_modules)

            dependencies_for_subproject[subproject] = dependencies_for_this.keys()

            for dependency in dependencies_for_this:
                for name in dependencies_for_this[dependency]:
                    dependencies[dependency].add(name)

    else: dependencies = introspection.get_all_dependencies()

# If a script name is given
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
        script_path = fs.join(introspection.pts_do_dir, matches[0][0], matches[0][1])

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

    try:
        from lxml import html
    except Exception: raise RuntimeError("You need the 'lxml' package to be able to specify the '--canopy' option")

    page_as_string = requests.get(url).content
    tree = html.fromstring(page_as_string)

    tables = [e for e in tree.iter() if e.tag == 'table']
    table = tables[-1]

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

data = dict()

present_dependencies = []
not_present_dependencies = []

# Loop over the packages and check their presence
for dependency in sorted(dependencies, key=str.lower):

    # Get the list of PTS scripts for this dependency
    script_list = dependencies[dependency]

    #print(script_list)

    # Skip packages from the standard library, unless the appropriate flag is enabled
    if introspection.is_std_lib(dependency) and not arguments.standard: continue

    # Check presency and version
    if dependency in packages:
        present = True
        version = packages[dependency] if arguments.version else None
    else:
        present = introspection.is_present_package(dependency)
        version = None

    if present: present_dependencies.append(dependency)
    else: not_present_dependencies.append(dependency)

    # Check whether the package is available through Canopy
    in_canopy = dependency.lower() in canopy_packages if canopy_packages is not None else None

    # Check whether the package is available through pip
    if arguments.pip:

        results = list(search(dependency))
        in_pip = False
        for result in results:
            if result["name"] == dependency:
                in_pip = True
                break
    else: in_pip = None

    # Get description
    if arguments.description:

        results = list(search(dependency))
        description = None
        for result in results:
            if result["name"] == dependency:
                description = result["summary"]
                if description.lower().startswith(dependency.lower() + ": "): description = description.split(": ")[1]
                break
    else: description = None

    data[dependency] = [script_list, present, version, in_canopy, in_pip, description]

if arguments.subprojects:

    #dependencies_for_subproject = defaultdict(set)
    #for dependency in data:
    #    script_list = data[dependency][0]
    #    for script_path in script_list:
    #        subproject = script_path.split("/pts/")[1].split("/")[0]
    #        if subproject == "do":
    #            subproject = script_path.split("do/")[1].split("/")[0]
    #        dependencies_for_subproject[subproject].add(dependency)

    print("")

    for subproject in dependencies_for_subproject:

        print(fmt.underlined + subproject + fmt.reset + ":")
        print("")

        for dependency in dependencies_for_subproject[subproject]:

            if dependency not in data: continue # dependencies that were not added to the data because they belong to the standard library

            script_list = data[dependency][0]
            present = data[dependency][1]
            version = data[dependency][2]
            in_canopy = data[dependency][3]
            in_pip = data[dependency][4]
            description = data[dependency][5]

            if present:
                print("   " + fmt.bold + fmt.green + dependency + fmt.reset)
            else: print("   " + fmt.bold + fmt.red + dependency + fmt.reset)

        print("")

else:

    print("")

    if len(present_dependencies) > 0:

        print(fmt.bold + fmt.green + "present:" + fmt.reset)
        print("")

        for dependency in present_dependencies:
            print_dependency(dependency, data)

        print("")

    if len(not_present_dependencies) > 0:

        print(fmt.bold + fmt.red + "not present:" + fmt.reset)
        print("")

        for dependency in not_present_dependencies:
            print_dependency(dependency, data)

        print("")

# -----------------------------------------------------------------
