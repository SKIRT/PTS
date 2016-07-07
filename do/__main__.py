#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# -----------------------------------------------------------------
# Main module for the do package
# -----------------------------------------------------------------

## \package pts.do.__main__ Execute one of the scripts in this package from any current directory.
#
# Before proceeding, ensure that your login script includes the extra lines as described in \ref InstallMacSetUp,
# and that you have logged in again after that change.
#
# To execute a python script named \c example.py residing in this directory, enter "pts example" in a Terminal window.
# This will work regardless of your current directory.
# You can include command line arguments as well, as in "pts example arg1 arg2".
# You can shorten the script name to the first few letters as long as there are no two scripts with matching names.
# For example "pts exa arg1 arg2" would still execute \c example.py, assuming that only one script has a name
# starting with the string "exa".
#
# Use "ipts" rather than "pts" to enter interactive python mode (with the >>> prompt) after executing the script.
# This is useful for testing or experimenting with pts functionality: the script imports the relevant
# pts module(s) and initializes some key objects that then can be used from the interactive python prompt.
#

# -----------------------------------------------------------------

# Import standard modules
import sys
import importlib
from operator import itemgetter

# Import the relevant PTS modules
from pts.core.tools import inspection
from pts.core.tools import filesystem as fs
from pts.core.tools import time

# -----------------------------------------------------------------

def show_all_available(scripts, tables=None):

    print("No match found. Available commands:")

    # Combine scripts and tables
    for subproject in tables:
        table = tables[subproject]
        for i in range(len(table["Command"])):
            scripts.append((subproject, table["Command"][i] + ".py", table["Description"][i]))

    # Sort on the 'do' subfolder name
    scripts = sorted(scripts, key=itemgetter(0))

    current_dir = None
    for script in scripts:

        description = "  " + script[2] if len(script) == 3 else ""

        if current_dir == script[0]:
            print(" " * len(current_dir) + "/" + script[1][:-3] + description)
        else:
            print(script[0] + "/" + script[1][:-3] + description)
            current_dir = script[0]

# -----------------------------------------------------------------

def show_possible_matches(matches, table_matches=None, tables=None):

    print("The command you provided is ambiguous. Possible matches:")

    # Combine script and table matches
    for subproject, index in table_matches:
        command = tables[subproject]["Command"][index]
        description = tables[subproject]["Description"][index]
        matches.append((subproject, command + ".py", description))

    # Sort on the 'do' subfolder name
    matches = sorted(matches, key=itemgetter(0))

    current_dir = None
    for script in matches:

        description = "  " + script[2] if len(script) == 3 else ""

        if current_dir == script[0]:
            print(" " * len(current_dir) + "/" + script[1][:-3] + description)
        else:
            print(script[0] + "/" + script[1][:-3] + description)
            current_dir = script[0]

# -----------------------------------------------------------------

# Get the name of the script to execute
script_name = sys.argv[1] if len(sys.argv) > 1 else None

# Find matches in existing do scripts
scripts = inspection.get_scripts()
tables = inspection.get_arguments_tables()
matches = inspection.find_matches_scripts(script_name, scripts)
table_matches = inspection.find_matches_tables(script_name, tables)

# No match
if len(matches) + len(table_matches) == 0: show_all_available(scripts, tables)

# If there is a unique match in an existing script, return it
elif len(matches) == 1 and len(table_matches) == 0:

    match = matches[0]

    # Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
    target = fs.join(inspection.pts_do_dir, match[0], match[1])
    sys.argv[0] = target
    del sys.argv[1]
    print "Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:])
    exec open(target)

# If there is an unique match in a table
elif len(table_matches) == 1 and len(matches) == 0:

    subproject, index = table_matches[0]
    command_name = tables[subproject]["Command"][index]
    class_path_relative = tables[subproject]["Path"][index]
    class_path = "pts." + subproject + "." + class_path_relative
    module_path, class_name = class_path.rsplit('.', 1)

    subproject_path = inspection.pts_subproject_dir(subproject)

    sys.argv[0] = fs.join(inspection.pts_root_dir, module_path.replace(".", "/") + ".py") # this is actually not necessary (and not really correct, it's not like we are calling the module where the class is..)
    del sys.argv[1] # but this is important

    module = importlib.import_module(module_path)

    cls = getattr(module, class_name)

    from pts.core.tools import logging

    # Import the configuration

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    configuration_module = importlib.import_module(configuration_module_path)

    config = getattr(configuration_module, "config")

    # Read the configuration settings from the command-line arguments
    config.read()

    # Determine the log file path
    logfile_path = fs.join(config.fixed["log_path"], time.unique_name("log") + ".txt") if config.arguments.report else None

    # Determine the log level
    level = "DEBUG" if config.arguments.debug else "INFO"

    # Initialize the logger
    log = logging.setup_log(level=level, path=logfile_path)
    log.start("Starting " + command_name + " ...")

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Run the instance
    inst.run()

# Show possible matches if there are more than just one
else: show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
