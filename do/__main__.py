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

# Import the relevant PTS modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import time

# -----------------------------------------------------------------

# Get the name of the script to execute
script_name = sys.argv[1] if len(sys.argv) > 1 else None

# Find matches in existing do scripts
scripts = introspection.get_scripts()
tables = introspection.get_arguments_tables()

if script_name is None:

    print("Welcome to PTS")
    introspection.show_all_available(scripts, tables)
    exit()

matches = introspection.find_matches_scripts(script_name, scripts)
table_matches = introspection.find_matches_tables(script_name, tables)

# No match
if len(matches) + len(table_matches) == 0: introspection.show_all_available(scripts, tables)

# If there is a unique match in an existing script, return it
elif len(matches) == 1 and len(table_matches) == 0:

    match = matches[0]

    # Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
    target = fs.join(introspection.pts_do_dir, match[0], match[1])
    sys.argv[0] = target
    del sys.argv[1]
    print "Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:])
    exec open(target)

# If there is an unique match in a table
elif len(table_matches) == 1 and len(matches) == 0:

    subproject, index = table_matches[0]
    command_name = tables[subproject]["Command"][index]
    description = tables[subproject]["Description"][index]
    class_path_relative = tables[subproject]["Path"][index]
    class_path = "pts." + subproject + "." + class_path_relative
    module_path, class_name = class_path.rsplit('.', 1)

    subproject_path = introspection.pts_subproject_dir(subproject)

    sys.argv[0] = fs.join(introspection.pts_root_dir, module_path.replace(".", "/") + ".py") # this is actually not necessary (and not really correct, it's not like we are calling the module where the class is..)
    del sys.argv[1] # but this is important

    module = importlib.import_module(module_path)

    cls = getattr(module, class_name)

    # Import things
    from pts.core.tools import logging
    from pts.core.basics.configuration import ConfigurationReader

    ## GET THE CONFIGURATION DEFINITION

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    configuration_module = importlib.import_module(configuration_module_path)

    definition = getattr(configuration_module, "definition")

    ## CREATE THE CONFIGURATION

    # Create the ConfigurationReader
    reader = ConfigurationReader(command_name, description, log_path="log")

    # Create the configuration from the definition and from reading the command line arguments
    config = reader.read(definition)

    ## SETUP LOGGER

    # Determine the log file path
    logfile_path = fs.join(config.log_path, time.unique_name("log") + ".txt") if config.report else None

    # Determine the log level
    level = "DEBUG" if config.debug else "INFO"

    # Initialize the logger
    log = logging.setup_log(level=level, path=logfile_path)
    log.start("Starting " + command_name + " ...")

    ## DO WHAT HAS TO BE DONE

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Run the instance
    inst.run()

# Show possible matches if there are more than just one
else: introspection.show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
