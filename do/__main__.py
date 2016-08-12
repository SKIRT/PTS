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
import argparse
import time as _time

# Import the relevant PTS modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import time, parsing
from pts.do.commandline import show_all_available, show_possible_matches

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser(prog="pts")
parser.add_argument("do_command", type=str, help="the name of the PTS do command (preceeded by the subproject name and a slash if ambigious; i.e. 'subproject/do_command')", default=None)
parser.add_argument("--interactive", action="store_true", help="use interactive mode for the configuration")
parser.add_argument("--arguments", action="store_true", help="use argument mode for the configuration")
parser.add_argument("--configfile", type=str, help="use a configuration file")
parser.add_argument("--remote", type=str, help="launch the PTS command remotely")
parser.add_argument("--keep", action="store_true", help="keep the remote output")
parser.add_argument("--input", type=parsing.directory_path, help="the name/path of the input directory")
parser.add_argument("--output", type=parsing.directory_path, help="the name/path of the output directory")
parser.add_argument("options", nargs=argparse.REMAINDER, help="options for the specific do command")

# -----------------------------------------------------------------

# Find possible PTS commands
scripts = introspection.get_scripts()
tables = introspection.get_arguments_tables()
if len(sys.argv) == 1: # nothing but 'pts' is provided

    print("")
    print("  ### Welcome to PTS ### ")
    print("")
    parser.print_help()
    print("")
    show_all_available(scripts, tables)
    parser.exit()

# -----------------------------------------------------------------

# Parse the command-line arguments
args = parser.parse_args()

# -----------------------------------------------------------------

# Get the name of the do script
script_name = args.do_command

# Determine the configuration method
configuration_method = None
if args.interactive: configuration_method = "interactive"
elif args.arguments: configuration_method = "arguments"
elif args.configfile is not None: configuration_method = "file:" + args.configfile

# Construct clean arguments list
sys.argv = ["pts", args.do_command] + args.options

# Find matches
matches = introspection.find_matches_scripts(script_name, scripts)
table_matches = introspection.find_matches_tables(script_name, tables)

# No match
if len(matches) + len(table_matches) == 0: show_all_available(scripts, tables)

# If there is a unique match in an existing script, return it
elif len(matches) == 1 and len(table_matches) == 0:

    if args.remote is not None: raise ValueError("This do command cannot be executed remotely")

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

    configuration_method_table = tables[subproject]["Configuration method"][index]

    subproject_path = introspection.pts_subproject_dir(subproject)

    sys.argv[0] = fs.join(introspection.pts_root_dir, module_path.replace(".", "/") + ".py") # this is actually not necessary (and not really correct, it's not like we are calling the module where the class is..)
    del sys.argv[1] # but this is important

    # Get the class of the configurable of which an instance has to be created
    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)

    # Import things
    from pts.core.tools import logging
    from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter, InteractiveConfigurationSetter, FileConfigurationSetter

    ## GET THE CONFIGURATION DEFINITION

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    try:
        configuration_module = importlib.import_module(configuration_module_path)
        #has_configuration = True
        definition = getattr(configuration_module, "definition")
    except ImportError:
        logging.log.warning("No configuration definition found for the " + class_name + " class")
        #has_configuration = False
        definition = ConfigurationDefinition() # Create new configuration definition

    ## CREATE THE CONFIGURATION

    # If not specified on the command line (before the command name), then use the default specified in the commands.dat file
    if configuration_method is None: configuration_method = configuration_method_table

    # Create the configuration setter
    if configuration_method == "interactive": setter = InteractiveConfigurationSetter(command_name, description)
    elif configuration_method == "arguments": setter = ArgumentConfigurationSetter(command_name, description)
    elif configuration_method.startswith("file"):
        configuration_filepath = configuration_method.split(":")[1]
        setter = FileConfigurationSetter(configuration_filepath, command_name, description)
    else: raise ValueError("Invalid configuration method: " + configuration_method)

    # Create the configuration from the definition and from reading the command line arguments
    config = setter.run(definition)

    ## SAVE THE CONFIG
    if config.write_config:
        config_file_path = fs.join(config.config_dir_path(), command_name + ".cfg")
        config.save(config_file_path)
    ##

    # If the PTS command has to be executed remotely
    if args.remote is not None:

        # Additional imports
        from pts.core.basics.remote import Remote

        # Exact command name
        exact_command_name = subproject + "/" + command_name

        # Determine the log level
        level = "DEBUG" if config.debug else "INFO"
        log = logging.setup_log(level=level)
        log.start("Starting " + exact_command_name + " on remote host " + args.remote + " ...")

        # Debugging
        log.debug("Initializing the remote ...")

        # Initialize the remote execution environment
        remote = Remote()
        remote.setup(args.remote)

        # Run PTS remotely
        task = remote.run_pts(exact_command_name, config, keep_remote_output=args.keep)

        # Succesfully submitted
        log.success("Succesfully submitted the PTS job to the remote host")

    # The PTS command has to be executed locally
    else:

        ## SETUP LOGGER

        # Determine the log file path
        logfile_path = fs.join(config.log_path, time.unique_name("log") + ".txt") if config.report else None

        # Determine the log level
        level = "DEBUG" if config.debug else "INFO"

        # Initialize the logger
        log = logging.setup_log(level=level, path=logfile_path)
        log.start("Starting " + command_name + " ...")

        # Record starting time
        start = _time.clock()

        ## DO WHAT HAS TO BE DONE

        # Create the class instance, configure it with the configuration settings
        inst = cls(config)

        # Run the instance
        inst.run()

        # Record end time
        end = _time.clock()
        seconds = end - start

        # Succesfully finished
        log.success("Finished " + command_name + " in " + str(seconds) + " seconds")

# Show possible matches if there are more than just one
else: show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
