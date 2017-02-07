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
import argparse
import time as _time

# Import the relevant PTS modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import time
from pts.do.commandline import show_all_available, show_possible_matches, start_target
from pts.modeling.welcome import welcome as welcome_modeling
from pts.magic.welcome import welcome as welcome_magic
from pts.dustpedia.welcome import welcome as welcome_dustpedia
from pts.core.basics.configuration import create_configuration

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser(prog="pts")
parser.add_argument("do_command", type=str, help="the name of the PTS do command (preceeded by the subproject name and a slash if ambigious; i.e. 'subproject/do_command')", default=None, nargs='?')
parser.add_argument("--version", "-v", action="store_true", help="show the PTS version")
parser.add_argument("--interactive", action="store_true", help="use interactive mode for the configuration")
parser.add_argument("--arguments", action="store_true", help="use argument mode for the configuration")
parser.add_argument("--configfile", type=str, help="use a configuration file")
parser.add_argument("--rerun", action="store_true", help="use the last used configuration")
parser.add_argument("--remote", type=str, help="launch the PTS command remotely")
parser.add_argument("--keep", action="store_true", help="keep the remote output")
parser.add_argument("--input", type=str, help="the name/path of the input directory")
parser.add_argument("--output", type=str, help="the name/path of the output directory")
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

if args.version:

    version = introspection.pts_version()
    print(version)
    exit()

# Check input and output options, should be directories
if args.input is not None and not fs.is_directory(args.input): raise ValueError("Input path should be an existing directory")
if args.output is not None and not fs.is_directory(args.output): raise ValueError("Output path should be an existing directory")

# -----------------------------------------------------------------

# Get the name of the do script
script_name = args.do_command

# Determine the configuration method
configuration_method = None
if args.interactive: configuration_method = "interactive"
elif args.arguments: configuration_method = "arguments"
elif args.configfile is not None: configuration_method = "file:" + args.configfile
elif args.rerun: configuration_method = "last"

# Construct clean arguments list
sys.argv = ["pts", args.do_command] + args.options

# Find matches
matches = introspection.find_matches_scripts(script_name, scripts)
table_matches = introspection.find_matches_tables(script_name, tables)

# No match
if len(matches) + len(table_matches) == 0:

    from pts.core.tools import logging
    log = logging.setup_log()
    log.error("Unknown command: " + script_name)
    show_all_available(scripts, tables)

# If there is a unique match in an existing script, return it
elif len(matches) == 1 and len(table_matches) == 0:

    if args.remote is not None: raise ValueError("This do command cannot be executed remotely")

    match = matches[0]

    # Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
    target = fs.join(introspection.pts_do_dir, match[0], match[1])
    sys.argv[0] = target
    del sys.argv[1]
    print "Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:])

    command_name = match[1]

    # Set target
    def start(): exec open(target)

    # Start
    start_target(command_name, start)

# If there is an unique match in a table
elif len(table_matches) == 1 and len(matches) == 0:

    # Resolve
    subproject, index = table_matches[0]
    resolved = introspection.resolve_from_match(subproject, tables[subproject], index)

    # Get properties
    command_name = resolved.command_name
    hidden = resolved.hidden
    description = resolved.description
    module_path = resolved.module_path
    class_name = resolved.class_name
    configuration_method_table = resolved.configuration_method
    configuration_module_path = resolved.configuration_module_path
    subproject_path = introspection.pts_subproject_dir(subproject)

    # Set
    sys.argv[0] = fs.join(introspection.pts_root_dir, module_path.replace(".", "/") + ".py") # this is actually not necessary (and not really correct, it's not like we are calling the module where the class is..)
    del sys.argv[1] # but this is important

    # Welcome message
    if subproject == "modeling": welcome_modeling()
    elif subproject == "magic": welcome_magic()
    elif subproject == "dustpedia": welcome_dustpedia()

    # Get the class
    cls = introspection.get_class(module_path, class_name)

    # Import things
    from pts.core.tools import logging

    # Get the configuration definition
    definition = introspection.get_configuration_definition(configuration_module_path)

    # If not specified on the command line (before the command name), then use the default specified in the commands.dat file
    if configuration_method is None: configuration_method = configuration_method_table

    # Create the configuration
    config = create_configuration(definition, command_name, description, configuration_method)

    ## SAVE THE CONFIG
    if config.write_config:
        config_file_path = fs.join(config.config_dir_path(), command_name + ".cfg")
        config.saveto(config_file_path)
    ##

    if not args.rerun:

        if not fs.is_directory(introspection.pts_user_config_dir): fs.create_directory(introspection.pts_user_config_dir)

        # NEW: CACHE THE CONFIG
        config_cache_path = fs.join(introspection.pts_user_config_dir, command_name + ".cfg")
        config.saveto(config_cache_path)

    ##

    # Mark begin of modeling command for history
    if subproject == "modeling" and command_name != "setup" and command_name != "model_galaxy" and command_name != "model_sed":

        from ..modeling.component.component import load_modeling_history
        history = load_modeling_history(fs.cwd())
        history.add_entry(command_name)
        history.save()

    # If the PTS command has to be executed remotely
    if args.remote is not None:

        # Additional imports
        from pts.core.remote.remote import Remote

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
        start = _time.time()

        ## DO WHAT HAS TO BE DONE

        # Create the class instance, configure it with the configuration settings
        inst = cls(config)

        # Run the instance
        inst.run()

        # Record end time
        end = _time.time()
        seconds = end - start

        # Succesfully finished
        log.success("Finished " + command_name + " in " + str(seconds) + " seconds")

    # Mark the end of this modeling script
    if subproject == "modeling" and command_name != "setup" and command_name != "model_galaxy" and command_name != "model_sed":
        history.mark_end()
        history.save()

# Show possible matches if there are more than just one
else:

    # Show error
    from pts.core.tools import logging
    log = logging.setup_log()
    log.error("The command you provided is ambigious. Possible matches:")
    show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
