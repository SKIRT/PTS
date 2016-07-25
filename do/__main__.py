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

# Import the relevant PTS modules
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs
from pts.core.tools import time
from pts.do.commandline import show_all_available, show_possible_matches

# -----------------------------------------------------------------

# Create the command-line parser
parser = argparse.ArgumentParser(prog="pts")
parser.add_argument("do_command", type=str, help="the name of the PTS do command (preceeded by the subproject name and a slash if ambigious; i.e. 'subproject/do_command')", default=None)
parser.add_argument("--interactive", action="store_true", help="use interactive mode for the configuration")
parser.add_argument("--arguments", action="store_true", help="use argument mode for the configuration")
parser.add_argument("--configfile", type=str, help="use a configuration file")
parser.add_argument("--remote", type=str, help="launch the PTS command remotely")
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

    module = importlib.import_module(module_path)

    cls = getattr(module, class_name)

    # Import things
    from pts.core.tools import logging
    from pts.core.basics.configuration import ArgumentConfigurationSetter, InteractiveConfigurationSetter, FileConfigurationSetter

    ## GET THE CONFIGURATION DEFINITION

    configuration_name = tables[subproject]["Configuration"][index]
    if configuration_name == "--": configuration_name = command_name
    configuration_module_path = "pts." + subproject + ".config." + configuration_name

    configuration_module = importlib.import_module(configuration_module_path)

    definition = getattr(configuration_module, "definition")

    ## CREATE THE CONFIGURATION

    # If not specified on the command line (before the command name), then use the default specified in the commands.dat file
    if configuration_method is None: configuration_method = configuration_method_table

    # Create the configuration setter
    if configuration_method == "interactive": setter = InteractiveConfigurationSetter(command_name, description, log_path="log")
    elif configuration_method == "arguments": setter = ArgumentConfigurationSetter(command_name, description, log_path="log")
    elif configuration_method.startswith("file"):
        configuration_filepath = configuration_method.split(":")[1]
        setter = FileConfigurationSetter(configuration_filepath, command_name, description, log_path="log")
    else: raise ValueError("Invalid configuration method: " + configuration_method)

    # Create the configuration from the definition and from reading the command line arguments
    config = setter.run(definition)

    # If the PTS command has to be executed remotely
    if args.remote is not None:

        # Additional imports
        import tempfile
        from pts.core.basics.remote import Remote

        unique_session_name = time.unique_name(command_name)

        ## SETUP LOGGER

        # Determine the log level
        level = "DEBUG" if config.debug else "INFO"
        log = logging.setup_log(level=level)
        log.start("Starting " + command_name + " on remote host " + args.remote + " ...")

        ##

        # Debugging
        log.debug("Initializing the remote ...")

        # Initialize the remote execution environment
        remote = Remote()
        remote.setup(args.remote)

        # Create a remote temporary directory
        remote_temp_path = remote.temp_directory

        ##

        ## CHANGE THE LOG PATH TO A REMOTE PATH

        # Always create a log file while executing remotely
        config.report = True
        config.log_path = fs.join(remote_temp_path, time.unique_name("log") + ".txt")

        ##

        # Debugging
        log.debug("Saving the configuration file locally ...")

        # Determine path to the temporarily saved local configuration file
        temp_path = tempfile.gettempdir()
        temp_conf_path = fs.join(temp_path, unique_session_name + ".cfg")

        # Save the configuration file to the temporary directory
        config.save(temp_conf_path)

        # Debugging
        log.debug("Uploading the configuration file to '" + remote_temp_path + "' ...")

        # Upload the config file
        remote_conf_path = fs.join(remote_temp_path, fs.name(temp_conf_path))
        remote.upload(temp_conf_path, remote_temp_path)

        # Remove the original config file
        fs.remove_file(temp_conf_path)

        # Debugging
        log.debug("Creating a script for remote execution ...")

        # Determine the path to the remote equivalent of this file
        remote_main_path = fs.join(remote.pts_package_path, "do", "__main__.py")

        # Create a bash script
        temp_script_path = fs.join(temp_path, unique_session_name + ".py")

        with open(temp_script_path, 'w') as script_file:

            script_file.write("#!/usr/bin/env python\n")
            script_file.write("# -*- coding: utf8 -*-\n")
            script_file.write("\n")
            script_file.write("python " + remote_main_path + " --configfile " + remote_conf_path + " " + command_name + "\n")

        #print(temp_path)
        #exit()

        # Execute the script
        remote.start_screen(unique_session_name, temp_script_path, remote_temp_path)

        # Remove the local script file
        fs.remove_file(temp_script_path)

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

        ## DO WHAT HAS TO BE DONE

        # Create the class instance, configure it with the configuration settings
        inst = cls(config)

        # Run the instance
        inst.run()

# Show possible matches if there are more than just one
else: show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------
