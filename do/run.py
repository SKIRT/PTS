#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.run Run a PTS command.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import sys
import importlib

# Import the relevant PTS classes and modules
from pts.core.tools import introspection
from pts.core.tools import time
from pts.core.tools import filesystem as fs
from .commandline import start_and_clear
from pts.modeling.welcome import welcome as welcome_modeling
from pts.magic.welcome import welcome as welcome_magic
from pts.evolve.welcome import welcome as welcome_evolve
from pts.dustpedia.welcome import welcome as welcome_dustpedia
from pts.core.basics.configuration import create_configuration
from pts.modeling.setup import setup as setup_modeling, finish as finish_modeling
from pts.modeling.setup import check_modeling_cwd
from pts.magic.setup import setup as setup_magic, finish as finish_magic
from pts.evolve.setup import setup as setup_evolve, finish as finish_evolve
from pts.dustpedia.setup import setup as setup_dustpedia, finish as finish_dustpedia
from pts.do.commandline import show_all_available, show_possible_matches, print_welcome
from pts.core.basics.filemonitor import FileMonitor

# -----------------------------------------------------------------

def no_input(parser, scripts, tables):

    """
    This function ...
    :param parser:
    :param scripts:
    :param tables:
    :return:
    """

    print_welcome()
    parser.print_help()
    print("")
    show_all_available(scripts, tables)
    print("")
    parser.exit()

# -----------------------------------------------------------------

def show_version():

    """
    This function ...
    :return:
    """

    version = introspection.pts_version()
    print(version)
    exit()

# -----------------------------------------------------------------

def no_match(script_name, scripts, tables):

    """
    This function ...
    :param script_name:
    :param scripts:
    :param tables:
    :return:
    """

    from pts.core.basics.log import setup_log
    log = setup_log()
    log.error("Unknown command: " + script_name)
    show_all_available(scripts, tables)

# -----------------------------------------------------------------

def multiple_matches(matches, table_matches, tables):

    """
    Thisfunction ...
    :param matches:
    :param table_matches:
    :param tables:
    :return:
    """

    # Show error
    from pts.core.basics.log import setup_log
    log = setup_log()
    log.error("The command you provided is ambigious. Possible matches:")
    show_possible_matches(matches, table_matches, tables)

# -----------------------------------------------------------------

def run_script(matches, args):

    """
    This function ...
    :param matches:
    :param args:
    :return:
    """

    if args.remote is not None: raise ValueError("This do command cannot be executed remotely")

    match = matches[0]

    # Execute the matching script, after adjusting the command line arguments so that it appears that the script was executed directly
    target = fs.join(introspection.pts_do_dir, match[0], match[1])
    sys.argv[0] = target
    del sys.argv[1]
    print("Executing: " + match[0] + "/" + match[1] + " " + " ".join(sys.argv[1:]))

    #command_name = match[1]

    # Set target
    #def start(): exec open(target)

    # Start # DOESN'T WORK WHEN THE SCRIPT FILE DEFINES A FUNCTION
    #start_target(command_name, start)
    exec open(target)

# -----------------------------------------------------------------

def run_configurable(table_matches, args, tables):

    """
    This function ...
    :param table_matches:
    :param args:
    :param tables:
    :return:
    """

    # Determine the configuration method
    configuration_method = None
    if args.interactive: configuration_method = "interactive"
    elif args.arguments: configuration_method = "arguments"
    elif args.configfile is not None: configuration_method = "file:" + args.configfile
    elif args.rerun: configuration_method = "last"

    # Regenerate the configuration method option
    if args.interactive: configuration_method_argument = "--interactive"
    elif args.arguments: configuration_method_argument = "--arguments"
    elif args.configfile is not None: configuration_method_argument = "--configfile '" + args.configfile + "'"
    elif args.rerun: configuration_method_argument = "--rerun"
    else: configuration_method_argument = ""

    # Resolve
    subproject, index = table_matches[0]
    resolved = introspection.resolve_from_match(subproject, tables[subproject], index)

    # Get properties
    title = resolved.title
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

    # Get a list of the leftover arguments
    leftover_arguments = sys.argv[1:]

    # Welcome message
    if subproject == "modeling": welcome_modeling()
    elif subproject == "magic": welcome_magic()
    elif subproject == "dustpedia": welcome_dustpedia()
    elif subproject == "evolve": welcome_evolve()

    # Special
    if subproject == "modeling": check_modeling_cwd(command_name, fs.cwd())

    # Get the configuration definition
    definition = introspection.get_configuration_definition_pts_not_yet_in_pythonpath(configuration_module_path)

    # If not specified on the command line (before the command name), then use the default specified in the commands.dat file
    if configuration_method is None: configuration_method = configuration_method_table

    # Check whether arguments are passed and the configuration method is interactive
    if configuration_method == "interactive" and len(leftover_arguments) > 0: raise ValueError("Arguments on the command-line are not supported by default for this command. Run with pts --arguments to change this behaviour.")

    # Create the configuration
    config = create_configuration(definition, command_name, description, configuration_method)

    ## SAVE THE CONFIG if requested
    if config.write_config:
        config_filepath = config.config_file_path(command_name)
        config.saveto(config_filepath)
    else: config_filepath = None

    # If this is not a re-run
    if not args.rerun:
        if not fs.is_directory(introspection.pts_user_config_dir): fs.create_directory(introspection.pts_user_config_dir)
        # CACHE THE CONFIG
        config_cache_path = fs.join(introspection.pts_user_config_dir, command_name + ".cfg")
        config.saveto(config_cache_path)

    # Setup function
    if subproject == "modeling": setup_modeling(command_name, fs.cwd(), configuration_method_argument)
    elif subproject == "magic": setup_magic(command_name, fs.cwd())
    elif subproject == "dustpedia": setup_dustpedia(command_name, fs.cwd())
    elif subproject == "evolve": setup_evolve(command_name, fs.cwd())

    # Initialize the logger
    log = initialize_pts(config, remote=args.remote, command_name=command_name)

    # Exact command name
    exact_command_name = subproject + "/" + command_name

    # If the PTS command has to be executed remotely
    if args.remote is not None: run_remotely(exact_command_name, config, args.keep, args.remote, log)

    # The PTS command has to be executed locally
    else: run_locally(exact_command_name, module_path, class_name, config, args.input_files, args.output_files, args.output, log)

    # Finish function
    if subproject == "modeling": finish_modeling(command_name, fs.cwd(), config_path=config_filepath)
    elif subproject == "magic": finish_magic(command_name, fs.cwd())
    elif subproject == "dustpedia": finish_dustpedia(command_name, fs.cwd())
    elif subproject == "evolve": finish_evolve(command_name, fs.cwd())

# -----------------------------------------------------------------

def initialize_pts(config, remote=None, command_name=None):

    """
    This function ...
    :param config:
    :param remote:
    :param command_name:
    :return:
    """

    # Initialize the logger
    log = initialize_log(config, remote=remote, command_name=command_name)

    # Initialize the file monitor
    if log.is_debug:
        monitor = FileMonitor(short=True)
        monitor.patch()

    # Return the log
    return log

# -----------------------------------------------------------------

def initialize_log(config, remote=None, command_name=None):

    """
    This function ...
    :param config:
    :parma remote:
    :param command_name:
    :return:
    """

    # Determine the log level
    level = "INFO"
    if config.debug: level = "DEBUG"
    if config.brief: level = "SUCCESS"

    # Determine memuse flag
    if hasattr(config, "memuse"): memuse = config.memuse
    else: memuse = False

    # Determine log path
    #if args.remote is None: logfile_path = fs.join(config.log_path, time.unique_name("log") + ".txt") if config.report else None
    #else: logfile_path = None

    # Determine name of log file
    if command_name is not None: filename = command_name
    else: filename = "log"

    # Determine the log file path
    if remote is None: logfile_path = fs.join(config.log_path, time.unique_name(filename, separator="__") + ".txt") if config.report else None
    else: logfile_path = None

    # Initialize the logger
    from ..core.basics.log import setup_log
    log = setup_log(level=level, path=logfile_path, memory=memuse)

    # Return the logger
    return log

# -----------------------------------------------------------------

def run_locally(command_name, module_path, class_name, config, input_files, output_files, output, log):

    """
    This function ...
    :param command_name: 
    :param module_path:
    :param class_name:
    :param config:
    :param input_files:
    :param output_files:
    :param output:
    :param log:
    :return: 
    """

    # Start message
    log.start("Starting " + command_name + " ...")

    # Get the class
    cls = introspection.get_class(module_path, class_name)

    # Create the class instance, configure it with the configuration settings
    inst = cls(config)

    # Set input files
    input_dict = {}
    if input_files is not None:

        # Loop over the names of the input variables
        for name in input_files:

            # Get class path
            classpath, filepath = input_files[name]
            modulepath, classname = classpath.rsplit(".", 1)

            # Get input class
            input_module = importlib.import_module(modulepath)
            input_class = getattr(input_module, classname)

            # Open the input file
            input_object = input_class(filepath)

            # Set to input dict
            input_dict[name] = input_object

    # Start
    start_and_clear(command_name, inst.run, **input_dict)

    # Write output files
    if output_files is not None:

        types = dict()

        # Loop over the names of the attributes for output
        for name in output_files:

            # Get filepath
            filepath = output_files[name]

            # Get the output object
            output_object = getattr(inst, name)

            # Set the type
            types[name] = type(output_object).__module__ + "." + type(output_object).__name__

            # Save the output object
            real_filepath = filepath + "." + type(output_object).default_extension
            output_object.saveto(real_filepath)

        # Save the types
        from pts.core.tools import serialization
        types_path = fs.join(output, "types.dat")
        serialization.write_dict(types, types_path)

# -----------------------------------------------------------------

def run_remotely(command_name, config, keep, host_id, log):

    """
    This function ...
    :param command_name:
    :param config:
    :param keep:
    :param host_id:
    :param log:
    :return: 
    """

    # Additional imports
    from pts.core.remote.remote import Remote

    # Start message
    log.start("Starting " + command_name + " on remote host " + host_id + " ...")

    # Debugging
    log.debug("Initializing the remote ...")

    # Initialize the remote execution environment
    remote = Remote(host_id=host_id)

    # Run PTS remotely
    task = remote.run_pts(command_name, config, keep_remote_output=keep)

    # Succesfully submitted
    log.success("Succesfully submitted the PTS job to the remote host")

# -----------------------------------------------------------------
