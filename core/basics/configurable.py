#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.configurable Contains the Configurable class, a class for representing classes that can be
#  configured with a configuration file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import traceback
from collections import OrderedDict
from abc import ABCMeta, abstractproperty, abstractmethod

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from .configuration import find_command
from .log import log
from ..tools.utils import lazyproperty

# -----------------------------------------------------------------

itemize_symbols = ["-", "*", ">"]

# -----------------------------------------------------------------

def write_input(input_dict, path, light=False):

    """
    This function ...
    :param input_dict:
    :param path:
    :param light:
    :return:
    """

    # Import things
    from ..tools.serialization import write_dict
    from ..tools import introspection

    # Dictionary for remainder input (not saved as file)
    remainder = dict()

    # Dictionary for the class paths
    classes = dict()

    # Loop over the input items
    for name in input_dict:

        # Get the class
        cls = input_dict[name].__class__

        # Should be saved as a file
        if hasattr(cls, "default_extension"):

            # Determine the filename
            filename = name + "." + input_dict[name].default_extension
            filepath = fs.join(path, filename)  # local temporary path

            # If 'light' is enabled, don't write out files that can't get stringifyied
            if not light:

                # Save the object, but don't change its internal path
                original_path = input_dict[name].path
                input_dict[name].saveto(filepath)
                input_dict[name].path = original_path

            subproject, relative_class_subproject = introspection.get_class_path(input_dict[name].__class__)
            classes[name] = subproject + "." + relative_class_subproject

        # Add to remainder dict
        else: remainder[name] = input_dict[name]

    # Write the remainder dictionary
    remainder_path = fs.join(path, "input.dat")
    write_dict(remainder, remainder_path)

    # Write the classes dictionary, if necessary
    if len(classes) > 0:
        classes_path = fs.join(path, "classes.dat")
        write_dict(classes, classes_path)

# -----------------------------------------------------------------

class NotSavedPlaceHolder(object):

    """
    This class ...
    """

    def __init__(self, filepath):

        """
        This function ...
        :param filepath:
        """

        self.filepath = filepath

# -----------------------------------------------------------------

def load_input(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Import things
    from ..tools.serialization import load_dict
    from ..tools import introspection

    input_dict = dict()

    # Add the input.dat input
    input_file_path = fs.join(path, "input.dat")
    if fs.is_file(input_file_path):
        remainder = load_dict(input_file_path)
        for key in remainder: input_dict[key] = remainder[key]

    # Load the classes data
    classes_path = fs.join(path, "classes.dat")
    #classes = load_dict(classes_path)

    if fs.is_file(classes_path):

        classes = load_dict(classes_path)

        # Loop over the names
        for name in classes:

            # Get the class path
            class_path = classes[name]

            # Get the class
            cls = introspection.get_class_from_path(class_path)

            # Get the default extension
            filename = name + "." + cls.default_extension

            # Determine the full file path
            filepath = fs.join(path, filename)

            # Check whether present, and load
            if fs.is_file(filepath): input_dict[name] = cls.from_file(filepath)
            else: input_dict[name] = NotSavedPlaceHolder(filepath) # print("Input file '" + filepath + "' not present")

    # Return the input
    return input_dict

# -----------------------------------------------------------------

class Configurable(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta
    _log_section = None

    # -----------------------------------------------------------------

    _command_name = None

    # -----------------------------------------------------------------

    @classmethod
    def command_name(cls):

        """
        This function ...
        :return:
        """

        if cls._command_name is not None: return cls._command_name

        # Find the corresponding command
        command_name, class_name, configuration_module_path, description = find_command(cls)

        # Check if found
        if command_name is None: raise ValueError("Could not find the command name for the " + cls.__name__ + " class")

        # Set the command name
        cls._command_name = command_name

        # Return the command
        return command_name

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set configuration
        if kwargs.pop("no_config", False): self.config = None
        else: self.config = self.get_config(*args, **kwargs)

        # Set the detached calculations flag
        self.detached = False

    # -----------------------------------------------------------------

    @abstractmethod
    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Call the implementation
        self._run(**kwargs)

        # 3. Finish
        self.finish()

    # -----------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :return:
        """

        if self._log_section is not None: log.remove_subsection()

    # -----------------------------------------------------------------

    @classmethod
    def get_config(cls, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get config from positional args
        if len(args) == 0: config = None
        elif len(args) == 1: config = args[0]
        else: raise ValueError("Can only one positional argument, which is the configuration (dictionary)")

        # Get config from kwargs
        if "config" in kwargs:
            if config is not None: raise ValueError("Config was passed as positional argument, so cannot be passed as keyword argument as well")
            else: config = kwargs.pop("config")

        # Get other settings from kwargs
        interactive = kwargs.pop("interactive", False)
        unlisted = kwargs.pop("unlisted", False)
        cwd = kwargs.pop("cwd", None)
        prompt_optional = kwargs.pop("prompt_optional", None)
        use_default = kwargs.pop("use_default", None)
        check_required = kwargs.pop("check_required", True)
        extend_config = kwargs.pop("extend_config", False)

        from .configuration import get_config_for_class

        # Get the config
        if unlisted:
            from .map import Map
            assert isinstance(config, Map)
            return config
        else: return get_config_for_class(cls, config, interactive=interactive, cwd=cwd, prompt_optional=prompt_optional,
                                          use_default=use_default, check_required=check_required, extend_config=extend_config)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Set the log section
        if self._log_section is not None: log.add_subsection(self._log_section)

        # NEW: WRITE THE CONFIGURATION
        if self.config.write_config: self.config.saveto(self.config.config_file_path(self.command_name()))

    # -----------------------------------------------------------------

    @property
    def class_name(self):

        """
        This function ...
        :return:
        """

        name = type(self).__name__
        return name

    # -----------------------------------------------------------------

    @property ## I THINK THIS FUNCTION CAN BE REMOVED (IT SHOULD) AND REPLACED BY CLASS_NAME
    def name(self):

        """
        This function ...
        :return:
        """

        name = type(self).__name__.lower()
        if "plotter" in name: return "plotter"
        else: return name

    # -----------------------------------------------------------------

    @property
    def input_path(self):

        """
        This function ...
        :return:
        """

        # If 'input' defined in the config
        if "input" in self.config:

            full_input_path = fs.absolute_or_in(self.config.input, self.config.path)
            if not fs.is_directory(full_input_path): raise ValueError("The input directory does not exist")
            return full_input_path

        # Else, use the working directory as input directory
        else: return self.config.path

    # -----------------------------------------------------------------

    def input_path_file(self, name, check=False):

        """
        This function ...
        :param name:
        :param check:
        :return:
        """

        path = fs.join(self.input_path, name)
        if check and not fs.is_file(path): raise ValueError("Input file '" + name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    def input_path_directory(self, name, check=False):

        """
        This function ...
        :param name:
        :param check:
        :return:
        """

        path = fs.join(self.input_path, name)
        if check and not fs.is_directory(path): raise ValueError("Input directory '" + name + "' does not exist")
        return path

    # -----------------------------------------------------------------

    @property
    def output_path(self):

        """
        This function ...
        :return:
        """

        # If 'output' is defined in the config
        if "output" in self.config and self.config.output is not None:

            full_output_path = fs.absolute_or_in(self.config.output, self.config.path)
            if not fs.is_directory(full_output_path): fs.create_directory(full_output_path)
            return full_output_path

        # Else, use the working directory as output directory
        else: return self.config.path

    # -----------------------------------------------------------------

    def output_path_file(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.output_path, name)

    # -----------------------------------------------------------------

    def output_path_directory(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        # Create and return path
        if create: return fs.create_directory_in(self.output_path, name)
        else: return fs.join(self.output_path, name)

# -----------------------------------------------------------------

class HierarchicConfigurable(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(HierarchicConfigurable, self).__init__(*args, **kwargs)

        # The children
        self.children = dict()

    # -----------------------------------------------------------------

    def __getattr__(self, attr):

        """
        This function ...
        Overriding __getattr__ should be fine (will not break the default behaviour) -- __getattr__ is only called
        as a last resort i.e. if there are no attributes in the instance that match the name.
        :param attr:
        :return:
        """

        if attr.startswith("__") and attr.endswith("__"): raise AttributeError("Can't delegate this attribute")
        return self.children[attr]

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Delete its children
        self.children = dict()

    # -----------------------------------------------------------------

    def add_child(self, name, type, config=None):

        """
        This function ...
        :param name:
        :param type:
        :param config:
        :return:
        """

        if name in self.children: raise ValueError("Child with this name already exists")

        # new ...
        if config is None: config = {}
        config["output_path"] = self.config.output_path
        config["input_path"] = self.config.input_path

        self.children[name] = type(config)

    # -----------------------------------------------------------------

    def setup_before(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def setup_after(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

class InvalidCommandError(Exception):

    """
    This class ...
    """

    def __init__(self, message, command):

        """
        Thisf unction ...
        :param message:
        :param command:
        """

        # Call the base class constructor with the parameters it needs
        super(InvalidCommandError, self).__init__(message)

        # The command
        self.command = command

# -----------------------------------------------------------------

class InteractiveConfigurable(Configurable):

    """
    This class ...
    """

    # Define class properties
    __metaclass__ = ABCMeta
    _commands = None

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(InteractiveConfigurable, self).__init__(*args, **kwargs)

        # The commands that have been executed
        self.commands = []

        # The commands that have failed to be executed
        self.failed = []

    # -----------------------------------------------------------------

    def run_commands(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Running commands ...")

        # Loop over the commands
        for command in self.config.commands:

            # Debugging
            log.debug("Running '" + command + "' ...")

            # Process command, give error if fails
            success = self.try_command(command)

            # Add command
            if success: self.commands.append(command)
            else: self.failed.append(command)

    # -----------------------------------------------------------------

    def interactive(self):

        """
        This function ...
        :return:
        """

        from .configuration import prompt_string

        # Inform the user
        log.info("Entering interactive mode ...")

        # Enter loop
        while True:

            # Get next command, break if no command is given
            command = prompt_string("command", "command to be executed")
            if not command: break

            # DEVELOPER COMMAND
            if command.startswith("$"):
                self._run_developer_command(command)
                continue

            # REPEAT LAST COMMAND
            if command == "last": command = self.last_command

            # Process command
            success = self.try_command(command)

            # Add command
            if success: self.commands.append(command)
            else: self.failed.append(command)

    # -----------------------------------------------------------------

    @property
    def nhistory_commands(self):

        """
        This function ...
        :return:
        """

        return len(self._history)

    # -----------------------------------------------------------------

    @property
    def has_history_commands(self):

        """
        This function ...
        :return:
        """

        return self.nhistory_commands > 0

    # -----------------------------------------------------------------

    @property
    def ncommands(self):

        """
        This function ...
        :return:
        """

        return len(self.commands)

    # -----------------------------------------------------------------

    @property
    def nfailed(self):

        """
        This function ...
        :return:
        """

        return len(self.failed)

    # -----------------------------------------------------------------

    @property
    def has_commands(self):

        """
        This function ...
        :return:
        """

        return self.ncommands > 0

    # -----------------------------------------------------------------

    @property
    def has_failed(self):

        """
        This function ...
        :return:
        """

        return self.nfailed > 0

    # -----------------------------------------------------------------

    def try_command(self, command):

        """
        This function ...
        :return:
        """

        success = True

        # Try
        try: self.process_command(command)

        # Invalid command
        except InvalidCommandError as e:

            if "Invalid command" in e.message: log.warning("Invalid command: '" + e.command + "'")
            else: log.warning("Invalid command: '" + e.command + "': " + e.message)

            # NO SUCCESS
            success = False

        # Other exception
        except Exception as e:

            message = str(e)

            # Too few arguments
            if "too few arguments" in message:

                log.error("Too few arguments")
                # Get first word == key
                key = command.split(" ")[0]
                if self.has_subcommands(key):

                    from ..tools import strings
                    splitted = strings.split_except_within_double_quotes(command, add_quotes=False)

                    if len(splitted) > 1: usage = self.get_usage_for_key(splitted[1], self._commands[key])
                    else:
                        # Set subcommands with descriptions
                        subcommands = self.get_subcommands(key)
                        subcommands_descriptions = OrderedDict()
                        for subkey in subcommands:
                            function_name, pass_command, description, subject = subcommands[subkey]
                            subcommands_descriptions[subkey] = description
                        from .configuration import ConfigurationDefinition, get_usage
                        definition = ConfigurationDefinition(write_config=False)
                        definition.add_required("subcommand", "string", "subcommand", choices=subcommands_descriptions)
                        usage = get_usage(key, definition, add_logging=False, add_cwd=False)

                    for line in usage: log.error(line)

                else:

                    usage = self.get_usage_for_key(key, self._commands)
                    for line in usage: log.error(line)

            # Invalid choice
            elif "invalid choice" in message:

                log.error("Invalid choice")
                log.error(str(e))

            # Shown unknown error
            else:

                traceback.print_exc()
                log.error(str(e))

            # NO SUCCESS
            success = False

        # Return the success flag
        return success

    # -----------------------------------------------------------------

    def process_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Run the command
        self._run_command(command, self._commands)

    # -----------------------------------------------------------------

    def _run_command(self, command, cmds):

        """
        This function ...
        :param command:
        :param cmds:
        :return:
        """

        # Get first word
        first = command.split(" ")[0]

        # Check whether interactive
        if first.startswith("*"):
            interactive = True
            command = command[1:]
            first = first[1:]
        else: interactive = False

        # Find key
        if first not in cmds: raise InvalidCommandError("Invalid command: '" + first + "'", command)
        key = first

        # Has subcommands
        if self.has_subcommands(key): self._run_subcommand(command, self._commands[key])

        # Regular command
        else: self._run_nosubcommand_impl(command, cmds, interactive=interactive)

    # -----------------------------------------------------------------

    def _run_developer_command(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        # Clean
        command = command[1:].strip()

        # Create evaluation command
        eval_command = "self." + command

        # Evaluate
        result = eval(eval_command)

        # Show the result, if applicable
        if result is not None: print(result)

    # -----------------------------------------------------------------

    def _run_subcommand(self, command, subcommands, interactive=False):

        """
        This function ...
        :param command:
        :param subcommands:
        :param interactive:
        :return:
        """

        from ..tools import strings
        from ..tools import formatting as fmt
        from .configuration import ConfigurationDefinition, get_help, prompt_string

        # Get first word == key
        key = command.split(" ")[0]

        # Get the possible subcommands
        #subcommands = self.get_subcommands(key)

        # Get command without main command
        subcommand = strings.split_at_first(command, key)[1].strip()
        #print("SUBCOMMAND", subcommand)
        #print("SUBCOMMANDS", subcommands.keys())
        #print(strings.startswith_any(subcommand, subcommands))

        # No subcommand?
        if not strings.startswith_any(subcommand, subcommands):

            # Show help for the main command
            #if "-h" in subcommand or "--help" in subcommand:
            if "--help" in subcommand: # WE DO NOT SUPPORT -h ANYMORE, E.G. IT CAN BE IN '--host'

                # Set subcommands with descriptions
                subcommands_descriptions = OrderedDict()
                for subkey in subcommands:
                    function_name, pass_command, description, subject = subcommands[subkey]
                    subcommands_descriptions[subkey] = description

                # Create definition for the subcommands
                definition = ConfigurationDefinition(write_config=False)
                definition.add_required("subcommand", "string", "subcommand", choices=subcommands_descriptions)
                help = get_help(key, definition, add_logging=False, add_cwd=False)

                # Show help
                for line in help: print(fmt.red + line + fmt.reset)

            # Nothing more than the main command is given as input
            elif subcommand == "":

                # Interactive mode: prompt between the different subcommands
                if interactive:

                    # Set subcommands with descriptions
                    subcommands_descriptions = OrderedDict()
                    for subkey in subcommands:
                        function_name, pass_command, description, subject = subcommands[subkey]
                        subcommands_descriptions[subkey] = description

                    # Prompt for the subcommand
                    subcommand = prompt_string("subcommand", "subcommand", choices=subcommands_descriptions)

                    # Run the subcommand in interactive mode
                    self._run_nosubcommand_impl(subcommand, subcommands, main_command=key, interactive=True)

                # Not enough input
                else: raise InvalidCommandError("Not enough input for '" + key + "' command", command)

            # Invalid command
            else: raise InvalidCommandError("Invalid command: '" + subcommand + "'", command)

        # Run the command
        else: self._run_command_impl(subcommand, subcommands, main_command=key, interactive=interactive)

    # -----------------------------------------------------------------

    def _run_command_impl(self, command, cmds, main_command=None, interactive=False):

        """
        This function ...
        :param command:
        :param cmds:
        :param main_command:
        :param interactive:
        :return:
        """

        #print("cmd", command)
        #print("cmds", cmds)

        # Get first word == key
        key = command.split(" ")[0]
        has_subcommands = isinstance(cmds[key], dict)

        # Has subcommands
        if has_subcommands: self._run_subcommand(command, cmds[key])

        # Regular command
        else: self._run_nosubcommand_impl(command, cmds, main_command=main_command, interactive=interactive)

    # -----------------------------------------------------------------

    def _run_nosubcommand_impl(self, command, cmds, main_command=None, interactive=False):

        """
        This function ...
        :param command:
        :param cmds:
        :param main_command:
        :param interactive:
        :return:
        """

        from ..tools import formatting as fmt

        # Get first word == key
        key = command.split(" ")[0]

        #print(key)
        #print(cmds[key])

        # Get function name and description
        function_name, pass_command, description, subject = cmds[key]
        #print(function_name, pass_command, description, subject)

        # Show help for the command
        #if "-h" in command or "--help" in command:
        if "--help" in command:  # WE DO NOT SUPPORT -h ANYMORE, E.G. IT CAN BE IN '--host'

            # Get the help info
            help = self.get_help_for_key(key, cmds, main_command=main_command)

            # Show help
            if help is None: print(fmt.red + "no input required" + fmt.reset)
            else:
                for line in help: print(fmt.red + line + fmt.reset)

        # Actually run
        else:

            # Get the function
            function = getattr(self, function_name)

            # Call the function
            if pass_command: function(command, interactive=interactive)
            else: function(interactive=interactive)

    # -----------------------------------------------------------------

    def has_subcommands(self, command):
        return isinstance(self._commands[command], dict)

    # -----------------------------------------------------------------

    def get_subcommands(self, command):

        """
        This function ...
        :param command:
        :return:
        """

        #if command not in self._subcommands: raise InvalidCommandError("Invalid command: '" + command + "'", command)
        if command not in self._commands: raise InvalidCommandError("Invalid command: '" + command + "'", command)
        #return self._subcommands[command]
        return self._commands[command]

    # -----------------------------------------------------------------

    def get_config_from_definition(self, name, definition, command, interactive=False):

        """
        This function ...
        :param name:
        :param definition:
        :param command:
        :param interactive:
        :return:
        """

        from .configuration import prompt_settings, parse_arguments

        # Get settings interactively
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=command, error="exception", exit_on_help=False, initialize=False, add_logging=False, add_cwd=False)

        # Remove path
        config.pop("_path")

        # Return the configuration
        return config

    # -----------------------------------------------------------------

    def get_definition_for_function_name(self, function_name):

        """
        This function ...
        :param function_name:
        :return:
        """

        from ..tools import strings

        # Get definition
        definition_property_name = strings.split_at_last(function_name, "_command")[0] + "_definition"
        definition = getattr(self, definition_property_name, None)

        # Return
        return definition

    # -----------------------------------------------------------------

    def get_kwargs_for_function_name(self, function_name):

        """
        This function ...
        :param function_name:
        :return:
        """

        # Get kwargs
        kwargs_property_name = function_name.split("_command")[0] + "_kwargs"
        kwargs = getattr(self, kwargs_property_name, {})

        # Return
        return kwargs

    # -----------------------------------------------------------------

    def get_help_for_key(self, key, cmds, main_command=None):

        """
        This function ...
        :param key:
        :param cmds:
        :param main_command:
        :return:
        """

        from .configuration import get_help

        # Get properties
        function_name, pass_command, description, subject = cmds[key]
        if subject is None: return None # no help for this command (simple command that doesn't need input)

        # Get the definition
        definition = self.get_definition_for_function_name(function_name)

        # Get the kwargs
        kwargs = self.get_kwargs_for_function_name(function_name)

        # Set command name
        if main_command is None: name = key
        else: name = main_command + " " + key

        # Get help lines
        if subject is None: help = self.get_help_command(definition, name=name)
        else:

            # Get the definition automatically
            get_definition_function_name = "get_" + subject + "_command_definition"
            get_definition_function = getattr(self, get_definition_function_name, None)
            if get_definition_function is None: raise ValueError("Invalid subject '" + subject + "'")
            definition = get_definition_function(definition, **kwargs)

            # Get the help
            help = get_help(name, definition, add_logging=False, add_cwd=False)

        # Return the help info
        return help

    # -----------------------------------------------------------------

    def get_usage_for_key(self, key, cmds, main_command=None):

        """
        This function ...
        :param key:
        :param cmds:
        :param main_command:
        :return:
        """

        from .configuration import get_usage

        # Get properties
        function_name, pass_command, description, subject = cmds[key]
        #if subject is None: return None # no usage for this command (simple command that doesn't need input)

        # Get the definition
        #print(function_name)
        definition = self.get_definition_for_function_name(function_name)

        # Get the kwargs
        kwargs = self.get_kwargs_for_function_name(function_name)

        # Set command name
        if main_command is None: name = key
        else: name = main_command + " " + key

        #print(key, definition)

        # Get usage lines
        if subject is None: usage = self.get_usage_command(definition, name=name)
        else:

            # Get the definition
            get_definition_function_name = "get_" + subject + "_command_definition"
            get_definition_function = getattr(self, get_definition_function_name, None)
            if get_definition_function is None: raise ValueError("Invalid subject '" + subject + "'")
            definition = get_definition_function(definition, **kwargs)

            # Get the usage
            usage = get_usage(name, definition, add_logging=False, add_cwd=False)

        # Return the usage info
        return usage

    # -----------------------------------------------------------------

    def get_usage_command(self, definition, name):

        """
        This function ...
        :param definition:
        :param name:
        :return:
        """

        from .configuration import get_usage

        # Return the usage
        return get_usage(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def get_help_command(self, definition, name):

        """
        This function ...
        :param definition:
        :param name:
        :return:
        """

        from .configuration import get_help

        # Return the usage
        return get_help(name, definition, add_logging=False, add_cwd=False)

    # -----------------------------------------------------------------

    def show_help(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Inform the user
        log.info("Showing help ...")

        print("")

        # Loop over the commands
        for key in self._commands:

            # Show help
            self._show_help_impl(self._commands[key], key)

        print("")

    # -----------------------------------------------------------------

    def _show_help_impl(self, spec, key):

        """
        This function ...
        :param spec:
        :param key:
        :return:
        """

        from .configuration import ConfigurationDefinition, get_usage
        from ..tools import formatting as fmt

        # Get properties
        if isinstance(spec, dict):
            function_name = pass_command = None
            if hasattr(spec, "subject"): subject = spec.subject
            else: subject = None
            if hasattr(spec, "description"): description = spec.description
            else: description = "no description"
        else: function_name, pass_command, description, subject = spec

        # Show
        print(" - " + fmt.bold + key + fmt.reset + ": " + description)

        #if function_name is None and not self.has_subcommands(key): raise RuntimeError("Something went wrong")

        # Show usage
        if pass_command:

            # Get usage
            usage = self.get_usage_for_key(key, self._commands)

            # Show
            for line in usage: print("    " + fmt.blue + line + fmt.reset)

        # Command with subcommands
        elif self.has_subcommands(key):

            # Set subcommands with descriptions
            subcommands = self.get_subcommands(key)
            subcommands_descriptions = OrderedDict()

            # Get description for each subcommand
            for subkey in subcommands:
                sspec = subcommands[subkey]
                if isinstance(sspec, dict): description = sspec.description if hasattr(sspec, "description") else "no description"
                else: function_name, pass_command, description, subject = subcommands[subkey]
                subcommands_descriptions[subkey] = description

            # Create definition for the subcommands
            definition = ConfigurationDefinition(write_config=False)
            definition.add_required("subcommand", "string", "subcommand", choices=subcommands_descriptions)
            usage = get_usage(key, definition, add_logging=False, add_cwd=False)

            # Show usage
            for line in usage: print("    " + fmt.blue + line + fmt.reset)

            # Show help for subcommands
            self._show_help_subcommands(subcommands, key)

        # No input needed
        else: print("    " + fmt.blue + "no input" + fmt.reset)

    # -----------------------------------------------------------------

    def _show_help_subcommands(self, subcommands, main_command, level=1):

        """
        This function ...
        :param subcommands:
        :param main_command:
        :param level:
        :return:
        """

        from ..tools import formatting as fmt

        prefix = "   " * level
        symbol = itemize_symbols[level]

        # Loop over the commands
        for key in subcommands:

            spec = subcommands[key]

            # Get description
            if isinstance(spec, dict):
                function_name = pass_command = None
                if hasattr(spec, "subject"): subject = spec.subject
                else: subject = None
                if hasattr(spec, "description"): description = spec.description
                else: description = "no description"
            else: function_name, pass_command, description, subject = subcommands[key]

            # Show
            print(prefix + " " + symbol + " " + fmt.bold + key + fmt.reset + ": " + description)

            # Show usage
            # if not pass_command or subject is None: continue
            if pass_command:

                # Get usage
                usage = self.get_usage_for_key(key, subcommands, main_command=main_command)

                # Show
                for line in usage: print(prefix + "   " + fmt.blue + line + fmt.reset)

            # Has subcommands: THIRD LEVEL
            elif isinstance(spec, dict):

                subsubcommands = spec
                #print(subsubcommands)
                self._show_help_subcommands(subsubcommands, main_command + " " + key, level=level+1)

            # No input expected
            else: print(prefix + "     " + fmt.blue + "no input" + fmt.reset)

    # -----------------------------------------------------------------

    def get_config_from_command(self, command, definition, name=None, index=1, interactive=False):

        """
        This function ...
        :param command:
        :param definition:
        :param name:
        :param index:
        :param interactive:
        :return:
        """

        from ..tools import strings
        from .configuration import prompt_settings, parse_arguments

        # Parse
        splitted = strings.split_except_within_double_quotes(command, add_quotes=False)
        if name is None: name = splitted[0]

        # Set parse command
        parse_command = splitted[index:]

        # Interactively get the settings
        if interactive: config = prompt_settings(name, definition, initialize=False, add_logging=False, add_cwd=False, add_config_path=False)

        # Parse arguments
        else: config = parse_arguments(name, definition, command=parse_command, error="exception", exit_on_help=False,
                                 initialize=False, add_logging=False, add_cwd=False)

        # Remove path
        config.pop("_path")

        # Return the configuration
        return config

    # -----------------------------------------------------------------

    @property
    def history_pts_user_path(self):

        """
        This function ...
        :return:
        """

        from ..tools import introspection
        return introspection.pts_user_history_dir

    # -----------------------------------------------------------------

    @abstractproperty
    def history_filename(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def history_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path
        return fs.join(self.history_pts_user_path, self.history_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def _history(self):

        """
        This function ...
        :return:
        """

        return fs.get_lines(self.history_path) if self.has_history else []

    # -----------------------------------------------------------------

    @property
    def has_history(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.history_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def show_history_definition(self):

        """
        This function ...
        :return:
        """

        from .configuration import ConfigurationDefinition
        definition = ConfigurationDefinition(write_config=False)
        definition.add_flag("all", "show all of the history (also from previous sessions)", False)
        definition.add_optional("path", "string", "write the history to a file")
        return definition

    # -----------------------------------------------------------------

    def show_history_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Parse command
        config = self.get_config_from_command(command, self.show_history_definition)

        # Show history
        self.show_history(all=config.all, path=config.path)

    # -----------------------------------------------------------------

    @property
    def last_command(self):

        """
        This function ...
        :return:
        """

        if not self.has_commands:
            if not self.has_history_commands: return None
            else: return self._history[-1]
        else: return self.commands[-1]

    # -----------------------------------------------------------------

    def show_history(self, all=False, path=None):

        """
        This function ...
        :param all:
        :param path:
        :return:
        """

        # Inform the user
        log.info("Showing history of commands ...")

        # Lines to write
        lines = []

        # Show all history?
        if all:
            print("")
            for command in self._history:
                print(command)
                lines.append(command)

        print("")
        for command in self.commands:
            print(command)
            lines.append(command)
        print("")

        # Write to file
        if path is not None: fs.write_lines(path, lines)

    # -----------------------------------------------------------------

    def write_history(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the history ...")

        # Write
        fs.add_lines(self.history_path, self.commands, create=True)

# -----------------------------------------------------------------
