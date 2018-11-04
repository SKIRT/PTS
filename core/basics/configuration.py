#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.configuration Contains the configuration class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
from types import MethodType
from abc import ABCMeta, abstractmethod
from types import NoneType
import sys
import argparse
from collections import OrderedDict
import StringIO
from functools import partial

# Import the relevant PTS classes and modules
from .map import Map
from ..tools import parsing, stringify
from ..tools import filesystem as fs
from ..basics.log import log
from .composite import SimplePropertyComposite
from ..tools import introspection
from ..tools import numbers, types
from ..tools.stringify import tostr
from ..tools import sequences

# -----------------------------------------------------------------

subtypes = dict()
subtypes["integer"] = ["positive_integer", "negative_integer"]
subtypes["real"] = ["fraction", "positive_real", "negative_real"]
subtypes["string"] = ["file_path", "directory_path"]
subtypes["quantity"] = ["photometric_quantity", "photometric_density_quantity", "photometric_brightness_quantity", "time_quantity", "length_quantity", "temperature_quantity", "mass_quantity", "mass_density_quantity", "mass_surface_density_quantity", "data_quantity"]
subtypes["unit"] = ["photometric_unit", "photometric_density_unit"]
subtypes["filter"] = ["broad_band_filter", "narrow_band_filter"]
subtypes["list"] = ["integer_list", "real_list", "quantity_list", "ascending_real_list", "descending_real_list", "ascending_integer_list", "descending_integer_list", "ascending_quantity_list", "descending_quantity_list"]

related_types = []
related_types.append(["integer", "positive_integer", "negative_integer", "even_integer", "even_positive_integer", "even_negative_integer", "odd_integer", "odd_positive_integer", "odd_negative_integer"])
related_types.append(["real", "fraction", "positive_real", "negative_real"])
related_types.append(["string", "file_path", "directory_path", "string_no_spaces"])
related_types.append(["quantity", "photometric_quantity", "photometric_density_quantity"])
related_types.append(["unit", "photometric_unit", "photometric_density_unit"])
related_types.append(["filter", "narrow_band_filter", "broad_band_filter"])
related_types.append(["broad_band_filter_list", "lazy_filter_list", "narrow_band_list"])

# -----------------------------------------------------------------

def parse_logging_arguments(name, description=None, add_cwd=True):

    """
    This function ...
    :param name:
    :return:
    """

    definition = ConfigurationDefinition(write_config=False)
    return parse_arguments(name, definition, description, add_logging=True, add_cwd=add_cwd)

# -----------------------------------------------------------------

def create_definition(**kwargs):

    """
    This function ...
    definition = create_definition(a, b, c, d=None, e="default", f=False) # CANNOT WORK
    # ONLY OPTIONAL AND FLAGS
    :param args:
    :param kwargs:
    :return:
    """

    definition = ConfigurationDefinition(write_config=False)
    for name in kwargs:

        default = kwargs[name]
        description = "value for '" + name + "'"

        if types.is_boolean_type(default): definition.add_flag(name, description, default)
        elif types.is_none(default): raise ValueError("Type cannot be defined when None is given as default")
        else:
            # Add optional argument
            ptype = stringify.get_parsing_type(default)
            definition.add_optional(name, ptype, description, default=default)

    # Return the definition
    return definition

# -----------------------------------------------------------------

def prompt_settings(name, definition, description=None, add_logging=True, add_cwd=True, initialize=True, add_config_path=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param initialize:
    :param add_config_path:
    :return:
    """

    # Create the configuration
    setter = InteractiveConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)
    config = setter.run(definition, prompt_optional=True)

    # Initialize PTS
    if initialize:
        from ...do.run import initialize_pts
        initialize_pts(config)

    # Return the configuration
    return config

# -----------------------------------------------------------------

def parse_arguments(name, definition, description=None, add_logging=True, add_cwd=True, command=None, error="exit", exit_on_help=True, initialize=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param command:
    :param error:
    :param exit_on_help:
    :param initialize:
    :return:
    """

    # Create the configuration
    setter = ArgumentConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd, error=error, exit_on_help=exit_on_help)
    config = setter.run(definition, command=command)

    # Initialize PTS
    if initialize:
        from ...do.run import initialize_pts
        initialize_pts(config)

    # Return the configuration
    return config

# -----------------------------------------------------------------

def get_usage(name, definition, description=None, add_logging=True, add_cwd=True, error="exit", exit_on_help=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param error:
    :param exit_on_help:
    :return:
    """

    # Create the configuration
    setter = ArgumentConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd,
                                         error=error, exit_on_help=exit_on_help)

    # Return the usage info lines
    return setter.get_usage(definition)

# -----------------------------------------------------------------

def print_usage(name, definition, description=None, add_logging=True, add_cwd=True, error="exit", exit_on_help=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param error:
    :param exit_on_help:
    :return:
    """

    # Create the configuration
    setter = ArgumentConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd,
                                         error=error, exit_on_help=exit_on_help)

    # Print the usage info
    return setter.print_usage(definition)

# -----------------------------------------------------------------

def get_help(name, definition, description=None, add_logging=True, add_cwd=True, error="exit", exit_on_help=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param error:
    :param exit_on_help:
    :return:
    """

    # Create the configuration
    setter = ArgumentConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd,
                                         error=error, exit_on_help=exit_on_help)

    # Return the help info lines
    return setter.get_help(definition)

# -----------------------------------------------------------------

def print_help(name, definition, description=None, add_logging=True, add_cwd=True, error="exit", exit_on_help=True):

    """
    This function ...
    :param name:
    :param definition:
    :param description:
    :param add_logging:
    :param add_cwd:
    :param error:
    :param exit_on_help:
    :return:
    """

    # Create the configuration
    setter = ArgumentConfigurationSetter(name, description=description, add_logging=add_logging, add_cwd=add_cwd,
                                         error=error, exit_on_help=exit_on_help)

    # Print the help info
    return setter.print_help(definition)

# -----------------------------------------------------------------

def create_configuration_passive(command_name, class_name, configuration_module_path, config_dict, description=None, cwd=None):

    """
    This function ...
    :param command_name:
    :param class_name:
    :param configuration_module_path:
    :param config_dict:
    :param description:
    :param cwd:
    :return:
    """

    # Find definition
    #try:
    #    configuration_module = importlib.import_module(configuration_module_path)
    #    definition = getattr(configuration_module, "definition")
    #except ImportError:
    #    log.warning("No configuration definition found for the " + class_name + " class")
    #    definition = ConfigurationDefinition()  # Create new configuration definition

    # Get the definition
    definition = get_definition(class_name, configuration_module_path, cwd=cwd)

    ## CREATE THE CONFIGURATION

    # Create the configuration setter
    if config_dict is None: config_dict = dict()  # no problem if all options are optional
    setter = DictConfigurationSetter(config_dict, command_name, description)

    # Create the configuration from the definition and from the provided configuration dictionary
    config = setter.run(definition)

    # Set the working directory
    if cwd is not None: config.path = cwd

    # Return the configuration
    return config

# -----------------------------------------------------------------

def prompt_yn(name, description, default=None):

    """
    This function ...
    :param name:
    :param description:
    :param default:
    :return:
    """

    # Create definition
    definition = ConfigurationDefinition(write_config=False)
    definition.add_flag(name, description, default=default)

    # Create setter
    setter = InteractiveConfigurationSetter("proceed", add_logging=False, add_cwd=False, add_config_path=False)

    # Get the answer
    while True:
        config = setter.run(definition, prompt_optional=True)
        if config[name] is None: log.warning("Answer with yes (y) or no (n)")
        else: return config[name]

# -----------------------------------------------------------------

def prompt_proceed(description=None, default=None):

    """
    This function ...
    :param description:
    :param default:
    :return:
    """

    if description is None: description = "proceed?"

    # Create definition
    definition = ConfigurationDefinition(write_config=False)
    definition.add_flag("proceed", description, default=default)

    # Create setter
    setter = InteractiveConfigurationSetter("proceed", add_logging=False, add_cwd=False, add_config_path=False)

    # Get the answer
    while True:
        config = setter.run(definition, prompt_optional=True)
        if config.proceed is None: log.warning("Answer with yes (y) or no (n)")
        else: return config.proceed

# -----------------------------------------------------------------

def prompt_finish(description=None):

    """
    This function ...
    :param description:
    :return:
    """

    if description is None: description = "press ENTER to finish"
    result = raw_input(description)

# -----------------------------------------------------------------

def prompt_automatic(name, description, default, choices=None, default_alias=None, suggestions=None):

    """
    This function ...
    :param name:
    :param description:
    :param default:
    :param choices:
    :param default_alias:
    :param suggestions:
    :return:
    """

    # None is not allowed
    if default is None: raise ValueError("Default value cannot be None")

    # Determine parsing type
    ptype, string = stringify.stringify(default)

    # Prompt
    return prompt_variable(name, ptype, description, choices=choices, default=default, default_alias=default_alias, suggestions=suggestions)

# -----------------------------------------------------------------

def prompt_choice(name, description, choices):

    """
    This function lets the user pick one option from a list of defined choices
    :param name:
    :param description:
    :param choices:
    :return:
    """

    # Determine parsing type
    list_ptype, string = stringify.stringify_list(choices)
    ptype = list_ptype.rsplit("_list", 1)[0]

    # Prompt
    return prompt_variable(name, ptype, description, choices=choices)

# -----------------------------------------------------------------

def prompt_choices(name, description, choices):

    """
    This function lets the user pick multiple options from a list of defined choices
    :param name:
    :param description:
    :param choices:
    :return:
    """

    # Determine parsing type
    list_ptype, string = stringify.stringify_list(choices)
    #ptype = list_ptype.rsplit("_list", 1)[0]

    # Prompt
    return prompt_variable(name, list_ptype, description, choices=choices)

# -----------------------------------------------------------------

def prompt_variable(name, parsing_type, description, choices=None, default=None, suggestions=None, required=True,
                    default_alias=None, convert_default=False, cue=None):

    """
    This function ....
    :param name: 
    :param parsing_type: 
    :param description: 
    :param choices:
    :param default:
    :param suggestions:
    :param required:
    :param default_alias:
    :param convert_default:
    :param cue:
    :return: 
    """

    # Create definition
    definition = ConfigurationDefinition(write_config=False)

    # Add setting
    if default is not None: definition.add_optional(name, parsing_type, description, choices=choices, suggestions=suggestions, default=default, default_alias=default_alias, convert_default=convert_default)
    elif required: definition.add_required(name, parsing_type, description, choices=choices, suggestions=suggestions)
    else: definition.add_optional(name, parsing_type, description, choices=choices, suggestions=suggestions)

    # Create setter
    setter = InteractiveConfigurationSetter(name, add_logging=False, add_cwd=False, add_config_path=False, cue=cue)

    # Get the answer
    config = setter.run(definition, prompt_optional=True)
    return config[name]

# -----------------------------------------------------------------

def prompt_fixed(name, description, value):

    """
    This function is just to show a fixed variable to the user
    :param name:
    :param description:
    :param value:
    :return:
    """

    # Create definition
    definition = ConfigurationDefinition(write_config=False)

    # Add setting
    definition.add_fixed(name, description, value)

    # Create setter
    setter = InteractiveConfigurationSetter(name, add_logging=False, add_cwd=False, add_config_path=False)

    # Get the answer
    config = setter.run(definition, prompt_optional=False)
    return config[name]

# -----------------------------------------------------------------

def prompt_string_list(name, description, choices=None, default=None, required=True, all_default=False, default_alias=None):

    """
    This function ...
    :param name: 
    :param description: 
    :param choices: 
    :param default: 
    :param required:
    :param all_default:
    :param default_alias:
    :return: 
    """

    # Check
    if default is not None and all_default: raise ValueError("Cannot specify 'default' with 'all_default' enabled")

    # All default
    if all_default:
        if choices is None: raise ValueError("Choices are not given")
        default = choices
        default_alias = "all"

    # Prompt
    return prompt_variable(name, "string_list", description, choices=choices, default=default, required=required, default_alias=default_alias)

# -----------------------------------------------------------------

def prompt_string(name, description, choices=None, default=None, required=True, cue=None):

    """
    This function ...
    :param name:
    :param description:
    :param choices:
    :param default:
    :param required:
    :param cue:
    :return:
    """

    return prompt_variable(name, "string", description, choices=choices, default=default, required=required, cue=cue)

# -----------------------------------------------------------------

def prompt_filepath(name, description):

    """
    Thisf unction ...
    :param name:
    :param description:
    :return:
    """

    return prompt_variable(name, "file_path", description, required=True)

# -----------------------------------------------------------------

def prompt_index(name, description, choices):

    """
    Thisf unction ...
    :param name:
    :param description:
    :param choices:
    :return:
    """

    return prompt_string(name, description, choices=choices, required=True)

# -----------------------------------------------------------------

def prompt_integer(name, description, choices=None, default=None, required=True):

    """
    Thisf unction ...
    :param name:
    :param description:
    :param choices:
    :param default:
    :param required:
    :return:
    """

    return prompt_variable(name, "integer", description, choices=choices, default=default, required=required)

# -----------------------------------------------------------------

def prompt_real(name, description, choices=None, default=None, required=True):

    """
    This function ...
    :param name:
    :param description:
    :param choices:
    :param default:
    :param required:
    :return:
    """

    return prompt_variable(name, "real", description, choices=choices, default=default, required=required)

# -----------------------------------------------------------------

def prompt_weights(name, description, choices=None, default=None, required=True):

    """
    This function ...
    :param name: 
    :param description: 
    :param choices: 
    :param default: 
    :param required: 
    :return: 
    """

    return prompt_variable(name, "weights", description, choices=choices, default=default, required=required)

# -----------------------------------------------------------------

def create_configuration_interactive(definition, command_name, description, **kwargs):

    """
    This function ...
    :param definition: 
    :param command_name: 
    :param description: 
    :return: 
    """

    return create_configuration(definition, command_name, description, "interactive", **kwargs)

# -----------------------------------------------------------------

def create_configuration(definition, command_name, description, configuration_method, **kwargs):

    """
    This function ...
    :param definition:
    :param command_name:
    :param description:
    :param configuration_method:
    :param kwargs:
    :return:
    """

    ## CREATE THE CONFIGURATION

    prompt_optional = kwargs.pop("prompt_optional", None)

    # Create the configuration setter
    if configuration_method == "interactive": setter = InteractiveConfigurationSetter(command_name, description, **kwargs)
    elif configuration_method == "arguments": setter = ArgumentConfigurationSetter(command_name, description, **kwargs)
    elif configuration_method.startswith("file"):
        configuration_filepath = configuration_method.split(":")[1]
        setter = FileConfigurationSetter(configuration_filepath, command_name, description, **kwargs)
    elif configuration_method == "last":
        configuration_filepath = fs.join(introspection.pts_user_config_dir, command_name + ".cfg")
        if not fs.is_directory(introspection.pts_user_config_dir): fs.create_directory(introspection.pts_user_config_dir)
        if not fs.is_file(configuration_filepath): raise RuntimeError("Cannot use rerun (config file not present)")
        setter = FileConfigurationSetter(configuration_filepath, command_name, description, **kwargs)
    else: raise ValueError("Invalid configuration method: " + configuration_method)

    # Create the configuration from the definition and from reading the command line arguments
    if configuration_method == "interactive": config = setter.run(definition, prompt_optional=prompt_optional)
    else: config = setter.run(definition)

    # Return the configuration
    return config

# -----------------------------------------------------------------

def create_configuration_flexible(name, definition, settings=None, default=False):

    """
    This function ...
    :param name:
    :param definition:
    :param settings:
    :param default:
    :return: 
    """

    ## A test settings dict is given
    if settings is not None:

        # Create the configuration
        setter = DictConfigurationSetter(settings, name, add_logging=False, add_cwd=False, add_config_path=False)
        config = setter.run(definition)

    # Settings are not given, default flag is added
    elif default:

        # Create the configuration
        setter = PassiveConfigurationSetter(name, add_cwd=False, add_logging=False, add_config_path=False)
        config = setter.run(definition)

    # No test configuration is given and default flag is not added
    else:

        # Create the configuration
        setter = InteractiveConfigurationSetter(name, add_cwd=False, add_logging=False, add_config_path=False)
        config = setter.run(definition, prompt_optional=True)

    # Return the configuration
    return config

# -----------------------------------------------------------------

def get_config_for_class(cls, config=None, interactive=False, cwd=None, prompt_optional=None, use_default=None,
                         check_required=True, extend_config=False):

    """
    This function ...
    :param cls:
    :param config:
    :param interactive:
    :param cwd:
    :param prompt_optional:
    :param use_default:
    :param check_required:
    :param extend_config: EXTEND CONFIGS TO ALL SETTINGS EVEN IF PASSED AS ACTUAL CONFIGURATION INSTANCE (don't assume config contains all)
    :return:
    """

    # If config is specified
    if config is not None:

        # Actual configuration object
        if isinstance(config, Configuration):

            # config, interactive=False, cwd=None, prompt_optional=None, use_default=None, check_required=True
            if extend_config: return extend_config_for_class(cls, config, interactive=interactive, cwd=cwd, prompt_optional=prompt_optional, use_default=use_default, check_required=check_required)
            else: return config

        # Other dictionary-like: always extend
        elif isinstance(config, dict): return extend_config_for_class(cls, config, interactive=interactive, cwd=cwd, prompt_optional=prompt_optional, use_default=use_default, check_required=check_required)

        # Not a valid config argument
        else: raise ValueError("Config should be Configuration, dictionary or None")

    # Look for the config
    else: return create_config_for_class(cls)

# -----------------------------------------------------------------

def create_config_for_class(cls, interactive=False, cwd=None, prompt_optional=None, use_default=None):

    """
    This function ...
    :param cls:
    :param interactive:
    :param cwd:
    :param prompt_optional:
    :param use_default:
    :return:
    """

    # Find the command
    command_name, class_name, configuration_module_path, description = find_command(cls)

    # If command is found
    if command_name is not None: return _create_config_from_module_path(class_name, configuration_module_path, interactive=interactive, cwd=cwd, prompt_optional=prompt_optional, use_default=use_default)

    # Command is not found
    else: return create_new_config(class_name, cwd=cwd, use_default=use_default)

# -----------------------------------------------------------------

def create_new_config(class_name, cwd=None, use_default=None):

    """
    This function ...
    :param cwd:
    :param use_default:
    :return:
    """

    # Cannot pass use_default
    if use_default is not None: raise ValueError("Cannot specifiy 'use_default': command is not recognized so configuration definition cannot be found")

    # Create an empty definition
    definition = ConfigurationDefinition(write_config=False)
    setter = InteractiveConfigurationSetter(class_name, add_logging=False)

    # Create new config
    config = setter.run(definition, prompt_optional=False)

    # Set the path
    if cwd is not None: config.path = cwd

    # Return the configuration
    return config

# -----------------------------------------------------------------

def _create_config_from_module_path(class_name, configuration_module_path, interactive=False, cwd=None, prompt_optional=None,
                                    use_default=None):

    """
    This function ...
    :param class_name:
    :param configuration_module_path:
    :param interactive:
    :param cwd:
    :param prompt_optional:
    :param use_default:
    :return:
    """

    # Get definition
    definition = get_definition(class_name, configuration_module_path, cwd=cwd)

    ## CREATE THE CONFIGURATION

    # Create configuration setter
    if interactive:

        setter = InteractiveConfigurationSetter(class_name, add_logging=False)
        config = setter.run(definition, prompt_optional=prompt_optional, use_default=use_default)

    # Not interactive
    else:

        setter = PassiveConfigurationSetter(class_name, add_logging=False)
        config = setter.run(definition)

    # Set the path
    if cwd is not None: config.path = cwd

    # Return the configuration
    return config

# -----------------------------------------------------------------

def extend_config_for_class(cls, config, interactive=False, cwd=None, prompt_optional=None, use_default=None,
                            check_required=True):

    """
    This function ...
    :param cls:
    :param config:
    :param interactive:
    :param cwd:
    :param prompt_optional:
    :param use_default:
    :param check_required:
    :return:
    """

    # Find the command
    command_name, class_name, configuration_module_path, description = find_command(cls)

    # Get configuration definition
    if command_name is not None: definition = get_definition(class_name, configuration_module_path, cwd=cwd)
    else: definition = ConfigurationDefinition(write_config=False)

    # Create the configuration interactively with pre-defined settings
    if interactive:

        # Create configuration with InteractiveConfigurationSetter
        setter = InteractiveConfigurationSetter(class_name, add_logging=False)
        config = setter.run(definition, settings=config, prompt_optional=prompt_optional, use_default=use_default)

    # Not interactive, use dict configuration setter
    else:

        # Create the DictConfigurationSetter
        setter = DictConfigurationSetter(config, command_name, description, check_required=check_required)
        config = setter.run(definition)

    # Set the path
    if cwd is not None: config.path = cwd

    # Return the configuration
    return config

# -----------------------------------------------------------------

def find_command(cls):

    """
    This function ...
    :return:
    """

    from ..tools import introspection

    tables = introspection.get_arguments_tables()
    # table_matches = introspection.find_matches_tables(script_name, tables)

    import inspect

    # Determine relative PTS path of the passed class
    class_name = cls.__name__
    class_path = inspect.getfile(cls).split(".py")[0]
    relative_class_path = class_path.rsplit("pts/")[1]
    relative_class_pts = relative_class_path.replace("/", ".") + "." + class_name

    # Determine subproject and relative class path
    subproject, relative_class_subproject = relative_class_pts.split(".", 1)

    # Get the correct table
    table = tables[subproject]

    command_name = None
    description = None
    configuration_name = None
    configuration_module_path = None

    # Loop over the entries in the table
    for i in range(len(table["Path"])):

        # print(table["Path"][i], relative_class_subproject)

        #print(table["Path"][i], relative_class_subproject)

        if table["Path"][i] == relative_class_subproject:

            command_name = table["Command"][i]
            description = table["Description"][i]

            if command_name.startswith("*"): command_name = command_name[1:]

            configuration_name = table["Configuration"][i]
            if configuration_name == "--": configuration_name = command_name
            configuration_module_path = "pts." + subproject + ".config." + configuration_name

            break

    # Return the command name
    return command_name, class_name, configuration_module_path, description

# -----------------------------------------------------------------

def get_definition(class_name, configuration_module_path, cwd=None):

    """
    This function ...
    :return:
    """

    import importlib

    ## GET THE CONFIGURATION DEFINITION
    try:
        #print(cwd)
        if cwd is not None: original_cwd = fs.change_cwd(cwd)
        else: original_cwd = None
        configuration_module = importlib.import_module(configuration_module_path)
        # has_configuration = True
        definition = getattr(configuration_module, "definition")
        if original_cwd is not None: fs.change_cwd(original_cwd)
    except ImportError:
        log.warning("No configuration definition found for the " + class_name + " class")
        # has_configuration = False
        definition = ConfigurationDefinition(write_config=False)  # Create new configuration definition

    # Return the configuration definition
    return definition

# -----------------------------------------------------------------

def combine_configs(*args):

    """
    This function ...
    :param args:
    :return:
    """

    # Initialize a new configuration
    config = Configuration()

    # Loop over the configurations
    for cfg in args:
        for label in cfg: config[label] = cfg[label]

    # Return the resulting configuration
    return config

# -----------------------------------------------------------------

def are_related_types(type_a, type_b):

    """
    This function ...
    :param type_a:
    :param type_b:
    :return:
    """

    for lst in related_types:
        if type_a in lst and type_b in lst: return True

    return False

# -----------------------------------------------------------------

def parent_type(type_name):

    """
    This function ...
    :param type_name:
    :return:
    """

    if type_name in subtypes.keys(): return type_name
    for label in subtypes:
        if type_name in subtypes[label]: return label
    return None

# -----------------------------------------------------------------

class Configuration(Map):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        """

        # Call the constructor of the base class
        super(Configuration, self).__init__(*args, **kwargs)

        # The path
        self._path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create the config
        config = cls()

        # Load the settings
        with open(path, 'r') as configfile: load_mapping(configfile, config)

        # Set the path
        config._path = path

        # Return the config
        return config

    # -----------------------------------------------------------------

    @classmethod
    def from_string(cls, string):

        """
        This function ...
        :param string:
        :return:
        """

        # Create the config
        config = cls()

        # Load the settings
        lines = string.split("\n")
        load_mapping(iter(lines), config)

        # Return the config
        return config

    # -----------------------------------------------------------------

    @classmethod
    def from_remote_file(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        # Create the config
        config = cls()

        # Load the settings
        lines_iterator = remote.read_lines(path)
        load_mapping(lines_iterator, config)

        # Set the path: NO
        #config._path = path

        # Return the config
        return config

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    def log_dir_path(self):

        """
        This function ...
        :return:
        """

        if "log_path" in self:
            if self["log_path"] is not None: return fs.absolute_or_in(self["log_path"], self["path"])
            else: return self.output_path()
        else: return None

    # -----------------------------------------------------------------

    def config_file_path(self, command_name=None):

        """
        This function ...
        :param command_name:
        :return:
        """

        # Config path specified
        if "config_path" in self and self["config_path"] is not None:

            # File or directory?
            name = fs.name(self["config_path"])

            # File, determine full path
            if "." in name: filepath = fs.absolute_or_in(self["config_path"], self["path"])

            # Directory, set file path from directory path
            else:

                dirpath = fs.absolute_or_in(self["config_path"], self["path"])

                # Create directory if necessary
                if not fs.is_directory(dirpath): fs.create_directory(dirpath)

                # Determine filepath
                filename = command_name + ".cfg" if command_name is not None else "config.cfg"
                filepath = fs.join(dirpath, filename)

            # Retrun the config file path
            return filepath

        # No config path specified
        else:

            # Determine filepath
            dirpath = self.output_path()

            # Determine the filename
            if command_name is not None: base_filename = command_name
            else: base_filename = "config"

            # Add timestamp?
            if "add_timestamp" in self and self["add_timestamp"]:
                from ..tools import time
                timestamp = time.filename_timestamp()
                filename = base_filename + "__" + timestamp + ".cfg"
            else: filename = base_filename + ".cfg"

            # Determine the full filepath
            filepath = fs.join(dirpath, filename)

            # Return the config file path
            return filepath

    # -----------------------------------------------------------------

    #@property
    def output_path(self):

        """
        This function ...
        :return:
        """

        # If 'output' is defined in the config
        if "output" in self and self["output"] is not None:

            full_output_path = fs.absolute_or_in(self["output"], self["path"])
            if not fs.is_directory(full_output_path): fs.create_directory(full_output_path)
            return full_output_path

        # Else, use the working directory as output directory
        else: return self["path"]

    # -----------------------------------------------------------------

    def to_string(self):

        """
        This function ...
        :return:
        """

        lines = []
        mapping_to_lines(lines, self)
        return "\n".join(lines)

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self._path is None: raise RuntimeError("Path is not defined")

        # Save
        self.saveto(self._path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write
        with open(path, 'w') as configfile: write_mapping(configfile, self)

        # Update the path
        self._path = path

# -----------------------------------------------------------------

def open_mapping(filepath):

    """
    This function ...
    :param filepath: 
    :return: 
    """

    parameters = Map()
    with open(filepath, "r") as fh: load_mapping(fh, parameters)
    return parameters

# -----------------------------------------------------------------

def open_box(filepath):

    """
    This function ...
    :param filepath:
    :return:
    """

    from box import Box
    parameters = Box(ordered_box=True)
    with open(filepath, "r") as fh: load_box(fh, parameters)
    return parameters

# -----------------------------------------------------------------

def load_mapping(fh, parameters, indent=""):
    return load_parameters(fh, parameters, Map, indent=indent)

# -----------------------------------------------------------------

def load_box(fh, parameters, indent=""):
    from box import Box
    return load_parameters(fh, parameters, Box, indent=indent, cls_kwargs={"ordered_dict": True})

# -----------------------------------------------------------------

def load_parameters(mappingfile, mapping, cls, indent="", cls_kwargs=None):

    """
    This function ...
    :param mappingfile
    :param mapping:
    :param cls:
    :param indent:
    :param cls_kwargs:
    :return:
    """

    if cls_kwargs is None: cls_kwargs = {}

    # STATES:
    # 0: empty line
    # 1: have description
    # 2: we got everything for a certain name
    # 3: section will begin (we expect a { now)
    # 4: inside section

    state = 0
    description = None

    # Loop over the lines in the file
    for line in mappingfile:

        # Strip end-of-line character
        line = line.rstrip("\n")

        # Remove the indent
        line = line.lstrip(indent)

        # Empty line
        if line == "":

            #print("empty", line)
            if state != 2 and state != 0: raise RuntimeError("At empty line and previous state was " + str(state))

            state = 0
            description = None
            continue

        # Description
        if line.startswith("#"):

            #print("comment", line)

            if description is not None:
                assert state == 1
                state = 1
                description += line.split("#")[1].strip()
            else:
                assert state == 0
                state = 1
                description = line.split("#")[1].strip()
            continue

        elif ":" in line:

            # Remove comment after the declaration
            if "#" in line: line = line.split("#")[0]

            try: before, after = line.split(":", 1)
            except ValueError: print("ERROR processing: ", line)

            #print("declaration", before, after)

            if after.strip() == "":

                # Check if there is a specification
                if before.endswith("]"):
                    name = before.split("[")[0].strip()
                    specification = before.split("[")[1].split("]")[0].strip()
                else:
                    name = before
                    specification = None

                if specification is None or specification == "section":
                    state = 3 # state 3: just before "{" is expected to start subdefinition
                    continue
                else:
                    if specification.endswith("list"): value = []
                    elif specification.endswith("tuple"): value = ()
                    elif specification.endswith("dictionary"): value = {}
                    else: raise RuntimeError("Encountered empty value for '" + name + "' of type '" + specification)

                    # Set the value
                    mapping[name] = value

            else:

                state = 2

                if before.endswith("]"):

                    name = before.split("[")[0].strip()
                    specification = before.split("[")[1].split("]")[0].strip()

                else:

                    name = before
                    specification = None

                # Should not happen
                if specification is None: raise ValueError("Invalid line (no specification): " + line)

                # Strip away leading or trailing spaces
                value = after.strip()

                #if specification is None:
                #    value = eval(value)

                if value == "None": value = None
                else:

                    #print(line)
                    parsing_function = getattr(parsing, specification)
                    #print(value)
                    value = parsing_function(value)

                #print("definition", value)

                # TODO: incorporate descriptions ??
                mapping[name] = value

                description = None

        # Start of section
        elif line.startswith("{"):

            assert state == 3
            state = 4

            if before.endswith("]"):

                name = before.split("[")[0].strip()
                specification = before.split("[")[1].split("]")[0].strip()

                assert specification == "section"

            else: name = before

            # Initialize the submapping
            mapping[name] = cls(**cls_kwargs)

            # Load the submapping
            load_mapping(mappingfile, mapping[name], indent=indent+"    ")

            state = 2
            description = None

        # End of section or complete definition
        elif line.startswith("}"): return

# -----------------------------------------------------------------

def save_mapping(path, mapping):

    """
    This function ...
    :param path:
    :param mapping:
    :return:
    """

    with open(path, 'w') as fh: write_mapping(fh, mapping)

# -----------------------------------------------------------------

def save_box(path, box):

    """
    This function ...
    :param path:
    :param box:
    :return:
    """

    with open(path, 'w') as fh: write_box(fh, box)

# -----------------------------------------------------------------

def write_mapping(mappingfile, mapping, indent=""):

    """
    This function ...
    :param mappingfile:
    :param mapping:
    :param indent:
    :return:
    """

    index = 0
    length = len(mapping)
    for name in mapping:

        # Skip internal stuff
        if name.startswith("_"): continue

        value = mapping[name]

        if isinstance(value, Map):

            print(indent + name + ":", file=mappingfile)
            print(indent + "{", file=mappingfile)
            write_mapping(mappingfile, value, indent=indent+"    ")
            print(indent + "}", file=mappingfile)

        else:
            ptype, string = stringify.stringify(mapping[name])
            #if string.startswith(" "): print(ptype, string)
            if ptype is None: ptype = "UNKNOWN"
            print(indent + name + " [" + ptype + "]: " + string, file=mappingfile)

        if index != length - 1: print("", file=mappingfile)
        index += 1

# -----------------------------------------------------------------

def write_box(boxfile, box, indent=""):

    """
    This function ...
    :param boxfile:
    :param box:
    :param indent:
    :return:
    """

    from box import Box

    index = 0
    length = len(box)
    for name in box:

        # Skip internal stuff
        if name.startswith("_"): continue

        value = box[name]

        if isinstance(value, Box):

            print(indent + name + ":", file=boxfile)
            print(indent + "{", file=boxfile)
            write_mapping(boxfile, value, indent=indent+"    ")
            print(indent + "}", file=boxfile)

        else:
            ptype, string = stringify.stringify(box[name])
            if ptype is None: ptype = "UNKNOWN"
            print(indent + name + " [" + ptype + "]: " + string, file=boxfile)

        if index != length - 1: print("", file=boxfile)
        index += 1

# -----------------------------------------------------------------

def mapping_to_lines(lines, mapping, indent=""):

    """
    This function ...
    :param lines:
    :param mapping:
    :param indent:
    :return:
    """

    index = 0
    length = len(mapping)
    for name in mapping:

        value = mapping[name]

        if isinstance(value, Map):

            lines.append(indent + name + ":")
            lines.append(indent + "{")
            mapping_to_lines(lines, value, indent=indent+"    ")
            lines.append(indent + "}")

        else:

            ptype, string = stringify.stringify(mapping[name])
            lines.append(indent + name + " [" + ptype + "]: " + string)

        if index != length - 1: lines.append("")
        index += 1

# -----------------------------------------------------------------

def print_mapping(mapping, indent="", empty_lines=True):

    """
    This function ...
    :param mapping:
    :param indent:
    :param empty_lines:
    :return:
    """

    # Create output string
    output = StringIO.StringIO()

    # Write mapping to string buffer
    print("")
    write_mapping(output, mapping, indent=indent)
    print("")

    # Show contents
    contents = output.getvalue()
    for line in contents.split("\n"):
        if not empty_lines and line.strip() == "": continue
        print(line)

    # Close the string buffer
    output.close()

# -----------------------------------------------------------------

def prompt_mapping(parameters, contains=None, not_contains=None, exact_name=None, exact_not_name=None, startswith=None,
                   endswith=None, label=None, descriptions=None, choices=None, fixed=None, suggestions=None,
                   required=True, add_suggestions=False):

    """
    This function ...
    :param parameters:
    :param contains:
    :param not_contains:
    :param exact_name:
    :param exact_not_name:
    :param startswith:
    :param endswith:
    :param label:
    :param descriptions:
    :param choices:
    :param fixed:
    :param suggestions:
    :param required:
    :param add_suggestions:
    :return:
    """

    has_changed = False
    used_suggestions = []

    # Adapt
    for name in parameters:

        # Skip properties related to configuration
        if name == "config_path": continue

        if suggestions is not None and name in suggestions: used_suggestions.append(name)

        # Checks
        if contains is not None and contains not in name: continue
        if not_contains is not None and not_contains in name: continue
        if exact_name is not None and name != exact_name: continue
        if exact_not_name is not None and name == exact_not_name: continue
        if startswith is not None and not name.startswith(startswith): continue
        if endswith is not None and not name.endswith(endswith): continue

        # Get properties
        description = descriptions[name] if descriptions is not None and name in descriptions else "no description"
        choics = choices[name] if choices is not None and name in choices else None
        suggestns = suggestions[name] if suggestions is not None and name in suggestions else None
        default = parameters[name]
        ptype, string = stringify.stringify(default)

        # No ptype: value was probably None: expect any kind of property
        if ptype is None or ptype == "None": ptype = "any"

        # Add label to description
        if label is not None: description = description + " [" + label + "]"

        # Fixed variable: show value and description
        if fixed is not None and name in fixed:
            value = prompt_fixed(name, description, default)
            continue

        # Ask for the new value
        value = prompt_variable(name, ptype, description, choices=choics, default=default, required=required, suggestions=suggestns)
        if default is None and value == "": continue

        # Set the property
        if value != default:

            # Debugging
            log.debug("Changing the value of '" + name + "' to '" + tostr(value) + "' ...")
            log.debug("Original value: '" + string + "'")

            # Set the new value
            parameters[name] = value

            # Set flag
            has_changed = True

    # Add suggested
    if suggestions is not None and add_suggestions:
        for name in suggestions:

            if name in used_suggestions: continue
            values = suggestions[name]
            if len(values) > 1: raise ValueError("Multiple suggestions")
            value = values[0]
            parameters[name] = value
            has_changed = True

            # Show the suggested value
            description = descriptions[name] if descriptions is not None and name in descriptions else "no description"
            value = prompt_fixed(name, description, value)

    # Return flag
    return has_changed

# -----------------------------------------------------------------

class ConfigurationDefinition(object):

    """
    This function ...
    """

    def __init__(self, prefix=None, log_path=None, config_path=None, write_config=True, add_timestamp=False):

        """
        This function ...
        :param prefix:
        :param log_path:
        :param config_path:
        :param write_config:
        :param add_timestamp:
        """

        # Prefix
        self.prefix = prefix

        # Log path and config path
        self.log_path = log_path
        self.config_path = config_path
        self.write_config = write_config
        self.add_timestamp = add_timestamp

        # Dictionary of sections
        self.sections = OrderedDict()
        self.section_descriptions = dict()

        # Dictionary of fixed parameters
        self.fixed = OrderedDict()

        # Dictionary of required parameters
        self.required = OrderedDict()

        # Dictionary of positional optional parameters
        self.pos_optional = OrderedDict()

        # Dictionary of optional parameters
        self.optional = OrderedDict()

        # Dictionary of flags
        self.flags = OrderedDict()

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, log_path=None, config_path=None, write_config=True, add_timestamp=False):

        """
        This function ...
        :param path:
        :param log_path:
        :param config_path:
        :param write_config:
        :param add_timestamp:
        :return:
        """

        # Create the definition
        definition = cls(log_path=log_path, config_path=config_path, write_config=write_config, add_timestamp=add_timestamp)

        # Load the definition
        with open(path, 'r') as configfile: load_definition(configfile, definition)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    @classmethod
    def from_defaults(cls, defaults, log_path=None, config_path=None, write_config=True, add_timestamp=False, positional=False, positional_flags=False):

        """
        This function ...
        :param defaults:
        :param log_path:
        :param config_path:
        :param write_config:
        :param add_timestamp:
        :param positional:
        :param positional_flags:
        :return:
        """

        # Create the definition
        definition = cls(log_path=log_path, config_path=config_path, write_config=write_config, add_timestamp=add_timestamp)

        # Add settings
        for name in defaults:

            # Get the value
            value = defaults[name]

            # Get the type
            ptype, string = stringify.stringify(value)

            # Add setting
            if ptype == "boolean":
                if positional and positional_flags: definition.add_positional_optional(name, ptype, "no description", default=value)
                else: definition.add_flag(name, "no description", default=value)
            else:
                if positional: definition.add_positional_optional(name, ptype, "no description", default=value)
                else: definition.add_optional(name, ptype, "no description", default=value)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def copy(self, sections=True, fixed=True, required=True, pos_optional=True, optional=True, flags=True):

        """
        This function ...
        :param sections:
        :param fixed:
        :param required:
        :param pos_optional:
        :param optional:
        :param flags:
        :return:
        """

        # Create new
        new = ConfigurationDefinition(prefix=self.prefix, log_path=self.log_path, config_path=self.config_path, write_config=self.write_config, add_timestamp=self.add_timestamp)

        # Add sections
        if sections:
            new.sections = copy.deepcopy(self.sections)
            new.section_descriptions = copy.deepcopy(self.section_descriptions)

        # Add fixed
        if fixed: new.fixed = copy.deepcopy(self.fixed)

        # Add required
        if required: new.required = copy.deepcopy(self.required)

        # Add pos optional
        if pos_optional: new.pos_optional = copy.deepcopy(self.pos_optional)

        # Add optional
        if optional: new.optional = copy.deepcopy(self.optional)

        # Add flags
        if flags: new.flags = copy.deepcopy(self.flags)

        # Return the copy
        return new

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.fixed) + len(self.required) + len(self.pos_optional) + len(self.optional) + len(self.flags)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Write the definition
        with open(path, 'w') as configfile: write_definition(self, configfile)

    # -----------------------------------------------------------------

    def set_arguments(self, parser, prefix=None):

        """
        This function ...
        :param parser:
        :param prefix:
        :return:
        """

        # Required
        index = 0
        for name in self.required:

            # Get properties
            real_type = self.required[name].type
            description = self.required[name].description
            choices = self.required[name].choices
            suggestions = self.required[name].suggestions
            min_value = self.required[name].min_value
            max_value = self.required[name].max_value
            forbidden = self.required[name].forbidden

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None: name = prefix + "/" + name

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list") or real_type.__name__.endswith("_tuple"): choices = None

            # Add argument to argument parser

            #print("real_type", real_type.__name__)
            #print(len(self.pos_optional))
            #print(index, len(self.required) - 1)

            # Construct type
            the_type = construct_type(real_type, min_value, max_value, forbidden)

            if suggestions is not None: description += " [suggestions: " + stringify.stringify(suggestions)[1] + "]"

            # If this is the last required argument, and it is of 'string' type, and there are no positional arguments,
            # we can allow the string to contain spaces without having to add quotation marks
            if real_type.__name__ == "string" and len(self.pos_optional) == 0 and index == len(self.required) - 1:
                parser.add_argument(name, nargs='+', help=description, choices=choices)
            else: parser.add_argument(name, type=the_type, help=description, choices=choices)

            index += 1

        # Positional optional
        for name in self.pos_optional:

            # Get properties
            real_type = self.pos_optional[name].type
            description = self.pos_optional[name].description
            default = self.pos_optional[name].default
            choices = self.pos_optional[name].choices
            suggestions = self.pos_optional[name].suggestions
            min_value = self.pos_optional[name].min_value
            max_value = self.pos_optional[name].max_value
            forbidden = self.pos_optional[name].forbidden

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None: name = prefix + "/" + name

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list") or real_type.__name__.endswith("_tuple"): choices = None

            # Construct type
            if min_value is not None or max_value is not None or forbidden is not None: the_type = construct_type(real_type, min_value, max_value, forbidden)
            else: the_type = real_type

            if suggestions is not None: description += " [suggestions: " + stringify.stringify(suggestions)[1] + "]"

            # Add argument to argument parser
            parser.add_argument(name, type=the_type, help=description, default=default, nargs='?', choices=choices)

        # Optional
        for name in self.optional:

            # Get properties
            real_type = self.optional[name].type
            description = self.optional[name].description
            default = self.optional[name].default
            choices = self.optional[name].choices
            letter = self.optional[name].letter
            suggestions = self.optional[name].suggestions
            min_value = self.optional[name].min_value
            max_value = self.optional[name].max_value
            forbidden = self.optional[name].forbidden

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None:
                name = prefix + "/" + name
                letter = None

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list") or real_type.__name__.endswith("_tuple"): choices = None

            # Construct type
            the_type = construct_type(real_type, min_value, max_value, forbidden)

            if suggestions is not None: description += " [suggestions: " + stringify.stringify(suggestions)[1] + "]"

            # Add the argument
            if letter is None: parser.add_argument("--" + name, type=the_type, help=description, default=default, choices=choices)
            else: parser.add_argument("-" + letter, "--" + name, type=the_type, help=description, default=default, choices=choices)

        # Flag
        for name in self.flags:

            # Get properties
            description = self.flags[name].description
            letter = self.flags[name].letter
            default = self.flags[name].default

            # Default is True
            if default is True:

                name = "not_" + name
                description = "don't " + description

                if prefix is not None:
                    name = prefix + "/" + name
                    letter = None

                # Add the argument
                if letter is None: parser.add_argument("--" + name, action="store_true", help=description)
                else: parser.add_argument("-" + letter, "--" + name, action="store_true", help=description)

            # Default is False
            elif default is False:

                if prefix is not None:
                    name = prefix + "/" + name
                    letter = None

                # Add the argument
                if letter is None: parser.add_argument("--" + name, action="store_true", help=description)
                else: parser.add_argument("-" + letter, "--" + name, action="store_true", help=description)

            # Default is None
            elif default is None:

                letter = None

                name1 = name
                description1 = description

                name2 = "not_" + name
                description2 = "don't " + description

                if prefix is not None:
                    name1 = prefix + "/" + name1
                    name2 = prefix + "/" + name2

                # Add 2 arguments to the parser
                parser.add_argument("--" + name1, action="store_true", help=description1)
                parser.add_argument("--" + name2, action="store_true", help=description2)

            # Invalid option
            else: raise ValueError("Invalid option for default of flag '" + name + "': " + tostr(default))

        # Add arguments of sections
        for section_name in self.sections:

            #if self.prefix is None: section_prefix = section_name
            #else: section_prefix

            if prefix is None: section_prefix = section_name
            else: section_prefix = prefix + "/" + section_name

            self.sections[section_name].set_arguments(parser, section_prefix)

    # -----------------------------------------------------------------

    def get_settings(self, settings, arguments, prefix=None):

        """
        This function ...
        :param settings:
        :param arguments:
        :param prefix:
        :return:
        """

        #print(arguments)

        # Add fixed
        for name in self.fixed: settings[name] = self.fixed[name].value

        # Add required
        for name in self.required:

            # Get properties
            real_type = self.required[name].type
            choices = self.required[name].choices

            if prefix is not None: argument_name = prefix + "/" + name
            else: argument_name = name

            # Get the value
            value = getattr(arguments, argument_name)

            # Check each value in the list, if applicable
            if real_type.__name__.endswith("_list") and choices is not None:
                # Check whether they are all in the choices
                for item in value:
                    if not item in choices: raise ValueError("Element '" + str(item) + "' not recognized. Options are: " + stringify.stringify(choices)[1])

            # Convert from list into string if string contains spaces and was the last required argument
            if real_type.__name__ == "string" and isinstance(value, list): value = " ".join(value)

            # Set the value
            settings[name] = value

        # Add positional optional
        for name in self.pos_optional:

            # Get properties
            real_type = self.pos_optional[name].type
            choices = self.pos_optional[name].choices

            if prefix is not None: argument_name = prefix + "/" + name
            else: argument_name = name

            # Get the value
            value = getattr(arguments, argument_name)

            # Check each value in the list, if applicable
            if real_type.__name__.endswith("_list") and choices is not None and value is not None:
                # Check whether they are all in the choices
                for item in value:
                    if not item in choices: raise ValueError("Element '" + str(item) + "' not recognized. Options are: " + stringify.stringify(choices)[1])

            # Set the value
            settings[name] = value

        # Add optional
        for name in self.optional:

            # Get properties
            real_type = self.optional[name].type
            choices = self.optional[name].choices

            if prefix is not None: argument_name = prefix + "/" + name
            else: argument_name = name

            # Get the value
            value = getattr(arguments, argument_name)

            # Check each value in the list, if applicable
            if real_type.__name__.endswith("_list") and choices is not None and value is not None:
                # Check whether they are all in the choices
                for item in value:
                    if not item in choices: raise ValueError("Element '" + str(item) + "' not recognized. Options are: " + stringify.stringify(choices)[1])

            # Set the value
            settings[name] = value

        # Add flags
        for name in self.flags:

            # Get properties
            default = self.flags[name].default

            # Default is True
            if default is True:

                command_name = "not_" + name

                if prefix is not None: argument_name = prefix + "/" + command_name
                else: argument_name = command_name

                # Get the value
                settings[name] = not getattr(arguments, argument_name)

            # Default is False
            elif default is False:

                command_name = name

                if prefix is not None: argument_name = prefix + "/" + command_name
                else: argument_name = command_name

                # Get the value
                settings[name] = getattr(arguments, argument_name)

            # Default is None
            elif default is None:

                command_name1 = name
                command_name2 = "not_" + name

                if prefix is not None:
                    argument_name1 = prefix + "/" + command_name1
                    argument_name2 = prefix + "/" + command_name2
                else:
                    argument_name1 = command_name1
                    argument_name2 = command_name2

                value1 = getattr(arguments, argument_name1)
                value2 = getattr(arguments, argument_name2)

                # Argument 1 is True
                if value1:

                    if value2: raise ValueError("Options '" + argument_name1 + "' and '" + argument_name2 + "' cannot both be enabled because they contradict each other")
                    settings[name] = True

                # Argument 1 is False
                else:

                    # Argument 2 is True (this is the negated one!)
                    if value2: settings[name] = False

                    # Argument 2 is False (so none of both options are specified: use the default of None)
                    else: settings[name] = None

            # Invalid option
            else: raise ValueError("Invalid value for default of flag '" + name + "': " + tostr(default))

        # Add the configuration settings of the various sections
        for name in self.sections:

            # Create a map for the settings
            settings[name] = Map()

            if prefix is None: section_prefix = name
            else: section_prefix = prefix + "/" + name

            # Recursively add the settings
            definition = self.sections[name]
            description = self.section_descriptions[name]
            definition.get_settings(settings[name], arguments, section_prefix)

    # -----------------------------------------------------------------

    def import_settings(self, definition, fixed=True, required=True, pos_optional=True, optional=True, flags=True,
                        sections=True, recursive=False, required_to="required", pos_optional_to="pos_optional",
                        optional_to="optional"):

        """
        This function ...
        :param definition:
        :param fixed:
        :param required:
        :param pos_optional:
        :param optional:
        :param flags:
        :param sections:
        :param recursive:
        :param required_to:
        :param pos_optional_to:
        :param optional_to:
        :return:
        """

        # Add fixed
        if fixed: self.import_fixed(definition)

        # Add required
        if required: self.import_required(definition, to=required_to)

        # Add positional optional
        if pos_optional: self.import_pos_optional(definition, to=pos_optional_to)

        # Add optional
        if optional: self.import_optional(definition, to=optional_to)

        # Add flags
        if flags: self.import_flags(definition)

        # Add sections
        if sections: self.import_sections(definition, recursive=recursive, fixed=fixed, required=required,
                                          pos_optional=pos_optional, optional=optional, flags=flags,
                                          required_to=required_to, pos_optional_to=pos_optional_to, optional_to=optional_to)

    # -----------------------------------------------------------------

    def import_properties(self, definition, names):

        """
        This function ...
        :param definition:
        :param names:
        :return:
        """

        # Loop over the property names
        for name in names:

            # Get the entry
            if definition.is_fixed(name): self.fixed[name] = definition.fixed[name].copy()
            elif definition.is_required(name): self.required[name] = definition.required[name].copy()
            elif definition.is_positional_optional(name): self.pos_optional[name] = definition.pos_optional[name].copy()
            elif definition.is_optional(name): self.optional[name] = definition.optional[name].copy()
            elif definition.is_flag(name): self.flags[name] = definition.flags[name].copy()
            else: raise ValueError("Property '" + name + "' does not exist in definition")

    # -----------------------------------------------------------------

    def import_fixed(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Loop over the settings
        for name in definition.fixed:

            if name in self.fixed: raise ValueError("Already fixed argument '" + name + "' in this definition")
            self.fixed[name] = definition.fixed[name].copy()

    # -----------------------------------------------------------------

    def import_required(self, definition, to="required"):

        """
        This function ...
        :param definition:
        :param to:
        :return:
        """

        # Loop over the settings
        for name in definition.required:

            # To required
            if to == "required":

                if name in self.required: raise ValueError("Already required argument '" + name + "' in this definition")
                self.required[name] = definition.required[name].copy()

            # To positional optional
            elif to == "pos_optional":

                if name in self.pos_optional: raise ValueError("Already positional optional argument '" + name + "' in this definition")
                self.pos_optional[name] = definition.required[name].copy()

            # To optional
            elif to == "optional":

                if name in self.optional: raise ValueError("Already optional argument '" + name + "' in this definition")
                self.optional[name] = definition.required[name].copy()

            # Invalid
            else: raise ValueError("Invalid")

    # -----------------------------------------------------------------

    def import_pos_optional(self, definition, to="pos_optional"):

        """
        This function ...
        :param definition:
        :param to:
        :return:
        """

        # Loop over the settings
        for name in definition.pos_optional:

            # To required
            if to == "required":

                if name in self.required: raise ValueError("Already required argument '" + name + "' in this definition")
                self.required[name] = definition.pos_optional[name].copy()

            # To positional optional
            elif to == "pos_optional":

                if name in self.pos_optional: raise ValueError("Already positional optional argument '" + name + "' in this definition")
                self.pos_optional[name] = definition.pos_optional[name].copy()

            # To optional
            elif to == "optional":

                if name in self.optional: raise ValueError("Already optional argument '" + name + "' in this definition")
                self.optional[name] = definition.pos_optional[name].copy()

            # Invalid
            else: raise ValueError("Invalid")

    # -----------------------------------------------------------------

    def import_optional(self, definition, to="optional"):

        """
        This function ...
        :param definition:
        :param to:
        :return:
        """

        # Loop over the settings
        for name in definition.optional:

            # To required
            if to == "required":

                if name in self.required: raise ValueError("Already required argument '" + name + "' in this definition")
                self.required[name] = definition.optional[name].copy()

            # Positional optional
            elif to == "pos_optional":

                if name in self.pos_optional: raise ValueError("Already positional optional argument '" + name + "' in this definition")
                self.pos_optional[name] = definition.optional[name].copy()

            # Optional
            elif to == "optional":

                if name in self.optional: raise ValueError("Already optional argument '" + name + "' in this definition")
                self.optional[name] = definition.optional[name].copy()

            # Invalid
            else: raise ValueError("Invalid")

    # -----------------------------------------------------------------

    def import_flags(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Loop over the settings
        for name in definition.flags:

            if name in self.flags: raise ValueError("Already flag '" + name + "' in this definition")
            self.flags[name] = definition.flags[name].copy()

    # -----------------------------------------------------------------

    def import_sections(self, definition, recursive=False, fixed=True, required=True, pos_optional=True, optional=True,
                        flags=True, required_to="required", pos_optional_to="pos_optional", optional_to="optional"):

        """
        This function ...
        :param definition: 
        :param recursive:
        :param fixed:
        :param required:
        :param pos_optional:
        :param optional:
        :param flags:
        :param required_to:
        :param pos_optional_to:
        :param optional_to:
        :return: 
        """

        # Loop over the to be imported sections
        for name in definition.sections:

            # Existing section
            if name in self.sections:

                if recursive:
                    self.sections[name].import_settings(definition.sections[name], fixed=fixed, required=required,
                                                        pos_optional=pos_optional, optional=optional, flags=flags,
                                                        sections=True, recursive=True, required_to=required_to,
                                                        pos_optional_to=pos_optional_to, optional_to=optional_to)
                else: raise ValueError("Already a section '" + name + "' in this definition")

            # New section
            else:
                self.sections[name] = definition.sections[name].copy()
                self.section_descriptions[name] = definition.section_descriptions[name]

    # -----------------------------------------------------------------

    def add_section(self, name, description):

        """
        This function ...
        :param name:
        :param description:
        :return:
        """

        # Determine prefix
        #if self.prefix is None: prefix = name
        #else: prefix = self.prefix + "/" + name
        prefix = name

        # Add the section
        self.sections[name] = ConfigurationDefinition(prefix=prefix)
        self.section_descriptions[name] = description

    # -----------------------------------------------------------------

    @property
    def section_names(self):

        """
        This function ...
        :return:
        """

        return self.sections.keys()

    # -----------------------------------------------------------------

    def is_section(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.section_names

    # -----------------------------------------------------------------

    def remove_section(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.section_names: raise ValueError("Not a section: '" + name + "'")
        del self.sections[name]
        del self.section_descriptions[name]

    # -----------------------------------------------------------------

    def remove_all_sections(self):

        """
        This function ...
        :return:
        """

        for name in self.section_names: self.remove_section(name)

    # -----------------------------------------------------------------

    def import_section(self, name, description, definition, sections=True, fixed=True, required=True, pos_optional=True,
                       optional=True, flags=True):

        """
        This function ...
        :param name:
        :param description:
        :param definition:
        :param sections:
        :param fixed:
        :param required:
        :param pos_optional:
        :param optional:
        :param flags:
        :return:
        """

        # Determine prefix
        prefix = name

        # Add the section
        self.sections[name] = definition.copy(sections=sections, fixed=fixed, required=required,
                                              pos_optional=pos_optional, optional=optional, flags=flags)
        self.sections[name].prefix = prefix
        self.section_descriptions[name] = description

    # -----------------------------------------------------------------

    def import_section_from_composite_class(self, name, description, cls):

        """
        This funtion ...
        """

        # First create 'default' instance
        instance = cls()
        self.import_section_from_properties(name, description, instance)

    # -----------------------------------------------------------------

    def import_section_from_properties(self, name, description, instance):

        """
        This function creates a configuration section from a SimplePropertyComposite-derived class
        :param name:
        :param description:
        :param instance:
        :return:
        """

        default_values = vars(instance)

        # Add empty section
        self.add_section(name, description)

        # Loop over the properties
        for pname in default_values:

            if pname.startswith("_"): continue

            default_value = default_values[pname]

            # Section
            if isinstance(default_value, SimplePropertyComposite):

                # Recursive call to this function
                pdescription = instance._descriptions[pname]
                self.sections[name].import_section_from_properties(pname, pdescription, default_value)

            else:

                # Get the type, description and the choices
                ptype = instance._ptypes[pname]
                pdescription = instance._descriptions[pname]
                choices = instance._choices[pname]

                # Add to the definition
                if ptype == "boolean": self.sections[name].add_flag(pname, pdescription, default_value)
                else: self.sections[name].add_optional(pname, ptype, pdescription, default_value, choices=choices)

    # -----------------------------------------------------------------

    @property
    def property_names(self):

        """
        This function ...
        :return:
        """

        return self.fixed_names + self.required_names + self.positional_optional_names + self.optional_names + self.flag_names

    # -----------------------------------------------------------------

    def remove_setting(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name in self.fixed_names: self.remove_fixed(name)
        elif name in self.required_names: self.remove_required(name)
        elif name in self.positional_optional_names: self.remove_positional_optional(name)
        elif name in self.optional_names: self.remove_optional(name)
        elif name in self.flag_names: self.remove_flag(name)
        else: raise ValueError("Not a setting: '" + name + "'")

    # -----------------------------------------------------------------

    def add_fixed(self, name, description, value):

        """
        This function ...
        :param name:
        :param description:
        :param value:
        :return:
        """

        self.fixed[name] = Map(description=description, value=value)

    # -----------------------------------------------------------------

    @property
    def fixed_names(self):

        """
        This function ...
        :return:
        """

        return self.fixed.keys()

    # -----------------------------------------------------------------

    def is_fixed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.fixed_names

    # -----------------------------------------------------------------

    def remove_fixed(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.fixed_names: raise ValueError("Not a fixed setting: '" + name + "'")
        del self.fixed[name]

    # -----------------------------------------------------------------

    def remove_all_fixed(self):

        """
        This function ...
        :return:
        """

        for name in self.fixed_names: self.remove_fixed(name)

    # -----------------------------------------------------------------

    def add_required(self, name, user_type, description, choices=None, dynamic_list=False, suggestions=None,
                     min_value=None, max_value=None, forbidden=None):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param choices:
        :param dynamic_list:
        :param suggestions:
        :param min_value:
        :param max_value:
        :param forbidden:
        :return:
        """

        # Check
        if choices is not None and suggestions is not None: raise ValueError("Cannot specify both choices and suggestions at the same time")

        # Get the real type
        real_type = get_real_type(user_type)

        # Add
        self.required[name] = Map(type=real_type, description=description, choices=choices, dynamic_list=dynamic_list,
                                  suggestions=suggestions, min_value=min_value, max_value=max_value, forbidden=forbidden)

    # -----------------------------------------------------------------

    @property
    def required_names(self):

        """
        This function ...
        :return:
        """

        return self.required.keys()

    # -----------------------------------------------------------------

    def is_required(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.required_names

    # -----------------------------------------------------------------

    def remove_required(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.required_names: raise ValueError("Not a required setting: '" + name + "'")
        del self.required[name]

    # -----------------------------------------------------------------

    def remove_all_required(self):

        """
        This function ...
        :return:
        """

        for name in self.required_names: self.remove_required(name)

    # -----------------------------------------------------------------

    def add_positional_optional(self, name, user_type, description, default=None, choices=None,
                                convert_default=False, dynamic_list=False, suggestions=None, min_value=None,
                                max_value=None, forbidden=None, default_alias=None):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param choices:
        :param convert_default:
        :param dynamic_list:
        :param suggestions:
        :param min_value:
        :param max_value:
        :param forbidden:
        :param default_alias:
        :return:
        """

        # Check
        if choices is not None and suggestions is not None: raise ValueError("Cannot specify both choices and suggestions at the same time")

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if default is not None:

            # Convert or check default value
            if convert_default: default = parse_default(default, user_type, real_type)
            else: default = check_default(default, user_type)

        # Check default
        if default is not None and choices is not None:

            # List-type default value
            if types.is_sequence(default) or types.is_tuple(default):
                if not sequences.is_subset(default, choices): raise ValueError("The default value '" + tostr(default, delimiter=", ") + "' does not contain a subset of the choices (" + tostr(choices, delimiter=", ") + ")")

            # Regular default value
            elif default not in choices: raise ValueError("The default value '" + tostr(default) + "' is not one of the choices (" + tostr(choices, delimiter=", ") + ")")

        # Add
        self.pos_optional[name] = Map(type=real_type, description=description, default=default, choices=choices,
                                      dynamic_list=dynamic_list, suggestions=suggestions, min_value=min_value,
                                      max_value=max_value, forbidden=forbidden, default_alias=default_alias)

    # -----------------------------------------------------------------

    @property
    def positional_optional_names(self):

        """
        This function ...
        :return:
        """

        return self.pos_optional.keys()

    # -----------------------------------------------------------------

    def is_positional_optional(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.positional_optional_names

    # -----------------------------------------------------------------

    def remove_positional_optional(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.positional_optional_names: raise ValueError("Not a positional optional setting: '" + name + "'")
        del self.pos_optional[name]

    # -----------------------------------------------------------------

    def remove_all_positional_optional(self):

        """
        This function ...
        :return:
        """

        for name in self.positional_optional_names: self.remove_positional_optional(name)

    # -----------------------------------------------------------------

    def add_optional(self, name, user_type, description, default=None, choices=None, letter=None, convert_default=False,
                     dynamic_list=False, suggestions=None, min_value=None, max_value=None, forbidden=None,
                     default_alias=None, convert_choices=False, convert_suggestions=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param choices:
        :param letter:
        :param convert_default:
        :param dynamic_list:
        :param suggestions:
        :param min_value:
        :param max_value:
        :param forbidden:
        :param default_alias:
        :param convert_choices:
        :param convert_suggestions:
        :return:
        """

        # Check
        if choices is not None and suggestions is not None: raise ValueError("Cannot specify both choices and suggestions at the same time")

        if self.prefix is not None and letter is not None: raise ValueError("Cannot assign letter argument for child configuration definition")

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if default is not None:

            # Convert or check default value
            if convert_default: default = parse_default(default, user_type, real_type)
            else: default = check_default(default, user_type)

        # Convert choices
        if choices is not None and convert_choices:
            choices = [parse_default(choice, user_type, real_type) for choice in choices]

        # Convert suggestions
        if suggestions is not None and convert_suggestions:
            suggestions = [parse_default(suggestion, user_type, real_type) for suggestion in suggestions]

        # Check default
        if default is not None and choices is not None:

            # List-type default value
            if types.is_sequence(default) or types.is_tuple(default):
                if not sequences.is_subset(default, choices): raise ValueError("The default value '" + tostr(default, delimiter=", ") + "' does not contain a subset of the choices (" + tostr(choices, delimiter=", ") + ")")

            # Regular default value
            elif default not in choices: raise ValueError("The default value '" + tostr(default) + "' is not one of the choices (" + tostr(choices, delimiter=", ") + ")")

        # Add
        self.optional[name] = Map(type=real_type, description=description, default=default, choices=choices,
                                  letter=letter, dynamic_list=dynamic_list, suggestions=suggestions,
                                  min_value=min_value, max_value=max_value, forbidden=forbidden, default_alias=default_alias)

    # -----------------------------------------------------------------

    @property
    def optional_names(self):

        """
        This function ...
        :return:
        """

        return self.optional.keys()

    # -----------------------------------------------------------------

    def is_optional(self, name):

        """
        This function ...
        :param name
        :return:
        """

        return name in self.optional_names

    # -----------------------------------------------------------------

    def remove_optional(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.optional_names: raise ValueError("Not an optional setting: '" + name + "'")
        del self.optional[name]

    # -----------------------------------------------------------------

    def remove_all_optional(self):

        """
        This function ...
        :return:
        """

        for name in self.optional_names: self.remove_optional(name)

    # -----------------------------------------------------------------

    def add_flag(self, name, description, default=False, letter=None, convert_default=False):

        """
        This function ...
        :param name:
        :param description:
        :param default:
        :param letter:
        :param convert_default:
        :return:
        """

        # Convert default
        if convert_default: default = get_real_value(default, parsing.boolean)

        # Add
        self.flags[name] = Map(description=description, letter=letter, default=default)

    # -----------------------------------------------------------------

    @property
    def flag_names(self):

        """
        This function ...
        :return:
        """

        return self.flags.keys()

    # -----------------------------------------------------------------

    def is_flag(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.flag_names

    # -----------------------------------------------------------------

    def remove_flag(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name not in self.flag_names: raise ValueError("Not a flag setting: '" + name + "'")
        del self.flags[name]

    # -----------------------------------------------------------------

    def remove_all_flags(self):

        """
        This function ...
        :return:
        """

        for name in self.flag_names: self.remove_flag(name)

# -----------------------------------------------------------------

class ConfigurationSetter(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, name, description=None, add_logging=True, add_cwd=True, add_config_path=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        """

        # Set the name and description
        self.name = name
        self.description = description

        # The configuration definition
        self.definition = None

        # Set options
        self.add_logging = add_logging
        self.add_cwd = add_cwd
        self.add_config_path = add_config_path

        # The configuration
        self.config = Configuration()

    # -----------------------------------------------------------------

    @abstractmethod
    def run(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def set_logging_and_cwd(self, definition=None):

        """
        This function ...
        :return:
        """

        if definition is None: definition = self.definition

        cwd_path = fs.cwd()

        # Add logging options
        if self.add_logging:

            # Log path to absolute path
            log_path = fs.absolute_or_in(definition.log_path, cwd_path) if definition.log_path is not None else cwd_path # set absolute log path

            definition.add_fixed("log_path", "the directory for the log file be written to", log_path)
            definition.add_flag("debug", "enable debug output", letter="d")
            definition.add_flag("brief", "brief output", letter="b")
            definition.add_flag("memuse", "log memory usage", letter="m")

            # Report?
            if definition.log_path is not None: definition.add_fixed("report", "write a report file", True) # if log path is defined in definition, always report
            else: definition.add_flag("report", "write a report file")  # otherwise, ask

        # Add config path
        if self.add_config_path and definition.write_config:

            # Set the path to the directory where the configuration file should be saved
            if definition.config_path is not None: definition.add_fixed("config_path", "directory for the configuration file to be written to", definition.config_path)
            else: definition.add_optional("config_path", "directory_path", "directory for the configuration file to be written to (relative to the working directory or absolute) (if None, the output directory is used)")

            # Write config?
            if definition.config_path is not None: definition.add_fixed("write_config", "write the configuration", True) # if config path is defined in the definition, always write
            else: definition.add_flag("write_config", "write the configuration") # otherwise, ask

        # Add timestamp?
        if definition.write_config and definition.add_timestamp: definition.add_fixed("add_timestamp", "add timestamp to configuration file", True)

        # Add the path to the current working directory
        if self.add_cwd: definition.add_fixed("path", "the working directory", cwd_path)

# -----------------------------------------------------------------

class InteractiveConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True, add_config_path=True, cue=None):

        """
        The constructor ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        :param cue:
        """

        # Call the constructor of the base class
        super(InteractiveConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

        # Set cue
        if cue is None: cue = "   : "
        self.cue = cue

    # -----------------------------------------------------------------

    def run(self, definition, prompt_optional=None, settings=None, use_default=None):

        """
        This function ...
        :param definition:
        :param prompt_optional:
        :param settings:
        :param use_default:
        :return:
        """

        # Set the definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Do interactive
        self.interactive(prompt_optional, settings=settings, use_default=use_default)

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def interactive(self, prompt_optional, settings=None, use_default=None):

        """
        This function ...
        :param prompt_optional:
        :param settings:
        :param use_default:
        :return:
        """

        if prompt_optional is None:

            # Ask whether optional parameters have to be shown, or to just use the default values
            log.info("Do you want to configure optional settings (y or n)?")
            log.info("Press ENTER to use the default (True)")

            while True:
                answer = raw_input(self.cue)
                if answer == "":
                    prompt_optional = True
                    break
                else:
                    try:
                        prompt_optional = parsing.boolean(answer)
                        break
                    except ValueError: log.warning("Invalid input. Try again.")

        # Get the settings from an interactive prompt
        add_settings_interactive(self.config, self.definition, prompt_optional=prompt_optional, settings=settings,
                                 use_default=use_default, cue=self.cue)

# -----------------------------------------------------------------

def error_to_exception(self, message):

    """
    This function ...
    :param message:
    :return:
    """

    raise ValueError(message)

# -----------------------------------------------------------------

class NoExitHelpAction(argparse.Action):

    """
    This class ...
    """

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):

        """
        The constructor ...
        :param option_strings:
        :param dest:
        :param default:
        :param help:
        """

        super(NoExitHelpAction, self).__init__( option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)

    # -----------------------------------------------------------------

    def __call__(self, parser, namespace, values, option_string=None):

        """
        This function ...
        :param parser:
        :param namespace:
        :param values:
        :param option_string:
        :return:
        """

        parser.print_help()

# -----------------------------------------------------------------

class ArgumentConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True, add_config_path=True, error="exit",
                 exit_on_help=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        :param error:
        :param exit_on_help:
        """

        # Call the constructor of the base class
        super(ArgumentConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

        # Create the command-line parser
        self.parser = argparse.ArgumentParser(prog=name, description=description)

        # Set error handling
        if error == "exit": pass
        elif error == "exception": self.parser.error = MethodType(error_to_exception, self.parser) # change the function
        else: raise ValueError("Invalid option for error handling")

        # Set help handling
        if not exit_on_help: self.parser.register('action', 'help', NoExitHelpAction) # DOESN'T WORK??

        # The parsed arguments
        self.arguments = None

        # The command as a string
        self.command = None

    # -----------------------------------------------------------------

    @staticmethod
    def get_arguments():

        """
        This function ...
        :return:
        """

        return sys.argv[1:]

    # -----------------------------------------------------------------

    def get_help(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set logging and cwd
        self.set_logging_and_cwd(definition)

        # Set arguments
        definition.set_arguments(self.parser)

        # Return the help lines
        lines = self.parser.format_help().split("\n")
        return [line for line in lines if line] # ignore empty lines

    # -----------------------------------------------------------------

    def get_usage(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set logging and cwd
        self.set_logging_and_cwd(definition)

        # Set arguments
        definition.set_arguments(self.parser)

        # Return the usage lines
        lines = self.parser.format_usage().split("\n")
        return [line for line in lines if line]  # ignore empty lines

    # -----------------------------------------------------------------

    def print_help(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set logging and cwd
        self.set_logging_and_cwd(definition)

        # Set arguments
        definition.set_arguments(self.parser)

        # Print the help
        return self.parser.print_help()

    # -----------------------------------------------------------------

    def print_usage(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set logging and cwd
        self.set_logging_and_cwd(definition)

        # Set arguments
        definition.set_arguments(self.parser)

        # Print the help
        return self.parser.print_usage()

    # -----------------------------------------------------------------

    def run(self, definition, command=None):

        """
        This function ...
        :param definition:
        :param command:
        :return:
        """

        # Set definition
        self.definition = definition
        self.command = command

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Parse
        self.parse()

        # Create config
        self.create_config()

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def parse(self):

        """
        This function ...
        :return:
        """

        # Set arguments
        self.definition.set_arguments(self.parser)

        # Let the parser parse
        self.arguments = self.parser.parse_args(args=self.command)

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Add the settings
        self.definition.get_settings(self.config, self.arguments)

# -----------------------------------------------------------------

class FileConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, path, name, description=None, add_logging=True, add_cwd=True, add_config_path=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        """

        # Call the constructor of the base class
        super(FileConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

        # Set the path to the specified configuration file
        self.path = path

        # The user-provided configuration
        self.user_config = None

    # -----------------------------------------------------------------

    def run(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Load the user provided configuration
        self.load_config()

        # Create the config (check it against the definition)
        self.create_config()

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def load_config(self):

        """
        This function ...
        :return:
        """

        self.user_config = Configuration.from_file(self.path)

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # TODO: TEMPORARY

        self.config = self.user_config

# -----------------------------------------------------------------

class DictConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, dictionary, name, description=None, add_logging=True, add_cwd=True, add_config_path=True,
                 check_required=True):

        """
        The constructor ...
        :param dictionary:
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        :param check_required:
        """

        # Call the constructor of the base class
        super(DictConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

        # Set the user-provided dictionary
        self.dictionary = dictionary

        # Flag
        self.check_required = check_required

    # -----------------------------------------------------------------

    def run(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Create the configuration
        self.create_config()

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Add the settings to the configuration
        #config, definition, dictionary
        add_settings_from_dict(self.config, self.definition, self.dictionary, check_required=self.check_required)

# -----------------------------------------------------------------

class PassiveConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True, add_config_path=True):

        """
        The constructor ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        """

        # Call the constructor of the base class
        super(PassiveConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

    # -----------------------------------------------------------------

    def run(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Create the configuration
        self.create_config()

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def create_config(self):

        """
        This function ...
        :return:
        """

        # Add the settings to the configuration
        add_settings_default(self.config, self.definition)

# -----------------------------------------------------------------

class GraphicalConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, path, name, description=None, add_logging=True, add_cwd=True, add_config_path=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        :param add_config_path:
        """

        # Call the constructor of the base class
        super(GraphicalConfigurationSetter, self).__init__(name, description, add_logging=add_logging, add_cwd=add_cwd, add_config_path=add_config_path)

        # ...

    # -----------------------------------------------------------------

    def run(self, definition):

        """
        This function ...
        :param definition:
        :return:
        """

        # Set definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

# -----------------------------------------------------------------

def get_real_type(user_type):

    """
    This function ...
    :param user_type:
    :return:
    """

    if types.is_string_type(user_type): return getattr(parsing, user_type)
    else: return user_type

# -----------------------------------------------------------------

def get_real_value(default, user_type):

    """
    This function ...
    :param default:
    :param user_type:
    :return:
    """

    return user_type(default)

# -----------------------------------------------------------------

def load_definition(configfile, definition):

    """
    This function ...
    :param configfile:
    :param definition:
    :return:
    """

    # STATES:
    # 0: empty line
    # 1: have description
    # 2: we got everything for a certain name
    # 3: section will begin (we expect a { now)

    state = 0
    description = None

    index = 0

    # Loop over the lines in the file
    for line in configfile:

        # Strip end-of-line character
        line = line.rstrip("\n").lstrip()

        # Empty line
        if line == "":
            if state != 2:
                raise RuntimeError("State is not 2 (but " + str(state) + " and empty line encountered")  # 2 means we got everything for a certain name
            state = 0
            description = None
            continue

        # Description
        if line.startswith("#"):

            if description is not None:
                assert state == 1
                state = 1
                description += line.split("#")[1].strip()
            else:
                assert state == 0
                state = 1
                description = line.split("#")[1].strip()
            continue

        elif ":" in line:

            before, after = line.split(":")

            #after_after = after.split(" #")[0] if len(after.split(" #")) > 1 else ""
            after_without_comment = after.split(" #")[0]
            #print("after after", after_after)
            if after_without_comment.strip() == "":

                state = 3  # state 3: just before "{" is expected to start subdefinition
                continue

            else:

                state = 2
                name = before.split("[")[0].strip()
                specification = before.split("[")[1].split("]")[0].strip().split(", ")

                choices_string = None
                if " #" in after:
                    value = after.split(" #")[0]
                    comment_string = after.split("#")[1].strip()
                    if "choices = " in comment_string:
                        choices_string = comment_string.split("choices = ")[1]
                else: value = after

                # Strip value for leading and trailing spaces
                value = value.strip()

                kind = specification[0]
                print(specification)
                user_type = specification[1].strip() if kind != "flag" else None

                # Set user type
                if user_type == "None":
                    user_type = NoneType
                    value = None
                elif kind != "flag": user_type = getattr(parsing, user_type)

                # Convert choices to real list of the user_type-'list' type
                choices = None
                if choices_string is not None:

                    user_type_list_name = user_type.__name__ + "_list"
                    user_type_list = getattr(parsing, user_type_list_name)
                    choices = user_type_list(choices_string)

                if kind == "fixed":

                    value = user_type(value) if value is not None else None # convert here because the add_fixed function doesn't bother with types now
                    definition.add_fixed(name, description, value)

                elif kind == "required": definition.add_required(name, user_type, description, choices=choices)
                elif kind == "pos_optional": definition.add_positional_optional(name, user_type, description, default=value, choices=choices, convert_default=True)
                elif kind == "optional": definition.add_optional(name, user_type, description, default=value, choices=choices, letter=None, convert_default=True)
                elif kind == "flag":

                    value = parsing.boolean(value) if value is not None else False
                    definition.add_flag(name, description, default=value, letter=None)

                else: raise ValueError("Invalid kind of argument: " + kind)

                description = None

        # Start of
        elif line.startswith("{"):

            assert state == 3
            state = 4

            name = before.split("[")[0].strip()
            specification = before.split("[")[1].split("]")[0].strip()

            assert specification == "section"

            # Initialize the configuration definition
            definition.add_section(name, description)

            # Load into the section
            load_definition(configfile, definition.sections[name])

            state = 2
            description = None

        # End of section or complete definition
        elif line.startswith("}"): return

# -----------------------------------------------------------------

def write_definition(definition, configfile, indent=""):

    """
    This function ...
    :param definition:
    :param configfile:
    :param indent:
    :return:
    """

    # Fixed
    for name in definition.fixed:

        # Get properties
        description = definition.fixed[name].description
        value = definition.fixed[name].value

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [fixed]: " + str(value), file=configfile)

    # Required
    for name in definition.required:

        # Get properties
        real_type = definition.required[name].type
        description = definition.required[name].description
        choices = definition.required[name].choices
        dynamic_list = definition.required[name].dynamic_list
        suggestions = definition.required[name].suggestions
        min_value = definition.required[name].min_value
        max_value = definition.required[name].max_value
        forbidden = definition.required[name].forbidden

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [required, " + real_type.__name__ + "]: None" + choices_string, file=configfile)

    # Positional optional
    for name in definition.pos_optional:

        # Get properties
        real_type = definition.pos_optional[name].type
        description = definition.pos_optional[name].description
        default = definition.pos_optional[name].default
        choices = definition.pos_optional[name].choices
        dynamic_list = definition.pos_optional[name].dynamic_list
        suggestions = definition.pos_optional[name].suggestions
        min_value = definition.pos_optional[name].min_value
        max_value = definition.pos_optional[name].max_value
        forbidden = definition.pos_optional[name].forbidden

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [pos_optional, " + real_type.__name__ + "]: " + str(default) + choices_string, file=configfile)

    # Optional
    for name in definition.optional:

        # Get properties
        real_type = definition.optional[name].type
        description = definition.optional[name].description
        default = definition.optional[name].default
        choices = definition.optional[name].choices
        letter = definition.optional[name].letter
        dynamic_list = definition.optional[name].dynamic_list
        suggestions = definition.optional[name].suggestions
        min_value = definition.optional[name].min_value
        max_value = definition.optional[name].max_value
        forbidden = definition.optional[name].forbidden

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [optional, " + real_type.__name__ + "]: " + str(default) + choices_string, file=configfile)

    # Flag
    for name in definition.flags:

        # Get properties
        description = definition.flags[name].description
        letter = definition.flags[name].letter
        default = definition.flags[name].default  # True or False

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [flag]: " + str(default), file=configfile)

    # Sections
    for section_name in definition.sections:

        # Get properties
        section_definition = definition.sections[section_name]
        section_description = definition.section_descriptions[section_name]

        print(indent + "# " + section_description, file=configfile)
        print(indent + section_name + "[section]:", file=configfile)
        print(indent + "{", file=configfile)

        # Write the section definition
        write_definition(section_definition, configfile, indent + "    ")

        print(indent + "}", file=configfile)

# -----------------------------------------------------------------

def add_settings_from_dict(config, definition, dictionary, check_required=True):

    """
    This function ...
    :param config:
    :param definition:
    :param dictionary:
    :param check_required:
    :return:
    """

    dict_that_is_emptied = copy.deepcopy(dictionary)

    removed_keys = []

    # Fixed
    for name in definition.fixed:

        # Get properties
        description = definition.fixed[name].description
        value = definition.fixed[name].value

        # Set the value in the config
        config[name] = value

    # Required
    for name in definition.required:

        if name not in dictionary:

            if check_required: raise ValueError("The option '" + name + "' is not specified in the configuration dictionary")
            else: value = None

        else:

            # Get properties
            real_type = definition.required[name].type
            description = definition.required[name].description
            choices = definition.required[name].choices
            dynamic_list = definition.required[name].dynamic_list
            suggestions = definition.required[name].suggestions
            min_value = definition.required[name].min_value
            max_value = definition.required[name].max_value
            forbidden = definition.required[name].forbidden

            # Get the value specified in the dictionary
            value = dict_that_is_emptied.pop(name)
            removed_keys.append(name)

            # Checks
            check_list(name, value, real_type)
            check_tuple(name, value, real_type)
            check_dictionary(name, value, real_type)
            check_type(name, value, real_type)
            check_forbidden(name, value, real_type, forbidden)
            check_min(name, value, real_type, min_value)
            check_max(name, value, real_type, max_value)
            check_choices(name, value, real_type, choices)

        # Set the value
        config[name] = value

    # Positional optional
    for name in definition.pos_optional:

        # Get properties
        real_type = definition.pos_optional[name].type
        description = definition.pos_optional[name].description
        default = definition.pos_optional[name].default
        choices = definition.pos_optional[name].choices
        dynamic_list = definition.pos_optional[name].dynamic_list
        suggestions = definition.pos_optional[name].suggestions
        min_value = definition.pos_optional[name].min_value
        max_value = definition.pos_optional[name].max_value
        forbidden = definition.pos_optional[name].forbidden

        # Check if this option is specified in the dictionary
        if name in dictionary:

            # Get the value
            value = dict_that_is_emptied.pop(name)
            removed_keys.append(name)

            if value is not None:

                # Checks
                check_list(name, value, real_type)
                check_tuple(name, value, real_type)
                check_dictionary(name, value, real_type)
                check_type(name, value, real_type)
                check_forbidden(name, value, real_type, forbidden)
                check_min(name, value, real_type, min_value)
                check_max(name, value, real_type, max_value)
                check_choices(name, value, real_type, choices)

        # Use the default value otherwise
        else: value = default

        # Set the value
        config[name] = value

    # Optional
    for name in definition.optional:

        # Get properties
        real_type = definition.optional[name].type
        description = definition.optional[name].description
        default = definition.optional[name].default
        choices = definition.optional[name].choices
        letter = definition.optional[name].letter
        dynamic_list = definition.optional[name].dynamic_list
        suggestions = definition.optional[name].suggestions
        min_value = definition.optional[name].min_value
        max_value = definition.optional[name].max_value
        forbidden = definition.optional[name].forbidden

        # Check if this option is specified in the dictionary
        if name in dictionary:

            # Get the value
            value = dict_that_is_emptied.pop(name)
            removed_keys.append(name)

            if value is not None:

                # Checks
                check_list(name, value, real_type)
                check_tuple(name, value, real_type)
                check_dictionary(name, value, real_type)
                check_type(name, value, real_type)
                check_forbidden(name, value, real_type, forbidden)
                check_min(name, value, real_type, min_value)
                check_max(name, value, real_type, max_value)
                check_choices(name, value, real_type, choices)

        # Use the default value otherwise
        else: value = default

        # Set the value
        config[name] = value

    # Flags
    for name in definition.flags:

        # Get properties
        default = definition.flags[name].default  # True or False

        # Check if this option is specified in the dictionary
        if name in dictionary:
            value = dict_that_is_emptied.pop(name)
            removed_keys.append(name)

        # Use the default value otherwise
        else: value = default

        # Set the boolean value
        config[name] = value

    # Add the configuration settings of the various sections
    for name in definition.sections:

        # Create a map for the settings
        config[name] = Map()

        # Recursively add the settings
        section_definition = definition.sections[name]

        if name in dictionary: section_dictionary = dict_that_is_emptied[name]
        else: section_dictionary = dict() # new empty dict

        # Add the settings
        removed_section_keys = add_settings_from_dict(config[name], section_definition, section_dictionary, check_required=check_required)

        if name in dictionary:
            #print(dict_that_is_emptied[name])
            for key in removed_section_keys:
                del dict_that_is_emptied[name][key]
            #print(dict_that_is_emptied[name])
            if len(dict_that_is_emptied[name]) == 0: del dict_that_is_emptied[name]

    # Add leftover settings
    add_nested_dict_values_to_map(config, dict_that_is_emptied)

    # Return the list of removed keys
    return removed_keys

# -----------------------------------------------------------------

def check_list(name, value, real_type):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :return:
    """

    # For lists: check here
    # List-type setting
    if real_type.__name__.endswith("_list") and not isinstance(value, list): raise ValueError("The option '" + name + "' must be a list")

# -----------------------------------------------------------------

def check_tuple(name, value, real_type):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :return:
    """

    # For tuples: check here
    if real_type.__name__.endswith("_tuple") and not isinstance(value, tuple): raise ValueError("The option '" + name + "' must be a tuple")

# -----------------------------------------------------------------

def check_dictionary(name, value, real_type):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :return:
    """

    # For dictionaries: check here
    if real_type.__name__.endswith("_dictionary") and not isinstance(value, dict): raise ValueError("The option '" + name + "' must be a dictionary")

# -----------------------------------------------------------------

def check_type(name, value, real_type):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :return:
    """

    # TODO: Use derived and related types to check the type of the given value

    # List-type setting
    #if real_type.__name__.endswith("_list"):
    #    pass

    # Single-value setting
    #else:
    #    pass

    return

# -----------------------------------------------------------------

def check_forbidden(name, value, real_type, forbidden):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :param forbidden:
    :return:
    """

    # Check forbidden
    if forbidden is None: return

    # Check whether forbidden is a list
    if not isinstance(forbidden, list): raise ValueError("Forbidden values for '" + name + "' must be a list")

    # List-type setting
    if real_type.__name__.endswith("_list"):

        # Loop over the entries
        for entry in value:
            if entry in forbidden: raise ValueError("The value '" + stringify.stringify(entry)[1] + "' in '" + name + "' is forbidden")

    # Single-value type setting
    else:

        if value in forbidden: raise ValueError("The value '" + stringify.stringify(value)[1] + " for '" + name + "' is forbidden")

# -----------------------------------------------------------------

def check_min(name, value, real_type, min_value):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :param min_value:
    :return:
    """

    # Check min
    if min_value is not None:

        # List type setting
        if real_type.__name__.endswith("_list"):

            # Loop over the entries
            for entry in value:

                if entry < min_value: raise ValueError("The value '" + stringify.stringify(entry)[1] + "' for '" + name + "' is too low: minimum value is '" + stringify.stringify(min_value)[1] + "'")

        # Single-value type setting
        else:

            if value < min_value: raise ValueError("The value '" + stringify.stringify(value)[1] + "' for '" + name + "' is too low: minimum value is '" + stringify.stringify(min_value)[1] + "'")

# -----------------------------------------------------------------

def check_max(name, value, real_type, max_value):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :param max_value:
    :return:
    """

    if max_value is not None:

        # List type setting
        if real_type.__name__.endswith("_list"):

            # Loop over the entries
            for entry in value:

                if entry > max_value: raise ValueError("The value '" + stringify.stringify(entry)[1] + "' for '" + name + "' is too high: maximum value is '" + stringify.stringify(max_value)[1] + "'")

        # Single-value type setting
        else:

            if value > max_value: raise ValueError("The value '" + stringify.stringify(value)[1] + "' for '" + name + "' is too high: maximum value is '" + stringify.stringify(max_value)[1] + "'")

# -----------------------------------------------------------------

def check_choices(name, value, real_type, choices):

    """
    This function ...
    :param name:
    :param value:
    :param real_type:
    :param choices:
    :return:
    """

    # No choices
    if choices is None: return

    # Check that choices is list
    if not isinstance(choices, list) and not isinstance(choices, dict): raise ValueError("Choices for '" + name + "' must be a list")

    # More than one choice
    if len(choices) > 1:

        # Convert the choices into a string
        choices_string = stringify.stringify(choices.keys())[1] if isinstance(choices, dict) else stringify.stringify(choices)[1]

        # Check whether the given value matches the allowed choices

        # List-type setting
        if real_type.__name__.endswith("_list"):  # list-type setting

            # Loop over the entries in the given list
            for entry in value:

                # Check whether each entry appears in the choices
                if entry not in choices: raise ValueError("The value '" + stringify.stringify(entry)[1] + "' in '" + name + "' is not valid, choices are: '" + choices_string + "'")

        # Single-value setting and not None
        else:
            if value not in choices: raise ValueError("The value '" + stringify.stringify(value)[1] + "' for the option '" + name + "' is not valid, choices are: '" + choices_string + "'")

    # Exactly one choice
    else:

        # List-type setting
        if real_type.__name__.endswith("_list"):  # list-type setting

            if isinstance(choices, dict): exact_value = [choices.keys()[0]]
            else: exact_value = [choices[0]]

        # Single-value setting
        else:

            if isinstance(choices, dict): exact_value = choices.keys()[0]
            else: exact_value = choices[0]

        # Check with the value that the only exact choice
        if value != exact_value: raise ValueError("The value '" + stringify.stringify(value)[1] + "' for '" + name + "' is not valid: only allowed option is '" + stringify.stringify(exact_value)[1] + "'")

# -----------------------------------------------------------------

def add_settings_default(config, definition):

    """
    This function ...
    :param config:
    :param definition:
    :return:
    """

    # Fixed
    for name in definition.fixed:
        value = definition.fixed[name].value
        config[name] = value

    # Required
    for name in definition.required:

        # Set the value
        config[name] = None

    # Positional optional
    for name in definition.pos_optional:

        # Get properties
        default = definition.pos_optional[name].default

        # Set the value
        config[name] = default

    # Optional
    for name in definition.optional:

        # Get properties
        default = definition.optional[name].default

        # Set the value
        config[name] = default

    # Flags
    for name in definition.flags:

        # Get properties
        default = definition.flags[name].default  # True or False

        # Set the boolean value
        config[name] = default

    # Add the configuration settings of the various sections
    for name in definition.sections:

        # Create a map for the settings
        config[name] = Map()

        # Recursively add the settings
        section_definition = definition.sections[name]

        # Add the settings
        add_settings_default(config[name], section_definition)

# -----------------------------------------------------------------

def add_settings_interactive(config, definition, prompt_optional=True, settings=None, use_default=None, cue="   : "):

    """
    This function ...
    :param config:
    :param definition:
    :param prompt_optional:
    :param settings:
    :param use_default:
    :param cue:
    :return:
    """

    from ..tools import formatting as fmt

    # Fixed
    for name in definition.fixed:

        # Check
        if use_default is not None and name in use_default: raise ValueError("Cannot use the default value for a fixed setting")

        # Check in settings
        if settings is not None and name in settings:
            config[name] = settings[name]
            continue

        # Get properties
        description = definition.fixed[name].description
        value = definition.fixed[name].value

        # Give name and description
        log.success(name + ": " + description)

        # Inform the user
        log.info("Using fixed value '" + tostr(value) + "' for " + name)

        # Set the value
        config[name] = value

    # Required
    for name in definition.required:

        # Check in settings
        if settings is not None and name in settings:
            config[name] = settings[name]
            continue

        # Check
        if use_default is not None and name in use_default: raise ValueError("Cannot use the default value for a required setting")

        # Get properties
        real_type = definition.required[name].type
        description = definition.required[name].description
        choices = definition.required[name].choices
        dynamic_list = definition.required[name].dynamic_list
        suggestions = definition.required[name].suggestions
        min_value = definition.required[name].min_value
        max_value = definition.required[name].max_value
        forbidden = definition.required[name].forbidden

        # Give name and description
        log.success(name + ": " + description)

        # Get list of choices and a dict of their descriptions
        if choices is not None:
            choices_list = choices.keys() if isinstance(choices, dict) else choices
            choice_descriptions = choices if isinstance(choices, dict) else None
        else: choices_list = choice_descriptions = None

        # No choices
        if choices is None:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                # Dynamic list
                if dynamic_list:

                    real_base_type = getattr(parsing, real_type.__name__.split("_list")[0])

                    # Construct type
                    the_type = construct_type(real_base_type, min_value, max_value, forbidden)

                    log.info("Provide a value for a list item and press ENTER. Leave blank and press ENTER to end the list")

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = []
                    while True:
                        answer = raw_input(cue)
                        if answer == "": break # end of the list
                        try:
                            #single_value = real_base_type(answer)
                            single_value = the_type(answer)
                            value.append(single_value)
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

                # Not a dynamic list
                else:

                    log.info("Provide the values, seperated by commas")

                    # Construct type
                    the_type = construct_type(real_type, min_value, max_value, forbidden)

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = []  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input(cue)
                        try:
                            #value = real_type(answer)
                            value = the_type
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            # Single-value setting
            else:

                log.info("Provide a value")

                # Construct type
                the_type = construct_type(real_type, min_value, max_value, forbidden)

                # Show suggestions
                if suggestions is not None:
                    log.info("Suggestions:")
                    for suggestion in suggestions:
                        if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                        else: suggestion_description = ""
                        log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                value = None # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    answer = raw_input(cue)
                    try:
                        #value = real_type(answer)
                        value = the_type(answer)
                        break
                    except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # More than one choice
        elif len(choices) > 1:

            # List-type setting
            if real_type.__name__.endswith("_list"): # list-type setting

                log.info("Choose one or more of the following options (separated only by commas)")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = None # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the numbers of the choice
                    answer = raw_input(cue)
                    try:
                        indices = parsing.integer_list(answer)
                        value = [choices_list[index] for index in indices] # value is a list
                        break
                    except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

            # Single-value setting
            else:

                log.info("Choose one of the following options")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = None  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the number of the choice
                    answer = raw_input(cue)
                    try:
                        index = parsing.integer(answer)
                        value = choices_list[index]
                        break
                    except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # Only one choice
        #elif len(choices) == 1:
        else:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                if isinstance(choices, dict):
                    log.info("Only one option: automatically using a list of this value '[" + str(choices.keys()[0]) + "]' for " + name)
                    value = [choices.keys()[0]]
                else:
                    # Inform the user
                    log.info("Only one option: automatically using a list of this value '[" + str(choices[0]) + "]' for " + name)
                    value = [choices[0]]

            # Single-value setting
            else:

                if isinstance(choices, dict):
                    # Inform the user
                    log.info("Only one option: automatically using value of '" + str(choices.keys()[0]) + "' for " + name)
                    value = choices.keys()[0]
                else:
                    # Inform the user
                    log.info("Only one option: automatically using value of '" + str(choices[0]) + "' for " + name)
                    value = choices[0]

        # Set the value
        config[name] = value

    # Positional optional
    for name in definition.pos_optional:

        # Check in settings
        if settings is not None and name in settings:
            config[name] = settings[name]
            continue

        # Get properties
        real_type = definition.pos_optional[name].type
        description = definition.pos_optional[name].description
        default = definition.pos_optional[name].default
        choices = definition.pos_optional[name].choices
        dynamic_list = definition.pos_optional[name].dynamic_list
        suggestions = definition.pos_optional[name].suggestions
        min_value = definition.pos_optional[name].min_value
        max_value = definition.pos_optional[name].max_value
        forbidden = definition.pos_optional[name].forbidden
        default_alias = definition.pos_optional[name].default_alias

        # Check
        if use_default is not None and name in use_default:
            config[name] = default
            continue

        # Get list of choices and a dict of their descriptions
        if choices is not None:
            choices_list = choices.keys() if isinstance(choices, dict) else choices
            choice_descriptions = choices if isinstance(choices, dict) else None
        else: choices_list = choice_descriptions = None

        if not prompt_optional:
            value = default
            # Set the value
            config[name] = value
            continue

        # Give name and description
        log.success(name + ": " + description)

        # Show default value
        if default_alias is not None: log.info("Press ENTER to use the default value (" + default_alias + ")")
        else: log.info("Press ENTER to use the default value (" + stringify.stringify(default)[1] + ")")

        # Choices are not given
        if choices_list is None:

            # List-typ setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                # Dynamic list
                if dynamic_list:

                    log.info(
                        "or provide other values. Enter a value and press ENTER. To end the list, leave blank and press ENTER.")

                    real_base_type = getattr(parsing, real_type.__name__.split("_list")[0])

                    # Construct type
                    the_type = construct_type(real_base_type, min_value, max_value, forbidden)

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = []  # to remove warning
                    while True:
                        answer = raw_input(cue)
                        if answer == "":
                            break  # end of list
                        else:
                            try:
                                # single_value = real_type(answer)
                                single_value = the_type(answer)
                                value.append(single_value)
                            except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

                # Not a dynamic list
                else:

                    log.info("or provide other values, separated by commas")

                    # Construct type
                    the_type = construct_type(real_type, min_value, max_value, forbidden)

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input(cue)
                        if answer == "":
                            value = default
                            break
                        else:
                            try:
                                # value = real_type(answer)
                                value = the_type(answer)
                                break
                            except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

            # Not a list
            else:

                log.info("or provide another value")

                # Construct type
                the_type = construct_type(real_type, min_value, max_value, forbidden)

                # Show suggestions
                if suggestions is not None:
                    log.info("Suggestions:")
                    for suggestion in suggestions:
                        if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                        else: suggestion_description = ""
                        log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            # value = real_type(answer)
                            value = the_type(answer)
                            break
                        except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # Choices are given
        #if choices_list is not None:
        # More than one choice or no default
        elif len(choices) > 1 or default is None:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                log.info("or choose one or more of the following options (separated only by commas)")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = default # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the numbers of the choice
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            indices = parsing.integer_list(answer)
                            value = [choices_list[index] for index in indices] # value is a list
                            break
                        except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("or choose one of the following options")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the number of the choice
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            index = parsing.integer(answer)
                            value = choices_list[index]
                            break
                        except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # Exactly one choice AND a default value
        else:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                if isinstance(choices, dict):
                    log.info("Only one option: automatically using a list of this value '[" + tostr(choices.keys()[0]) + "]' for " + name)
                    value = [choices.keys()[0]]
                    assert value == default
                else:
                    # Inform the user
                    log.info("Only one option: automatically using a list of this value '[" + tostr(choices[0]) + "]' for " + name)
                    value = [choices[0]]
                    assert value == default

            # Single-value setting
            else:

                if isinstance(choices, dict):
                    # Inform the user
                    log.info(
                        "Only one option: automatically using value of '" + tostr(choices.keys()[0]) + "' for " + name)
                    value = choices.keys()[0]
                    assert value == default
                else:
                    # Inform the user
                    log.info("Only one option: automatically using value of '" + tostr(choices[0]) + "' for " + name)
                    value = choices[0]
                    assert value == default

        # Set the value
        config[name] = value

    # Optional
    for name in definition.optional:

        # Check settings
        if settings is not None and name in settings:
            config[name] = settings[name]
            continue

        # Get properties
        real_type = definition.optional[name].type
        description = definition.optional[name].description
        default = definition.optional[name].default
        choices = definition.optional[name].choices
        letter = definition.optional[name].letter
        dynamic_list = definition.optional[name].dynamic_list
        suggestions = definition.optional[name].suggestions
        min_value = definition.optional[name].min_value
        max_value = definition.optional[name].max_value
        forbidden = definition.optional[name].forbidden
        default_alias = definition.optional[name].default_alias

        # Check
        if use_default is not None and name in use_default:
            config[name] = default
            continue

        # Get list of choices and a dict of their descriptions
        if choices is not None:
            choices_list = choices.keys() if isinstance(choices, dict) else choices
            choice_descriptions = choices if isinstance(choices, dict) else None
        else: choices_list = choice_descriptions = None

        if not prompt_optional:
            value = default
            # Set the value
            config[name] = value
            continue

        # Give name and description
        log.success(name + ": " + description)

        # Show default value
        if default_alias is not None: log.info("Press ENTER to use the default value (" + default_alias + ")")
        else: log.info("Press ENTER to use the default value (" + stringify.stringify(default)[1] + ")")

        # Choices are not given
        if choices_list is None:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                # Dynamic list
                if dynamic_list:

                    log.info("or provide other values. Enter a value and press ENTER. To end the list, leave blank and press ENTER.")

                    real_base_type = getattr(parsing, real_type.__name__.split("_list")[0])

                    # Construct type
                    the_type = construct_type(real_base_type, min_value, max_value, forbidden)

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = []  # to remove warning
                    while True:
                        answer = raw_input(cue)
                        if answer == "":
                            break  # end of list
                        else:
                            try:
                                # single_value = real_type(answer)
                                single_value = the_type(answer)
                                value.append(single_value)
                            except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

                # Not dynamic list
                else:

                    log.info("or provide other values, separated by commas")

                    # Construct type
                    the_type = construct_type(real_type, min_value, max_value, forbidden)

                    # Show suggestions
                    if suggestions is not None:
                        log.info("Suggestions:")
                        for suggestion in suggestions:
                            if isinstance(suggestion, dict): suggestion_description = ": " + suggestions[suggestion]
                            else: suggestion_description = ""
                            log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold + suggestion_description)

                    value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input(cue)
                        #print(answer)
                        if answer == "":
                            value = default
                            break
                        else:
                            try:
                                # value = real_type(answer)
                                value = the_type(answer)
                                break
                            except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

            # Single-value setting
            else:

                log.info("or provide another value")

                # Construct type
                the_type = construct_type(real_type, min_value, max_value, forbidden)

                # Show suggestions
                if suggestions is not None:
                    log.info("Suggestions:")
                    for suggestion in suggestions:
                        if isinstance(suggestions, dict): suggestion_description = ": " + suggestions[suggestion]
                        else: suggestion_description = ""
                        log.info(" - " + fmt.bold + stringify.stringify(suggestion)[1] + fmt.reset_bold +  suggestion_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the input
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            # value = real_type(answer)
                            value = the_type(answer)
                            break
                        except ValueError as e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # Choices are given
        #if choices_list is not None:
        # If more choices are given or no default is given
        elif len(choices_list) > 1 or default is None:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                log.info("or choose one or more of the following options (separated only by commas)")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:

                    # Get the numbers of the choice
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            indices = parsing.integer_list(answer)
                            value = [choices_list[index] for index in indices] # value is a list here
                            break
                        except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

            # Not a list
            else:

                log.info("or choose one of the following options")

                for index, choice_value in enumerate(choices_list):
                    label = tostr(choice_value)
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[choice_value]
                    log.info(" - [" + str(index) + "] " + fmt.bold + label + fmt.reset_bold + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:

                    # Get the number of the choice
                    answer = raw_input(cue)
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            index = parsing.integer(answer)
                            value = choices_list[index] # if we are here, no error was raised
                            break
                        except (ValueError, IndexError) as e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # No choices
        else:

            # List-type setting
            if real_type.__name__.endswith("_list"):  # list-type setting

                if isinstance(choices, dict):
                    log.info("Only one option: automatically using a list of this value '[" + tostr(choices.keys()[0]) + "]' for " + name)
                    value = [choices.keys()[0]]
                    assert value == default
                else:
                    # Inform the user
                    log.info("Only one option: automatically using a list of this value '[" + tostr(choices[0]) + "]' for " + name)
                    value = [choices[0]]
                    assert value == default

            # Single-value setting
            else:

                if isinstance(choices, dict):
                    # Inform the user
                    log.info(
                        "Only one option: automatically using value of '" + tostr(choices.keys()[0]) + "' for " + name)
                    value = choices.keys()[0]
                    assert value == default
                else:
                    # Inform the user
                    log.info("Only one option: automatically using value of '" + tostr(choices[0]) + "' for " + name)
                    value = choices[0]
                    assert value == default

        # Set the value
        config[name] = value

    # Flags
    for name in definition.flags:

        # Check settings
        if settings is not None and name in settings:
            config[name] = settings[name]
            continue

        # Get properties
        description = definition.flags[name].description
        letter = definition.flags[name].letter
        default = definition.flags[name].default  # True or False

        # Check
        if use_default is not None and name in use_default:
            config[name] = default
            continue

        if not prompt_optional:
            value = default
            # Set the value
            config[name] = value
            continue

        # Give name and description
        log.success(name + ": " + description)

        # Ask the question
        log.info("Do you want '" + name + "' to be enabled or not (y or n) or press ENTER for the default (" + tostr(default) + ")")

        value = default  # to remove warning from IDE that value could be referenced (below) without assignment
        while True:
            answer = raw_input(cue)
            if answer == "":
                value = default
                break
            else:
                try:
                    value = parsing.boolean(answer) # if this passes without error, we have valid input
                    break
                except ValueError as e:
                    # Give warning and go to the next iteration
                    log.warning("Invalid input: " + str(e) + ". Try again.")

        # Set the value
        config[name] = value

    # Add the configuration settings of the various sections
    for name in definition.sections:

        # Create a map for the settings
        config[name] = Map()

        # Recursively add the settings
        section_definition = definition.sections[name]
        section_description = definition.section_descriptions[name]

        # Give name and description
        log.success(name + ": " + section_description + " (section)")

        # Check in settings
        settings_section = settings[name] if settings is not None and name in settings else None

        # Add the settings
        add_settings_interactive(config[name], section_definition, prompt_optional=prompt_optional, settings=settings_section)

# -----------------------------------------------------------------

def add_nested_dict_values_to_map(mapping, dictionary):

    """
    This function ...
    :param mapping:
    :param dictionary:
    :return:
    """

    # Loop over the dictionary items
    for name in dictionary:

        # Get the value
        value = dictionary[name]

        # Recursively call this function if the value is an other dictionary
        if isinstance(value, dict):

            # Already existing submapping
            if name in mapping: add_nested_dict_values_to_map(mapping[name], value)
            else:
                mapping[name] = Map()
                add_nested_dict_values_to_map(mapping[name], value)

        # Else, add the value to the mapping
        else: mapping[name] = value

# -----------------------------------------------------------------

def smart_type(argument, real_type, min_value, max_value, forbidden):

    """
    This function ...
    :param argument:
    :param real_type:
    :param min_value:
    :param max_value:
    :param forbidden:
    :return:
    """

    parsed = real_type(argument)
    if min_value is not None and parsed < min_value: raise ValueError("Value should be higher than " + stringify.stringify_not_list(min_value)[1])
    if max_value is not None and parsed > max_value: raise ValueError("Value should be lower than " + stringify.stringify_not_list(max_value)[1])
    if forbidden is not None and parsed in forbidden: raise ValueError("Value " + stringify.stringify_not_list(parsed) + " is forbidden")

    # All checks passed, return the parsed value
    return parsed

# -----------------------------------------------------------------

def smart_list_type(argument, real_base_type, min_value, max_value, forbidden):

    """
    This function ...
    :param argument:
    :param real_base_type:
    :param min_value:
    :param max_value:
    :param forbidden:
    :return:
    """

    arguments = [argument.strip() for argument in argument.split(",")]
    return [smart_type(argument, real_base_type, min_value, max_value, forbidden) for argument in arguments]

# -----------------------------------------------------------------

def construct_type(real_type, min_value, max_value, forbidden):

    """
    This function ...
    :param real_type:
    :param min_value:
    :param max_value:
    :param forbidden:
    :return:
    """

    # List type
    if real_type.__name__.endswith("list") and not hasattr(parsing, real_type.__name__):

        # Determine base type parsing function
        base_type = getattr(parsing, real_type.__name__.split("_list")[0])
        the_type = partial(smart_list_type, **{"real_base_type": base_type, "min_value": min_value, "max_value": max_value, "forbidden": forbidden})

        # Return the new function
        return the_type

    else:

        # Construct the actual type
        the_type = partial(smart_type, **{"real_type": real_type, "min_value": min_value, "max_value": max_value, "forbidden": forbidden})

        # Return the new function
        return the_type

# -----------------------------------------------------------------

def parse_default(default, user_type, real_type):

    """
    This function ...
    :param default:
    :param user_type:
    :param real_type:
    :return:
    """

    user_type_name = user_type if types.is_string_type(user_type) else user_type.__name__

    if user_type_name.endswith("list") and types.is_sequence(default):

        real_base_type = getattr(parsing, user_type_name.split("_list")[0])
        default = [get_real_value(arg, real_base_type) for arg in default]

    elif user_type_name.endswith("tuple") and types.is_tuple(default):

        real_base_type = getattr(parsing, user_type_name.split("_tuple")[0])
        default = [get_real_value(arg, real_base_type) for arg in default]

    else: default = get_real_value(default, real_type)

    # Return the default value
    return default

# -----------------------------------------------------------------

def check_default(default, user_type):

    """
    This function ...
    :param default:
    :param user_type:
    :return:
    """

    default_type, default_string = stringify.stringify(default)
    #print(default_type, default_string)
    if default_type != user_type and not are_related_types(default_type, user_type):

        # List-like property
        if user_type.endswith("list") or user_type.endswith("tuple"):

            if user_type.endswith("list"): base_type = user_type.split("_list")[0]
            else: base_type = user_type.split("_tuple")[0]

            # Ascending and descending lists
            if base_type.startswith("ascending"):
                base_type = base_type.split("ascending_")[1]
                direction = "ascending"
            elif base_type.startswith("descending"):
                base_type = base_type.split("descending_")[1]
                direction = "descending"
            else: direction = None

            # e.g. for 'integer_and_string_list' > integer_or_string
            if "_and_" in base_type: base_type = base_type.replace("_and_", "_or_")

            new_default = []
            for value in default:
                try:
                    #print(value, type(value), base_type)
                    value = check_default_single_value(value, base_type)
                    new_default.append(value)
                except ValueError: raise ValueError("Default value '" + tostr(default) + "' is not of the right type '" + user_type + "'")
            default = new_default

            # Check if ascending or descending
            if direction is not None:
                if direction == "ascending" and not sequences.is_ascending(default): raise ValueError("Default list is not ascending")
                if direction == "descending" and not sequences.is_descending(default): raise ValueError("Default list is not descending")

        # Single-value property
        else: default = try_to_convert_to_type(default, user_type)

    # Return the check default value
    return default

# -----------------------------------------------------------------

def check_default_single_value(default, user_type):

    """
    Thi function ...
    :param value:
    :param user_type:
    :return:
    """

    default_type, default_string = stringify.stringify(default)
    if default_type != user_type and not are_related_types(default_type, user_type):

        default = try_to_convert_to_type(default, user_type)

    # Return the checked default value
    return default

# -----------------------------------------------------------------

def try_to_convert_to_type(default, user_type):

    """
    This function ...
    :param default:
    :param user_type:
    :return:
    """

    if user_type == "mixed": return default
    elif parent_type(user_type) == "integer" and numbers.is_integer(default): return int(default)
    elif parent_type(user_type) == "real" and numbers.is_integer(default): return float(default)
    elif types.is_string_type(default): return try_to_convert_from_string(default, user_type)
    else:
        #raise ValueError("Default value '" + str(default) + "' could not be converted to the right type '" + user_type + "'")
        try: return try_to_convert_from_string(tostr(default), user_type)
        except ValueError: raise ValueError("Default value '" + tostr(default) + "' could not be converted to the right type '" + user_type + "'")

# -----------------------------------------------------------------

def try_to_convert_from_string(string, user_type):

    """
    This function ...
    :param string: 
    :param user_type: 
    :return: 
    """

    # Get the parsing function
    parsing_function = getattr(parsing, user_type)

    try:
        value = parsing_function(string)
        return value
    except ValueError: raise ValueError("String '" + tostr(string) + "' could not be converted to the right type '" + user_type + "'")

# -----------------------------------------------------------------
