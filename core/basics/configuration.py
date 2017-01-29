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
import inspect
from abc import ABCMeta, abstractmethod
from types import NoneType
import sys
import argparse
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .map import Map
from ..tools import parsing, stringify
from ..tools import filesystem as fs
from ..tools.logging import log
from .composite import SimplePropertyComposite

# -----------------------------------------------------------------

subtypes = dict()
subtypes["integer"] = ["positive_integer", "negative_integer"]
subtypes["real"] = ["fraction", "positive_real", "negative_real"]
subtypes["string"] = ["file_path"]
subtypes["quantity"] = ["photometric_quantity", "photometric_density_quantity"]
subtypes["unit"] = ["photometric_unit", "photometric_density_unit"]

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

    #@property
    def config_dir_path(self):

        """
        The directory where the config file should be saved
        :return:
        """

        if "config_path" in self:
            if self["config_path"] is not None: return fs.absolute_or_in(self["config_path"], self["path"]) # absolute path or relative to the working directory
            else: return self.output_path()
        else: return None

    # -----------------------------------------------------------------

    #@property
    def output_path(self):

        """
        This function ...
        :return:
        """

        # If 'output' is defined in the config
        if "output" in self:

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

def load_mapping(mappingfile, mapping, indent=""):

    """
    This function ...
    :param mappingfile
    :param mapping:
    :param indent:
    :return:
    """

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

            try:
                before, after = line.split(":")
            except ValueError: print("ERROR processing: ", line)

            #print("declaration", before, after)

            if after.strip() == "":

                state = 3  # state 3: just before "{" is expected to start subdefinition
                continue

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

                    parsing_function = getattr(parsing, specification)
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
            mapping[name] = Map()

            # Load the submapping
            load_mapping(mappingfile, mapping[name], indent=indent+"    ")

            state = 2
            description = None

        # End of section or complete definition
        elif line.startswith("}"): return

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
            print(indent + name + " [" + ptype + "]: " + string, file=mappingfile)

        if index != length - 1: print("", file=mappingfile)
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

class ConfigurationDefinition(object):

    """
    This function ...
    """

    def __init__(self, prefix=None, log_path=None, config_path=None, write_config=True):

        """
        This function ...
        :param prefix:
        :param log_path:
        :param config_path:
        :param write_config:
        """

        # Prefix
        self.prefix = prefix

        # Log path and config path
        self.log_path = log_path
        self.config_path = config_path
        self.write_config = write_config

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
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Create the definition
        definition = cls()

        # Load the definition
        with open(path, 'r') as configfile: load_definition(configfile, definition)

        # Return the definition
        return definition

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
        for name in self.required:

            # Get properties
            real_type = self.required[name].type
            description = self.required[name].description
            choices = self.required[name].choices

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None: name = prefix + "/" + name

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list"): choices = None

            # Add argument to argument parser
            parser.add_argument(name, type=real_type, help=description, choices=choices)

        # Positional optional
        for name in self.pos_optional:

            # Get properties
            real_type = self.pos_optional[name].type
            description = self.pos_optional[name].description
            default = self.pos_optional[name].default
            choices = self.pos_optional[name].choices

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None: name = prefix + "/" + name

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list"): choices = None

            # Add argument to argument parser
            parser.add_argument(name, type=real_type, help=description, default=default, nargs='?', choices=choices)

        # Optional
        for name in self.optional:

            # Get properties
            real_type = self.optional[name].type
            description = self.optional[name].description
            default = self.optional[name].default
            choices = self.optional[name].choices
            letter = self.optional[name].letter

            # Add prefix
            #if self.prefix is not None: name = self.prefix + "/" + name

            if prefix is not None:
                name = prefix + "/" + name
                letter = None

            # Don't set choices for 'list'-type argument values, the choices here are allowed to be entered in any combination. Not just one of the choices is expected.
            if real_type.__name__.endswith("_list"): choices = None

            # Add the argument
            if letter is None: parser.add_argument("--" + name, type=real_type, help=description, default=default, choices=choices)
            else: parser.add_argument("-" + letter, "--" + name, type=real_type, help=description, default=default, choices=choices)

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
            else: raise ValueError("Invalid option for default of flag '" + name + "': " + str(default))

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
                    if not item in choices: raise ValueError("Element '" + item + "' not recognized. Options are: " + stringify.stringify(choices)[1])

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
                    if not item in choices: raise ValueError("Element '" + item + "' not recognized. Options are: " + stringify.stringify(choices)[1])

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
                    if not item in choices: raise ValueError("Element '" + item + "' not recognized. Options are: " + stringify.stringify(choices)[1])

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
            else: raise ValueError("Invalid value for default of flag '" + name + "': " + str(default))

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

    def import_section(self, name, description, definition):

        """
        This function ...
        :param name:
        :param description:
        :param definition:
        :return:
        """

        # Determine prefix
        #if self.prefix is None: prefix = name
        #else: prefix = self.prefix + "/" + name
        prefix = name

        # Add the section
        self.sections[name] = copy.deepcopy(definition)
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

            default_value = default_values[pname]

            # Section
            if isinstance(default_value, SimplePropertyComposite):

                # Recursive call to this function
                pdescription = instance._descriptions[pname]
                self.sections[name].import_section_from_properties(pname, pdescription, default_value)

            else:

                # Get the type, description and the choices
                ptype = instance._types[pname]
                pdescription = instance._descriptions[pname]
                choices = instance._choices[pname]

                # Add to the definition
                if ptype == "boolean": self.sections[name].add_flag(pname, pdescription, default_value)
                else: self.sections[name].add_optional(pname, ptype, pdescription, default_value, choices=choices)

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

    def add_required(self, name, user_type, description, choices=None, dynamic_list=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param choices:
        :param dynamic_list:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Add
        self.required[name] = Map(type=real_type, description=description, choices=choices, dynamic_list=dynamic_list)

    # -----------------------------------------------------------------

    def add_positional_optional(self, name, user_type, description, default=None, choices=None,
                                convert_default=False, dynamic_list=False):

        """
        This function ...
        :param name:
        :param user_type:
        :param description:
        :param default:
        :param choices:
        :param convert_default:
        :param dynamic_list:
        :return:
        """

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if default is not None:
            if convert_default: default = get_real_value(default, real_type)
            else: # check default
                default_type, default_string = stringify.stringify(default)
                if default_type != user_type and user_type not in subtypes[default_type]:
                    raise ValueError("Default value is not of the right type: " + default_type + " instead of " + user_type)

        # Add
        self.pos_optional[name] = Map(type=real_type, description=description, default=default, choices=choices, dynamic_list=dynamic_list)

    # -----------------------------------------------------------------

    def add_optional(self, name, user_type, description, default=None, choices=None, letter=None, convert_default=False, dynamic_list=False):

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
        :return:
        """

        if self.prefix is not None and letter is not None: raise ValueError("Cannot assign letter argument for child configuration definition")

        # Get the real type
        real_type = get_real_type(user_type)

        # Get the real default value
        if default is not None:
            if convert_default: default = get_real_value(default, real_type)
            else: # check default
                default_type, default_string = stringify.stringify(default)
                if default_type != user_type and user_type not in subtypes[default_type]:
                    raise ValueError("Default value is not of the right type: " + default_type + " instead of " + user_type)

        # Add
        self.optional[name] = Map(type=real_type, description=description, default=default, choices=choices, letter=letter, dynamic_list=dynamic_list)

    # -----------------------------------------------------------------

    def add_flag(self, name, description, default=False, letter=None):

        """
        This function ...
        :param name:
        :param description:
        :param default:
        :param letter:
        :return:
        """

        # Add
        self.flags[name] = Map(description=description, letter=letter, default=default)

# -----------------------------------------------------------------

class ConfigurationSetter(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, name, description=None, add_logging=True, add_cwd=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Set the name and description
        self.name = name
        self.description = description

        # The configuration definition
        self.definition = None

        # Set options
        self.add_logging = add_logging
        self.add_cwd = add_cwd

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

    def set_logging_and_cwd(self):

        """
        This function ...
        :return:
        """

        cwd_path = fs.cwd()

        # Add logging options
        if self.add_logging:

            # Log path to absolute path
            log_path = fs.absolute_or_in(self.definition.log_path, cwd_path) if self.definition.log_path is not None else cwd_path # set absolute log path

            self.definition.add_fixed("log_path", "the directory for the log file be written to", log_path)
            self.definition.add_flag("debug", "enable debug output")
            self.definition.add_flag("report", "write a report file")

        # Add config path
        if self.definition.write_config:

            # Set the path to the directory where the configuration file should be saved
            if self.definition.config_path is not None: self.definition.add_fixed("config_path", "the directory for the configuration file to be written to", self.definition.config_path)
            else: self.definition.add_optional("config_path", "directory_path", "the directory for the configuration file to be written to (relative to the working directory or absolute) (if None, the output directory is used)")

        # Add the path to the current working directory
        if self.add_cwd: self.definition.add_fixed("path", "the working directory", cwd_path)

# -----------------------------------------------------------------

class InteractiveConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True):

        """
        The constructor ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(InteractiveConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

    # -----------------------------------------------------------------

    def run(self, definition, prompt_optional=None):

        """
        This function ...
        :param definition:
        :param prompt_optional:
        :return:
        """

        # Set the definition
        self.definition = definition

        # Set logging and cwd
        self.set_logging_and_cwd()

        # Do interactive
        self.interactive(prompt_optional)

        # Return the config
        return self.config

    # -----------------------------------------------------------------

    def interactive(self, prompt_optional):

        """
        This function ...
        :param prompt_optional:
        :return:
        """

        if prompt_optional is None:

            # Ask whether optional parameters have to be shown, or to just use the default values
            log.info("Do you want to configure optional settings (y or n)?")
            log.info("Press ENTER to use the default (True)")

            answer = raw_input("   : ")
            if answer == "": prompt_optional = True
            else: prompt_optional = parsing.boolean(answer)

        add_settings_interactive(self.config, self.definition, prompt_optional=prompt_optional)

# -----------------------------------------------------------------

class ArgumentConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(ArgumentConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

        # Create the command-line parser
        self.parser = argparse.ArgumentParser(prog=name, description=description)

        # The parsed arguments
        self.arguments = None

    # -----------------------------------------------------------------

    @staticmethod
    def get_arguments():

        """
        This function ...
        :return:
        """

        return sys.argv[1:]

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
        self.arguments = self.parser.parse_args()

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

    def __init__(self, path, name, description=None, add_logging=True, add_cwd=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(FileConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

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

    def __init__(self, dictionary, name, description=None, add_logging=True, add_cwd=True):

        """
        The constructor ...
        :param dictionary:
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(DictConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

        # Set the user-provided dictionary
        self.dictionary = dictionary

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
        add_settings_from_dict(self.config, self.definition, self.dictionary)

# -----------------------------------------------------------------

class PassiveConfigurationSetter(ConfigurationSetter):

    """
    This class ...
    """

    def __init__(self, name, description=None, add_logging=True, add_cwd=True):

        """
        The constructor ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(PassiveConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

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

    def __init__(self, path, name, description=None, add_logging=True, add_cwd=True):

        """
        This function ...
        :param name:
        :param description:
        :param add_logging:
        :param add_cwd:
        """

        # Call the constructor of the base class
        super(GraphicalConfigurationSetter, self).__init__(name, description, add_logging, add_cwd)

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

    if isinstance(user_type, basestring): return getattr(parsing, user_type)
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

                kind = specification[0]
                user_type = specification[1].strip() if kind != "flag" else None

                #if user_type == "str": user_type = str
                #elif user_type == "int": user_type = int
                #elif user_type == "float": user_type = float
                #elif user_type == "bool": user_type = bool # for fixed flags ...
                #elif user_type == "None": # for fixed things that should be None
                #    user_type = NoneType
                #    value = None

                if user_type == "None":
                    user_type = NoneType
                    value = None
                elif kind != "flag":
                    user_type = getattr(parsing, user_type)

                # Convert choices to real list of the user_type-'list' type
                choices = None
                if choices_string is not None:

                    user_type_list_name = user_type.__name__ + "_list"
                    user_type_list = getattr(parsing, user_type_list_name)
                    choices = user_type_list(choices_string)

                if kind == "fixed":

                    value = user_type(value) if value is not None else None # convert here because the add_fixed function doesn't bother with types now
                    definition.add_fixed(name, description, value)

                elif kind == "required":

                    definition.add_required(name, user_type, description, choices=choices)

                elif kind == "pos_optional":

                    definition.add_positional_optional(name, user_type, description, default=value, choices=choices, convert_default=True)

                elif kind == "optional":

                    definition.add_optional(name, user_type, description, default=value, choices=choices, letter=None, convert_default=True)

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

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [required, " + str(real_type) + "]: None" + choices_string, file=configfile)

    # Positional optional
    for name in definition.pos_optional:

        # Get properties
        real_type = definition.pos_optional[name].type
        description = definition.pos_optional[name].description
        default = definition.pos_optional[name].default
        choices = definition.pos_optional[name].choices

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [pos_optional, " + str(real_type) + "]: " + str(default) + choices_string, file=configfile)

    # Optional
    for name in definition.optional:

        # Get properties
        real_type = definition.optional[name].type
        description = definition.optional[name].description
        default = definition.optional[name].default
        choices = definition.optional[name].choices
        letter = definition.optional[name].letter

        choices_string = ""
        if isinstance(choices, dict): choices_string = " # choices = " + stringify.stringify(choices.keys())[1]
        elif choices is not None: choices_string = " # choices = " + stringify.stringify(choices)[1]

        print(indent + "# " + description, file=configfile)
        print(indent + name + " [optional, " + str(real_type) + "]: " + str(default) + choices_string, file=configfile)

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

def add_settings_from_dict(config, definition, dictionary):

    """
    This function ...
    :param config:
    :param definition:
    :param dictionary:
    :return:
    """

    dict_that_is_emptied = copy.deepcopy(dictionary)

    # Fixed
    for name in definition.fixed:

        value = definition.fixed[name].value
        config[name] = value

    # Required
    for name in definition.required:

        if not name in dictionary: raise ValueError("The option '" + name + "' is not specified in the configuration dictionary")

        choices = definition.required[name].choices

        # Get the value specified in the dictionary
        #value = dictionary[name]
        value = dict_that_is_emptied.pop(name)

        # TODO: check type?

        # Check with choices
        if choices is not None:
            choices_string = str(choices.keys()) if isinstance(choices, dict) else str(choices)
            if value not in choices: raise ValueError("The value of '" + str(value) + "' for the option '" + name + "' is not valid: choices are: " + choices_string)

        # Set the value
        config[name] = value

    # Positional optional
    for name in definition.pos_optional:

        # Get properties
        default = definition.pos_optional[name].default
        choices = definition.pos_optional[name].choices

        # Check if this option is specified in the dictionary
        if name in dictionary:

            #value = dictionary[name]
            value = dict_that_is_emptied.pop(name)

            # TODO: check type?

            # Check with choices
            if choices is not None:
                choices_string = str(choices.keys()) if isinstance(choices, dict) else str(choices)
                if value not in choices: raise ValueError("The value of '" + str(value) + "' for the option '" + name + "' is not valid: choices are: " + choices_string)

        # Use the default value otherwise
        else: value = default

        # Set the value
        config[name] = value

    # Optional
    for name in definition.optional:

        # Get properties
        default = definition.optional[name].default
        choices = definition.optional[name].choices

        # Check if this option is specified in the dictionary
        if name in dictionary:

            #value = dictionary[name]
            value = dict_that_is_emptied.pop(name)

            # TODO: check type?

            # Check with choices
            if choices is not None:
                choices_string = str(choices.keys()) if isinstance(choices, dict) else str(choices)
                if value not in choices: raise ValueError("The value of '" + str(value) + "' for the option '" + name + "' is not valid: choices are: " + choices_string)

        # Use the default value otherwise
        else: value = default

        # Set the value
        config[name] = value

    # Flags
    for name in definition.flags:

        # Get properties
        #letter = definition.flags[name].letter
        default = definition.flags[name].default  # True or False

        # Check if this option is specified in the dictionary
        if name in dictionary: value = dict_that_is_emptied.pop(name)

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
        add_settings_from_dict(config[name], section_definition, section_dictionary)

    # Add leftover settings
    add_nested_dict_values_to_map(config, dict_that_is_emptied)

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

def add_settings_interactive(config, definition, prompt_optional=True):

    """
    This function ...
    :return:
    """

    # Fixed
    for name in definition.fixed:

        # Get properties
        description = definition.fixed[name].description
        value = definition.fixed[name].value

        # Give name and description
        log.success(name + ": " + description)

        # Inform the user
        log.info("Using fixed value '" + str(value) + "' for " + name)

        # Set the value
        config[name] = value

    # Required
    for name in definition.required:

        # Get properties
        real_type = definition.required[name].type
        description = definition.required[name].description
        choices = definition.required[name].choices
        dynamic_list = definition.required[name].dynamic_list

        # Give name and description
        log.success(name + ": " + description)

        # Get list of choices and a dict of their descriptions
        if choices is not None:
            choices_list = choices.keys() if isinstance(choices, dict) else choices
            choice_descriptions = choices if isinstance(choices, dict) else None
        else: choices_list = choice_descriptions = None

        if choices is not None:

            if real_type.__name__.endswith("_list"): # list-type setting

                log.info("Choose one or more of the following options (separated only by commas)")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = None # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the numbers of the choice
                    answer = raw_input("   : ")
                    try:
                        indices = parsing.integer_list(answer)
                        value = [choices_list[index] for index in indices] # value is a list
                        break
                    except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("Choose one of the following options")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = None  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:

                    # Get the number of the choice
                    answer = raw_input("   : ")
                    try:
                        index = parsing.integer(answer)
                        value = choices_list[index]
                        break
                    except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

        else:

            if real_type.__name__.endswith("_list"):  # list-type setting

                if dynamic_list:

                    real_base_type = getattr(parsing, real_type.__name__.split("_list")[0])

                    log.info("Provide a value for a list item and press ENTER. Leave blank and press ENTER to end the list")

                    value = []
                    while True:
                        answer = raw_input("   : ")
                        if answer == "": break # end of the list
                        try:
                            single_value = real_base_type(answer)
                            value.append(single_value)
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

                else:

                    log.info("Provide the values, seperated by commas")

                    value = []  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input("   : ")
                        try:
                            value = real_type(answer)
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("Provide a value")

                value = None # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    answer = raw_input("   : ")
                    try:
                        value = real_type(answer)
                        break
                    except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

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

        #
        log.info("Press ENTER to use the default value (" + stringify.stringify(default)[1] + ")")

        if choices_list is not None:

            if real_type.__name__.endswith("_list"):  # list-type setting

                log.info("or choose one or more of the following options (separated only by commas)")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = default # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the numbers of the choice
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            indices = parsing.integer_list(answer)
                            value = [choices_list[index] for index in indices] # value is a list
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("or choose one of the following options")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the number of the choice
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            index = parsing.integer(answer)
                            value = choices_list[index]
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

        else:

            if real_type.__name__.endswith("_list"):  # list-type setting

                if dynamic_list:

                    log.info("or provide other values. Enter a value and press ENTER. To end the list, leave blank and press ENTER.")

                    value = [] # to remove warning
                    while True:
                        answer = raw_input("   : ")
                        if answer == "": break # end of list
                        else:
                            try:
                                single_value = real_type(answer)
                                value.append(single_value)
                            except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

                else:

                    log.info("or provide other values, separated by commas")

                    value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input("   : ")
                        if answer == "":
                            value = default
                            break
                        else:
                            try:
                                value = real_type(answer)
                                break
                            except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("or provide another value")

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            value = real_type(answer)
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

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

        #
        log.info("Press ENTER to use the default value (" + stringify.stringify(default)[1] + ")")

        if choices_list is not None:

            if real_type.__name__.endswith("_list"):  # list-type setting

                log.info("or choose one or more of the following options (separated only by commas)")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:

                    # Get the numbers of the choice
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            indices = parsing.integer_list(answer)
                            value = [choices_list[index] for index in indices] # value is a list here
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("or choose one of the following options")

                for index, label in enumerate(choices_list):
                    choice_description = ""
                    if choice_descriptions is not None: choice_description = ": " + choice_descriptions[label]
                    log.info(" - [" + str(index) + "] " + label + choice_description)

                value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                while True:

                    # Get the number of the choice
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            index = parsing.integer(answer)
                            value = choices_list[index] # if we are here, no error was raised
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

        else:

            if real_type.__name__.endswith("_list"):  # list-type setting

                if dynamic_list:

                    log.info("or provide other values. Enter a value and press ENTER. To end the list, leave blank and press ENTER.")

                    value = [] # to remove warning
                    while True:
                        answer = raw_input("   : ")
                        if answer == "": break # end of list
                        else:
                            try:
                                single_value = real_type(answer)
                                value.append(single_value)
                            except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

                else:

                    log.info("or provide other values, separated by commas")

                    value = default  # to remove warning from IDE that value could be referenced (below) without assignment
                    while True:
                        answer = raw_input("   : ")
                        if answer == "":
                            value = default
                            break
                        else:
                            try:
                                value = real_type(answer)
                                break
                            except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

            else:

                log.info("or provide another value")

                value = default # to remove warning from IDE that value could be referenced (below) without assignment
                while True:
                    # Get the input
                    answer = raw_input("   : ")
                    if answer == "":
                        value = default
                        break
                    else:
                        try:
                            value = real_type(answer)
                            break
                        except ValueError, e: log.warning("Invalid input: " + str(e) + ". Try again.")

        # Set the value
        config[name] = value

    # Flags
    for name in definition.flags:

        # Get properties
        description = definition.flags[name].description
        letter = definition.flags[name].letter
        default = definition.flags[name].default  # True or False

        if not prompt_optional:
            value = default
            # Set the value
            config[name] = value
            continue

        # Give name and description
        log.success(name + ": " + description)

        # Ask the question
        log.info("Do you want '" + name + "' to be enabled or not (y or n) or press ENTER for the default (" + str(default) + ")")

        value = default  # to remove warning from IDE that value could be referenced (below) without assignment
        while True:
            answer = raw_input("   : ")
            if answer == "":
                value = default
                break
            else:
                try:
                    value = parsing.boolean(answer) # if this passes without error, we have valid input
                    break
                except ValueError, e:
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

        # Add the settings
        add_settings_interactive(config[name], section_definition, prompt_optional=prompt_optional)

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

            mapping[name] = Map()
            add_nested_dict_values_to_map(mapping[name], value)

        # Else, add the value to the mapping
        else: mapping[name] = value

# -----------------------------------------------------------------
