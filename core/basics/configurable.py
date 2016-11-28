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
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..tools.logging import log

# -----------------------------------------------------------------

class Configurable(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

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

        # Set the command name
        cls._command_name = command_name

        # Return the command
        return command_name

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # If config is specified
        if config is not None:

            from .configuration import Configuration, ConfigurationDefinition, DictConfigurationSetter

            if isinstance(config, Configuration): self.config = config
            elif isinstance(config, dict):

                # Find the command
                command_name, class_name, configuration_module_path, description = find_command(self.__class__)

                # Get configuration definition
                if command_name is not None:
                    # Get definition
                    definition = get_definition(class_name, configuration_module_path)
                else: definition = ConfigurationDefinition(write_config=False)

                # Create the DictConfigurationSetter
                setter = DictConfigurationSetter(config, command_name, description)

                # Set the configuration
                self.config = setter.run(definition)

            # Not a valid config argument
            else: raise ValueError("Config should be Configuration, dictionary or None")

        # Look for the config
        else:

            # Find the command
            command_name, class_name, configuration_module_path, description = find_command(self.__class__)

            if command_name is not None:

                # Get definition
                definition = get_definition(class_name, configuration_module_path)

                ## CREATE THE CONFIGURATION

                from pts.core.basics.configuration import ConfigurationDefinition, PassiveConfigurationSetter

                # If not specified on the command line (before the command name), then use the default specified in the commands.dat file
                #if configuration_method is None: configuration_method = configuration_method_table

                setter = PassiveConfigurationSetter(class_name, add_logging=False)

                # Create the configuration from the definition
                self.config = setter.run(definition)

                #log.warning("The object has not been configured yet")

            else:

                from .configuration import ConfigurationDefinition
                from .configuration import InteractiveConfigurationSetter

                definition = ConfigurationDefinition(write_config=False)
                setter = InteractiveConfigurationSetter(class_name, add_logging=False)

                # Create new config
                self.config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # nothing is (yet) required here
        pass

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

    def input_path_file(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.input_path, name)

    # -----------------------------------------------------------------

    @property
    def output_path(self):

        """
        This function ...
        :return:
        """

        # If 'output' is defined in the config
        if "output" in self.config:

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

class HierarchicConfigurable(Configurable):

    """
    This class ...
    """

    def __init__(self, config):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(HierarchicConfigurable, self).__init__(config)

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

def find_command(cls):

    """
    This function ...
    :return:
    """

    from ..tools import introspection

    tables = introspection.get_arguments_tables()
    # table_matches = introspection.find_matches_tables(script_name, tables)

    import inspect

    class_name = cls.__name__

    class_path = inspect.getfile(cls).split(".py")[0]

    relative_class_path = class_path.rsplit("pts/")[1]

    relative_class_pts = relative_class_path.replace("/", ".") + "." + class_name

    subproject, relative_class_subproject = relative_class_pts.split(".", 1)

    # print(subproject, relative_class_subproject)

    # exit()

    # Get the correct table
    table = tables[subproject]

    command_name = None
    description = None
    configuration_name = None
    configuration_module_path = None

    # print(table)

    for i in range(len(table["Path"])):

        # print(table["Path"][i], relative_class_subproject)

        if table["Path"][i] == relative_class_subproject:

            command_name = table["Command"][i]
            description = table["Description"][i]

            configuration_name = table["Configuration"][i]
            if configuration_name == "--": configuration_name = command_name
            configuration_module_path = "pts." + subproject + ".config." + configuration_name

            break

    # Return the command name
    return command_name, class_name, configuration_module_path, description

# -----------------------------------------------------------------

def get_definition(class_name, configuration_module_path):

    """
    This function ...
    :return:
    """

    import importlib

    # Import things
    from pts.core.basics.configuration import ConfigurationDefinition

    ## GET THE CONFIGURATION DEFINITION
    try:
        configuration_module = importlib.import_module(configuration_module_path)
        # has_configuration = True
        definition = getattr(configuration_module, "definition")
    except ImportError:
        log.warning("No configuration definition found for the " + class_name + " class")
        # has_configuration = False
        definition = ConfigurationDefinition(write_config=False)  # Create new configuration definition

    return definition

# -----------------------------------------------------------------
