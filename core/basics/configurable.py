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
import os
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..tools import configuration
from ..tools import filesystem as fs
from ..tools.logging import log

# -----------------------------------------------------------------

class Configurable(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config):

        """
        The constructor ...
        :param config:
        """

        if config is not None: self.config = config

        # Look for the config
        else:

            from ..tools import introspection

            tables = introspection.get_arguments_tables()
            #table_matches = introspection.find_matches_tables(script_name, tables)

            import inspect

            class_name = self.__class__.__name__

            class_path = inspect.getfile(self.__class__).split(".py")[0]

            relative_class_path = class_path.rsplit("pts/")[1]

            relative_class_pts = relative_class_path.replace("/", ".") + "." + class_name

            subproject, relative_class_subproject = relative_class_pts.split(".", 1)

            #print(subproject, relative_class_subproject)

            #exit()

            # Get the correct table
            table = tables[subproject]

            command_name = None
            description = None
            configuration_name = None
            configuration_module_path = None

            #print(table)

            for i in range(len(table["Path"])):

                #print(table["Path"][i], relative_class_subproject)

                if table["Path"][i] == relative_class_subproject:

                    command_name = table["Command"][i]
                    description = table["Description"][i]

                    configuration_name = table["Configuration"][i]
                    if configuration_name == "--": configuration_name = command_name
                    configuration_module_path = "pts." + subproject + ".config." + configuration_name

                    break

            if command_name is not None:

                #print(configuration_module_path)

                import importlib

                # Import things
                #from pts.core.tools import logging
                from pts.core.basics.configuration import ConfigurationDefinition, PassiveConfigurationSetter

                ## GET THE CONFIGURATION DEFINITION
                try:
                    configuration_module = importlib.import_module(configuration_module_path)
                    # has_configuration = True
                    definition = getattr(configuration_module, "definition")
                except ImportError:
                    log.warning("No configuration definition found for the " + class_name + " class")
                    # has_configuration = False
                    definition = ConfigurationDefinition(write_config=False)  # Create new configuration definition

                #print(definition.sections)

                ## CREATE THE CONFIGURATION

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

class OldConfigurable(object):

    """
    This class ...
    """

    def __init__(self, config, subpackage):

        """
        The constructor ...
        :param config:
        :param subpackage:
        :return:
        """

        # Call the constructor of the base class
        super(OldConfigurable, self).__init__()

        # -- Attributes --

        # Set the configuration object
        self.config = configuration.set(subpackage, self.name, config)

        # The children of the object
        self.children = dict()

    # -----------------------------------------------------------------

    def __getattr__(self, attr):

        """
        This function ...
        Overriding __getattr__ should be fine (will not break the default behaviour) -- __getattr__ is only called
        as a last resort i.e. if there are no attributes in the instance that match the name.
        :param name:
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

    def setup(self):
        
        """
        This function ...
        """

        pass

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

    def full_input_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name is None: return None

        if os.path.isabs(name): return name
        elif "input_path" in self.config and self.config.input_path is not None: return os.path.join(self.config.input_path, name)
        else: return os.path.join(os.getcwd(), name)

    # -----------------------------------------------------------------

    def full_output_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if name is None: return None

        if os.path.isabs(name): return name
        elif "output_path" in self.config and self.config.output_path is not None: return os.path.join(self.config.output_path, name)
        else: return os.path.join(os.getcwd(), name)

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        name = type(self).__name__.lower()
        if "plotter" in name: return "plotter"
        else: return name

# -----------------------------------------------------------------
