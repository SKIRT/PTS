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
from .configuration import find_command

# -----------------------------------------------------------------

def write_input(input_dict, path):

    """
    This function ...
    :param input_dict:
    :param path:
    :return:
    """

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

    # Write the classes dictionary
    classes_path = fs.join(path, "classes.dat")
    write_dict(classes, classes_path)

# -----------------------------------------------------------------

def load_input(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ..tools.serialization import load_dict

    input_dict = dict()

    # Add the input.dat input
    input_file_path = fs.join(path, "input.dat")
    if fs.is_file(input_file_path):
        remainder = load_dict(input_file_path)
        for key in remainder: input_dict[key] = remainder[key]

    # Load the classes data
    classes_path = fs.join(path, "classes.dat")
    #classes = load_dict(classes_path)

    # TODO: write this!

    # Loop over the files in the directory
    #for filepath, filename in fs.files_in_path(path, exact_not_name="input"):

    # Return the input
    return input_dict

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

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Set configuration
        self.config = self.get_config(*args, **kwargs)

        # Set the detached calculations flag
        self.detached = False

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

        from .configuration import get_config_for_class

        # Get the config
        if unlisted:
            from .map import Map
            assert isinstance(config, Map)
            return config
        else: return get_config_for_class(cls, config, interactive=interactive, cwd=cwd)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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

    def output_path_directory(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Create and return path
        return fs.create_directory_in(self.output_path, name)

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
