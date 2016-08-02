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
        #else:

        else:

            from .configuration import ConfigurationDefinition
            from .configuration import InteractiveConfigurationSetter

            definition = ConfigurationDefinition()
            setter = InteractiveConfigurationSetter(self.class_name, add_logging=False)

            # Create new config
            self.config = setter.run(definition, prompt_optional=False)

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
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
