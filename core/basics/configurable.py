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

# Import the relevant PTS classes and modules
from .loggable import Loggable
from ..tools import configuration

# -----------------------------------------------------------------

class Configurable(Loggable):

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
        super(Configurable, self).__init__()

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

        # Determine the full path to the log file
        path = self.full_output_path(self.config.logging.path) if self.config.logging.path is not None else None

        # Call the setup function of the Loggable base class
        super(Configurable, self).setup(self.config.logging.level, path)

        # Call the setup functions of the children
        for child_name in self.children:

            child = self.children[child_name]

            # Set the input and output path for the child
            child.config.input_path = self.config.input_path
            child.config.output_path = self.config.output_path

            # Options for logging
            child.config.logging.level = "WARNING"
            child.config.logging.path = self.config.logging.path

            # Set the log level and path for the different children of this object, if cascading is enabled
            if self.config.logging.cascade:

                # Galaxy extractor
                child.config.logging.cascade = True
                child.config.logging.level = self.config.logging.level

    # -----------------------------------------------------------------

    def add_child(self, name, type, config=None):

        """
        This function ...
        :param name:
        :param instance:
        :return:
        """

        if name in self.children: raise ValueError("Child with this name already exists")
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

        if os.path.isabs(name): return name
        elif "output_path" in self.config and self.config.output_path is not None: return os.path.join(self.config.output_path, name)
        else: return os.path.join(os.getcwd(), name)

# -----------------------------------------------------------------
