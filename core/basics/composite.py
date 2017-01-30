#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.composite Contains the SimplePropertyComposite class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy
import warnings
from abc import ABCMeta
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools import parsing
from ..tools.logging import log
from ..tools import formatting as fmt
from ..tools.stringify import stringify

# -----------------------------------------------------------------

class SimplePropertyComposite(object):

    """
    This class ...
    @DynamicAttrs
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        """

        # The path
        self._path = None

        # The descriptions
        self._descriptions = dict()

        # The parsing types
        self._ptypes = dict()

        # The choices
        self._choices = dict()

        # The sections
        self._sections = OrderedDict()

    # -----------------------------------------------------------------

    def add_property(self, name, ptype, description, default_value=None, choices=None):

        """
        This function ...
        :param name:
        :param ptype:
        :param description:
        :param default_value:
        :param choices:
        :return:
        """

        # Check
        if hasattr(self, name): raise ValueError("A property with the name '" + name + "' already exists")

        # Set the ptype
        self._ptypes[name] = ptype

        # Set the description
        self._descriptions[name] = description

        # Set the choices
        self._choices[name] = choices

        # Set the attribute with the default value
        setattr(self, name, default_value)

    # -----------------------------------------------------------------

    def add_section(self, name, description):

        """
        This function ...
        :param name:
        :param description:
        :return:
        """

        # Set the description
        self._descriptions[name] = description

        # Set an attribute that is a nested SimplePropertyComposite
        setattr(self, name, SimplePropertyComposite())

    # -----------------------------------------------------------------

    def __setattr__(self, name, value):

        """
        This function ...
        :param name:
        :param value:
        :return:
        """

        if name.startswith("_"):
            super(SimplePropertyComposite, self).__setattr__(name, value)
            return

        if value is None: pass
        elif isinstance(value, SimplePropertyComposite): assert name in self._descriptions
        else:

            # Check the type
            ptype, string = stringify(value)

            # Try converting the string back to the type it actually needs to be
            the_type = self._types[name]
            parsing_function = getattr(parsing, the_type)
            try: value = parsing_function(string)
            except ValueError: raise ValueError("The value given is of the wrong type: '" + ptype + "', must be '" + the_type + "'")

        # Actually set the attribute
        super(SimplePropertyComposite, self).__setattr__(name, value)

    # -----------------------------------------------------------------

    def set_properties(self, properties):

        """
        This function ...
        :param properties:
        :return:
        """

        # Loop over all the options defined in the 'options' dictionary
        for name in properties:

            # Check whether an option with this name exists in this class
            if hasattr(self, name):

                # Check if the option is composed of other options (a Map), or if it is just a simple variable
                #if isinstance(getattr(self, name), Map):
                if isinstance(getattr(self, name), SimplePropertyComposite):

                    assert isinstance(properties[name], dict)  # dict, or Map (is derived from dict)
                    getattr(self, name).set_items(properties[name])

                # If it is a simple variable, just use setattr to set the attribute of this class
                else: setattr(self, name, properties[name])

            # If the option does not exist, ignore it but give a warning
            else: warnings.warn("The property '" + name + "' does not exist")

    # -----------------------------------------------------------------

    @classmethod
    def from_dict(cls, dictionary):

        """
        This function ...
        :param dictionary
        :return:
        """

        # Create a new instance
        composite = cls()

        # Set the properties from the dictionary
        composite.set_properties(dictionary)

        # Return the new property composite
        return composite

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        lines = []

        # Loop over the variables
        for name in vars(self):

            dtype, value = stringify(getattr(self, name))
            line = " - " + fmt.bold +  name + fmt.reset + ": " + value
            lines.append(line)

        return "\n".join(lines)

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading " + cls.__name__ + " from " + path + " ...")

        properties = dict()

        with open(path, 'r') as f:

            for line in f:

                if "Type:" in line: continue

                line = line[:-1]
                if not line: continue

                name, rest = line.split(": ")
                value, dtype = rest.split(" [")
                dtype = dtype.split("]")[0]

                # Set the property value
                if dtype == "None" or value.strip() == "None": properties[name] = None
                else: properties[name] = getattr(parsing, dtype)(value)

        # Create the class instance
        composite = cls(**properties)

        # Set the path
        composite._path = path

        # Return
        return composite

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether the path is valid
        if self._path is None: raise RuntimeError("The path is not defined")

        # Save
        self.saveto(self._path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the " + self.__class__.__name__ + " to " + path + " ...")

        # Write the properties
        with open(path, 'w') as fh:

            # Print the type
            print("Type:", self.__class__.__name__, file=fh)

            # Loop over the variables
            for name in vars(self):

                # Skip internal variables
                if name.startswith("_"): continue

                dtype, value = stringify(getattr(self, name))
                actual_dtype = self._types[name]
                print(name + ":", value + " [" + actual_dtype + "]", file=fh)

        # Update the path
        self._path = None

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------
