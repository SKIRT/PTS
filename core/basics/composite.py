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
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ..tools import parsing
from ..tools.logging import log
from .configuration import stringify_not_list

# -----------------------------------------------------------------

class SimplePropertyComposite(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

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

        with open(path, 'r') as instrumentfile:

            for line in instrumentfile:

                if "Type:" in line: continue

                name, rest = line.split(": ")
                value, dtype = rest.split("[")
                dtype = dtype.split("]")[0]

                # Set the property value
                if dtype == "None": properties[name] = None
                else: properties[name] = getattr(parsing, dtype)(value)

        # Create the class instance
        return cls(**properties)

    # -----------------------------------------------------------------

    def save(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Inform the user
        log.info("Saving the " + self.__class__.__name__ + " to " + path + " ...")

        # Write the properties
        with open(path, 'w') as instrumentfile:

            # Print the type
            print("Type:", self.__class__.__name__, file=instrumentfile)

            # Loop over the variables
            for name in vars(self):

                dtype, value = stringify_not_list(getattr(self, name))
                print(name + ":", value + " [" + dtype + "]", file=instrumentfile)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

# -----------------------------------------------------------------
