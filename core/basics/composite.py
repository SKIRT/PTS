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
from ..tools import formatting as fmt
from ..tools.stringify import stringify

# -----------------------------------------------------------------

class SimplePropertyComposite(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self):

        """
        This function ...
        """

        # The path
        self._path = None

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
                if dtype == "None": properties[name] = None
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
        with open(path, 'w') as instrumentfile:

            # Print the type
            print("Type:", self.__class__.__name__, file=instrumentfile)

            # Loop over the variables
            for name in vars(self):

                # Skip internal variables
                if name.startswith("_"): continue

                dtype, value = stringify(getattr(self, name))
                print(name + ":", value + " [" + dtype + "]", file=instrumentfile)

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
