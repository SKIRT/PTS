#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.containers Contains container classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools import types

# -----------------------------------------------------------------

class KeyList(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.contents = OrderedDict()

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return: 
        """

        return len(self.contents)

    # -----------------------------------------------------------------

    def __iter__(self):

        """
        THis function ...
        :return: 
        """

        for name in self.contents: yield self.contents[name]

    # -----------------------------------------------------------------

    @property
    def keys(self):

        """
        This function ...
        :return: 
        """

        return self.contents.keys()

    # -----------------------------------------------------------------

    @property
    def values(self):

        """
        This function ...
        :return: 
        """

        return self.contents.values()

    # -----------------------------------------------------------------

    @property
    def items(self):

        """
        This function ...
        :return: 
        """

        return self.contents.items()

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_key):

        """
        This function ...
        :param index_or_key:
        :return: 
        """

        # Get the key
        if types.is_integer_type(index_or_key): key = self.keys[index_or_key]
        else: key = index_or_key

        # Return the element
        return self.contents[key]

    # -----------------------------------------------------------------

    def append(self, key, element):

        """
        This function ...
        :param key: 
        :param element: 
        :return: 
        """

        # Check the key
        if key in self.keys: raise ValueError("Already an element with the key '" + str(key) + "'")
        self.contents[key] = element

    # -----------------------------------------------------------------

    def replace(self, key, element):

        """
        This function ...
        :param key: 
        :param element: 
        :return: 
        """

        # Check the key
        if key not in self.keys: raise ValueError("The element with key '" + str(key) + "' does not exist")

        # Get the old value
        old = self.contents[key]

        # Replace
        self.contents[key] = element

        # Return the old value
        return old

    # -----------------------------------------------------------------

    def __contains__(self, key):

        """
        This function ...
        :param key: 
        :return: 
        """

        return key in self.keys

    # -----------------------------------------------------------------

    def sort(self, cmp=None, key=None, reverse=False):

        """
        This function ...
        :param cmp: 
        :param key: 
        :param reverse: 
        :return: 
        """

        # Create key function
        if key is None: key_function = lambda x: x[1]
        else: key_function = lambda x: key(x[1])

        # Create cmp function
        if cmp is None: cmp_function = None
        else: cmp_function = lambda x,y: cmp_function(x[1], y[1]) # cmp_function(a,b)

        # Sort the contents
        self.contents = OrderedDict(sorted(self.contents.iteritems(), cmp=cmp_function, key=key_function, reverse=reverse))

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return: 
        """

        # Create new
        new = self.__class__()

        # Add elements from this into new
        for key, element in self.items: new.append(key, element)

        # Return the new copy
        return new

# -----------------------------------------------------------------

class NamedList(KeyList):

    """
    This class ...
    """

    @property
    def names(self):

        """
        This function ...
        :return: 
        """

        return self.contents.keys()

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_name):

        """
        THis function ...
        :param index_or_name: 
        :return: 
        """

        # Get the filter
        if types.is_string_type(index_or_name): name = index_or_name
        elif types.is_integer_type(index_or_name): name = self.names[index_or_name]
        else: raise ValueError("Invalid index or name: must be integer or string")

        # Return the element
        return self.contents[name]

    # -----------------------------------------------------------------

    def append(self, name, element):

        """
        This function ...
        :param name: 
        :param element: 
        :return: 
        """

        if name in self.names: raise ValueError("Already an element with the name '" + name + "'")
        self.contents[name] = element

    # -----------------------------------------------------------------

    def replace(self, name, element):

        """
        This fucntion ...
        :param name: 
        :param element: 
        :return: 
        """

        # Check the name
        if name not in self.names: raise ValueError("The element with name '" + name + "' does not exist")

        # Get the old value
        old = self.contents[name]

        # Replace
        self.contents[name] = element

        # Return the old value
        return old

# -----------------------------------------------------------------
