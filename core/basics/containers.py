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
import copy
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..tools import types
from ..filter.filter import parse_filter, Filter

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

    def to_dictionary(self):

        """
        This function ...
        :return: 
        """

        return copy.copy(self.contents)

    # -----------------------------------------------------------------

    def remove_all(self):

        """
        This function ...
        :return: 
        """

        self.contents = OrderedDict()

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

        # Slice 'list[a:b]'
        if isinstance(index_or_key, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_key.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                key = self.keys[index]
                value = self[index]
                new.contents[key] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_key, tuple):

            new = self.__class__()
            for i_or_k in index_or_key:
                # Get the actual key and value
                if types.is_integer_type(i_or_k): key = self.keys[i_or_k]
                else: key = i_or_k
                value = self[key]
                new.contents[key] = value
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_key):

            # Get the key
            key = self.keys[index_or_key]
            return self.contents[key]

        # Assume it is a proper key
        else: return self.contents[index_or_key]

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_key, value):

        """
        Set an element
        :param index_or_key:
        :param value:
        """

        # Get the name
        if types.is_integer_type(index_or_key): key = self.keys[index_or_key]
        else: key = index_or_key

        # Replace
        self.replace(key, value)

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
        #if cmp is None: cmp_function = None
        #else: cmp_function = lambda x,y: cmp_function(x[1], y[1]) # cmp_function(a,b)
        cmp_function = cmp

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

class FilterBasedList(KeyList):

    """
    This function ...
    """

    @property
    def filters(self):

        """
        THis function ...
        :return: 
        """

        return self.keys

    # -----------------------------------------------------------------

    @property
    def filter_names(self):

        """
        This function ...
        :return: 
        """

        return [str(fltr) for fltr in self.filters]

    # -----------------------------------------------------------------

    def __contains__(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        return fltr in self.keys

    # -----------------------------------------------------------------

    def __getitem__(self, index_or_filter):

        """
        This function ...
        :param index_or_filter:
        :return: 
        """

        # Slice 'list[a:b]'
        if isinstance(index_or_filter, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_filter.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                fltr = self.filters[index]
                value = self[fltr]
                new.contents[fltr] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_filter, tuple):

            # Initialize new
            new = self.__class__()
            for i_or_f in index_or_filter:
                if types.is_integer_type(i_or_f):
                    # Get the key and vlaue
                    fltr = self.filters[i_or_f]
                    value = self[fltr]
                    new.contents[fltr] = value
                elif types.is_string_type(i_or_f):
                    # Convert into filter
                    fltr = parse_filter(i_or_f)
                    value = self[fltr]
                    new.contents[fltr] = value
                elif isinstance(i_or_f, Filter):
                    value = self[i_or_f]
                    new.contents[i_or_f] = value
                else: raise ValueError("Invalid value in tuple: must be index, filter string or filter")
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_filter):

            # Get the key
            fltr = self.filters[index_or_filter]
            return self.contents[fltr]

        # Get the filter
        elif types.is_string_type(index_or_filter):

            # Convert into filter
            fltr = parse_filter(index_or_filter)
            return self.contents[fltr]

        # Filter
        elif isinstance(index_or_filter, Filter): return self.contents[index_or_filter]

        # Invalid index or filter
        else: raise ValueError("Invalid input: must be slice of indices, tuple of indices and/or filter (strings), single index, filter string or filter")

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_filter, value):

        """
        This functio n...
        :param index_or_filter: 
        :param value: 
        :return: 
        """

        # Get the filters
        if types.is_string_type(index_or_filter): index_or_filter = parse_filter(index_or_filter)

        # Call the function of the base class
        self.replace(index_or_filter, value)

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

        # Slice 'list[a:b]'
        if isinstance(index_or_name, slice):

            # Create new class instance
            new = self.__class__()
            indices = index_or_name.indices(len(self))
            for index in range(indices[0], indices[1], indices[2]):
                # Get the key and value
                name = self.names[index]
                value = self[name]
                new.contents[name] = value
            # Return the new list
            return new

        # Tuple 'list[a,b,c]'
        elif isinstance(index_or_name, tuple):

            new = self.__class__()
            for i_or_n in index_or_name:
                if types.is_integer_type(i_or_n): name = self.names[i_or_n]
                elif types.is_string_type(i_or_n): name = i_or_n
                else: raise ValueError("Invalid item in tuple: must be index or name")
                value = self.contents[name]
                new[name] = value
            # Return the new list
            return new

        # Integer type 'list[i]'
        elif types.is_integer_type(index_or_name):

            # Get the name
            name = self.names[index_or_name]
            return self.contents[name]

        # String: name
        elif types.is_string_type(index_or_name):
            name = index_or_name
            return self.contents[name]

        # Invalid
        else: raise ValueError("Key must be slice, tuple, index or name")

    # -----------------------------------------------------------------

    def __setitem__(self, index_or_name, value):

        """
        Set an element
        :param index_or_name:
        :param value:
        """

        # Get the name
        if types.is_string_type(index_or_name): name = index_or_name
        elif types.is_integer_type(index_or_name): name = self.names[index_or_name]
        else: raise ValueError("Invalid index or name: must be integer or string")

        # Replace
        self.replace(name, value)

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
