#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.map Contains the Map class, a dictionary-like object that provides
#  access to its values by using the 'dot' notation.

from __future__ import print_function

# Import standard modules
import warnings
from collections import OrderedDict

# -----------------------------------------------------------------

class Map(dict):

    """
    With this class you can use the Map object like another dictionary (including json serialization) or with the dot notation.
    Credit: 'epool' (http://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary)
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.iteritems():
                self[k] = v

    # -----------------------------------------------------------------

    def __getattr__(self, attr):
        if attr.startswith("__") and attr.endswith("__"): raise AttributeError
        return self.get(attr)

    # -----------------------------------------------------------------

    def __setattr__(self, key, value):
        if key.startswith("__") and key.endswith("__"): super(Map, self).__setattr__(key, value)
        self.__setitem__(key, value)

    # -----------------------------------------------------------------

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    # -----------------------------------------------------------------

    def __delattr__(self, item):
        self.__delitem__(item)

    # -----------------------------------------------------------------

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]

    # -----------------------------------------------------------------

    def __len__(self):

        count = 0
        for label in self:
            if label.startswith("_"): continue
            count += 1
        return count

    # -----------------------------------------------------------------

    def __eq__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        if len(self) != len(other):
            print("lengths not equal:", self, other)
            return False

        for label in self:

            if label.startswith("_"): continue

            value = self[label]
            other_value = other[label]

            if isinstance(value, list):

                for a in value:
                    if a not in other_value:
                        print("value " + str(a) + " not in list " + str(other_value))
                        return False

            elif value != other_value:

                #try: keys = value.keys() # if mapping-like
                #except AttributeError: print("value " + str(value) + " and " + str(other_value) + " are not equal")
                return False

        return True

    # -----------------------------------------------------------------

    def __ne__(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        return not self.__eq__(other)

    # -----------------------------------------------------------------

    def set_items(self, items):

        """
        This function allows setting multiple items in the Map from a dictionary (or dictionary-like)
        :param items:
        :return:
        """

        # Loop over all the items in the 'items' dictionary
        for key in items:

            if key.startswith("_"): continue

            # Check whether an item with this key exists in this Map
            if key in self:

                # Check if the item is composed of other options (i.e. this is a nested Map), or if it is just a simple variable
                if isinstance(self[key], Map): self[key].set_items(items[key])

                # If it is a simple variable, just set the corresponding item of this Map
                else: self[key] = items[key]

            # If the key is not present, show a warning
            else: warnings.warn("An item with the key '" + key + "' is not present in the Map")

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function makes a recursive shallow copy
        :return:
        """

        # Initialize map
        new = Map()

        # Loop over all items
        for key in self:

            # Copy nested Map
            if isinstance(self[key], Map): new[key] = self[key].copy()

            # Just set the item, don't copy
            else: new[key] = self[key]

        # Return the copy
        return new

# -----------------------------------------------------------------

class ConfigurationTest(Map):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        super(ConfigurationTest, self).__init__(*args, **kwargs)

        # Types and descriptions
        #self.types = dict()
        self.descriptions = dict()

    # -----------------------------------------------------------------

    def set_description(self, key, value):

        """
        This function ...
        :param key:
        :param value:
        :return:
        """

        self.descriptions[key] = value

    # -----------------------------------------------------------------

    def items_with_descriptions(self):

        """
        This function ...
        :return:
        """

        # The list of items
        items = []

        for key in self:

            if key == "descriptions": continue

            value = self[key]
            description = self.descriptions[key] if key in self.descriptions else None

            items.append((key, value, description))

        # Return the items
        return items

    # -----------------------------------------------------------------

    def __getattr__(self, item):

        """
        This function ...
        :param item:
        :return:
        """

        try:
            #print(item in self)
            #value = self.get(item)
            value = self.__getitem__(item)
            #print("heeere", item in self)
            return value
        #except AttributeError:
        except KeyError:
            self[item] = ConfigurationTest()
            #print("here")
            return self[item]

# -----------------------------------------------------------------
