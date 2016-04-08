#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.map Contains the Map class, a dictionary-like object that provides
#  access to its values by using the 'dot' notation.

# Import standard modules
import warnings

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

    def set_items(self, items):

        """
        This function allows setting multiple items in the Map from a dictionary (or dictionary-like)
        :param items:
        :return:
        """

        # Loop over all the items in the 'items' dictionary
        for key in items:

            # Check whether an item with this key exists in this Map
            if key in self:

                # Check if the item is composed of other options (i.e. this is a nested Map), or if it is just a simple variable
                if isinstance(self[key], Map): self[key].set_items(items[key])

                # If it is a simple variable, just set the corresponding item of this Map
                else: self[key] = items[key]

            # If the key is not present, show a warning
            else: warnings.warn("An item with the key '" + key + "' is not present in the Map")

# -----------------------------------------------------------------
