#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

"""
With the class in this module you can use the Map object like another dictionary(including json serialization) or with the dot notation.
Credit: 'epool' (http://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary)
"""

# -----------------------------------------------------------------

class Map(dict):
    
    """
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
        return self.get(attr)

    # -----------------------------------------------------------------

    def __setattr__(self, key, value):
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