#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.filesystem Provides useful functions for manipulating the local file system.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import json
import pickle

# -----------------------------------------------------------------

def to_json(object):

    """
    This function ...
    :param self:
    :return:
    """

    return json.dumps(object, default=lambda o: o.__dict__, sort_keys=True, indent=4)

# -----------------------------------------------------------------

def dump(object, path, method="pickle"):

    """
    This function ...
    :param object:
    :param path:
    :param method:
    :return:
    """

    # Serialize using pickle
    if method == "pickle":

        pickle.dump(object, open(path, 'wb'))

    # Serialize to the json format
    elif method == "json":

        with open(path, 'w') as out_file:
            json.dump(object, out_file, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    # Not a valid serialization method
    else: raise ValueError("Not a valid method")

# -----------------------------------------------------------------

def load(path):

    """
    This function ...
    :param path:
    :return:
    """

    return pickle.load(open(path, 'r'))

# -----------------------------------------------------------------
