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

# Import the relevant PTS classes and modules
from . import parsing, stringify
from ..basics.map import Map

# -----------------------------------------------------------------

def write_data_tuple(data, path):

    """
    This function ...
    :param data:
    :param path:
    :return:
    """

    with open(path, 'w') as fh: write_data_tuple_impl(fh, data)

# -----------------------------------------------------------------

def write_data_tuple_impl(fh, data, indent=""):

    """
    This function ...
    :param fh:
    :param data:
    :param indent:
    :return:
    """

    for index in range(len(data)):

        value = data[index]

        if isinstance(value, list):

            #print(indent + "[0] [" + ptype + "]: " + string, file=fh)
            print(indent + "[" + str(index) + "] [list]:", file=fh)
            write_list_impl(fh, value, indent=indent + "    ", use_serialization=True)

        else:

            ptype, string = stringify.stringify(data[index])
            print(indent + "[" + str(index) + "] [" + ptype + "]: " + string, file=fh)

# -----------------------------------------------------------------

def write_dict(dct, path):

    """
    This function ...
    :param dct:
    :param path:
    :return:
    """

    with open(path, 'w') as fh: write_dict_impl(fh, dct)

# -----------------------------------------------------------------

def write_dict_impl(dictfile, dct, indent=""):

    """
    This function ...
    :param dictfile:
    :param dct:
    :param indent:
    :return:
    """

    index = 0
    length = len(dct)
    for name in dct:

        name_ptype, name_string = stringify.stringify_not_list(name)

        value = dct[name]

        if isinstance(value, Map):

            print(indent + "[" + name_ptype + "] " + name_string + " [Map]:", file=dictfile)
            print(indent + "{", file=dictfile)
            write_dict_impl(dictfile, value, indent=indent + "    ")
            print(indent + "}", file=dictfile)

        elif isinstance(value, dict):

            print(indent + "[" + name_ptype + "] " + name_string + " [dict]:", file=dictfile)
            print(indent + "{", file=dictfile)
            write_dict_impl(dictfile, value, indent=indent+"    ")
            print(indent + "}", file=dictfile)

        elif isinstance(value, tuple):

            print(indent + "[" + name_ptype + "] " + name_string + " [tuple]:", file=dictfile)
            print(indent + "{", file=dictfile)
            write_data_tuple_impl(dictfile, value, indent=indent+"    ")
            print(indent + "}", file=dictfile)

        else:

            ptype, string = stringify.stringify(dct[name])
            print(indent + "[" + name_ptype + "] " + name_string + " [" + ptype + "]: " + string, file=dictfile)

        if index != length - 1: print("", file=dictfile)
        index += 1

# -----------------------------------------------------------------

def load_dict(path):

    """
    This function ...
    :param path:
    :return:
    """

    dct = dict()
    with open(path, 'r') as fh: load_dict_impl(fh, dct)
    return dct

# -----------------------------------------------------------------

def load_dict_impl(dictfile, dct, indent=""):

    """
    This function ...
    :param dictfile:
    :param dct:
    :param indent:
    :return:
    """

    # Loop over the lines in the file
    for line in dictfile:

        # Strip end-of-line character
        line = line.rstrip("\n")

        # Empty line
        if not line: continue

        #end = indent + "}"
        #print(list(line))
        #print(list(line))
        #print(list(end))

        #if line.startswith(end): return

        if line.strip() == "}": return

        if ":" in line:

            name_and_specification = line.split(":")[0].strip()

            #stripped = name_and_specification.strip()

            #if name_and_specification.endswith("]"):
            if line.split(":")[1].strip() != "":

            #if "[" in name_and_specification and "]" in name_and_specification:

                name_ptype = name_and_specification.split("]")[0][1:]
                name_string = name_and_specification.split("] ")[1].split(" [")[0]
                name_value = getattr(parsing, name_ptype)(name_string)

                ptype = line.split(name_string + " [")[1].split("]")[0]
                string = line.split(":")[1].strip()

                #print(string)
                #print(list(name))
                #print(list(indent))
                #name = name.split(indent)[1]
                #string = string.split(indent)[1]

                parsing_function = getattr(parsing, ptype)
                value = parsing_function(string)

                dct[name_value] = value

            else:

                name_ptype = name_and_specification.split("[")[1].split("]")[0]
                name_string = name_and_specification.split("] ")[1].split(" [")[0]

                new_indent = indent + "   "

                name_value = getattr(parsing, name_ptype)(name_string)

                #print(name_and_specification)

                map_or_dict = name_and_specification.split(name_string + " [")[1].split("]")[0]

                if map_or_dict == "dict": dct[name_value] = dict()
                elif map_or_dict == "Map": dct[name_value] = Map()

                load_dict_impl(dictfile, dct[name_value], new_indent)

# -----------------------------------------------------------------

def write_list(lst, path):

    """
    This function ...
    :param lst:
    :param path:
    :return:
    """

    with open(path, 'w') as fh: write_list_impl(fh, lst)

# -----------------------------------------------------------------

def write_list_impl(listfile, lst, indent="", use_serialization=False):

    """
    This function ...
    :param listfile:
    :param lst:
    :param indent:
    :return:
    """

    for element in lst:

        if isinstance(element, list) and use_serialization:
            write_list_impl(indent + "[list]:")
            write_list_impl(listfile, element, indent=indent+"    ", use_serialization=True)
        else:
            ptype, string = stringify.stringify(element)
            listfile.write(indent + "[" + ptype + "] " + string + "\n")

# -----------------------------------------------------------------

def load_list(path):

    """
    This function ...
    :param path:
    :return:
    """

    lst = []
    with open(path, 'r') as fh: load_list_impl(fh, lst)
    return lst

# -----------------------------------------------------------------

def load_list_impl(listfile, lst):

    """
    This function ...
    :param listfile:
    :param lst:
    :return:
    """

    for line in listfile:

        line = line[:-1]
        if not line: continue

        ptype = line.split("]")[0].split("[")[1]
        string = line.split("]")[1].strip()

        parsing_function = getattr(parsing, ptype)
        value = parsing_function(string)

        lst.append(value)

# -----------------------------------------------------------------

def to_json(object):

    """
    This function ...
    :param self:
    :return:
    """

    return json.dumps(object, default=lambda o: o.__dict__, sort_keys=True, indent=4)

# -----------------------------------------------------------------

def dump(object, path, method="pickle", protocol=None):

    """
    This function ...
    :param object:
    :param path:
    :param method:
    :param protocol:
    :return:
    """

    # Serialize using pickle
    if method == "pickle":

        pickle.dump(object, open(path, 'wb'), protocol=protocol)

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
