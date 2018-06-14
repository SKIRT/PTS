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
from collections import OrderedDict

# Import the relevant PTS classes and modules
from . import parsing, stringify
from . import types
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
            #print(indent + "[" + str(index) + "] [list]:", file=fh)
            print(indent + "x [list]:", file=fh)
            write_list_impl(fh, value, indent=indent + "    ")

        else:

            ptype, string = stringify.stringify(data[index])
            #print(indent + "[" + str(index) + "] [" + ptype + "]: " + string, file=fh)
            print(indent + "x [" + ptype + "]: " + string, file=fh)

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

        #name_ptype, name_string = stringify.stringify_not_list(name)
        name_ptype, name_string = stringify.stringify(name)

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

        #elif isinstance(value, tuple):
        #    print(indent + "[" + name_ptype + "] " + name_string + " [tuple]:", file=dictfile)
        #    print(indent + "{", file=dictfile)
        #    write_data_tuple_impl(dictfile, value, indent=indent+"    ")
        #    print(indent + "}", file=dictfile)

        else:

            ptype, string = stringify.stringify(dct[name])
            if ":" in string: string = string.replace(":", "*colon*")
            print(indent + "[" + name_ptype + "] " + name_string + " [" + ptype + "]: " + string, file=dictfile)

        if index != length - 1: print("", file=dictfile)
        index += 1

# -----------------------------------------------------------------

def load_dict(path, ordered=False):

    """
    This function ...
    :param path:
    :param ordered:
    :return:
    """

    # Initialize dictionary
    if ordered: dct = OrderedDict()
    else: dct = dict()

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

        #print("LINE", line)

        # Line actually defines something
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

                # Replace
                string = string.replace("*colon*", ":")

                #print(string)
                #print(list(name))
                #print(list(indent))
                #name = name.split(indent)[1]
                #string = string.split(indent)[1]

                #print("PTYPE", ptype)
                #print("STRING", string)

                if ptype == "None":
                    value = None
                else:

                    try:

                        parsing_function = getattr(parsing, ptype)
                        value = parsing_function(string)

                        #print("VALUE", value)

                    except AttributeError:

                        # Special case: list of lists of something
                        if ptype.endswith("_list_list"):

                            single_parsing_function = getattr(parsing, ptype.split("_list_list")[0])
                            list_parsing_function = getattr(parsing, ptype.split("_list_list")[0] + "_list")
                            result = []
                            for item in string.split(", "): result.append(list_parsing_function(item))
                            value = result

                            #print("LIST_LIST value", value)

                        # List
                        elif ptype.endswith("_list"):

                            base_type = ptype.split("_list")[0]
                            single_parsing_function = getattr(parsing, base_type)
                            result = []
                            for item in string.split(","): result.append(single_parsing_function(item))
                            value = result

                            #print("LIST VALUE", value)

                        # Tuple
                        elif ptype.endswith("_tuple"):

                            single_parsing_function = getattr(parsing, ptype.split("_tuple")[0])
                            result = []
                            for item in string.split(","): result.append(single_parsing_function(item))
                            value = tuple(result)

                            #print("TUPLE VALUE", value)

                        # Don't know what to do
                        else: raise AttributeError("Could not find the type '" + ptype + "'")

                dct[name_value] = value

            # EMPTY LIST IS POSSIBLE, IT WILL HAVE EMPTINESS AS VALUE
            elif _is_empty_list_specification(name_and_specification):

                name_string = name_and_specification.split("] ")[1].split(" [")[0]

                # Get type and value of the key
                name_ptype = name_and_specification.split("[")[1].split("]")[0]
                name_value = getattr(parsing, name_ptype)(name_string)

                # Set empty list
                dct[name_value] = []

            # OTHER CASES SHOULD BE DICTS OR MAPPINGS OF WHICH THE VALUES ARE DEFINED WITHIN A NEW INDENTATION CLAUSE
            else:

                name_ptype = name_and_specification.split("[")[1].split("]")[0]
                name_string = name_and_specification.split("] ")[1].split(" [")[0]

                new_indent = indent + "   "

                name_value = getattr(parsing, name_ptype)(name_string)

                #print(name_and_specification)

                map_or_dict = name_and_specification.split(name_string + " [")[1].split("]")[0]

                #print(dct)

                if map_or_dict == "dict": dct[name_value] = dict()
                elif map_or_dict == "Map": dct[name_value] = Map()
                else: raise ValueError("Don't know how to proceed with '" + map_or_dict + "'")

                load_dict_impl(dictfile, dct[name_value], new_indent)

# -----------------------------------------------------------------

def _is_empty_list_specification(name_and_specification):

    """
    This function ...
    :param name_and_specification:
    :return:
    """

    #name_ptype = name_and_specification.split("[")[1].split("]")[0]
    name_string = name_and_specification.split("] ")[1].split(" [")[0]

    return name_and_specification.split(name_string + " [")[1].split("]")[0] == "list"

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

def contains_other_lists(lst):

    """
    This function ...
    :param lst:
    :return:
    """

    for element in lst:
        if isinstance(element, list): return True
    return False

# -----------------------------------------------------------------

def write_list_impl(listfile, lst, indent=""):

    """
    This function ...
    :param listfile:
    :param lst:
    :param indent:
    :return:
    """

    for element in lst:

        if isinstance(element, list) and contains_other_lists(element):

            print(indent + "[list]:", file=listfile)
            write_list_impl(listfile, element, indent=indent + "    ")

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

        with open(path, 'wb') as out_file:
            pickle.dump(object, out_file, protocol=protocol)

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

    with open(path, 'r') as in_file:
        obj = pickle.load(in_file)
    return obj

# -----------------------------------------------------------------
