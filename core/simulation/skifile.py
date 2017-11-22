#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.skifile Reading and updating a SKIRT parameter file.
#
# An instance of the SkiFile class in this module allows reading from and updating an existing ski file.

# -----------------------------------------------------------------

# Import standard modules
from lxml import etree
from ..tools import filesystem as fs
from ..tools import types
from .skifile7 import SkiFile7
from .skifile8 import SkiFile8
from .skifile7 import LabeledSkiFile7
from .skifile8 import LabeledSkiFile8
from ..tools.introspection import skirt_main_version, has_skirt

# -----------------------------------------------------------------

if not has_skirt(): version_number = 8
else: version_number = skirt_main_version()

# -----------------------------------------------------------------

if version_number == 8:
    SkiFile = SkiFile8
    LabeledSkiFile = LabeledSkiFile8
elif version_number == 7:
    SkiFile = SkiFile7
    LabeledSkiFile = LabeledSkiFile7
else: raise ValueError("Invalid SKIRT version: " + str(version_number))

# -----------------------------------------------------------------

def fix_ski_file(path):

    """
    This function ...
    :param path:
    :return:
    """

    from ...core.prep.smile import expected_types

    # Set replacement dictionary
    replacements = dict()
    for a, b in expected_types.items():
        from_string = a + ' type=""'
        to_string = a + ' type="' + b + '"'
        replacements[from_string] = to_string

    # Replace lines
    replace_ski_file_lines(path, replacements)

# -----------------------------------------------------------------

def replace_ski_file_lines(path, replacement_dict):

    """
    This function ...
    :param path:
    :param replacement_dict:
    :return:
    """

    new_lines = []

    which_system = None

    # Read the lines
    for line in fs.read_lines(path):

        if "<dustSystem" in line: which_system = "dust"
        elif "/dustSystem" in line: which_system = None

        if "<stellarSystem" in line: which_system = "stellar"
        elif "/stellarSystem" in line: which_system = None

        # Loop over the replacements
        for from_string in replacement_dict:
            to_string = replacement_dict[from_string]

            # Determine the new line
            if from_string in line: line = line.replace(from_string, to_string)

            if 'components type=""' in line:
                if which_system == "dust": line = line.replace('components type=""', 'components type="DustComp"')
                elif which_system == "stellar": line = line.replace('components type=""', 'components type="StellarComp"')
                else: raise RuntimeError("Something went wrong")

        # Add the line
        new_lines.append(line)

    # Remove the file
    fs.remove_file(path)

    # Write lines
    fs.write_lines(path, new_lines)

# -----------------------------------------------------------------

def get_element_path(element):

    """
    This function ...
    :param element:
    :return:
    """

    from ..tools import sequences

    name = element.tag

    parent = element.getparent()
    if parent is not None:

        # Get the children of the parent
        children = parent.getchildren()

        # The parent has multiple children
        nchildren = len(children)
        if nchildren > 1:

            # Get the child labels (names)
            labels = [child.tag for child in children if child.tag is not etree.Comment]
            nactual_children = len(labels)

            # Loop over the components (also includes the comments)
            number_of_components = 0
            i = 0
            ids = []
            while i < len(children):
                if children[i].tag is etree.Comment:
                    ids.append(children[i].text.strip())
                    i += 2  # skip the next component -> it is the component corresponding to this comment
                # No name -> add the index of this component as the ID
                else:
                    ids.append(number_of_components)
                    i += 1
                # Increment the number of components
                number_of_components += 1

            #print(ids)

            # Check the index of this element
            index = 0
            actual_index = 0
            while True:
                if children[actual_index].tag is etree.Comment:
                    actual_index += 1
                    continue
                if children[actual_index] == element: break
                index += 1
                actual_index += 1
                if actual_index == nchildren:
                    index = None
                    break
            if index is None: raise RuntimeError("Something went wrong")

            #print(labels)
            #print(index)

            # Get the label and the ID
            element_label = labels[index]
            element_id = ids[index]

            if types.is_integer_type(element_id):
                if element_id != index: raise RuntimeError("Something went wrong")
                identifier = str(index)
            elif types.is_string_type(element_id): identifier = "'" + element_id + "'"
            else: raise RuntimeError("Something went wrong")

            if sequences.is_unique(labels, element_label): path = get_element_path(parent) + "/" + name
            else: path = get_element_path(parent) + "/[" + identifier + "]"

        # No multiple children
        else: path = get_element_path(parent) + "/" + name

    # No parent
    else: path = name

    # Return the total path
    return path

# -----------------------------------------------------------------

def is_labeled(value):

    """
    This function ...
    :param value:
    :return:
    """

    return value.startswith("[") and value.endswith("]")

# -----------------------------------------------------------------

def get_label(value):

    """
    This function ...
    :param value:
    :return:
    """

    label = value.split("[")[1].split(":")[0]
    return label

#-----------------------------------------------------------------
