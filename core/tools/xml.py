#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.tools.xml Provides useful functions for dealing with XML files (SkiFile, SmileSchema, ...)

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

def get_unique_element(element, name):

    """
    This function returns the xml tree element with the specified name that is a child of the specified element
    :param element:
    :param name:
    :return:
    """

    # Get child element of the given element
    parents = element.xpath(name)

    # Check if only one child element is present
    if len(parents) == 0: raise ValueError("Invalid ski file: no '" + name + "' elements within '" + element.tag + "'")
    elif len(parents) > 1: raise ValueError("Invalid ski file: multiple '" + name + "' elements within '" + element.tag + "'")
    parent = parents[0]

    # Get the children of the parent
    children = parent.getchildren()

    # Check if only one child object is present
    if len(children) == 0: raise ValueError("Invalid ski file: no '" + name + "' elements within '" + element.tag + "'")
    elif len(children) > 1: raise ValueError("Invalid ski file: multiple '" + name + "' elements within '" + element.tag + "'")
    child = children[0]

    # Return the child element
    return child

# -----------------------------------------------------------------

def get_unique_element_direct(element, name):

    """
    This function ...
    :param element:
    :param name:
    :return:
    """

    # Get child element of the given element
    children = element.xpath(name)

    # Check if only one child element is present
    if len(children) == 0: raise ValueError("Invalid ski file: no '" + name + "' elements within '" + element.tag + "'")
    elif len(children) > 1: raise ValueError("Invalid ski file: multiple '" + name + "' elements within '" + element.tag + "'")

    # Return the unique child
    return children[0]

# -----------------------------------------------------------------

def get_list(element, name):

    """
    This function ...
    :param element:
    :param name:
    :return:
    """

    # Get the private properties
    properties = element.xpath(name)

    if len(properties) == 0: return []
    elif len(properties) == 1: return properties[0].getchildren()
    else: raise ValueError("Ambigious result: multiple children with name '" + name + "' for element '" + element + "'")

# -----------------------------------------------------------------

def recursive_dict(element):

    """
    This function ...
    :param element:
    :return:
    """

    return element.tag, dict(map(recursive_dict, element)) or element.text

# -----------------------------------------------------------------

def add_properties(element, dictionary):

    """
    This function ...
    :param element:
    :param dictionary:
    :return:
    """

    for key, value in element.items(): dictionary[key] = value

# -----------------------------------------------------------------

def add_children(element, dictionary):

    """
    This function ...
    :param element:
    :param dictionary:
    :return:
    """

    dictionary["children"] = dict()

    for child in element.getchildren():

        dictionary["children"][child.tag] = dict()

        add_properties(child, dictionary["children"][child.tag])
        add_children(child, dictionary["children"][child.tag])

# -----------------------------------------------------------------

def get_properties(element):

    """
    This function ...
    :param element:
    :return:
    """

    properties = dict()
    add_properties(element, properties)
    add_children(element, properties)
    return properties

# -----------------------------------------------------------------

def get_all_elements(root):

    """
    This function ...
    :param root:
    :return:
    """

    elements = []
    add_all_elements(root, elements)
    return elements

# -----------------------------------------------------------------

def add_all_elements(root, elements):

    """
    This function ...
    :param root:
    :param elements:
    :return:
    """

    for child in root.getchildren():
        elements.append(child)
        add_all_elements(child, elements)

# -----------------------------------------------------------------
