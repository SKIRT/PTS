#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.textfile Contains functions to load and write text in and out files for SKIRT.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..units.parsing import parse_unit as u
from ..tools.parsing import integer_or_real_or_string

# -----------------------------------------------------------------

def get_title(filepath):

    """
    This function ...
    :param filepath: 
    :return: 
    """

    title = None

    for line in fs.read_lines(filepath):

        # We are no longer at the header
        if not line.startswith("#"): break

        # We are not longer at the title
        if line.startswith("# column"): break
        if ":" in line: break

        line = line.split("#")[1].strip()

        if title is None: title = line
        else: title += "\n" + line

    # Return the title
    return title

# -----------------------------------------------------------------

def get_meta(filepath):

    """
    This function ...
    :param filepath: 
    :return: 
    """

    meta = dict()

    # Read the SED file header
    for line in fs.read_lines(filepath):

        # We are no longer at the header
        if not line.startswith("#"): break

        # Not a line saying something about columns (meta information)
        if line.startswith("# column"): continue

        if ":" not in line: continue

        first, second = line.split(":")

        second = integer_or_real_or_string(second)

        name = first.split("#")[1].strip()

        # Add to the meta information
        meta[name] = second

    # Return the meta information
    return meta

# -----------------------------------------------------------------

def get_descriptions_and_units(filepath, remote=None):

    """
    This function ...
    :param filepath:
    :param remote:
    :return: 
    """

    descriptions = []
    units = []

    if remote is not None: lines_generator = remote.read_lines(filepath)
    else: lines_generator = fs.read_lines(filepath)

    # Read the SED file header
    for line in lines_generator:

        # We are no longer at the header
        if not line.startswith("#"): break

        # Not a line saying something about columns (meta information)
        if not line.startswith("# column"): continue

        # Split the line to get the column index and the unit
        first, second = line.split(":")

        # Get basic column info
        index = int(first.split("column ")[1]) - 1
        description = second.split(";")[0].split("(")[0]

        assert index == len(units)

        # Add description
        descriptions.append(description)

        # Get unit info
        if second.endswith(")") and "(" in second:

            mathematical = second.split(";")[1].split("(")[0].strip() if len(second.split(";")) > 1 else None
            unit_string = second.split("(")[1].split(")")[0]

            # Check whether spectral density unit
            if mathematical is not None:
                density = mathematical.startswith("lambda*") or mathematical.startswith("nu*") if mathematical is not None else False
                density_strict = True
            else:
                density=False
                density_strict = False

            # Parse the unit
            unit = u(unit_string, density=density, density_strict=density_strict)

        # No unit
        else: unit = None

        # Add the unit
        units.append(unit)

    # Garbage-collect the generator (to close the file handles)
    del lines_generator

    # Return
    return descriptions, units

# -----------------------------------------------------------------

def get_units(filepath, remote=None):

    """
    This function ...
    :param filepath:
    :param remote:
    :return: 
    """

    descriptions, units = get_descriptions_and_units(filepath, remote=remote)
    return units

# -----------------------------------------------------------------
