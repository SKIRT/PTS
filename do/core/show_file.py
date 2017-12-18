#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_dictionary Show a file created from a PTS structure (composite, table, ...).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.tools.serialization import load_dict
from pts.core.basics.table import SmartTable
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.core.data.sed import load_sed
from pts.core.basics.composite import load_composite

# -----------------------------------------------------------------

composite = "composite"
table = "table"
dictionary = "dictionary"
sed = "sed"
filetypes = [composite, table, dictionary, sed]

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filetype", "string", "type of file", choices=filetypes)
definition.add_required("filename", "file_path", "path of the dictionary file")
definition.add_flag("latex", "print as latex")
config = parse_arguments("show_file", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

def load_structure(path, filetype):

    """
    This function ...
    :param path:
    :param filetype:
    :return:
    """

    # Composite
    if filetype == composite:

        # Load composite
        structure = load_composite(path)

        # Create table
        tab = SmartTable.from_composite(composite)

    # Table
    elif filetype == table:

        # Load table
        structure = SmartTable.from_file(path)
        tab = structure

    # Dictionary
    elif filetype == dictionary:

        # Load dictionary
        structure = load_dict(path)
        tab = SmartTable.from_dictionary(structure)

    # SED
    elif filetype == sed:

        # Load SED
        structure = load_sed(path)
        tab = structure

    # Invalid
    else: raise ValueError("Unrecognized filetype")

    # Return
    return structure, tab

# -----------------------------------------------------------------

def show_structure(structure, filetype):

    """
    This function ...
    :param structure:
    :param filetype:
    :return:
    """

    # Composite
    if filetype == composite: print(structure)

    # Table
    if filetype == table: print(structure)

    # Dictionary
    elif filetype == dictionary:

        bullet = "-"
        for label in structure:
            line = " " + bullet + " " + fmt.bold + label + fmt.reset + ": " + tostr(structure[label])
            print(line)

    # SED
    elif filetype == sed: print(structure)

    # Not recognized
    else: raise ValueError("Unrecognized filetype")

# -----------------------------------------------------------------

# Load
structure, tab = load_structure(config.filename, config.filetype)

# -----------------------------------------------------------------

# Latex representation
if config.latex:
    if tab is None: raise ValueError("Not supported")
    tab.print_latex()

# Regular representation
else: show_structure(structure, config.filetype)

# -----------------------------------------------------------------
