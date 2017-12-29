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
from pts.core.basics.distribution import newDistribution
from pts.core.plot.sed import plot_sed

# -----------------------------------------------------------------

composite = "composite"
table = "table"
dictionary = "dictionary"
sed = "sed"
distribution = "distribution"
filetypes = [composite, table, dictionary, sed, distribution]

# -----------------------------------------------------------------

# Create configuration
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filetype", "string", "type of file", choices=filetypes)
definition.add_required("filename", "file_path", "path of the dictionary file")
definition.add_flag("latex", "print as latex")
definition.add_optional("columns", "string_list", "only show these columns")
definition.add_flag("interactive", "display tables interactively", False)
definition.add_optional("sort", "string", "sort the entries on this column")
definition.add_flag("plot", "make a plot")
definition.add_optional("plot_path", "string", "plot output path")
config = parse_arguments("show_file", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

def load_structure(path, filetype, columns=None):

    """
    This function ...
    :param path:
    :param filetype:
    :param columns:
    :return:
    """

    # Composite
    if filetype == composite:

        # Load composite
        structure = load_composite(path)

        # Create table
        tab = SmartTable.from_composite(composite)

        # Remove columns?
        if columns is not None: tab.remove_other_columns(columns)

    # Table
    elif filetype == table:

        # Load table
        structure = SmartTable.from_file(path)
        tab = structure

        # Remove columns?
        if columns is not None: tab.remove_other_columns(columns)

    # Dictionary
    elif filetype == dictionary:

        # Load dictionary
        structure = load_dict(path)
        tab = SmartTable.from_dictionary(structure)

        # Remove columns?
        if columns is not None: tab.remove_other_columns(columns)

    # SED
    elif filetype == sed:

        # Load SED
        structure = load_sed(path)
        tab = structure

        # Remove columns?
        if columns is not None: tab.remove_other_columns(columns)

    # Distribution
    elif filetype == distribution:

        # Load distribution
        #structure = Distribution.from_file(path)
        #tab = structure.as_table()
        structure = newDistribution.from_old_file(path)
        tab = structure

        # Remove columns?
        #if columns is not None: tab.remove_other_columns(columns) doesn't really make sense for distribution

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
    if filetype == table:
        if config.interactive: structure.more()
        else: print(structure)

    # Dictionary
    elif filetype == dictionary:

        bullet = "-"
        for label in structure:
            line = " " + bullet + " " + fmt.bold + label + fmt.reset + ": " + tostr(structure[label])
            print(line)

    # SED
    elif filetype == sed:
        if config.interactive: structure.more()
        else: print(structure)

    # Distribution
    elif filetype == distribution:
        #if config.interactive: structure.as_table().more()
        #else: print(structure.as_table())
        if config.interactive: structure.more()
        else: print(structure)

    # Not recognized
    else: raise ValueError("Unrecognized filetype")

# -----------------------------------------------------------------

def plot_structure(structure, filetype, filepath=None):

    """
    This function ...
    :param structure:
    :param filetype:
    :param filepath:
    :return:
    """

    # Composite
    if filetype == composite: raise NotImplementedError("Not implemented")

    # Table
    if filetype == table: raise NotImplementedError("Not implemented")

    # Dictionary
    elif filetype == dictionary: raise NotImplementedError("Not implemented")

    # SED
    elif filetype == sed: plot_sed(structure, path=filepath)

    # Distribution
    elif filetype == distribution: structure.plot(path=filepath)

    # Not recognized
    else: raise ValueError("Unrecognized filetype")

# -----------------------------------------------------------------

# Load
structure, tab = load_structure(config.filename, config.filetype, columns=config.columns)

# Sort?
if config.sort is not None: tab.sort(config.sort)

# -----------------------------------------------------------------

# Latex representation
if config.latex:
    if tab is None: raise ValueError("Not supported")
    tab.print_latex()

# Regular representation
else: show_structure(structure, config.filetype)

# -----------------------------------------------------------------

# Plot
if config.plot: plot_structure(structure, config.filetype)

# -----------------------------------------------------------------
