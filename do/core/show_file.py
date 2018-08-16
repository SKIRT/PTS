#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.show_file Show a file created from a PTS structure (composite, table, ...).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.plot.sed import plot_sed
from pts.core.plot.distribution import plot_distribution
from pts.core.tools import formatting as fmt
from pts.core.tools import sequences
from pts.core.tools.stringify import tostr
from pts.core.basics.structure import load_structure, filetypes, composite, table, dictionary, sed, distribution, regions

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filetype", "string", "type of file", choices=filetypes)
definition.add_required("filename", "file_path", "path of the dictionary file")

# As table?
definition.add_flag("table", "show as table")

# Displaying and formatting
definition.add_flag("interactive", "display tables interactively", False)
definition.add_flag("latex", "print as latex")

# Modyfying the table
definition.add_optional("columns", "string_list", "only show these columns")
definition.add_optional("sort", "string", "sort the entries on this column")

# Extra options
definition.add_flag("plot", "make a plot")
definition.add_optional("plot_path", "string", "plot output path")
definition.add_optional("plotting", "dictionary", "plotting options", dict())

# Formatting of the values
definition.add_flag("round", "round the table values")
definition.add_optional("ndecimal_places", "positive_integer", "number of decimal places when rounding")

# Additional options
definition.add_flag("unique", "show only the unique in each column")
definition.add_flag("sorted_unique", "sort the unique values")

# -----------------------------------------------------------------

# Create the configuration
config = parse_arguments("show_file", definition, add_logging=False, add_cwd=False)

# -----------------------------------------------------------------

def show_table(tab):

    """
    This function ...
    :param table:
    :return:
    """

    # Only unique values per column
    if config.unique:

        # Loop over the columns
        for colname in tab.column_names:

            # Get the values
            values = tab.get_column_values(colname, add_unit=False)
            unique_values = sequences.unique_values(values, ignore_none=True)
            if config.sorted_unique: unique_values = list(sorted(unique_values))
            nunique = len(unique_values)

            # Show the values
            print(" - " + fmt.bold + colname + fmt.reset_bold + ": " + tostr(unique_values, decimal_places=config.ndecimal_places, round=config.round) + " (" + str(nunique) + ")")

    # Interactive view
    elif config.interactive: tab.more()

    # Regular view
    else: fmt.print_table(tab, ndecimal_places=config.ndecimal_places, round=config.round)

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
    if filetype == table: show_table(structure)

    # Dictionary
    elif filetype == dictionary: fmt.print_dictionary(structure, bullet="-")

    # SED
    elif filetype == sed: show_table(structure)

    # Distribution
    elif filetype == distribution: show_table(structure)

    # Regions
    elif filetype == regions:
        for region in structure: print(region)

    # Not recognized
    else: raise ValueError("Unrecognized filetype")

# -----------------------------------------------------------------

def plot_structure(structure, filetype, filepath=None, **kwargs):

    """
    This function ...
    :param structure:
    :param filetype:
    :param filepath:
    :param kwargs:
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
    elif filetype == distribution: plot_distribution(structure, path=filepath, **kwargs)

    # Not recognized
    else: raise ValueError("Unrecognized filetype")

# -----------------------------------------------------------------

# Load
structure, tab = load_structure(config.filename, config.filetype)

# -----------------------------------------------------------------

# Modify table?
if tab is not None:

    # Remove columns?
    if config.columns is not None:
        if sequences.contains_more(config.columns, tab.column_names): raise ValueError("There are invalid column names: '" + tostr(sequences.get_other(config.columns, tab.column_names)) + "'")
        tab.remove_other_columns(config.columns)

    # Sort?
    if config.sort is not None: tab.sort(config.sort)

# -----------------------------------------------------------------

# Latex representation
if config.latex:
    if tab is None: raise ValueError("Not supported: not a table")
    tab.print_latex(round=config.round, ndecimal_places=config.ndecimal_places)

# Table
elif config.table:
    if tab is None: raise ValueError("Not supported: not a table")
    show_table(tab)

# Regular representation
else: show_structure(structure, config.filetype)

# -----------------------------------------------------------------
