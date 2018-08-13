#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.plot_file Plot a file created from a PTS structure (composite, table, ...).

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.plot.sed import plot_sed
from pts.core.plot.distribution import plot_distribution
from pts.core.basics.structure import load_structure, filetypes, composite, table, dictionary, sed, distribution, regions
from pts.core.tools import arrays
from pts.core.basics.distribution import Distribution

# -----------------------------------------------------------------

# Create configuration definition
definition = ConfigurationDefinition(write_config=False)
definition.add_required("filetype", "string", "type of file", choices=filetypes)
definition.add_required("filename", "file_path", "path of the dictionary file")

# Plot to path
definition.add_optional("path", "string", "plot output path")
#definition.add_optional("plotting", "dictionary", "plotting options", dict())

# Column
definition.add_optional("column", "string", "column name")
definition.add_flag("logarithmic", "logarithmic")

# -----------------------------------------------------------------

# Create the configuration
config = parse_arguments("plot_file", definition, add_logging=False, add_cwd=False)

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
    if filetype == table: #raise NotImplementedError("Not implemented")

        unit = structure.column_unit(config.column)
        column = arrays.plain_array(structure[config.column], unit=unit, array_unit=unit)
        #print(column)
        column = column[np.isfinite(column)]
        column = column[column != 0]
        distr = Distribution.from_values(config.column, column, logarithmic=config.logarithmic)
        plot_distribution(distr, logscale=config.logarithmic)

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
structure, tab = load_structure(config.filename, config.filetype, table_method="pandas")

# -----------------------------------------------------------------

plot_structure(structure, config.filetype)

# -----------------------------------------------------------------
