#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.basics.structure

# -----------------------------------------------------------------

# Import the relevant PTS classes and modules
from pts.core.tools.serialization import load_dict
from pts.core.basics.table import SmartTable
from pts.core.data.sed import load_sed
from pts.core.basics.composite import load_composite
from pts.core.basics.distribution import Distribution
from pts.magic.region.list import load_region_list

# -----------------------------------------------------------------

composite = "composite"
table = "table"
dictionary = "dictionary"
sed = "sed"
distribution = "distribution"
regions = "regions"
filetypes = [composite, table, dictionary, sed, distribution, regions]

# -----------------------------------------------------------------

def load_structure(path, filetype, table_method="lines"):

    """
    This function ...
    :param path:
    :param filetype:
    :param table_method:
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
        structure = SmartTable.from_file(path, method=table_method)
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

    # Distribution
    elif filetype == distribution:

        # Load distribution
        structure = Distribution.from_file(path)
        tab = structure

    # Regions
    elif filetype == regions:

        # Load regions
        structure = load_region_list(path)
        dictionaries = [region.__dict__ for region in structure]
        tab = SmartTable.from_dictionaries(*dictionaries)

    # Invalid
    else: raise ValueError("Unrecognized filetype")

    # Return
    return structure, tab

# -----------------------------------------------------------------
