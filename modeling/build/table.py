#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.table Contains the ModelsTable class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable

# -----------------------------------------------------------------

class ModelsTable(SmartTable):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ModelsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Name", str, None, "name for the model")
        self.add_column_info("Description", str, None, "description of the model")
        self.add_column_info("Old stars path", str, None, "source for the map of old stars")
        self.add_column_info("Young stars path", str, None, "source for the map of young stars")
        self.add_column_info("Ionizing stars path", str, None, "source for the map of ionizing stars")
        self.add_column_info("Dust path", str, None, "source for the map of the dust distribution")

    # -----------------------------------------------------------------

    def add_model(self, name, description, old_stars_path, young_stars_path, ionizing_stars_path, dust_path):

        """
        This function ...
        :param name:
        :param description:
        :param old_stars_path:
        :param young_stars_path:
        :param ionizing_stars_path:
        :param dust_path:
        :return:
        """

        # Set the values
        values = [name, description, old_stars_path, young_stars_path, ionizing_stars_path, dust_path]

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------
