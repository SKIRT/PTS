#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.tables Contains the ModelsTable and RepresentationsTable classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import arrays, tables

# -----------------------------------------------------------------

class RepresentationsTable(SmartTable):
        
    """
    This class ..."""

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RepresentationsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Name", str, None, "name for the representation")
        self.add_column_info("Model name", str, None, "name of the model")

    # -----------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return arrays.array_as_list(self["Name"])

    # -----------------------------------------------------------------

    def representations_for_model(self, model_name):

        """
        This function ...
        :param model_name:
        :return:
        """

        # Get the indices
        indices = tables.find_indices(self, model_name, column_name="Model name")

        # Get the representations
        representations = [self["Name"][index] for index in indices]

        # Return the representations
        return representations

    # -----------------------------------------------------------------

    def add_entry(self, name, model_name):

        """
        This function ...
        :param name:
        :param model_name:
        :return:
        """

        # Add row
        values = [name, model_name]
        self.add_row(values)

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

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return arrays.array_as_list(self["Name"])

    # -----------------------------------------------------------------

    def description_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Description"][index] if not self["Description"].mask[index] else None

    # -----------------------------------------------------------------

    def old_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Old stars path"][index]

    # -----------------------------------------------------------------

    def young_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Young stars path"][index]

    # -----------------------------------------------------------------

    def ionizing_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Ionizing stars path"][index]

    # -----------------------------------------------------------------

    def dust_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Dust path"][index]

# -----------------------------------------------------------------
