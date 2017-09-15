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

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.tools import arrays, tables
from ...core.tools import filesystem as fs

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

    def model_for_representation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = tables.find_index(self, name)
        return self["Model name"][index]

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
        self.add_column_info("Bulge path", str, None, "source for the stellar bulge component")
        self.add_column_info("Old stars path", str, None, "source for old stellar disk component")
        self.add_column_info("Young stars path", str, None, "source for young stellar disk component")
        self.add_column_info("Ionizing stars path", str, None, "source for ionizing stellar disk component")
        self.add_column_info("Additional stars paths", str, None, "source(s) for additional stellar component(s)")
        self.add_column_info("Dust path", str, None, "source for dust disk component")
        self.add_column_info("Additional dust paths", str, None, "source(s) for additional dust component(s)")

    # -----------------------------------------------------------------

    def add_model(self, name, description, bulge_path, old_stars_path, young_stars_path, ionizing_stars_path, dust_path,
                  additional_stellar_paths=None, additional_dust_paths=None):

        """
        This function ...
        :param name:
        :param description:
        :param bulge_path:
        :param old_stars_path:
        :param young_stars_path:
        :param ionizing_stars_path:
        :param dust_path:
        :param additional_stellar_paths:
        :param additional_dust_paths:
        :return:
        """

        if additional_stellar_paths is not None: additional_stars_string = ",".join(additional_stellar_paths)
        else: additional_stars_string = None
        if additional_dust_paths is not None: additional_dust_string = ",".join(additional_dust_paths)
        else: additional_dust_string = None

        # Set the values
        values = [name, description, bulge_path, old_stars_path, young_stars_path, ionizing_stars_path,
                  additional_stars_string, dust_path, additional_dust_string]

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

    def index_for_model(self, name):

        """
        This funtion ...
        :param name:
        :return:
        """

        return tables.find_index(self, name)

    # -----------------------------------------------------------------

    def description_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Description"][index] if not self["Description"].mask[index] else None

    # -----------------------------------------------------------------

    def stellar_component_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # INitialize dictionary for the paths
        paths = OrderedDict()

        # Get all component paths
        bulge_path = self.bulge_path_for_model(name)
        old_path = self.old_stars_path_for_model(name)
        young_path = self.young_stars_path_for_model(name)
        ionizing_path = self.ionizing_stars_path_for_model(name)
        additional_paths = self.additional_stars_paths_for_model(name)

        # Determine standard component names
        bulge_name = fs.name(bulge_path)
        old_name = fs.name(old_path)
        young_name = fs.name(young_path)
        ionizing_name = fs.name(ionizing_path)

        # Set standard paths
        paths[bulge_name] = bulge_path
        paths[old_name] = old_path
        paths[young_name] = young_path
        paths[ionizing_name] = ionizing_path

        # Set additional paths
        for path in additional_paths:
            name = fs.name(path)
            paths[name] = path

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    def dust_component_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Initialize dictionary for the paths
        paths = OrderedDict()

        # Get all component paths
        disk_path = self.dust_path_for_model(name)
        additional_paths = self.additional_dust_paths_for_model(name)

        # Determine standard component names
        disk_name = fs.name(disk_path)

        # Set standard paths
        paths[disk_name] = disk_path

        # Set additional paths
        for path in additional_paths:
            name = fs.name(path)
            paths[name] = path

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    def bulge_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Bulge path"][index]

    # -----------------------------------------------------------------

    def old_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Old stars path"][index]

    # -----------------------------------------------------------------

    def young_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Young stars path"][index]

    # -----------------------------------------------------------------

    def ionizing_stars_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Ionizing stars path"][index]

    # -----------------------------------------------------------------

    def additional_stars_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        if self["Additional stars paths"].mask[index]: return []
        else: return self["Additional stars paths"][index].split(",")

    # -----------------------------------------------------------------

    def has_additional_stars_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return not self["Additional stars paths"].mask[index] and self["Additional stars paths"][index] != ""

    # -----------------------------------------------------------------

    def dust_path_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return self["Dust path"][index]

    # -----------------------------------------------------------------

    def additional_dust_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        if self["Additional dust paths"].mask[index]: return []
        else: return self["Additional dust paths"][index].split(",")

    # -----------------------------------------------------------------

    def has_additional_dust_paths_for_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        index = self.index_for_model(name)
        return not self["Additional dust paths"].mask[index] and self["Additional dust paths"][index] != ""

# -----------------------------------------------------------------
