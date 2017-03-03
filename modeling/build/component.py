#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.component Contains the BuildComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools import filesystem as fs
from .table import ModelsTable

# -----------------------------------------------------------------

class BuildComponent(GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(BuildComponent, self).__init__(config, interactive)

        # Path to the models table
        self.models_table_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BuildComponent, self).setup(**kwargs)

        # Determine the path to the models table
        self.models_table_path = fs.join(self.models_path, "models.dat")

        # Initialize the models table if necessary
        if not fs.is_file(self.models_table_path):
            table = ModelsTable()
            table.saveto(self.models_table_path)

    # -----------------------------------------------------------------

    @property
    def models_table(self):

        """
        This function ...
        :return:
        """

        # Open the table
        return ModelsTable.from_file(self.models_table_path)

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return:
        """

        return self.models_table.names

# -----------------------------------------------------------------

def get_models_table_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return fs.join(modeling_path, "models", "models.dat")

# -----------------------------------------------------------------

def get_models_table(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return ModelsTable.from_file(get_models_table_path(modeling_path))

# -----------------------------------------------------------------

def get_model_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return get_models_table(modeling_path).names

# -----------------------------------------------------------------
