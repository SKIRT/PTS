#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.component Contains the FittingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..component.component import ModelingComponent
from ...core.tools import filesystem as fs
from .tables import RunsTable

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(config)

        # -- Attributes --

        # Runs table path
        self.runs_table_path = None

        # The database path
        self.database_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup(**kwargs)

        # Set the path to the runs table
        self.runs_table_path = fs.join(self.fit_path, "runs.dat")

        # Set the path to the database
        self.database_path = fs.join(self.fit_path, "database.db")

    # -----------------------------------------------------------------

    @property
    def runs_table(self):

        """
        This function ...
        :return:
        """

        return RunsTable.from_file(self.runs_table_path)

    # -----------------------------------------------------------------

    @property
    def run_names(self):

        """
        This function ...
        :return:
        """

        return self.runs_table.run_names

    # -----------------------------------------------------------------

    def model_for_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        return self.runs_table.model_for_run(run_name)

# -----------------------------------------------------------------

def get_run_names(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    fit_path = fs.join(modeling_path, "fit")
    return fs.directories_in_path(fit_path, returns="name")

# -----------------------------------------------------------------
