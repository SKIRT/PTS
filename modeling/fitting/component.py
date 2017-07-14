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
from astropy.table import Table

# Import the relevant PTS classes and modules
from ..component.component import ModelingComponent
from .tables import RunsTable
from .run import FittingRun
from .context import FittingContext

# -----------------------------------------------------------------

class FittingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FittingComponent, self).__init__(*args, **kwargs)

        # -- Attributes --

        self.context = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(FittingComponent, self).setup(**kwargs)

        # Load the fitting context
        self.context = FittingContext(self.fit_path)

    # -----------------------------------------------------------------

    @property
    def runs_table_path(self):

        """
        This function ...
        :return:
        """

        return self.context.runs_table_path

    # -----------------------------------------------------------------

    @property
    def database_path(self):

        """
        This function ...
        :return:
        """

        return self.context.database_path

    # -----------------------------------------------------------------

    @property
    def statistics_path(self):

        """
        This function ...
        :return:
        """

        return self.context.statistics_path

    # -----------------------------------------------------------------

    @property
    def populations_path(self):

        """
        This function ...
        :return:
        """

        return self.context.populations_path

    # -----------------------------------------------------------------

    @property
    def earth_instrument_name(self):

        """
        This function ...
        :return:
        """

        return self.context.earth_instrument_name

    # -----------------------------------------------------------------

    def load_fitting_run(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        model_name = self.model_for_run(name)
        return FittingRun(self.config.path, name, model_name)

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

    @property
    def statistics(self):

        """
        This function ...
        :return:
        """

        return Table.read(self.statistics_path)

# -----------------------------------------------------------------
