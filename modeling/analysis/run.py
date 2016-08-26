#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.run Contains the AnalysisRun class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from .info import AnalysisRunInfo
from ...core.simulation.skifile import SkiFile

# -----------------------------------------------------------------

class AnalysisRun(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        self.galaxy_name = None
        self.info = None

    # -----------------------------------------------------------------

    @classmethod
    def from_info(cls, info_path):

        """
        This function ...
        :param info_path:
        :return:
        """

        # Create the instance
        run = cls()

        # Set the analysis run info
        run.info = AnalysisRunInfo.from_file(info_path)

        # Set galaxy name
        modeling_path = fs.directory_of(fs.directory_of(run.info.path))
        run.galaxy_name = fs.name(modeling_path)

        # Return the analysis run object
        return run

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

    # -----------------------------------------------------------------

    @property
    def path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @property
    def out_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "out")

    # -----------------------------------------------------------------

    @property
    def extr_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "extr")

    # -----------------------------------------------------------------

    @property
    def plot_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "plot")

    # -----------------------------------------------------------------

    @property
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "misc")

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "attenuation")

    # -----------------------------------------------------------------

    @property
    def colours_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "colours")

    # -----------------------------------------------------------------

    @property
    def residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "residuals")

    # -----------------------------------------------------------------

    @property
    def heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "heating")

    # -----------------------------------------------------------------

    @property
    def heating_wavelength_grid_path(self):

        """
        This fucntion ...
        :return:
        """

        return fs.join(self.heating_path, "wavelength_grid.dat")

    # -----------------------------------------------------------------

    @property
    def heating_instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.heating_path, "instruments")

    # -----------------------------------------------------------------

    def heating_simulation_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_path, contribution)

    # -----------------------------------------------------------------

    def heating_ski_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_simulation_path_for_contribution(contribution), self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    def heating_output_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_simulation_path_for_contribution(contribution), "out")

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @property
    def ski_file_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.analysis_run_path, )

    # -----------------------------------------------------------------

    @property
    def ski_file(self):

        """
        This function ...
        :return:
        """

        return SkiFile(self.ski_file_path)

    # -----------------------------------------------------------------

    @property
    def dust_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "dust_grid.dg")

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "wavelength_grid.dat")

    # -----------------------------------------------------------------

    @property
    def instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "instruments")

# -----------------------------------------------------------------
