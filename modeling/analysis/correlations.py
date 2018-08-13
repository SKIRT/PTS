#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.correlations Contains the CorrelationsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools import introspection
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ...core.basics.scatter import Scatter2D
from ..core.data import Data3D
from ...magic.core.frame import Frame
from ...magic.tools.plotting import plot_scatters

# -----------------------------------------------------------------

# Define paths to reference sSFR to Funev data files
correlations_data_path = fs.join(introspection.pts_dat_dir("modeling"), "HeatingCorrelations")
m51_data_path = fs.join(correlations_data_path, "DeLooze2014sSFRHeating.dat")
m31_data_path = fs.join(correlations_data_path, "M31_sSFR_Fyoung.dat")
m51_name = "M51"
m31_name = "M31"

# -----------------------------------------------------------------

# Names for correlations
ssfr_funev_name = "sSFR-Funev"

# -----------------------------------------------------------------

class CorrelationsAnalyser(AnalysisRunComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(CorrelationsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Writing
        self.write()

        # Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CorrelationsAnalyser, self).setup()

    # -----------------------------------------------------------------

    @property
    def correlations_path(self):
        return self.analysis_run.correlations_path

    # -----------------------------------------------------------------

    @property
    def heating_path(self):
        return self.analysis_run.heating_path

    # -----------------------------------------------------------------

    @property
    def projected_heating_path(self):
        return fs.join(self.heating_path, "projected")

    # -----------------------------------------------------------------

    @property
    def cell_heating_path(self):
        return fs.join(self.heating_path, "cell")

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):
        return self.analysis_run.sfr_path

    # -----------------------------------------------------------------

    @property
    def projected_sfr_path(self):
        return fs.join(self.sfr_path, "projected")

    # -----------------------------------------------------------------

    @property
    def cell_sfr_path(self):
        return fs.join(self.sfr_path, "cell")

    # -----------------------------------------------------------------
    # M51 DATA
    # -----------------------------------------------------------------

    @lazyproperty
    def m51_column_names(self):

        """
        This function ...
        :return:
        """

        return fs.get_column_names(m51_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_data(self):

        """
        This function ...
        :return:
        """

        return np.loadtxt(m51_data_path, unpack=True)

    # -----------------------------------------------------------------

    @property
    def m51_log10_ssfr(self):

        """
        This function ...
        :return:
        """

        return self.m51_data[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_ssfr(self):

        """
        This function ...
        :return:
        """

        return 10**self.m51_log10_ssfr

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_fractions(self):

        """
        This function ...
        :return:
        """

        return 10**self.m51_data[1]

    # -----------------------------------------------------------------

    @lazyproperty
    def m51_scatter(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.m51_ssfr, self.m51_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # M31 data
    # -----------------------------------------------------------------

    @lazyproperty
    def m31_column_names(self):

        """
        Thisn function ...
        :return:
        """

        return fs.get_column_names(m31_data_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_data(self):

        """
        This function ...
        :return:
        """

        return np.loadtxt(m31_data_path, unpack=True)

    # -----------------------------------------------------------------

    @property
    def m31_log10_ssfr(self):

        """
        This function ...
        :return:
        """

        return self.m31_data[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_ssfr(self):

        """
        This function ...
        :return:
        """

        return 10**self.m31_log10_ssfr

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_fractions(self):

        """
        This function ...
        :return:
        """

        return self.m31_data[1] * 100.

    # -----------------------------------------------------------------

    @lazyproperty
    def m31_scatter(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.m31_ssfr, self.m31_fractions, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit,
                                 x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # sSFR-Funev
    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_funev_path(self):
        return fs.create_directory_in(self.correlations_path, ssfr_funev_name)

    # -----------------------------------------------------------------
    # sSFR-Funev cell scatter data
    # -----------------------------------------------------------------

    @property
    def ssfr_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_cells(self):
        return fs.is_file(self.ssfr_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def cell_funev_path(self):
        return fs.join(self.cell_heating_path, "total_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_path(self):
        return fs.join(self.cell_sfr_path, "ssfr.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_funev(self):
        return Data3D.from_file(self.cell_funev_path)

    # -----------------------------------------------------------------

    @property
    def cell_funev_values(self):
        return self.cell_funev.values

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr(self):
        return Data3D.from_file(self.cell_ssfr_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_values(self):
        return self.cell_ssfr.values

    # -----------------------------------------------------------------

    @property
    def ssfr_name(self):
        return "sSFR"

    # -----------------------------------------------------------------

    @property
    def funev_name(self):
        return "Funev"

    # -----------------------------------------------------------------

    @property
    def ssfr_unit(self):
        return self.cell_ssfr.unit

    # -----------------------------------------------------------------

    @property
    def ssfr_description(self):
        return "specific star formation rate"

    # -----------------------------------------------------------------

    @property
    def funev_description(self):
        return "fraction of dust heating by unevolved stars"

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_funev_cells_path", True, write=False)
    def ssfr_funev_cells(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.cell_ssfr_values, self.cell_funev_values, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # sSFR-Funev pixel scatter data
    # -----------------------------------------------------------------

    @property
    def ssfr_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_pixels(self):
        return fs.is_file(self.ssfr_funev_pixels)

    # -----------------------------------------------------------------

    @property
    def projected_heating_maps_path(self):
        return fs.join(self.projected_heating_path, "maps")

    # -----------------------------------------------------------------

    @property
    def pixel_funev_path(self):
        return fs.join(self.projected_heating_maps_path, "earth_diffuse.fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_funev(self):
        return Frame.from_file(self.pixel_funev_path)

    # -----------------------------------------------------------------

    @property
    def pixel_funev_values(self):
        return self.pixel_funev.values

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_earth.fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr(self):
        return Frame.from_file(self.pixel_ssfr_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_values(self):
        return self.pixel_ssfr.values

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_funev_pixels_path", True, write=False)
    def ssfr_funev_pixels(self):

        """
        This function ...
        :return:
        """

        return Scatter2D.from_xy(self.pixel_ssfr_values, self.pixel_funev_values, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # sSFR Funev scatter data
        self.write_ssfr_funev()

    # -----------------------------------------------------------------

    def write_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev scatter data ...")

        # Write the sSFR to Funev cell scatter data
        if self.do_write_ssfr_funev_cells: self.write_ssfr_funev_cells()

        # Write the sSFR to Funev pixel scatter data
        if self.do_write_ssfr_funev_pixels: self.write_ssfr_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_funev_cells(self):
        return not self.has_ssfr_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev dust cell scatter data ...")

        # Write
        self.ssfr_funev_cells.saveto(self.ssfr_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_funev_pixels(self):
        return not self.has_ssfr_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_funev_pixels(self):

        """
        This fnction ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev pixel scatter data ...")

        # Write
        self.ssfr_funev_pixels.saveto(self.ssfr_funev_pixels_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_ssfr_funev()

    # -----------------------------------------------------------------

    def plot_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR to Funev scatter data ...")

        # Cells
        if self.do_plot_ssfr_funev_cells: self.plot_ssfr_funev_cells()

        # Pixels
        if self.do_plot_ssfr_funev_pixels: self.plot_ssfr_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_cells_plot(self):
        return fs.is_file(self.ssfr_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_funev_cells(self):
        return not self.has_ssfr_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name] = self.ssfr_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_cells_title(self):
        return "Correlation between sSFR and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR to Funev dust cell scatter data ...")

        # Plot
        plot_scatters(self.ssfr_funev_cells_scatters, title=self.ssfr_funev_cells_title, x_scale="log", path=self.ssfr_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_funev_pixels(self):
        return not self.has_ssfr_funev_pixels_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name] = self.ssfr_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_funev_pixels_title(self):
        return "Correlation between sSFR and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plottin the sSFR to Funev pixel scatter data ...")

        # Plot
        plot_scatters(self.ssfr_funev_pixels_scatters, title=self.ssfr_funev_pixels_title, x_scale="log", path=self.ssfr_funev_pixels_plot_path)

# -----------------------------------------------------------------
