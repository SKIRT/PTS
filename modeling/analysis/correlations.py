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

# Import the relevant PTS classes and modules
from .component import AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools import introspection
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ...core.basics.scatter import Scatter2D

# -----------------------------------------------------------------

correlations_data_path = fs.join(introspection.pts_dat_dir("modeling"), "HeatingCorrelations")
m51_data_path = fs.join(correlations_data_path, "DeLooze2014sSFRHeating.dat")
m31_data_path = fs.join(correlations_data_path, "M31_sSFR_Fyoung.dat")

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
    # sSFR-Funev cell scatter data
    # -----------------------------------------------------------------

    @property
    def ssfr_funev_cells_path(self):
        return fs.join()

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_cells(self):
        return fs.is_file(self.ssfr_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_funev_cells_path", True, write=False)
    def ssfr_funev_cells(self):

        """
        This function ...
        :return:
        """

        #ssfr = self.model.total_ssfr_map_earth
        ##funev =

    # -----------------------------------------------------------------
    # sSFR-Funev pixel scatter data
    # -----------------------------------------------------------------

    @property
    def ssfr_funev_pixels_path(self):
        return fs.join()

    # -----------------------------------------------------------------

    @property
    def has_ssfr_funev_pixels(self):
        return fs.is_file(self.ssfr_funev_pixels)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_funev_pixels_path", True, write=False)
    def ssfr_funev_pixels(self):

        """
        This function ...
        :return:
        """

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
        log.info("Writing the sSFR Funev scatter data ...")

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

# -----------------------------------------------------------------
