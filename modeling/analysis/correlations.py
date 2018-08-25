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
from ...magic.tools.plotting import plot_scatters, plot_stilts

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

        # Replot?
        if self.config.replot:
            self.config.replot_ssfr_funev = True

        # Replot sSFR-Funev
        if self.config.replot_ssfr_funev:
            self.config.replot_ssfr_salim_funev = True
            self.config.replot_ssfr_ke_funev = True
            self.config.replot_ssfr_mappings_ke_funev = True

        # Salim
        if self.config.replot_ssfr_salim_funev:
            fs.remove_file(self.ssfr_salim_funev_cells_plot_path)
            fs.remove_file(self.ssfr_salim_funev_pixels_plot_path)

        # KE
        if self.config.replot_ssfr_ke_funev:
            fs.remove_file(self.ssfr_ke_funev_cells_plot_path)
            fs.remove_file(self.ssfr_ke_funev_pixels_plot_path)

        # MAPPINGS + KE
        if self.config.replot_ssfr_mappings_ke_funev:
            fs.remove_file(self.ssfr_mappings_ke_funev_cells_plot_path)
            fs.remove_file(self.ssfr_mappings_ke_funev_pixels_plot_path)

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

    @property
    def m51_ssfr_funev_path(self):
        return fs.join(self.ssfr_funev_path, "m51.dat")

    # -----------------------------------------------------------------

    @property
    def has_m51_ssfr_funev(self):
        return fs.is_file(self.m51_ssfr_funev_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_unit(self):
        return "Msun/yr"

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "m51_ssfr_funev_path", True, write=False)
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

        return self.m31_data[1] / 100.

    # -----------------------------------------------------------------

    @property
    def m31_ssfr_funev_path(self):
        return fs.join(self.ssfr_funev_path, "m31.dat")

    # -----------------------------------------------------------------

    @property
    def has_m31_ssfr_funev(self):
        return fs.is_file(self.m31_ssfr_funev_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "m31_ssfr_funev_path", True, write=False)
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
    def cell_funev_path(self):
        return fs.join(self.cell_heating_path, "total_fractions.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_funev(self):
        return fs.is_file(self.cell_funev_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_salim(self):
        return fs.is_file(self.cell_ssfr_salim_path)

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
    def cell_ssfr_salim(self):
        return Data3D.from_file(self.cell_ssfr_salim_path)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_salim_values(self):
        return self.cell_ssfr_salim.values

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
    def ssfr_salim_unit(self):
        return self.cell_ssfr_salim.unit

    # -----------------------------------------------------------------

    @property
    def ssfr_description(self):
        return "specific star formation rate"

    # -----------------------------------------------------------------

    @property
    def funev_description(self):
        return "fraction of dust heating by unevolved stars"

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_salim_mask(self):
        return np.isfinite(self.cell_ssfr_salim_values) * (self.cell_ssfr_salim_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_mask(self):
        return np.isfinite(self.cell_funev_values) * (self.cell_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_salim(self):
        return self.valid_cell_ssfr_salim_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_salim(self):
        return self.cell_ssfr_salim_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_salim(self):
        return self.cell_funev_values[self.valid_cell_mask_salim]

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_cells(self):
        return fs.is_file(self.ssfr_salim_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_salim_funev_cells_path", True, write=False)
    def ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_salim: raise IOError("The cell sSFR (Salim) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        #return Scatter2D.from_xy(self.cell_ssfr_values, self.cell_funev_values, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_salim, self.valid_cell_funev_values_salim, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_salim_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_ke(self):
        return fs.is_file(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_ke(self):
        return Data3D.from_file(self.cell_ssfr_ke_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_unit(self):
        return self.cell_ssfr_ke.unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_ke_values(self):
        return self.cell_ssfr_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_ke_mask(self):
        return np.isfinite(self.cell_ssfr_ke_values) * (self.cell_ssfr_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_mask(self):
        return np.isfinite(self.cell_funev_values) * (self.cell_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_ke(self):
        return self.valid_cell_ssfr_ke_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_ke(self):
        return self.cell_ssfr_ke_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_ke(self):
        return self.cell_funev_values[self.valid_cell_mask_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_cells(self):
        return fs.is_file(self.ssfr_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_ke_funev_cells_path", True, write=False)
    def ssfr_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_ke: raise IOError("The cell sSFR (K&E) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_ke, self.valid_cell_funev_values_ke, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_ke_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_path(self):
        return fs.join(self.cell_sfr_path, "ssfr_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_cell_ssfr_mappings_ke(self):
        return fs.is_file(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_ssfr_mappings_ke(self):
        return Data3D.from_file(self.cell_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_unit(self):
        return self.cell_ssfr_mappings_ke.unit

    # -----------------------------------------------------------------

    @property
    def cell_ssfr_mappings_ke_values(self):
        return self.cell_ssfr_mappings_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_mappings_ke_mask(self):
        return np.isfinite(self.cell_ssfr_mappings_ke_values) * (self.cell_ssfr_mappings_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_mask_mappings_ke(self):
        return self.valid_cell_ssfr_mappings_ke_mask * self.valid_cell_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_ssfr_values_mappings_ke(self):
        return self.cell_ssfr_mappings_ke_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_cell_funev_values_mappings_ke(self):
        return self.cell_funev_values[self.valid_cell_mask_mappings_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_cells(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_ke_funev_cells_path", True, write=False)
    def ssfr_mappings_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_cell_ssfr_mappings_ke: raise IOError("The cell sSFR (MAPPINGS + K&E) data is not present: run the SFR analysis first")
        if not self.has_cell_funev: raise IOError("The cell Funev data is not present: run the cell heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_cell_ssfr_values_mappings_ke, self.valid_cell_funev_values_mappings_ke, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_ke_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # sSFR-Funev pixel scatter data
    # -----------------------------------------------------------------

    @property
    def projected_heating_maps_path(self):
        return fs.join(self.projected_heating_path, "maps")

    # -----------------------------------------------------------------

    @property
    def pixel_funev_path(self):
        return fs.join(self.projected_heating_maps_path, "earth_diffuse.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_funev(self):
        return fs.is_file(self.pixel_funev_path)

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
    def pixel_ssfr_salim_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_salim_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_salim(self):
        return fs.is_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_salim(self):
        return Frame.from_file(self.pixel_ssfr_salim_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_salim_values(self):
        return self.pixel_ssfr_salim.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_salim_mask(self):
        return np.isfinite(self.pixel_ssfr_salim_values) * (self.pixel_ssfr_salim_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_mask(self):
        return np.isfinite(self.pixel_funev_values) * (self.pixel_funev_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_salim(self):
        return self.valid_pixel_ssfr_salim_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_salim(self):
        return self.pixel_ssfr_salim_values[self.valid_pixel_mask_salim]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_salim(self):
        return self.pixel_funev_values[self.valid_pixel_mask_salim]

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_salim.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_pixels(self):
        return fs.is_file(self.ssfr_salim_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_salim_funev_pixels_path", True, write=False)
    def ssfr_salim_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_salim: raise IOError("The sSFR (Salim) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        #return Scatter2D.from_xy(self.pixel_ssfr_values, self.pixel_funev_values, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_unit, x_description=self.ssfr_description, y_description=self.funev_description)
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_salim, self.valid_pixel_funev_values_salim, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_salim_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_ke_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_ke(self):
        return fs.is_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_ke(self):
        return Frame.from_file(self.pixel_ssfr_ke_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_ke_values(self):
        return self.pixel_ssfr_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_ke_mask(self):
        return np.isfinite(self.pixel_ssfr_ke_values) * (self.pixel_ssfr_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_ke(self):
        return self.valid_pixel_ssfr_ke_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_ke(self):
        return self.pixel_ssfr_ke_values[self.valid_pixel_mask_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_ke(self):
        return self.pixel_funev_values[self.valid_pixel_mask_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_pixels(self):
        return fs.is_file(self.ssfr_ke_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_ke_funev_pixels_path", True, write=False)
    def ssfr_ke_funev_pixels(self):

        """
        This functino ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_ke: raise IOError("The sSFR (K&E) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_ke, self.valid_pixel_funev_values_ke, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_ke_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_ke_path(self):
        return fs.join(self.projected_sfr_path, "ssfr_mappings_ke_earth.fits")

    # -----------------------------------------------------------------

    @property
    def has_pixel_ssfr_mappings_ke(self):
        return fs.is_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixel_ssfr_mappings_ke(self):
        return Frame.from_file(self.pixel_ssfr_mappings_ke_path)

    # -----------------------------------------------------------------

    @property
    def pixel_ssfr_mappings_ke_values(self):
        return self.pixel_ssfr_mappings_ke.values

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_mappings_ke_mask(self):
        return np.isfinite(self.pixel_ssfr_mappings_ke_values) * (self.pixel_ssfr_mappings_ke_values != 0)

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_mask_mappings_ke(self):
        return self.valid_pixel_ssfr_mappings_ke_mask * self.valid_pixel_funev_mask

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_ssfr_values_mappings_ke(self):
        return self.pixel_ssfr_mappings_ke_values[self.valid_pixel_mask_mappings_ke]

    # -----------------------------------------------------------------

    @lazyproperty
    def valid_pixel_funev_values_mappings_ke(self):
        return self.pixel_funev_values[self.valid_pixel_mask_mappings_ke]

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings_ke.dat")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_pixels(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_pixels_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(Scatter2D, "ssfr_mappings_ke_funev_pixels_path", True, write=False)
    def ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Checks
        if not self.has_pixel_ssfr_mappings_ke: raise IOError("The sSFR (MAPPINGS + K&E) frame is not present: run the SFR analysis first")
        if not self.has_pixel_funev: raise IOError("The Funev frame is not present: run the projected heating analysis first")

        # Create and return
        return Scatter2D.from_xy(self.valid_pixel_ssfr_values_mappings_ke, self.valid_pixel_funev_values_mappings_ke, x_name=self.ssfr_name, y_name=self.funev_name, x_unit=self.ssfr_mappings_ke_unit, x_description=self.ssfr_description, y_description=self.funev_description)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_ssfr_salim_funev_cells(self):
        return self.has_cell_ssfr_salim and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_salim_funev_pixels(self):
        return self.has_pixel_ssfr_salim and self.has_pixel_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_ke_funev_cells(self):
        return self.has_cell_ssfr_ke and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_ke_funev_pixels(self):
        return self.has_pixel_ssfr_ke and self.has_pixel_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_ke_funev_cells(self):
        return self.has_cell_ssfr_mappings_ke and self.has_cell_funev

    # -----------------------------------------------------------------

    @property
    def do_ssfr_mappings_ke_funev_pixels(self):
        return self.has_pixel_ssfr_mappings_ke and self.has_pixel_funev

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

        # Salim
        self.write_ssfr_salim_funev()

        # K&E
        self.write_ssfr_ke_funev()

        # MAPPINGS + K&E
        self.write_ssfr_mappings_ke_funev()

        # Write the M51 scatter data
        if self.do_write_m51_ssfr_funev: self.write_m51_ssfr_funev()

        # Write the M31 scatter data
        if self.do_write_m31_ssfr_funev: self.write_m31_ssfr_funev()

    # -----------------------------------------------------------------
    # SALIM
    # -----------------------------------------------------------------

    def write_ssfr_salim_funev(self):

        """
        This function ...
        :return:
        """

        # Write the sSFR to Funev cell scatter data
        if self.do_write_ssfr_salim_funev_cells: self.write_ssfr_salim_funev_cells()

        # Write the sSFR to Funev pixel scatter data
        if self.do_write_ssfr_salim_funev_pixels: self.write_ssfr_salim_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_salim_funev_cells(self):
        return self.do_ssfr_salim_funev_cells and not self.has_ssfr_salim_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR to Funev dust cell scatter data ...")

        # Write
        self.ssfr_salim_funev_cells.saveto(self.ssfr_salim_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_salim_funev_pixels(self):
        return self.do_ssfr_salim_funev_pixels and not self.has_ssfr_salim_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_salim_funev_pixels(self):

        """
        This fnction ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (Salim) to Funev pixel scatter data ...")

        # Write
        self.ssfr_salim_funev_pixels.saveto(self.ssfr_salim_funev_pixels_path)

    # -----------------------------------------------------------------
    # KENNICUTT & EVANS
    # -----------------------------------------------------------------

    def write_ssfr_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        self.write_ssfr_ke_funev_cells()

        # Pixels
        self.write_ssfr_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_ke_funev_cells(self):
        return self.do_ssfr_ke_funev_cells and not self.has_ssfr_ke_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_ke_funev_cells(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (K&E) to Funev cell scatter data ...")

        # Write
        self.ssfr_ke_funev_cells.saveto(self.ssfr_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_ke_funev_pixels(self):
        return self.do_ssfr_ke_funev_pixels and not self.has_ssfr_ke_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (K&E) to Funev pixel scatter data ...")

        # Write
        self.ssfr_ke_funev_pixels.saveto(self.ssfr_ke_funev_pixels_path)

    # -----------------------------------------------------------------
    # MAPPINGS + KENNICUTT & EVANS
    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        self.write_ssfr_mappings_ke_funev_cells()

        # Pixels
        self.write_ssfr_mappings_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_ke_funev_cells(self):
        return self.do_ssfr_mappings_ke_funev_cells and not self.has_ssfr_mappings_ke_funev_cells

    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev_cells(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS + K&E) to Funev cell scatter data ...")

        # Write
        self.ssfr_mappings_ke_funev_cells.saveto(self.ssfr_mappings_ke_funev_cells_path)

    # -----------------------------------------------------------------

    @property
    def do_write_ssfr_mappings_ke_funev_pixels(self):
        return self.do_ssfr_mappings_ke_funev_pixels and not self.has_ssfr_mappings_ke_funev_pixels

    # -----------------------------------------------------------------

    def write_ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the sSFR (MAPPINGS + K&E) to Funev pixel scatter data ...")

        # Write
        self.ssfr_mappings_ke_funev_pixels.saveto(self.ssfr_mappings_ke_funev_pixels_path)

    # -----------------------------------------------------------------
    # M51
    # -----------------------------------------------------------------

    @property
    def do_write_m51_ssfr_funev(self):
        return not self.has_m51_ssfr_funev

    # -----------------------------------------------------------------

    def write_m51_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the M51 sSFR to Funev scatter data ...")

        # Write
        self.m51_scatter.saveto(self.m51_ssfr_funev_path)

    # -----------------------------------------------------------------
    # M31
    # -----------------------------------------------------------------

    @property
    def do_write_m31_ssfr_funev(self):
        return not self.has_m31_ssfr_funev

    # -----------------------------------------------------------------

    def write_m31_ssfr_funev(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the M31 sSFR to Funev scatter data ...")

        # Write
        self.m31_scatter.saveto(self.m31_ssfr_funev_path)

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

        # Salim
        self.plot_ssfr_salim_funev()

        # K&E
        self.plot_ssfr_ke_funev()

        # MAPPINGS + K&E
        self.plot_ssfr_mappings_ke_funev()

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_salim_funev_cells: self.plot_ssfr_salim_funev_cells()

        # Pixels
        if self.do_plot_ssfr_salim_funev_pixels: self.plot_ssfr_salim_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_salim.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_cells_plot(self):
        return fs.is_file(self.ssfr_salim_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_salim_funev_cells(self):
        return self.do_ssfr_salim_funev_cells and not self.has_ssfr_salim_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_salim_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_salim_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_cells_title(self):
        return "Correlation between sSFR (Salim) and Funev of dust cells"

    # -----------------------------------------------------------------

    @property
    def ssfr_limits(self):
        return (1e-13,1e-9,)

    # -----------------------------------------------------------------

    @property
    def funev_limits(self):
        return (0,1,)

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (Salim) to Funev dust cell scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_salim_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_salim_funev_cells_title, path=self.ssfr_salim_funev_cells_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_salim_funev_cells_scatters, title=self.ssfr_salim_funev_cells_title, xlog=True,
                      path=self.ssfr_salim_funev_cells_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits,
                      density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_salim.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_salim_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_salim_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_salim_funev_pixels(self):
        return self.do_ssfr_salim_funev_pixels and not self.has_ssfr_salim_funev_pixels_plot

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_salim_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_salim_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_salim_funev_pixels_title(self):
        return "Correlation between sSFR (Salim) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_salim_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (Salim) to Funev pixel scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_salim_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_salim_funev_pixels_title, path=self.ssfr_salim_funev_pixels_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_salim_funev_pixels_scatters, title=self.ssfr_salim_funev_pixels_title, xlog=True,
                      path=self.ssfr_salim_funev_pixels_plot_path, xlimits=self.ssfr_limits, ylimits=self.funev_limits,
                      density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_ke_funev_cells(self):
        return self.do_ssfr_ke_funev_cells and not self.has_ssfr_ke_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_ke_funev_pixels(self):
        return self.do_ssfr_ke_funev_pixels and not self.has_ssfr_ke_funev_pixels_plot

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_ke_funev_cells: self.plot_ssfr_ke_funev_cells()

        # Pixels
        if self.do_plot_ssfr_ke_funev_pixels: self.plot_ssfr_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_cells_plot(self):
        return fs.is_file(self.ssfr_ke_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_ke_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_ke_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_cells_title(self):
        return "Correlation between sSFR (K&E) and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the sSFR (K&E) to Funev cell scatter data ...")

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_ke_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                            self.ssfr_description, self.funev_description,
                                            title=self.ssfr_ke_funev_cells_title, path=self.ssfr_ke_funev_cells_plot_path,
                                            ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_ke_funev_cells_scatters, title=self.ssfr_ke_funev_cells_title, xlog=True,
                          path=self.ssfr_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits,
                          density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_ke_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_ke_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_ke_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_ke_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_ke_funev_pixels_title(self):
        return "Correlation between sSFR (K&E) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_ke_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                           self.ssfr_description, self.funev_description,
                                           title=self.ssfr_ke_funev_pixels_title, path=self.ssfr_ke_funev_pixels_plot_path,
                                           ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_ke_funev_pixels_scatters, title=self.ssfr_ke_funev_pixels_title, xlog=True,
                          path=self.ssfr_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                          ylimits=self.funev_limits,
                          density=True)

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_ke_funev_cells(self):
        return self.do_ssfr_mappings_ke_funev_cells and not self.has_ssfr_mappings_ke_funev_cells_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_ssfr_mappings_ke_funev_pixels(self):
        return self.do_ssfr_mappings_ke_funev_pixels and not self.has_ssfr_mappings_ke_funev_pixels_plot

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev(self):

        """
        This function ...
        :return:
        """

        # Cells
        if self.do_plot_ssfr_mappings_ke_funev_cells: self.plot_ssfr_mappings_ke_funev_cells()

        # Pixels
        if self.do_plot_ssfr_mappings_ke_funev_pixels: self.plot_ssfr_mappings_ke_funev_pixels()

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_plot_path(self):
        return fs.join(self.ssfr_funev_path, "cells_mappings_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_cells_plot(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_cells_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_ke_funev_cells
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (cells)"] = self.ssfr_mappings_ke_funev_cells_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_cells_title(self):
        return "Correlation between sSFR (MAPPINGS + K&E) and Funev of dust cells"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev_cells(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_mappings_ke_funev_cells_scatter_paths, self.ssfr_name, self.funev_name,
                                            self.ssfr_description, self.funev_description,
                                            title=self.ssfr_mappings_ke_funev_cells_title, path=self.ssfr_mappings_ke_funev_cells_plot_path,
                                            ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_mappings_ke_funev_cells_scatters, title=self.ssfr_mappings_ke_funev_cells_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_cells_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits,
                              density=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_plot_path(self):
        return fs.join(self.ssfr_funev_path, "pixels_mappings_ke.pdf")

    # -----------------------------------------------------------------

    @property
    def has_ssfr_mappings_ke_funev_pixels_plot(self):
        return fs.is_file(self.ssfr_mappings_ke_funev_pixels_plot_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_scatters(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_ke_funev_pixels
        scatters[m51_name] = self.m51_scatter
        scatters[m31_name] = self.m31_scatter
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_scatter_paths(self):

        """
        This function ...
        :return:
        """

        scatters = OrderedDict()
        scatters[self.galaxy_name + " (pixels)"] = self.ssfr_mappings_ke_funev_pixels_path
        scatters[m51_name] = self.m51_ssfr_funev_path
        scatters[m31_name] = self.m31_ssfr_funev_path
        return scatters

    # -----------------------------------------------------------------

    @property
    def ssfr_mappings_ke_funev_pixels_title(self):
        return "Correlation between sSFR (MAPPINGS + K&E) and Funev of model pixels"

    # -----------------------------------------------------------------

    def plot_ssfr_mappings_ke_funev_pixels(self):

        """
        This function ...
        :return:
        """

        # Plot using TOPCAT's STILTS
        if self.config.topcat: plot_stilts(self.ssfr_mappings_ke_funev_pixels_scatter_paths, self.ssfr_name, self.funev_name,
                                        self.ssfr_description, self.funev_description,
                                        title=self.ssfr_mappings_ke_funev_pixels_title, path=self.ssfr_mappings_ke_funev_pixels_plot_path,
                                        ylimits=self.funev_limits, xlog=True, xlimits=self.ssfr_limits)

        # Plot using Matplotlib
        else: plot_scatters(self.ssfr_mappings_ke_funev_pixels_scatters, title=self.ssfr_mappings_ke_funev_pixels_title, xlog=True,
                              path=self.ssfr_mappings_ke_funev_pixels_plot_path, xlimits=self.ssfr_limits,
                              ylimits=self.funev_limits,
                              density=True)

# -----------------------------------------------------------------
