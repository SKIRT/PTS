#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.spire Contains the SPIRE class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from scipy.interpolate import interp1d

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.filter.broad import BroadBandFilter

# -----------------------------------------------------------------

# Online SPIRE handbook url
html_handbook_url = "http://herschel.esac.esa.int/Docs/SPIRE/html/spire_om.html"

# Local table paths
magic_dat_path = introspection.pts_dat_dir("magic")
kcol_temperature_beta1_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kcol_temperature_beta1.csv")
kcol_temperature_beta2_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kcol_temperature_beta2.csv")
kcol_spectral_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kcol_spectral.csv")

kbeam_temperature_beta1_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kbeam_temperature_beta1.csv")
kbeam_temperature_beta2_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kbeam_temperature_beta2.csv")
kbeam_spectral_table_path = fs.join(magic_dat_path, "SPIRE", "Colour_Corrections_Kbeam_spectral.csv")

# -----------------------------------------------------------------

class SPIRE(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Create the SPIRE filters
        self.psw = BroadBandFilter("SPIRE PSW")
        self.pmw = BroadBandFilter("SPIRE PMW")
        self.plw = BroadBandFilter("SPIRE PLW")

        # BASED ON TEMPERATURE
        # --------------------

        ## BETA = 1.5

        # Load the Kcol corrections table for beta = 1.5
        kcol_temperature_beta1_data = np.genfromtxt(kcol_temperature_beta1_table_path, delimiter=',')

        # Load point Kcol data
        self.kcol_temperature_beta1_point_f250 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 1], kind='cubic')
        self.kcol_temperature_beta1_point_f350 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 2], kind='cubic')
        self.kcol_temperature_beta1_point_f500 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 3], kind='cubic')

        # Load extended Kcol data
        self.kcol_temperature_beta1_extended_f250 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 4], kind='cubic')
        self.kcol_temperature_beta1_extended_f350 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 5], kind='cubic')
        self.kcol_temperature_beta1_extended_f500 = interp1d(kcol_temperature_beta1_data[:, 0], kcol_temperature_beta1_data[:, 6], kind='cubic')

        # Load the Kbeam corrections table
        kbeam_temperature_beta1_data = np.genfromtxt(kbeam_temperature_beta1_table_path, delimiter=',')

        # Load point Kbeam data
        #self.kbeam_temperature_beta1_point_f250 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 1], kind='cubic')
        #self.kbeam_temperature_beta1_point_f350 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 2], kind='cubic')
        #self.kbeam_temperature_beta1_point_f500 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 3], kind='cubic')

        # Load extended Kbeam data
        #self.kbeam_temperature_beta1_extended_f250 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 4], kind='cubic')
        #self.kbeam_temperature_beta1_extended_f350 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 5], kind='cubic')
        #self.kbeam_temperature_beta1_extended_f500 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 6], kind='cubic')

        # Load Kbeam data
        self.kbeam_temperature_beta1_f250 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 1], kind='cubic')
        self.kbeam_temperature_beta1_f350 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 2], kind='cubic')
        self.kbeam_temperature_beta1_f500 = interp1d(kbeam_temperature_beta1_data[:, 0], kbeam_temperature_beta1_data[:, 3], kind='cubic')

        ## BETA = 2.0

        # Load the Kcol corrections table for beta = 2.0
        kcol_temperature_beta2_data = np.genfromtxt(kcol_temperature_beta2_table_path, delimiter=',')

        # Load point Kcol data
        self.kcol_temperature_beta2_point_f250 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 1], kind='cubic')
        self.kcol_temperature_beta2_point_f350 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 2], kind='cubic')
        self.kcol_temperature_beta2_point_f500 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 3], kind='cubic')

        # Load extended Kcol data
        self.kcol_temperature_beta2_extended_f250 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 4], kind='cubic')
        self.kcol_temperature_beta2_extended_f350 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 5], kind='cubic')
        self.kcol_temperature_beta2_extended_f500 = interp1d(kcol_temperature_beta2_data[:, 0], kcol_temperature_beta2_data[:, 6], kind='cubic')

        # Load the Kbeam corrections table
        kbeam_temperature_beta2_data = np.genfromtxt(kbeam_temperature_beta2_table_path, delimiter=',')

        # Load point Kbeam data
        #self.kbeam_temperature_beta2_point_f250 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 1], kind='cubic')
        #self.kbeam_temperature_beta2_point_f350 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 2], kind='cubic')
        #self.kbeam_temperature_beta2_point_f500 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 3], kind='cubic')

        # Load extended Kbeam data
        #self.kbeam_temperature_beta2_extended_f250 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 4], kind='cubic')
        #self.kbeam_temperature_beta2_extended_f350 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 5], kind='cubic')
        #self.kbeam_temperature_beta2_extended_f500 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 6], kind='cubic')

        # Load Kbeam data
        self.kbeam_temperature_beta2_f250 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 1], kind='cubic')
        self.kbeam_temperature_beta2_f350 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 2], kind='cubic')
        self.kbeam_temperature_beta2_f500 = interp1d(kbeam_temperature_beta2_data[:, 0], kbeam_temperature_beta2_data[:, 3], kind='cubic')

        # BASED ON SPECTRAL INDEX
        # -----------------------

        ## Load the Kcol data
        kcol_spectral_data = np.genfromtxt(kcol_spectral_table_path, delimiter=',')

        # Load point Kcol data
        self.kcol_spectral_point_f250 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 1], kind='cubic')
        self.kcol_spectral_point_f350 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 2], kind='cubic')
        self.kcol_spectral_point_f500 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 3], kind='cubic')

        # Load extended Kcol data
        self.kcol_spectral_extended_f250 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 4], kind='cubic')
        self.kcol_spectral_extended_f350 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 4], kind='cubic')
        self.kcol_spectral_extended_f500 = interp1d(kcol_spectral_data[:, 0], kcol_spectral_data[:, 6], kind='cubic')

        ## Load the Kbeam data
        kbeam_spectral_data = np.genfromtxt(kbeam_spectral_table_path, delimiter=',')

        # Load point Kbeam data
        #self.kbeam_spectral_point_f250 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 1], kind='cubic')
        #self.kbeam_spectral_point_f350 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 2], kind='cubic')
        #self.kbeam_spectral_point_f500 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 3], kind='cubic')

        # Load extended Kbeam data
        #self.kbeam_spectral_extended_f250 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 4], kind='cubic')
        #self.kbeam_spectral_extended_f350 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 4], kind='cubic')
        #self.kbeam_spectral_extended_f500 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 6], kind='cubic')

        # Load Kbeam data
        self.kbeam_spectral_f250 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 4], kind='cubic')
        self.kbeam_spectral_f350 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 5], kind='cubic')
        self.kbeam_spectral_f500 = interp1d(kbeam_spectral_data[:, 0], kbeam_spectral_data[:, 6], kind='cubic')

    # -----------------------------------------------------------------

    def get_kcol_temperature(self, fltr, temperature, beta, extended=True):

        """
        This function ...
        :param fltr:
        :param temperature:
        :param beta:
        :param extended:
        :return:
        """

        if fltr == self.psw: return self.get_kcol_temperature_psw(temperature, beta, extended=extended)
        elif fltr == self.pmw: return self.get_kcol_temperature_pmw(temperature, beta, extended=extended)
        elif fltr == self.plw: return self.get_kcol_temperature_plw(temperature, beta, extended=extended)
        else: raise ValueError("Not a SPIRE filter")

    # -----------------------------------------------------------------

    def get_kcol_spectral(self, fltr, spectral_index, extended=True):

        """
        This function ...
        :param fltr:
        :param spectral_index:
        :param extended:
        :return:
        """

        if fltr == self.psw: return self.get_kcol_spectral_psw(spectral_index, extended=extended)
        elif fltr == self.pmw: return self.get_kcol_spectral_pmw(spectral_index, extended=extended)
        elif fltr == self.plw: return self.get_kcol_spectral_plw(spectral_index, extended=extended)
        else: raise ValueError("Not a SPIRE filter")

    # -----------------------------------------------------------------

    def get_kcol_temperature_psw(self, temperature, beta, extended=True):

        """
        This function ...
        :param temperature:
        :param beta:
        :param extended:
        :return:
        """

        if extended:
            if beta == 1.5: return self.kcol_temperature_beta1_extended_f250(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_extended_f250(temperature.to("K").value)
        else:
            if beta == 1.5: return self.kcol_temperature_beta1_point_f250(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_point_f250(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kcol_temperature_pmw(self, temperature, beta, extended=True):

        """
        This function ...
        :param temperature:
        :param beta:
        :param extended:
        :return:
        """

        if extended:
            if beta == 1.5: return self.kcol_temperature_beta1_extended_f350(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_extended_f350(temperature.to("K").value)
        else:
            if beta == 1.5: return self.kcol_temperature_beta1_point_f350(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_point_f350(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kcol_temperature_plw(self, temperature, beta, extended=True):

        """
        This function ...
        :param temperature:
        :param beta:
        :param extended:
        :return:
        """

        if extended:
            if beta == 1.5: return self.kcol_temperature_beta1_extended_f500(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_extended_f500(temperature.to("K").value)
        else:
            if beta == 1.5: return self.kcol_temperature_beta1_point_f500(temperature.to("K").value)
            elif beta == 2.0: return self.kcol_temperature_beta2_point_f500(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kcol_spectral_psw(self, spectral_index, extended=True):

        """
        This function ...
        :param spectral_index:
        :param extended:
        :return:
        """

        if extended: return self.kcol_spectral_extended_f250(spectral_index)
        else: return self.kcol_spectral_point_f250(spectral_index)

    # -----------------------------------------------------------------

    def get_kcol_spectral_pmw(self, spectral_index, extended=True):

        """
        This function ...
        :param spectral_index:
        :param extended:
        :return:
        """

        if extended: return self.kcol_spectral_extended_f350(spectral_index)
        else: return self.kcol_spectral_point_f350(spectral_index)

    # -----------------------------------------------------------------

    def get_kcol_spectral_plw(self, spectral_index, extended=True):

        """
        This function ...
        :param spectral_index:
        :param extended:
        :return:
        """

        if extended: return self.kcol_spectral_extended_f500(spectral_index)
        else: return self.kcol_spectral_point_f500(spectral_index)

    # -----------------------------------------------------------------

    def get_kbeam_temperature(self, fltr, temperature, beta):

        """
        This function ...
        :param fltr:
        :param temperature:
        :param beta:
        :return:
        """

        if fltr == self.psw: return self.get_kbeam_temperature_psw(temperature, beta)
        elif fltr == self.pmw: return self.get_kbeam_temperature_pmw(temperature, beta)
        elif fltr == self.plw: return self.get_kbeam_temperature_plw(temperature, beta)
        else: raise ValueError("Not a SPIRE filter")

    # -----------------------------------------------------------------

    def get_kbeam_spectral(self, fltr, spectral_index):

        """
        This function ...
        :param fltr:
        :param spectral_index:
        :return:
        """

        if fltr == self.psw: return self.get_kbeam_spectral_psw(spectral_index)
        elif fltr == self.pmw: return self.get_kbeam_spectral_pmw(spectral_index)
        elif fltr == self.plw: return self.get_kbeam_spectral_plw(spectral_index)
        else: raise ValueError("Not a SPIRE filter")

    # -----------------------------------------------------------------

    def get_kbeam_temperature_psw(self, temperature, beta):

        """
        This function ...
        :param temperature:
        :param beta:
        :param extended:
        :return:
        """

        if beta == 1.5: return self.kbeam_temperature_beta1_f250(temperature.to("K").value)
        elif beta == 2.0: return self.kbeam_temperature_beta2_f250(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kbeam_temperature_pmw(self, temperature, beta):

        """
        This function ...
        :param temperature:
        :param beta:
        :return:
        """

        if beta == 1.5: return self.kbeam_temperature_beta1_f350(temperature.to("K").value)
        elif beta == 2.0: return self.kbeam_temperature_beta2_f350(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kbeam_temperature_plw(self, temperature, beta):

        """
        This function ...
        :param temperature:
        :param beta:
        :return:
        """

        if beta == 1.5: return self.kbeam_temperature_beta1_f500(temperature.to("K").value)
        elif beta == 2.0: return self.kbeam_temperature_beta2_f500(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_kbeam_spectral_psw(self, spectral_index):

        """
        This function ...
        :param spectral_index:
        :return:
        """

        return self.kbeam_spectral_f250(spectral_index)

    # -----------------------------------------------------------------

    def get_kbeam_spectral_pmw(self, spectral_index):

        """
        This function ...
        :param spectral_index:
        :return:
        """

        return self.kbeam_spectral_f350(spectral_index)

    # -----------------------------------------------------------------

    def get_kbeam_spectral_plw(self, spectral_index):

        """
        This function ...
        :param spectral_index:
        :return:
        """

        return self.kbeam_spectral_f500(spectral_index)

# -----------------------------------------------------------------
