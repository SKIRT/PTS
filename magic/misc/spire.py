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
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.basics.filter import Filter

# -----------------------------------------------------------------

# Online SPIRE handbook url
html_handbook_url = "http://herschel.esac.esa.int/Docs/SPIRE/html/spire_om.html"

# Local table path
beam_correction_table_path = fs.join(introspection.pts_dat_dir("magic"), "SPIRE", "Colour_Corrections_KcolP_KcolE_SPIRE.csv")

# -----------------------------------------------------------------

class SPIRE(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Load the beam correction table
        Kbeam = np.genfromtxt("./Colour_Corrections_KcolP_KcolE_SPIRE.csv", delimiter=',')  # KcolE colour corections for SPIRE taken from the SPIRE handbook

        # Set SPIRE filters
        self.psw = Filter.from_string("SPIRE PSW")
        self.pmw = Filter.from_string("SPIRE PMW")
        self.plw = Filter.from_string("SPIRE PLW")

        # Interpolate
        self.f250 = interp1d(Kbeam[:, 0], Kbeam[:, 4], kind='cubic')
        self.f350 = interp1d(Kbeam[:, 0], Kbeam[:, 5], kind='cubic')
        self.f500 = interp1d(Kbeam[:, 0], Kbeam[:, 6], kind='cubic')

    # -----------------------------------------------------------------

    def get_ebeam_temperature(self, fltr, temperature):

        """
        This function ...
        :param fltr:
        :param temperature:
        :return:
        """

        # KEcorr factor dependant on model SED temperature and on SPIRE filter
        if fltr == self.psw: return self.get_ebeam_temperature_psw(temperature)
        elif fltr == self.pmw: return self.get_ebeam_temperature_pmw(temperature)
        elif fltr == self.plw: return self.get_ebeam_temperature_plw(temperature)
        else: raise ValueError("Not a SPIRE filter")

    # -----------------------------------------------------------------

    def get_ebeam_temperature_psw(self, temperature):

        """
        This function ...
        :param temperature:
        :return:
        """

        return self.f250(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_ebeam_temperature_pmw(self, temperature):

        """
        This function ...
        :param temperature:
        :return:
        """

        return self.f350(temperature.to("K").value)

    # -----------------------------------------------------------------

    def get_ebeam_temperature_plw(self, temperature):

        """
        This function ...
        :param temperature:
        :return:
        """

        return self.f500(temperature.to("K").value)

# -----------------------------------------------------------------
