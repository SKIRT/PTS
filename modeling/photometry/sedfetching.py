#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.sedfetching Contains the SEDFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np

# Import astronomical modules
from astroquery.vizier import Vizier

# Import the relevant PTS classes and modules
from ..core import ObservedSED
from ...core.basics.filter import Filter
from ...core.tools.logging import log
from ...core.basics.configurable import Configurable
from ..preparation import unitconversion

# -----------------------------------------------------------------

# SIMBAD: http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%40434630&Name=M++81&submit=display+all+measurements#lab_meas

# -----------------------------------------------------------------

class SEDFetcher(Configurable):

    """
    This class ...
    """
    
    def __init__(self, config=None):
    
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(SEDFetcher, self).__init__(config, "modeling")

        # -- Attributes --

        # The name of the galaxy
        self.galaxy_name = None

        # The Vizier querying object
        self.vizier = Vizier(keywords=["galaxies"])

        # The observed SED
        self.sed = ObservedSED()

        # The filters
        self.filters = dict()

    # -----------------------------------------------------------------

    def run(self, galaxy_name):

        """
        This function ...
        :param galaxy_name
        """

        # 1. Call the setup function
        self.setup(galaxy_name)

        # 2. If requested, query the GALEX ultraviolet atlas of nearby galaxies catalog (Gil de Paz+, 2007)
        if "galex" in self.config.catalogs: self.get_galex()

        # 3. If requested, query the 2MASS Extended sources catalog (IPAC/UMass, 2003-2006)
        if "2mass" in self.config.catalogs: self.get_2mass()

        # 4. If requested, query the SDSS Photometric Catalog, Release 7 (Adelman-McCarthy+, 2009)
        if "sdss" in self.config.catalogs: self.get_sdss()

        # 5. If requested, query the Radial distribution in SINGS galaxies (I.) catalogs (Munoz-Mateos+, 2009)
        if "sings" in self.config.catalogs: self.get_sings()

    # -----------------------------------------------------------------

    def setup(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Call the setup function of the base class
        super(SEDFetcher, self).setup()

        # Set the galaxy name
        self.galaxy_name = galaxy_name

        # Create a dictionary of filters
        keys = ["FUV", "NUV", "U", "B", "V", "J", "H", "K", "F12", "F25", "F60", "F100"]
        for key in keys: self.filters[key] = Filter.from_string(key)

    # -----------------------------------------------------------------

    def get_sdss(self):

        """
        This function ...
        :return:
        """

        # The SDSS Photometric Catalog, Release 7: "II/294/sdss7": ugriz (2009yCat.2294....0A)
        # Interesting columns:
        #
        # - "zsp": Spectroscopic redshift (when SpObjID>0)
        # - "e_zsp": Mean error on zsp
        # - "umag": Model magnitude in u filter (u) [mag]
        #   Note (5)  : The response curves of the SDSS filters (u g r i z) can be found on the SDSS pages: http://www.sdss.org/dr5/instruments/imager/index.html The central wavelength and FWHM are:
        #   -----------------------------------------------
        #   (nm)      u'      g'      r'      i'      z'
        #   -----------------------------------------------
        #   lambda   354.3    477.0   623.1   762.5   913.4
        #   FWHM      56.7    138.7   137.3   152.6    95.0
        #   -----------------------------------------------
        # - "e_umag": Mean error on umag (err_u) [mag]
        # - "gmag": Model magnitude in g filter (g) [mag]
        # - "e_gmag": Mean error on gmag (err_g) [mag]
        # - "rmag": Model magnitude in r filter (r) [mag]
        # - "e_rmag": Mean error on rmag (err_r) [mag]
        # - "imag": Model magnitude in i filter (i) [mag]
        # - "e_imag": Mean error on imag (err_i) [mag]
        # - "zmag": Model magnitude in z filter (z) [mag]
        # - "e_zmag": Mean error on zmag (err_z) [mag]
        # - "Q": Quality of the observation: 1=bad 2=acceptable 3=good 4=missing 5=hole
        sdss_result = self.vizier.query_object(self.galaxy_name, catalog="II/294/sdss7")

    # -----------------------------------------------------------------

    def get_galex(self):

        """
        This function ...
        :return:
        """

        # GALEX: "J/ApJS/173/185" of "J/ApJS/173/185/galex": B and V (2007ApJS..173..185G)
        # Interesting columns:
        #
        # - "MajAxis": Major-axis diameter of the D25 ellipse [arcmin]
        # - "MinAxis": Minor-axis diameter of the D25 ellipse [arcmin]
        # - "PA": Position angle of the D25 ellipse [deg]
        # - "Dist": Distance to the galaxy [Mpc]
        # - "Morph": Morphological type
        # - "T": Morphological type T
        # - "FUV": FUV (120-177nm) image mean sky background [ct/s]
        # - "sigFUV": Mean standard deviation of the FUV sky. Note (2): Measured by averaging the standard deviation within several regions around the position of the object. [ct/s]
        # - "e_FUV": The standard deviation of the mean FUV sky [ct/s]
        # - "NUV": Mean NUV (177-300nm) sky background [ct/s]
        # - "sigNUV": Mean standard deviation of the NUV sky [ct/s]
        # - "e_NUV": The standard deviation of the mean NUV sky [ct/s]
        # - "D25FUV": Observed D25 ellipse FUV band AB magnitude [mag]
        # - "e_D25FUV": Uncertainty in D25FUV [mag]
        # - "D25NUV": Observed D25 ellipse NUV band AB magnitude [mag]
        # - "e_D25NUV": Uncertainty in D25NUV [mag]
        # - "AFUV": Foreground FUV extinction [mag]
        # - "ANUV": Foreground NUV extinction [mag]
        # - "asyFUV": Observed asymptotic FUV (120-177nm) AB magnitude [mag]
        # - "e_asyFUV": Uncertainty in asyFUV [mag]
        # - "asyNUV": Observed asymptotic NUV (177-300nm) AB magnitude [mag]
        # - "e_asyNUV": Uncertainty in asyNUV [mag]
        # - "logFUV": Log of the FUV (120-177nm) luminosity [W]
        # - "logNUV": Log of the NUV (177-300nm) luminosity [W]
        # - "Umag": Johnson U band integrated magnitude [mag] Note (1): In the Vega scale. Published as part of the RC3 catalog, Cat. VII/155)
        # - "e_Umag": Uncertainty in Umag [mag]
        # - "Bmag": Johnson B band integrated magnitude [mag]
        # - "e_Bmag": Uncertainty in Bmag [mag]
        # - "Vmag": Johnson V band integrated magnitude [mag]
        # - "e_Vmag": Uncertainty in Vmag [mag]
        # - "Jmag": 2MASS J band total magnitude [mag]
        # - "e_Jmag": Uncertainty in Jmag [mag]
        # - "Hmag": 2MASS H band total magnitude [mag]
        # - "e_Hmag": Uncertainty in Hmag [mag]
        # - "Kmag": 2MASS K band total magnitude [mag]
        # - "e_Kmag": Uncertainty in Kmag [mag]
        # - "F12um": IRAS 12 micron flux density [Jy]
        # - "e_F12um": Uncertainty in F12um [Jy]
        # - "F25um": IRAS 25 micron flux density [Jy]
        # - "e_F25um": Uncertainty in F25um [Jy]
        # - "F60um": IRAS 60 micron flux density [Jy]
        # - "e_F60um": Uncertainty in F60um [Jy]
        # - "F100um": IRAS 100 micron flux density [Jy]
        # - "e_F100um": Uncertainty in F100um [Jy]

        result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJS/173/185/galex")
        # Result is a list of tables, we only have one table with one entry

        # All AB magnitudes
        fuv_mag = result[0][0]["asyFUV"]
        fuv_mag_error = result[0][0]["e_asyFUV"]
        nuv_mag = result[0][0]["asyNUV"]
        nuv_mag_error = result[0][0]["e_asyNUV"]

        # From Vega magnitude system to AB magnitudes
        u_mag = unitconversion.vega_to_ab(result[0][0]["Umag"], "U")
        u_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Umag"], "U")
        b_mag = unitconversion.vega_to_ab(result[0][0]["Bmag"], "B")
        b_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Bmag"], "B")
        v_mag = unitconversion.vega_to_ab(result[0][0]["Vmag"], "V")
        v_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Vmag"], "V")

        # From Vega magnitude system to AB magnitudes
        j_mag = unitconversion.vega_to_ab(result[0][0]["Jmag"], "J")
        j_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Jmag"], "J")
        h_mag = unitconversion.vega_to_ab(result[0][0]["Hmag"], "H")
        h_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Hmag"], "H")
        k_mag = unitconversion.vega_to_ab(result[0][0]["Kmag"], "Ks")
        k_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Kmag"], "Ks")

        # -- Calculate the fluxes in Jy and add them to the SED --

        # FUV

        fuv = unitconversion.ab_to_jansky(fuv_mag)
        fuv_error = unitconversion.ab_to_jansky(fuv_mag_error)
        self.sed.add_entry(self.filters["FUV"], fuv, fuv_error)

        # NUV flux
        nuv = unitconversion.ab_to_jansky(nuv_mag)
        nuv_error = unitconversion.ab_to_jansky(nuv_mag_error)
        self.sed.add_entry(self.filters["NUV"], nuv, nuv_error)

        # U band flux
        u = unitconversion.ab_to_jansky(u_mag)
        u_error = unitconversion.ab_to_jansky(u_mag_error)
        self.sed.add_entry(self.filters["U"], u, u_error)

        # B band flux
        b = unitconversion.ab_to_jansky(b_mag)
        b_error = unitconversion.ab_to_jansky(b_mag_error)
        self.sed.add_entry(self.filters["B"], b, b_error)

        # V band flux
        v = unitconversion.ab_to_jansky(v_mag)
        v_error = unitconversion.ab_to_jansky(v_mag_error)
        self.sed.add_entry(self.filters["V"], v, v_error)

        # J band flux
        j = unitconversion.ab_to_jansky(j_mag)
        j_error = unitconversion.ab_to_jansky(j_mag_error)
        self.sed.add_entry(self.filters["J"], j, j_error)

        # H band flux
        h = unitconversion.ab_to_jansky(h_mag)
        h_error = unitconversion.ab_to_jansky(h_mag_error)
        self.sed.add_entry(self.filters["H"], h, h_error)

        # K band flux
        k = unitconversion.ab_to_jansky(k_mag)
        k_error = unitconversion.ab_to_jansky(k_mag_error)
        self.sed.add_entry(self.filters["K"], k, k_error)

        # F12 band flux
        f12 = result[0][0]["F12um"]
        f12_error = result[0][0]["e_F12um"]
        self.sed.add_entry(self.filters["F12"], f12, f12_error)

        # F25 band flux
        f25 = result[0][0]["F25um"]
        f25_error = result[0][0]["e_F25um"]
        self.sed.add_entry(self.filters["F25"], f25, f25_error)

        # F60 band flux
        f60 = result[0][0]["F60um"]
        f60_error = result[0][0]["e_F60um"]
        self.sed.add_entry(self.filters["F60"], f60, f60_error)

        # F100 band flux
        f100 = result[0][0]["F100um"]
        f100_error = result[0][0]["e_F100um"]
        self.sed.add_entry(self.filters["F100"], f100, f100_error)

    # -----------------------------------------------------------------

    def get_2mass(self):

        """
        This function ...
        :return:
        """

        # 2MASS Extended Catalog: "VII/233/xsc": J, H, K (2006AJ....131.1163S)
        # Interesting columns:
        #
        # - "J.ext": J magnitude in r.ext (j_m_ext) [mag]
        # - "e_J.ext": σ(J.ext) (j_msig_ext) [mag]
        # - "H.ext": H magnitude in r.ext (h_m_ext) [mag]
        # - "e_H.ext": σ(H.ext) (h_msig_ext) [mag]
        # - "K.ext": J magnitude in r.ext (k_m_ext) [mag]
        # - "e_K.ext": σ(K.ext) (k_msig_ext) [mag]
        result = self.vizier.query_object(self.galaxy_name, catalog="VII/233/xsc")

    # -----------------------------------------------------------------

    def get_sings(self):

        """
        This function ...
        :return:
        """

        # Radial distribution in SINGS galaxies. I.
        # J/ApJ/703/1569/table1	(c)Sample (75 rows)
    	# J/ApJ/703/1569/table2	UV, optical and NIR surface photometry profiles (2206 rows)
    	# J/ApJ/703/1569/table3	UV, optical (SDSS) and NIR photometry profiles (2161 rows)
    	# J/ApJ/703/1569/table4	IRAC and MIPS surface photometry (4319 rows)
    	# J/ApJ/703/1569/table5	UV, optical, and near-IR asymptotic magnitudes for SINGS galaxies lacking SDSS data (43 rows)
    	# J/ApJ/703/1569/table6	UV, optical (SDSS), and near-IR asymptotic magnitudes (32 rows)
    	# J/ApJ/703/1569/table7	IRAC and MIPS asymptotic magnitudes (75 rows)
    	# J/ApJ/703/1569/table8	Non-parametrical morphological estimators (300 rows)


        result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/703/1569")

    # -----------------------------------------------------------------