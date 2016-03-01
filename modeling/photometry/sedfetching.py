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

# Import astronomical modules
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

# Import the relevant PTS classes and modules
from ..core import ObservedSED
from ...core.tools import tables, filesystem
from ...core.basics.filter import Filter
from ...core.tools.logging import log
from ...core.basics.errorbar import ErrorBar
from ...core.basics.configurable import Configurable
from ..preparation import unitconversion

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

        # The NGC ID of the galaxy
        self.ngc_id = None

        # The Vizier querying object
        self.vizier = Vizier(keywords=["galaxies"])
        self.vizier.ROW_LIMIT = -1

        # The Simbad querying object
        self.simbad = Simbad()
        self.simbad.ROW_LIMIT = -1

        # The observed SED
        self.seds = dict()

        # The filters
        self.filters = dict()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new SEDFetcher instance
        if arguments.config is not None: fetcher = cls(arguments.config)
        elif arguments.settings is not None: fetcher = cls(arguments.settings)
        else: fetcher = cls()

        # Set the output path
        if arguments.output_path is not None: fetcher.config.output_path = arguments.output_path

        # Run from the command line, so always write out the SEDs
        fetcher.config.write_seds = True
        fetcher.config.writing.seds_path = "SEDs"

        # Return the new instance
        return fetcher

    # -----------------------------------------------------------------

    def run(self, galaxy_name):

        """
        This function ...
        :param galaxy_name
        """

        # 1. Call the setup function
        self.setup(galaxy_name)

        # 2. If requested, query the GALEX ultraviolet atlas of nearby galaxies catalog (Gil de Paz+, 2007)
        if "GALEX" in self.config.catalogs: self.get_galex()

        # 3. If requested, query the 2MASS Extended sources catalog (IPAC/UMass, 2003-2006)
        if "2MASS" in self.config.catalogs: self.get_2mass()

        # 5. If requested, query the Radial distribution in SINGS galaxies (I.) catalogs (Munoz-Mateos+, 2009)
        if "SINGS" in self.config.catalogs: self.get_sings()

        # If requested, query the LVL global optical photometry (Cook+, 2014) catalog
        if "LVL" in self.config.catalogs: self.get_lvl()

        # If requested, query the Spitzer Local Volume Legacy: IR photometry (Dale+, 2009)
        if "Spitzer" in self.config.catalogs: self.get_spitzer()

        # If requested, query the Spitzer/IRS ATLAS project source (Hernan-Caballero+, 2011)
        if "Spitzer/IRS" in self.config.catalogs: self.get_spitzer_irs()

        # If requested, query the Compendium of ISO far-IR extragalactic data (Brauher+, 2008)
        if "IRAS" in self.config.catalogs: self.get_iras()

        # If requested, query the Imperial IRAS-FSC redshift catalogue (IIFSCz) (Wang+, 2009)
        if "IRAS-FSC" in self.config.catalogs: self.get_iras_fsc()

        # SPECIFIC for M81: not enabled, no time to figure out the unit conversion now
        #if self.ngc_id == "NGC 3031": self.get_m81()

        # Other interesting catalogs:
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/199/22
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/212/18/sample&-c=NGC%203031&-c.u=arcmin&-c.r=2&-c.eq=J2000&-c.geom=r&-out.max=50&-out.form=HTML%20Table&-oc.form=sexa
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/220/6

        # Writing
        self.write()

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

        # Get the NGC ID of the galaxy
        result = self.simbad.query_objectids(self.galaxy_name)
        for name in result["ID"]:
            if "NGC" in name:
                splitted = name.split("NGC")
                if splitted[0] == "":
                    number = int(splitted[1])
                    self.ngc_id = "NGC " + str(number)
                    break

        # Create a dictionary of filters
        keys = ["FUV", "NUV", "U", "B", "V", "R", "J", "H", "K", "IRAS 12", "IRAS 25", "IRAS 60", "IRAS 100", "I1", "I2", "I3", "I4", "MIPS 24", "MIPS 70", "MIPS 160", "SDSS u", "SDSS g", "SDSS r", "SDSS i", "SDSS z"]
        for key in keys: self.filters[key] = Filter.from_string(key)

    # -----------------------------------------------------------------

    def get_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the GALEX ultraviolet atlas of nearby galaxies ...")

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

        # Create an SED
        sed = ObservedSED()

        # All AB magnitudes

        # FUV --

        if not result[0]["asyFUV"].mask[0]:

            fuv_mag = result[0]["asyFUV"][0]
            fuv_mag_error = result[0][0]["e_asyFUV"]
            fuv_mag_lower = fuv_mag - fuv_mag_error
            fuv_mag_upper = fuv_mag + fuv_mag_error

            # flux
            fuv = unitconversion.ab_to_jansky(fuv_mag)
            fuv_lower = unitconversion.ab_to_jansky(fuv_mag_upper)
            fuv_upper = unitconversion.ab_to_jansky(fuv_mag_lower)
            fuv_error = ErrorBar(fuv_lower, fuv_upper, at=fuv)
            sed.add_entry(self.filters["FUV"], fuv, fuv_error)

        # NUV --

        if not result[0]["asyNUV"].mask[0]:

            nuv_mag = result[0][0]["asyNUV"]
            nuv_mag_error = result[0][0]["e_asyNUV"]
            nuv_mag_lower = nuv_mag - nuv_mag_error
            nuv_mag_upper = nuv_mag + nuv_mag_error

            # flux
            nuv = unitconversion.ab_to_jansky(nuv_mag)
            nuv_lower = unitconversion.ab_to_jansky(nuv_mag_upper)
            nuv_upper = unitconversion.ab_to_jansky(nuv_mag_lower)
            nuv_error = ErrorBar(nuv_lower, nuv_upper, at=nuv)
            sed.add_entry(self.filters["NUV"], nuv, nuv_error)

        # U band --

        if not result[0]["Umag"].mask[0]:

            # From Vega magnitude system to AB magnitudes
            u_mag = unitconversion.vega_to_ab(result[0][0]["Umag"], "U")
            u_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Umag"], "U")
            u_mag_lower = u_mag - u_mag_error
            u_mag_upper = u_mag + u_mag_error

            # U band flux
            u = unitconversion.ab_to_jansky(u_mag)
            u_lower = unitconversion.ab_to_jansky(u_mag_upper)
            u_upper = unitconversion.ab_to_jansky(u_mag_lower)
            u_error = ErrorBar(u_lower, u_upper, at=u)
            sed.add_entry(self.filters["U"], u, u_error)

        # B band --

        if not result[0]["Bmag"].mask[0]:

            b_mag = unitconversion.vega_to_ab(result[0][0]["Bmag"], "B")
            b_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Bmag"], "B")
            b_mag_lower = b_mag - abs(b_mag_error)
            b_mag_upper = b_mag + abs(b_mag_error)

            #print("bmag", b_mag)
            #print("bmagerror", b_mag_error)
            #print("bmaglower", b_mag_lower)
            #print("bmagupper", b_mag_upper)

            # B band flux
            b = unitconversion.ab_to_jansky(b_mag)
            b_lower = unitconversion.ab_to_jansky(b_mag_upper)
            b_upper = unitconversion.ab_to_jansky(b_mag_lower)
            b_error = ErrorBar(b_lower, b_upper, at=b)
            sed.add_entry(self.filters["B"], b, b_error)

        # V band --

        if not result[0]["Vmag"].mask[0]:

            v_mag = unitconversion.vega_to_ab(result[0][0]["Vmag"], "V")
            v_mag_error = unitconversion.vega_to_ab(result[0][0]["e_Vmag"], "V")
            v_mag_lower = v_mag - v_mag_error
            v_mag_upper = v_mag + v_mag_error

            # V band flux
            v = unitconversion.ab_to_jansky(v_mag)
            v_lower = unitconversion.ab_to_jansky(v_mag_upper)
            v_upper = unitconversion.ab_to_jansky(v_mag_lower)
            v_error = ErrorBar(v_lower, v_upper, at=v)
            sed.add_entry(self.filters["V"], v, v_error)


        # In 2MASS magnitude system -> can be converted directly into Jy (see below)

        # J band --

        if not result[0]["Jmag"].mask[0]:

            j_mag = result[0][0]["Jmag"]
            j_mag_error = result[0][0]["e_Jmag"]
            j_mag_lower = j_mag - j_mag_error
            j_mag_upper = j_mag + j_mag_error

            # J band flux
            j = unitconversion.photometry_2mass_mag_to_jy(j_mag, "J")
            j_lower = unitconversion.photometry_2mass_mag_to_jy(j_mag_upper, "J")
            j_upper = unitconversion.photometry_2mass_mag_to_jy(j_mag_lower, "J")
            j_error = ErrorBar(j_lower, j_upper, at=j)
            sed.add_entry(self.filters["J"], j, j_error)

        # H band --

        if not result[0]["Hmag"].mask[0]:

            h_mag = result[0][0]["Hmag"]
            h_mag_error = result[0][0]["e_Hmag"]
            h_mag_lower = h_mag - h_mag_error
            h_mag_upper = h_mag + h_mag_error

            # H band flux
            h = unitconversion.photometry_2mass_mag_to_jy(h_mag, "H")
            h_lower = unitconversion.photometry_2mass_mag_to_jy(h_mag_upper, "H")
            h_upper = unitconversion.photometry_2mass_mag_to_jy(h_mag_lower, "H")
            h_error = ErrorBar(h_lower, h_upper, at=h)
            sed.add_entry(self.filters["H"], h, h_error)

        # K band --

        if not result[0]["Kmag"].mask[0]:

            k_mag = result[0][0]["Kmag"]
            k_mag_error = result[0][0]["e_Kmag"]
            k_mag_lower = k_mag - k_mag_error
            k_mag_upper = k_mag + k_mag_error

            # K band flux
            k = unitconversion.photometry_2mass_mag_to_jy(k_mag, "Ks")
            k_lower = unitconversion.photometry_2mass_mag_to_jy(k_mag_upper, "Ks")
            k_upper = unitconversion.photometry_2mass_mag_to_jy(k_mag_lower, "Ks")
            k_error = ErrorBar(k_lower, k_upper, at=k)
            sed.add_entry(self.filters["K"], k, k_error)


        # F12 band flux

        if not result[0]["F12um"].mask[0]:

            f12 = result[0][0]["F12um"]
            f12_error = ErrorBar(result[0][0]["e_F12um"])
            sed.add_entry(self.filters["IRAS 12"], f12, f12_error)

        # F25 band flux

        if not result[0]["F25um"].mask[0]:

            f25 = result[0][0]["F25um"]
            f25_error = ErrorBar(result[0][0]["e_F25um"])
            sed.add_entry(self.filters["IRAS 25"], f25, f25_error)

        # F60 band flux

        if not result[0]["F60um"].mask[0]:

            f60 = result[0][0]["F60um"]
            f60_error = ErrorBar(result[0][0]["e_F60um"])
            sed.add_entry(self.filters["IRAS 60"], f60, f60_error)

        # F100 band flux

        if not result[0]["F100um"].mask[0]:

            f100 = result[0][0]["F100um"]
            f100_error = ErrorBar(result[0][0]["e_F100um"])
            sed.add_entry(self.filters["IRAS 100"], f100, f100_error)

        # Add the SED to the dictionary
        self.seds["GALEX"] = sed

    # -----------------------------------------------------------------

    def get_2mass(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the 2MASS Extended Catalog ...")

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

        # Create an SED
        sed = ObservedSED()

        # In 2MASS magnitude system -> can be converted directly into Jy (see below)
        j_mag = result[0][0]["J.ext"]
        j_mag_error = result[0][0]["e_J.ext"]
        j_mag_lower = j_mag - j_mag_error
        j_mag_upper = j_mag + j_mag_error

        h_mag = result[0][0]["H.ext"]
        h_mag_error = result[0][0]["e_H.ext"]
        h_mag_lower = h_mag - h_mag_error
        h_mag_upper = h_mag + h_mag_error

        k_mag = result[0][0]["K.ext"]
        k_mag_error = result[0][0]["e_K.ext"]
        k_mag_lower = k_mag - k_mag_error
        k_mag_upper = k_mag + k_mag_error

        # J band flux
        j = unitconversion.photometry_2mass_mag_to_jy(j_mag, "J")
        j_lower = unitconversion.photometry_2mass_mag_to_jy(j_mag_upper, "J")
        j_upper = unitconversion.photometry_2mass_mag_to_jy(j_mag_lower, "J")
        j_error = ErrorBar(j_lower, j_upper, at=j)
        sed.add_entry(self.filters["J"], j, j_error)

        # H band flux
        h = unitconversion.photometry_2mass_mag_to_jy(h_mag, "H")
        h_lower = unitconversion.photometry_2mass_mag_to_jy(h_mag_upper, "H")
        h_upper = unitconversion.photometry_2mass_mag_to_jy(h_mag_lower, "H")
        h_error = ErrorBar(h_lower, h_upper, at=h)
        sed.add_entry(self.filters["H"], h, h_error)

        # K band flux
        k = unitconversion.photometry_2mass_mag_to_jy(k_mag, "Ks")
        k_lower = unitconversion.photometry_2mass_mag_to_jy(k_mag_upper, "Ks")
        k_upper = unitconversion.photometry_2mass_mag_to_jy(k_mag_lower, "Ks")
        k_error = ErrorBar(k_lower, k_upper, at=k)
        sed.add_entry(self.filters["K"], k, k_error)

        # Add the SED to the dictionary
        self.seds["2MASS"] = sed

    # -----------------------------------------------------------------

    def get_sings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the SINGS catalog ...")

        # Radial distribution in SINGS galaxies. I.
        # J/ApJ/703/1569/table1	(c)Sample (75 rows)
    	# J/ApJ/703/1569/table2	UV, optical and NIR surface photometry profiles (2206 rows)
    	# J/ApJ/703/1569/table3	UV, optical (SDSS) and NIR photometry profiles (2161 rows)
    	# J/ApJ/703/1569/table4	IRAC and MIPS surface photometry (4319 rows)
    	# J/ApJ/703/1569/table5	UV, optical, and near-IR asymptotic magnitudes for SINGS galaxies lacking SDSS data (43 rows)
    	# J/ApJ/703/1569/table6	UV, optical (SDSS), and near-IR asymptotic magnitudes (32 rows)
    	# J/ApJ/703/1569/table7	IRAC and MIPS asymptotic magnitudes (75 rows)
    	# J/ApJ/703/1569/table8	Non-parametrical morphological estimators (300 rows)

        result = self.vizier.get_catalogs("J/ApJ/703/1569")

        # Result is a TableList with 8 tables (0 to 7)
        # We need:
        # - Table6 -> index 5
        # - Table7 -> index 6

        # Table6
        # - "Name": Galaxy name (NGC ...)
        # - "FUV": GALEX FUV 0.153um-band AB magnitude [mag]
        # - "e_FUV": FUV uncertainty [mag]
        # - "NUV": GALEX NUV 0.227um-band AB magnitude [mag]
        # - "e_NUV": NUV uncertainty [mag]
        # - "umag": SDSS u-band 0.354um AB magnitude [mag]
        # - "e_umag": umag uncertainty [mag]
        # - "gmag": SDSS g-band 0.477um AB magnitude [mag]
        # - "e_gmag": gmag uncertainty [mag]
        # - "rmag": SDSS r-band 0.623um AB magnitude [mag]
        # - "e_rmag": rmag uncertainty [mag]
        # - "imag": SDSS i-band 0.762um AB magnitude [mag]
        # - "e_imag": imag uncertainty [mag]
        # - "zmag": SDSS z-band 0.913um AB magnitude [mag]
        # - "e_zmag": zmag uncertainty [mag]
        # - "Jmag": 2MASS J-band 1.25um AB magnitude [mag]
        # - "e_Jmag": Jmag uncertainty [mag]
        # - "Hmag": 2MASS H-band 1.65um AB magnitude [mag]
        # - "e_Hmag": Hmag uncertainty [mag]
        # - "Ksmag": 2MASS Ks-band 2.17um AB magnitude [mag]
        # - "e_Ksmag": Ksmag uncertainty [mag]

        # Find the row index that corresponds with the specified galaxy

        galaxy_index = tables.find_index(result[5], self.ngc_id)

        # FUV
        fuv_mag = result[5][galaxy_index]["FUV"]
        fuv_mag_error = result[5][galaxy_index]["e_FUV"]
        fuv_mag_lower = fuv_mag - fuv_mag_error
        fuv_mag_upper = fuv_mag + fuv_mag_error

        # NUV
        nuv_mag = result[5][galaxy_index]["NUV"]
        nuv_mag_error = result[5][galaxy_index]["e_NUV"]
        nuv_mag_lower = nuv_mag - nuv_mag_error
        nuv_mag_upper = nuv_mag + nuv_mag_error

        # u
        u_mag = result[5][galaxy_index]["umag"]
        u_mag_error = result[5][galaxy_index]["e_umag"]
        u_mag_lower = u_mag - u_mag_error
        u_mag_upper = u_mag + u_mag_error

        # g
        g_mag = result[5][galaxy_index]["gmag"]
        g_mag_error = result[5][galaxy_index]["e_gmag"]
        g_mag_lower = g_mag - abs(g_mag_error)
        g_mag_upper = g_mag + abs(g_mag_error)

        #print("gmag", g_mag)
        #print("gmagerror", g_mag_error)
        #print("gmaglower", g_mag_lower)
        #print("gmagupper", g_mag_upper)

        # r
        r_mag = result[5][galaxy_index]["rmag"]
        r_mag_error = result[5][galaxy_index]["e_rmag"]
        r_mag_lower = r_mag - r_mag_error
        r_mag_upper = r_mag + r_mag_error

        # i
        i_mag = result[5][galaxy_index]["imag"]
        i_mag_error = result[5][galaxy_index]["e_imag"]
        i_mag_lower = i_mag - i_mag_error
        i_mag_upper = i_mag + i_mag_error

        # z
        z_mag = result[5][galaxy_index]["zmag"]
        z_mag_error = result[5][galaxy_index]["e_zmag"]
        z_mag_lower = z_mag - z_mag_error
        z_mag_upper = z_mag + z_mag_error

        # J
        j_mag = result[5][galaxy_index]["Jmag"]
        j_mag_error = result[5][galaxy_index]["e_Jmag"]
        j_mag_lower = j_mag - j_mag_error
        j_mag_upper = j_mag + j_mag_error

        # H
        h_mag = result[5][galaxy_index]["Hmag"]
        h_mag_error = result[5][galaxy_index]["e_Hmag"]
        h_mag_lower = h_mag - h_mag_error
        h_mag_upper = h_mag + h_mag_error

        # Ks
        k_mag = result[5][galaxy_index]["Ksmag"]
        k_mag_error = result[5][galaxy_index]["e_Ksmag"]
        k_mag_lower = k_mag - k_mag_error
        k_mag_upper = k_mag + k_mag_error


        # Create an SED
        sed = ObservedSED()

        # FUV
        fuv = unitconversion.ab_to_jansky(fuv_mag)
        fuv_lower = unitconversion.ab_to_jansky(fuv_mag_upper)
        fuv_upper = unitconversion.ab_to_jansky(fuv_mag_lower)
        fuv_error = ErrorBar(fuv_lower, fuv_upper, at=fuv)
        sed.add_entry(self.filters["FUV"], fuv, fuv_error)

        # NUV
        nuv = unitconversion.ab_to_jansky(nuv_mag)
        nuv_lower = unitconversion.ab_to_jansky(nuv_mag_upper)
        nuv_upper = unitconversion.ab_to_jansky(nuv_mag_lower)
        nuv_error = ErrorBar(nuv_lower, nuv_upper, at=nuv)
        sed.add_entry(self.filters["NUV"], nuv, nuv_error)

        # u
        u = unitconversion.ab_to_jansky(u_mag)
        u_lower = unitconversion.ab_to_jansky(u_mag_upper)
        u_upper = unitconversion.ab_to_jansky(u_mag_lower)
        u_error = ErrorBar(u_lower, u_upper, at=u)
        sed.add_entry(self.filters["SDSS u"], u, u_error)

        # g
        g = unitconversion.ab_to_jansky(g_mag)
        g_lower = unitconversion.ab_to_jansky(g_mag_upper)
        g_upper = unitconversion.ab_to_jansky(g_mag_lower)
        g_error = ErrorBar(g_lower, g_upper, at=g)
        sed.add_entry(self.filters["SDSS g"], g, g_error)

        # r
        r = unitconversion.ab_to_jansky(r_mag)
        r_lower = unitconversion.ab_to_jansky(r_mag_upper)
        r_upper = unitconversion.ab_to_jansky(r_mag_lower)
        r_error = ErrorBar(r_lower, r_upper, at=r)
        sed.add_entry(self.filters["SDSS r"], r, r_error)

        # i
        i = unitconversion.ab_to_jansky(i_mag)
        i_lower = unitconversion.ab_to_jansky(i_mag_upper)
        i_upper = unitconversion.ab_to_jansky(i_mag_lower)
        i_error = ErrorBar(i_lower, i_upper, at=i)
        sed.add_entry(self.filters["SDSS i"], i, i_error)

        # z
        z = unitconversion.ab_to_jansky(z_mag)
        z_lower = unitconversion.ab_to_jansky(z_mag_upper)
        z_upper = unitconversion.ab_to_jansky(z_mag_lower)
        z_error = ErrorBar(z_lower, z_upper, at=z)
        sed.add_entry(self.filters["SDSS z"], z, z_error)

        # J
        j = unitconversion.ab_to_jansky(j_mag)
        j_lower = unitconversion.ab_to_jansky(j_mag_upper)
        j_upper = unitconversion.ab_to_jansky(j_mag_lower)
        j_error = ErrorBar(j_lower, j_upper, at=j)
        sed.add_entry(self.filters["J"], j, j_error)

        # H
        h = unitconversion.ab_to_jansky(h_mag)
        h_lower = unitconversion.ab_to_jansky(h_mag_upper)
        h_upper = unitconversion.ab_to_jansky(h_mag_lower)
        h_error = ErrorBar(h_lower, h_upper, at=h)
        sed.add_entry(self.filters["H"], h, h_error)

        # Ks
        k = unitconversion.ab_to_jansky(k_mag)
        k_lower = unitconversion.ab_to_jansky(k_mag_upper)
        k_upper = unitconversion.ab_to_jansky(k_mag_lower)
        k_error = ErrorBar(k_lower, k_upper, at=k)
        sed.add_entry(self.filters["K"], k, k_error)

        # Table7: IRAC and MIPS asymptotic magnitudes
        # - "logF3.6": Spitzer/IRAC 3.6um flux density [logJy]
        # - "e_logF3.6": logF3.6 uncertainty [logJy]
        # - "logF4.5": Spitzer/IRAC 4.5um flux density [logJy]
        # - "e_logF4.5": logF4.5 uncertainty [logJy]
        # - "logF5.8": Spitzer/IRAC 5.8um flux density [logJy]
        # - "e_logF5.8": logF5.8 uncertainty [logJy]
        # - "logF8.0": Spiter/IRAC 8.0um flux density [logJy]
        # - "e_logF8.0": logF8.0 uncertainty [logJy]
        # - "logF24": Spiter/MIPS 24um flux density [logJy]
        # - "e_logF24": logF24 uncertainty [logJy]
        # - "logF70": Spiter/MIPS 70um flux density [logJy]
        # - "e_logF70": logF70 uncertainty [logJy]
        # - "logF160": Spiter/MIPS 160um flux density [logJy]
        # - "e_logF160": logF160 uncertainty [logJy]

        # Table7 -> index 6

        galaxy_index = tables.find_index(result[6], self.ngc_id)

        # 3.6 micron
        i1_log = result[6][galaxy_index]["logF3.6"]
        i1_log_error = result[6][galaxy_index]["e_logF3.6"]
        i1_log_lower = i1_log - i1_log_error
        i1_log_upper = i1_log + i1_log_error

        # 4.5 micron
        i2_log = result[6][galaxy_index]["logF4.5"]
        i2_log_error = result[6][galaxy_index]["e_logF4.5"]
        i2_log_lower = i2_log - i2_log_error
        i2_log_upper = i2_log + i2_log_error

        # 5.8 micron
        i3_log = result[6][galaxy_index]["logF5.8"]
        i3_log_error = result[6][galaxy_index]["e_logF5.8"]
        i3_log_lower = i3_log - i3_log_error
        i3_log_upper = i3_log + i3_log_error

        # 8.0 micron
        i4_log = result[6][galaxy_index]["logF8.0"]
        i4_log_error = result[6][galaxy_index]["e_logF8.0"]
        i4_log_lower = i4_log - i4_log_error
        i4_log_upper = i4_log + i4_log_error

        #print("i4log", i4_log)
        #print("i4_log_error", i4_log_error)
        #print("i4_log_lower", i4_log_lower)
        #print("i4_log_upper", i4_log_upper)

        # 24 micron
        mips24_log = result[6][galaxy_index]["logF24"]
        mips24_log_error = result[6][galaxy_index]["e_logF24"]
        mips24_log_lower = mips24_log - mips24_log_error
        mips24_log_upper = mips24_log + mips24_log_error

        # 70 micron
        mips70_log = result[6][galaxy_index]["logF70"]
        mips70_log_error = result[6][galaxy_index]["e_logF70"]
        mips70_log_lower = mips70_log - mips70_log_error
        mips70_log_upper = mips70_log + mips70_log_error

        # 160 micron
        mips160_log = result[6][galaxy_index]["logF160"]
        mips160_log_error = result[6][galaxy_index]["e_logF160"]
        mips160_log_lower = mips160_log - mips160_log_error
        mips160_log_upper = mips160_log + mips160_log_error

        # Calculate data points and errobars in Janskys, add to the SED

        # 3.6 micron
        i1 = 10.**i1_log
        i1_lower = 10.**i1_log_lower
        i1_upper = 10.**i1_log_upper
        i1_error = ErrorBar(i1_lower, i1_upper, at=i1)
        sed.add_entry(self.filters["I1"], i1, i1_error)

        # 4.5 micron
        i2 = 10.**i2_log
        i2_lower = 10.**i2_log_lower
        i2_upper = 10.**i2_log_upper
        i2_error = ErrorBar(i2_lower, i2_upper, at=i2)
        sed.add_entry(self.filters["I2"], i2, i2_error)

        # 5.8 micron
        i3 = 10.**i3_log
        i3_lower = 10.**i3_log_lower
        i3_upper = 10.**i3_log_upper
        i3_error = ErrorBar(i3_lower, i3_upper, at=i3)
        sed.add_entry(self.filters["I3"], i3, i3_error)

        # 8.0 micron
        i4 = 10.**i4_log
        i4_lower = 10.**i4_log_lower
        i4_upper = 10.**i4_log_upper
        i4_error = ErrorBar(i4_lower, i4_upper, at=i4)
        sed.add_entry(self.filters["I4"], i4, i4_error)

        # 24 micron
        mips24 = 10.**mips24_log
        mips24_lower = 10.**mips24_log_lower
        mips24_upper = 10.**mips24_log_upper
        mips24_error = ErrorBar(mips24_lower, mips24_upper, at=mips24)
        sed.add_entry(self.filters["MIPS 24"], mips24, mips24_error)

        # 70 micron
        mips70 = 10.**mips70_log
        mips70_lower = 10.**mips70_log_lower
        mips70_upper = 10.**mips70_log_upper
        mips70_error = ErrorBar(mips70_lower, mips70_upper, at=mips70)
        sed.add_entry(self.filters["MIPS 70"], mips70, mips70_error)

        # 160 micron
        mips160 = 10.**mips160_log
        mips160_lower = 10.**mips160_log_lower
        mips160_upper = 10.**mips160_log_upper
        mips160_error = ErrorBar(mips160_lower, mips160_upper, at=mips160)
        sed.add_entry(self.filters["MIPS 160"], mips160, mips160_error)

        # Add the SED to the dictionary
        self.seds["SINGS"] = sed

    # -----------------------------------------------------------------

    def get_lvl(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the LVL catalog ...")

        # Create an SED
        sed = ObservedSED()

        # "J/MNRAS/445/881": LVL global optical photometry (Cook+, 2014)
        #  - "J/MNRAS/445/881/sample": Galaxies of the Spitzer Local Volume Legacy (LVL): properties (table1) and R25 photometry
        #  - "J/MNRAS/445/881/table3": Photometry within the IR apertures of Dale et al. (2009, Cat. J/ApJ/703/517) (258 rows)
        #  - "J/MNRAS/445/881/table4": Photometry within the UV apertures of Lee et al. (2011, Cat. J/ApJS/192/6) (258 rows)

        result = self.vizier.query_object(self.galaxy_name, catalog="J/MNRAS/445/881/sample")

        # ALL IN AB MAGNITUDE SYSTEM
        # Umag
        # Bmag
        # Vmag
        # Rmag
        # umag
        # gmag
        # rmag
        # imag
        # zmag

        # On SimBad, only sdss bands are used, from /sample ...

        relevant_bands = [("U", "U"), ("B", "B"), ("V", "V"), ("R", "R"), ("u", "SDSS u"), ("g", "SDSS g"), ("r", "SDSS r"), ("i", "SDSS i"), ("z", "SDSS z")]
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = band_prefix_catalog + "mag"
            error_column_name = "e_" + column_name

            # Skip masked values
            if result[0][column_name].mask[0]: continue

            # AB magnitude
            magnitude = result[0][0][column_name]
            magnitude_error = result[0][0][error_column_name]
            magnitude_lower = magnitude - magnitude_error
            magnitude_upper = magnitude + magnitude_error

            # Convert to Jy
            fluxdensity = unitconversion.ab_to_jansky(magnitude)
            fluxdensity_lower = unitconversion.ab_to_jansky(magnitude_upper)
            fluxdensity_upper = unitconversion.ab_to_jansky(magnitude_lower)
            fluxdensity_error = ErrorBar(fluxdensity_lower, fluxdensity_upper, at=fluxdensity)

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["LVL"] = sed

    # -----------------------------------------------------------------

    def get_spitzer(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the Spitzer catalog ...")

        # "J/ApJ/703/517": The Spitzer Local Volume Legacy: IR photometry (Dale+, 2009)
        # "J/ApJ/703/517/sample": Galaxy sample (table 1) and infrared flux densities (table 2) (258 rows)

        # Create an SED
        sed = ObservedSED()

        result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/703/517/sample")

        # F1.25: 2MASS J band (1.25 micron) flux density [Jy]
        # e_F1.25: Uncertainty in F1.25 [Jy]
        # F1.65: 2MASS H band (1.65 micron) flux density [Jy]
        # e_F1.65: Uncertainty in F1.65 [Jy]
        # F2.17: 2MASS Ks band (2.17 micron) flux density [Jy]
        # e_F2.17: Uncertainty in F2.17 [Jy]
        # F3.6: Spitzer/IRAC 3.5 micron band flux density [Jy]
        # e_F3.6: Uncertainty in F3.6 [Jy]
        # F4.5: Spitzer/IRAC 4.5 micron band flux density [Jy]
        # e_F4.5: Uncertainty in F4.5 [Jy]
        # F5.8: Spitzer/IRAC 5.8 micron band flux density [Jy]
        # e_F5.8: Uncertainty in F5.8 [Jy]
        # F8.0: Spitzer/IRAC 8.0 micron band flux density [Jy]
        # e_F8.0: Uncertainty in F8.0 [Jy]
        # F24: Spitzer/MIPS 24 micron flux density [Jy]
        # e_F24: Uncertainty in F24 [Jy]
        # F70: Spitzer/MIPS 70 micron band flux density [Jy]
        # e_F70: Uncertainty in F70 [Jy]
        # F160: Spitzer/MIPS 160 micron band flux density [Jy]
        # e_F160: Uncertainty in F160 [Jy]

        relevant_bands = [("1.25", "J"), ("1.65", "H"), ("2.17", "K"), ("3.6", "I1"), ("4.5", "I2"), ("5.8", "I3"), ("8.0", "I4"), ("24", "MIPS 24"), ("70", "MIPS 70"), ("160", "MIPS 160")]
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F" + band_prefix_catalog
            error_column_name = "e_" + column_name

            # Skip masked values
            if result[0][column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = result[0][0][column_name]
            fluxdensity_error = ErrorBar(result[0][0][error_column_name])

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["Spitzer"] = sed

    # -----------------------------------------------------------------

    def get_spitzer_irs(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the Spitzer/IRS catalog ...")

        # "J/MNRAS/414/500": Spitzer/IRS ATLAS project source (Hernan-Caballero+, 2011)
        # - "J/MNRAS/414/500/catalog": Spitzer/IRS ATLAS project source catalog, version 1.0 (739 rows)

        # !! Parentheses () are converted into underscores _ in the resulting Astropy tables!!

        # F(3.6): IRAC1 (3.6um) flux density [Jy]
        # e_F(3.6): rms uncertainty on F(3.6) [Jy]
        # F(8.0): IRAC4 (8.0um) flux density [Jy]
        # e_F(8.0): rms uncertainty on F(8.0) [Jy]
        # F(24): 24um flux density [Jy]
        # e_F(24): rms uncertainty on F(24) [Jy]
        # F(70): 70um or similar band flux density [Jy]
        # e_F(70): rms uncertainty on F(70) [Jy]

        # Create an SED
        sed = ObservedSED()

        result = self.vizier.query_object(self.galaxy_name, catalog="J/MNRAS/414/500/catalog")

        relevant_bands = [("3.6", "I1"), ("8.0", "I4"), ("24", "MIPS 24"), ("70", "MIPS 70")]
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F_" + band_prefix_catalog + "_"
            error_column_name = "e_" + column_name

            # Skip masked values
            if result[0][column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = result[0][0][column_name]
            fluxdensity_error = ErrorBar(result[0][0][error_column_name])

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["Spitzer-IRS"] = sed

    # -----------------------------------------------------------------

    def get_iras(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the IRAS catalog ...")

        # "J/ApJS/178/280": Compendium of ISO far-IR extragalactic data (Brauher+, 2008)
        # - "J/ApJS/178/280/table1": *Galaxies and properties

        # F12: IRAS 12um band flux density [Jy]
        # F25: IRAS 25um band flux density [Jy]
        # F60: IRAS 60um band flux density [Jy]
        # F100: IRAS 100um band flux density [Jy]

        # No errors ...

        # Create an SED
        sed = ObservedSED()

        result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJS/178/280/table1")

        relevant_bands = [("12", "IRAS 12"), ("25", "IRAS 25"), ("60", "IRAS 60"), ("100", "IRAS 100")]
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F" + band_prefix_catalog

            # Skip masked values
            if result[0][column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = result[0][0][column_name]
            fluxdensity_error = ErrorBar(0.0)

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["IRAS"] = sed

    # -----------------------------------------------------------------

    def get_iras_fsc(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the IRAS-FSC catalog ...")

        # "J/MNRAS/398/109": Imperial IRAS-FSC redshift catalogue (IIFSCz) (Wang+, 2009)
        # - "J/MNRAS/398/109/iifsczv4": IIFSCz Catalogue (MRR+LW 18/04/09)[spectrum/SED] (60303 rows)

        # S12um: IRAS-FSC flux at 12um [Jy]
        # S25um: IRAS-FSC flux at 25um [Jy]
        # S60um: IRAS-FSC flux at 60um [Jy]
        # S100um: IRAS-FSC flux at 100um [Jy]

        # not used in Simbad:
        # umag: SDSS u magnitude [mag: which system??]
        # e_umag:
        # gmag:
        # e_gmag:
        # rmag:
        # e_rmag:
        # imag:
        # e_imag:
        # zmag:
        # e_zmag:

        # Jmag: 2MASS J magnitude [mag: which system??]
        # e_Jmag: rms uncertainty on Jmag [mag: which system??]
        # Hmag: 2MASS H magnitude [mag]
        # e_Hmag: rms uncertainty on Hmag [mag]
        # Kmag: 2MASS K magnitude [mag]
        # e_Kmag rms uncertainty on Kmag [mag]

        # Create an SED
        sed = ObservedSED()

        result = self.vizier.query_object(self.galaxy_name, catalog="J/MNRAS/398/109/iifsczv4")

        relevant_bands = [("12", "IRAS 12"), ("25", "IRAS 25"), ("60", "IRAS 60"), ("100", "IRAS 100")]
        for band_prefix_catalog, filter_name in relevant_bands:

            # Flux and error already in Jy
            fluxdensity = result[0][0]["S" + band_prefix_catalog + "um"]
            fluxdensity_error = ErrorBar(0.0)

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["IRAS-FSC"] = sed

    # -----------------------------------------------------------------

    def get_m81(self):

        """
        This function ...
        :return:
        """

        # UV through far-IR analysis of M81 (Perez-Gonzalez+, 2006)
        # "J/ApJ/648/987"
        # Two tables:
        # - "J/ApJ/648/987/table1": Positions and photometry (GLOBALBKG and LOCALBKG cases) for the regions selected at 160um resolution
        # - "J/ApJ/648/987/table2": Positions and photometry for the regions selected at 24um resolution

        # Table 1

        #result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/648/987/table1")

        # Looks like this (only 3 matches)
        #1	M81	09 55 32.2	+69 03 59.0	680.0	40.78	43.18	42.84	42.63	42.61	41.98	42.67	42.99
        #2	Reg02	09 55 32.2	+69 03 59.0	64.0	39.67	42.63	42.28	41.98	41.76	41.25	41.83	41.78
        #3	Reg03	09 55 32.2	+69 03 59.0	104.0	39.40	42.38	42.03	41.76	41.58	40.90	41.62	41.78

        #galaxy_index = tables.find_index(result, self.galaxy_name)

        # Table 2

        # Interesting rows:
        # - logLFUV: Log of the FUV luminosity [1e-7 W] or [erg/s]
        # - logLNUV: Log of the NUV luminosity [1e-7 W] or [erg/s]
        # - logLHa: Log of the H-alpha luminosity [1e-7 W] or [erg/s]
        # - logL8: Log of the 8um luminosity [1e-7 W] or [erg/s]
        # - logL24: Log of the 24um luminosity [1e-7 W] or [erg/s]

        result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/648/987/table2")

        galaxy_index = tables.find_index(result[0], self.galaxy_name)

        # Create an SED
        sed = ObservedSED()

        relevant_bands = [("FUV", "FUV"), ("NUV", "NUV"), ("Ha", "Ha"), ("8", "I4"), ("24", "MIPS 24")]
        for band_prefix_catalog, filter_name in relevant_bands:

            # Flux and error already in Jy
            log = result[0][galaxy_index]["logL" + band_prefix_catalog]
            flux = 10.**log # In erg/s

            frequency = 0.0 # TODO: calculate this!

            # Also: calculate the amount of flux per unit area ! : divide also by 4 pi d_galaxy**2 !

            fluxdensity = flux / frequency
            fluxdensity_error = ErrorBar(0.0)

            # Add data point to SED
            sed.add_entry(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Add the SED to the dictionary
        self.seds["IRAS-FSC"] = sed

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # If requested, write out the SEDs
        if self.config.write_seds: self.write_seds()

    # -----------------------------------------------------------------

    def write_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SEDs ...")

        # Determine the full path to the SEDs directory
        path = self.full_output_path(self.config.writing.seds_path)

        # Create the SEDs directory if necessary
        filesystem.create_directory(path)

        # Loop over the different SEDs
        for label in self.seds:

            # Debugging info
            log.debug("Writing SED from " + label)

            # Determine the path to the new SED file
            sed_path = filesystem.join(path, label + ".dat")

            # Save the SED at the specified location
            self.seds[label].save(sed_path)

# -----------------------------------------------------------------
