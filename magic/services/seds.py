#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.sedfetching Contains the SEDFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astroquery.vizier import Vizier
from astroquery.vizier.core import TableParseError
from astropy.units import spectral

# Import the relevant PTS classes and modules
from ...core.data.sed import ObservedSED
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ...core.filter.filter import parse_filter
from ...core.basics.log import log
from ...core.basics.errorbar import ErrorBar
from ...dustpedia.data.seds import SEDFetcher as DustPediaSEDFetcher
from ...core.basics.configurable import Configurable
from ..tools import catalogs
from ...modeling.preparation import unitconversion
from ...core.tools import formatting as fmt
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# INTERESTING LINK FOR HALPHA IMAGES AND FLUXES for CERTAIN GALAXIES: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+AS/137/495&-to=3

# -----------------------------------------------------------------

catalog_names = ["DustPedia", "GALEX", "2MASS", "SINGS", "LVL", "Spitzer", "Spitzer/IRS", "IRAS", "IRAS-FSC", "S4G",
                 "SINGS Spectroscopy", "Brown", "Planck"]

# -----------------------------------------------------------------

class SEDFetcher(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SEDFetcher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Determine the NGC name of the galaxy
        self.ngc_name = None

        # The Vizier querying object
        self.vizier = Vizier(columns=["**"])
        self.vizier.ROW_LIMIT = -1

        # The observed SED
        self.seds = dict()

        # The filters
        self.filters = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        """

        # 2. Get the SEDs
        self.get()

        # List the SEDs
        if self.config.list: self.list()

        # Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDFetcher, self).setup(**kwargs)

        # Create a dictionary of filters
        keys = ["Ha", "FUV", "NUV", "U", "B", "V", "R", "J", "H", "K", "IRAS 12", "IRAS 25", "IRAS 60", "IRAS 100",
                "I1", "I2", "I3", "I4", "MIPS 24", "MIPS 70", "MIPS 160", "SDSS u", "SDSS g", "SDSS r", "SDSS i",
                "SDSS z"]
        for key in keys: self.filters[key] = parse_filter(key)

        # Get the NGC name
        if "ngc_name" in kwargs:
            self.ngc_name = kwargs.pop("ngc_name")
        else:
            self.ngc_name = catalogs.get_ngc_name(self.config.galaxy_name)

    # -----------------------------------------------------------------

    def get(self):

        """
        This function ...
        :return:
        """

        # 2. Get the dustpedia SED
        if "DustPedia" in self.config.catalogs: self.get_dustpedia()

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

        # If requested, query the S4G catalog for IRAC fluxes
        if "S4G" in self.config.catalogs: self.get_s4g()

        # If requested, query the Spectroscopy and abundances of SINGS galaxies (Moustakas+, 2010) catalog
        if "SINGS Spectroscopy" in self.config.catalogs: self.get_emission_lines()

        # If requested, query the Atlas of UV-to-MIR galaxy SEDs (Brown+, 2014)
        if "Brown" in self.config.catalogs: self.get_brown()

        # If requested, query the Planck Catalog of Compact Sources Release 1 (Planck, 2013)
        if "Planck" in self.config.catalogs: self.get_planck()

        # Other interesting catalogs:
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/199/22
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/212/18/sample&-c=NGC%203031&-c.u=arcmin&-c.r=2&-c.eq=J2000&-c.geom=r&-out.max=50&-out.form=HTML%20Table&-oc.form=sexa
        # http://vizier.cfa.harvard.edu/viz-bin/VizieR-3?-source=J/ApJS/220/6

    # -----------------------------------------------------------------

    def get_dustpedia(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the DustPedia SED for this galaxy ...")

        # Create the DustPedia photometry object
        dustpedia = DustPediaSEDFetcher()

        # Configure
        dustpedia.config.galaxy_name = self.config.galaxy_name
        dustpedia.config.write = False

        # Run
        dustpedia.run(ngc_name=self.ngc_name)

        # Get the SED
        sed = dustpedia.sed

        # Add the SED
        self.seds["DustPedia"] = sed

    # -----------------------------------------------------------------

    def fetch_table(self, catalog, index=0, object_name=None):

        """
        This function ...
        :param catalog:
        :param index:
        :param object_name:
        :return:
        """

        # Try fetching the catalog
        try:
            if object_name is not None: result = self.vizier.query_object(object_name, catalog=catalog)
            else: result = self.vizier.get_catalogs(catalog)
        except TableParseError:
            log.warning("Could not fetch data from the '" + catalog + "' catalog")
            return None

        # Check the result
        if len(result) == 0:
            log.warning("The result from the '" + catalog + "' is empty")
            return None

        # Get the table
        table = result[index]

        # Return the table
        return table

    # -----------------------------------------------------------------

    def fetch_tables(self, catalog, indices, object_name=None):

        """
        This function ...
        :param catalog:
        :param indices:
        :return:
        """

        # Try fetching the catalog
        try:
            if object_name is not None: result = self.vizier.query_object(object_name, catalog=catalog)
            else: result = self.vizier.get_catalogs(catalog)
        except TableParseError:
            log.warning("Could not fetch data from the '" + catalog + "' catalog")
            #return None
            return [None] * len(indices)

        # Check the result
        if len(result) == 0:
            log.warning("The result from the '" + catalog + "' is empty")
            #return None
            return [None] * len(indices)

        # Get the tables
        tables = []
        for index in indices: tables.append(result[index])
        return tables

    # -----------------------------------------------------------------

    def get_galex(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the GALEX ultraviolet atlas of nearby galaxies ...")

        # Fetch
        table = self.fetch_table("J/ApJS/173/185/galex", object_name=self.config.galaxy_name)
        if table is None: return

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # All AB magnitudes

        # FUV --

        if "asyFUV" in table.colnames and not table["asyFUV"].mask[0]:

            fuv_mag = table["asyFUV"][0]
            fuv_mag_error = table[0]["e_asyFUV"]
            fuv_mag_lower = fuv_mag - fuv_mag_error
            fuv_mag_upper = fuv_mag + fuv_mag_error

            # flux
            fuv = unitconversion.ab_to_jansky(fuv_mag)
            fuv_lower = unitconversion.ab_to_jansky(fuv_mag_upper)
            fuv_upper = unitconversion.ab_to_jansky(fuv_mag_lower)
            fuv_error = ErrorBar(fuv_lower, fuv_upper, at=fuv)
            sed.add_point(self.filters["FUV"], fuv, fuv_error)

        # NUV --

        if "asyNUV" in table.colnames and not table["asyNUV"].mask[0]:

            nuv_mag = table[0]["asyNUV"]
            nuv_mag_error = table[0]["e_asyNUV"]
            nuv_mag_lower = nuv_mag - nuv_mag_error
            nuv_mag_upper = nuv_mag + nuv_mag_error

            # flux
            nuv = unitconversion.ab_to_jansky(nuv_mag)
            nuv_lower = unitconversion.ab_to_jansky(nuv_mag_upper)
            nuv_upper = unitconversion.ab_to_jansky(nuv_mag_lower)
            nuv_error = ErrorBar(nuv_lower, nuv_upper, at=nuv)
            sed.add_point(self.filters["NUV"], nuv, nuv_error)

        # U band --

        if "Umag" in table.colnames and not table["Umag"].mask[0]:

            # From Vega magnitude system to AB magnitudes
            u_mag = unitconversion.vega_to_ab(table[0]["Umag"], "U")
            u_mag_error = unitconversion.vega_to_ab(table[0]["e_Umag"], "U")
            u_mag_lower = u_mag - u_mag_error
            u_mag_upper = u_mag + u_mag_error

            # U band flux
            u = unitconversion.ab_to_jansky(u_mag)
            u_lower = unitconversion.ab_to_jansky(u_mag_upper)
            u_upper = unitconversion.ab_to_jansky(u_mag_lower)
            u_error = ErrorBar(u_lower, u_upper, at=u)
            sed.add_point(self.filters["U"], u, u_error)

        # B band --

        if "Bmag" in table.colnames and not table["Bmag"].mask[0]:

            b_mag = unitconversion.vega_to_ab(table[0]["Bmag"], "B")
            b_mag_error = unitconversion.vega_to_ab(table[0]["e_Bmag"], "B")
            b_mag_lower = b_mag - abs(b_mag_error)
            b_mag_upper = b_mag + abs(b_mag_error)

            # print("bmag", b_mag)
            # print("bmagerror", b_mag_error)
            # print("bmaglower", b_mag_lower)
            # print("bmagupper", b_mag_upper)

            # B band flux
            b = unitconversion.ab_to_jansky(b_mag)
            b_lower = unitconversion.ab_to_jansky(b_mag_upper)
            b_upper = unitconversion.ab_to_jansky(b_mag_lower)
            b_error = ErrorBar(b_lower, b_upper, at=b)
            sed.add_point(self.filters["B"], b, b_error)

        # V band --

        if "Vmag" in table.colnames and not table["Vmag"].mask[0]:

            v_mag = unitconversion.vega_to_ab(table[0]["Vmag"], "V")
            v_mag_error = unitconversion.vega_to_ab(table[0]["e_Vmag"], "V")
            v_mag_lower = v_mag - v_mag_error
            v_mag_upper = v_mag + v_mag_error

            # V band flux
            v = unitconversion.ab_to_jansky(v_mag)
            v_lower = unitconversion.ab_to_jansky(v_mag_upper)
            v_upper = unitconversion.ab_to_jansky(v_mag_lower)
            v_error = ErrorBar(v_lower, v_upper, at=v)
            sed.add_point(self.filters["V"], v, v_error)

        # In 2MASS magnitude system -> can be converted directly into Jy (see below)

        # J band --

        if "Jmag" in table.colnames and not table["Jmag"].mask[0]:

            j_mag = table[0]["Jmag"]
            j_mag_error = table[0]["e_Jmag"]
            j_mag_lower = j_mag - j_mag_error
            j_mag_upper = j_mag + j_mag_error

            # J band flux
            j = unitconversion.photometry_2mass_mag_to_jy(j_mag, "J")
            j_lower = unitconversion.photometry_2mass_mag_to_jy(j_mag_upper, "J")
            j_upper = unitconversion.photometry_2mass_mag_to_jy(j_mag_lower, "J")
            j_error = ErrorBar(j_lower, j_upper, at=j)
            sed.add_point(self.filters["J"], j, j_error)

        # H band --

        if "Hmag" in table.colnames and not table["Hmag"].mask[0]:

            h_mag = table[0]["Hmag"]
            h_mag_error = table[0]["e_Hmag"]
            h_mag_lower = h_mag - h_mag_error
            h_mag_upper = h_mag + h_mag_error

            # H band flux
            h = unitconversion.photometry_2mass_mag_to_jy(h_mag, "H")
            h_lower = unitconversion.photometry_2mass_mag_to_jy(h_mag_upper, "H")
            h_upper = unitconversion.photometry_2mass_mag_to_jy(h_mag_lower, "H")
            h_error = ErrorBar(h_lower, h_upper, at=h)
            sed.add_point(self.filters["H"], h, h_error)

        # K band --

        if "Kmag" in table.colnames and not table["Kmag"].mask[0]:

            k_mag = table[0]["Kmag"]
            k_mag_error = table[0]["e_Kmag"]
            k_mag_lower = k_mag - k_mag_error
            k_mag_upper = k_mag + k_mag_error

            # K band flux
            k = unitconversion.photometry_2mass_mag_to_jy(k_mag, "Ks")
            k_lower = unitconversion.photometry_2mass_mag_to_jy(k_mag_upper, "Ks")
            k_upper = unitconversion.photometry_2mass_mag_to_jy(k_mag_lower, "Ks")
            k_error = ErrorBar(k_lower, k_upper, at=k)
            sed.add_point(self.filters["K"], k, k_error)

        # F12 band flux

        if "F12um" in table.colnames and not table["F12um"].mask[0]:

            f12 = table[0]["F12um"]
            f12_error = ErrorBar(table[0]["e_F12um"])
            sed.add_point(self.filters["IRAS 12"], f12, f12_error)

        # F25 band flux

        if "F25um" in table.colnames and not table["F25um"].mask[0]:

            f25 = table[0]["F25um"]
            f25_error = ErrorBar(table[0]["e_F25um"])
            sed.add_point(self.filters["IRAS 25"], f25, f25_error)

        # F60 band flux

        if "F60um" in table.colnames and not table["F60um"].mask[0]:

            f60 = table[0]["F60um"]
            f60_error = ErrorBar(table[0]["e_F60um"])
            sed.add_point(self.filters["IRAS 60"], f60, f60_error)

        # F100 band flux

        if "F100um" in table.colnames and not table["F100um"].mask[0]:

            f100 = table[0]["F100um"]
            f100_error = ErrorBar(table[0]["e_F100um"])
            sed.add_point(self.filters["IRAS 100"], f100, f100_error)

        # Check the number of points
        if len(sed) == 0:
            log.warning("No photometry found in the GALEX catalog")
            return

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

        # Get the table
        table = self.fetch_table("VII/233/xsc", object_name=self.config.galaxy_name)
        if table is None: return

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # In 2MASS magnitude system -> can be converted directly into Jy (see below)
        if "J.ext" in table.colnames:

            j_mag = table[0]["J.ext"]
            j_mag_error = table[0]["e_J.ext"]
            j_mag_lower = j_mag - j_mag_error
            j_mag_upper = j_mag + j_mag_error

            # J band flux
            j = unitconversion.photometry_2mass_mag_to_jy(j_mag, "J")
            j_lower = unitconversion.photometry_2mass_mag_to_jy(j_mag_upper, "J")
            j_upper = unitconversion.photometry_2mass_mag_to_jy(j_mag_lower, "J")
            j_error = ErrorBar(j_lower, j_upper, at=j)
            sed.add_point(self.filters["J"], j, j_error)

        if "H.ext" in table.colnames:

            h_mag = table[0]["H.ext"]
            h_mag_error = table[0]["e_H.ext"]
            h_mag_lower = h_mag - h_mag_error
            h_mag_upper = h_mag + h_mag_error

            # H band flux
            h = unitconversion.photometry_2mass_mag_to_jy(h_mag, "H")
            h_lower = unitconversion.photometry_2mass_mag_to_jy(h_mag_upper, "H")
            h_upper = unitconversion.photometry_2mass_mag_to_jy(h_mag_lower, "H")
            h_error = ErrorBar(h_lower, h_upper, at=h)
            sed.add_point(self.filters["H"], h, h_error)

        if "K.ext" in table.colnames:
            k_mag = table[0]["K.ext"]
            k_mag_error = table[0]["e_K.ext"]
            k_mag_lower = k_mag - k_mag_error
            k_mag_upper = k_mag + k_mag_error

            # K band flux
            k = unitconversion.photometry_2mass_mag_to_jy(k_mag, "Ks")
            k_lower = unitconversion.photometry_2mass_mag_to_jy(k_mag_upper, "Ks")
            k_upper = unitconversion.photometry_2mass_mag_to_jy(k_mag_lower, "Ks")
            k_error = ErrorBar(k_lower, k_upper, at=k)
            sed.add_point(self.filters["K"], k, k_error)

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the 2MASS catalog")
            return

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

        # Try getting the tables
        table1, table2 = self.fetch_tables("J/ApJ/703/1569", (5,6))

        # Result is a TableList with 8 tables (0 to 7)
        # We need:
        # - Table6 -> index 5
        # - Table7 -> index 6

        # Find the row index that corresponds with the specified galaxy

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        if table1 is not None:

            # Find the index for the galaxy
            galaxy_index = tables.find_index(table1, self.ngc_name)

            # Check
            if galaxy_index is not None:

                # FUV
                fuv_mag = table1[galaxy_index]["FUV"]
                fuv_mag_error = table1[galaxy_index]["e_FUV"]
                fuv_mag_lower = fuv_mag - fuv_mag_error
                fuv_mag_upper = fuv_mag + fuv_mag_error

                # NUV
                nuv_mag = table1[galaxy_index]["NUV"]
                nuv_mag_error = table1[galaxy_index]["e_NUV"]
                nuv_mag_lower = nuv_mag - nuv_mag_error
                nuv_mag_upper = nuv_mag + nuv_mag_error

                # u
                u_mag = table1[galaxy_index]["umag"]
                u_mag_error = table1[galaxy_index]["e_umag"]
                u_mag_lower = u_mag - u_mag_error
                u_mag_upper = u_mag + u_mag_error

                # g
                g_mag = table1[galaxy_index]["gmag"]
                g_mag_error = table1[galaxy_index]["e_gmag"]
                g_mag_lower = g_mag - abs(g_mag_error)
                g_mag_upper = g_mag + abs(g_mag_error)

                # print("gmag", g_mag)
                # print("gmagerror", g_mag_error)
                # print("gmaglower", g_mag_lower)
                # print("gmagupper", g_mag_upper)

                # r
                r_mag = table1[galaxy_index]["rmag"]
                r_mag_error = table1[galaxy_index]["e_rmag"]
                r_mag_lower = r_mag - r_mag_error
                r_mag_upper = r_mag + r_mag_error

                # i
                i_mag = table1[galaxy_index]["imag"]
                i_mag_error = table1[galaxy_index]["e_imag"]
                i_mag_lower = i_mag - i_mag_error
                i_mag_upper = i_mag + i_mag_error

                # z
                z_mag = table1[galaxy_index]["zmag"]
                z_mag_error = table1[galaxy_index]["e_zmag"]
                z_mag_lower = z_mag - z_mag_error
                z_mag_upper = z_mag + z_mag_error

                # J
                j_mag = table1[galaxy_index]["Jmag"]
                j_mag_error = table1[galaxy_index]["e_Jmag"]
                j_mag_lower = j_mag - j_mag_error
                j_mag_upper = j_mag + j_mag_error

                # H
                h_mag = table1[galaxy_index]["Hmag"]
                h_mag_error = table1[galaxy_index]["e_Hmag"]
                h_mag_lower = h_mag - h_mag_error
                h_mag_upper = h_mag + h_mag_error

                # Ks
                k_mag = table1[galaxy_index]["Ksmag"]
                k_mag_error = table1[galaxy_index]["e_Ksmag"]
                k_mag_lower = k_mag - k_mag_error
                k_mag_upper = k_mag + k_mag_error

                # FUV
                fuv = unitconversion.ab_to_jansky(fuv_mag)
                fuv_lower = unitconversion.ab_to_jansky(fuv_mag_upper)
                fuv_upper = unitconversion.ab_to_jansky(fuv_mag_lower)
                fuv_error = ErrorBar(fuv_lower, fuv_upper, at=fuv)
                sed.add_point(self.filters["FUV"], fuv, fuv_error)

                # NUV
                nuv = unitconversion.ab_to_jansky(nuv_mag)
                nuv_lower = unitconversion.ab_to_jansky(nuv_mag_upper)
                nuv_upper = unitconversion.ab_to_jansky(nuv_mag_lower)
                nuv_error = ErrorBar(nuv_lower, nuv_upper, at=nuv)
                sed.add_point(self.filters["NUV"], nuv, nuv_error)

                # u
                u = unitconversion.ab_to_jansky(u_mag)
                u_lower = unitconversion.ab_to_jansky(u_mag_upper)
                u_upper = unitconversion.ab_to_jansky(u_mag_lower)
                u_error = ErrorBar(u_lower, u_upper, at=u)
                sed.add_point(self.filters["SDSS u"], u, u_error)

                # g
                g = unitconversion.ab_to_jansky(g_mag)
                g_lower = unitconversion.ab_to_jansky(g_mag_upper)
                g_upper = unitconversion.ab_to_jansky(g_mag_lower)
                g_error = ErrorBar(g_lower, g_upper, at=g)
                sed.add_point(self.filters["SDSS g"], g, g_error)

                # r
                r = unitconversion.ab_to_jansky(r_mag)
                r_lower = unitconversion.ab_to_jansky(r_mag_upper)
                r_upper = unitconversion.ab_to_jansky(r_mag_lower)
                r_error = ErrorBar(r_lower, r_upper, at=r)
                sed.add_point(self.filters["SDSS r"], r, r_error)

                # i
                i = unitconversion.ab_to_jansky(i_mag)
                i_lower = unitconversion.ab_to_jansky(i_mag_upper)
                i_upper = unitconversion.ab_to_jansky(i_mag_lower)
                i_error = ErrorBar(i_lower, i_upper, at=i)
                sed.add_point(self.filters["SDSS i"], i, i_error)

                # z
                z = unitconversion.ab_to_jansky(z_mag)
                z_lower = unitconversion.ab_to_jansky(z_mag_upper)
                z_upper = unitconversion.ab_to_jansky(z_mag_lower)
                z_error = ErrorBar(z_lower, z_upper, at=z)
                sed.add_point(self.filters["SDSS z"], z, z_error)

                # J
                j = unitconversion.ab_to_jansky(j_mag)
                j_lower = unitconversion.ab_to_jansky(j_mag_upper)
                j_upper = unitconversion.ab_to_jansky(j_mag_lower)
                j_error = ErrorBar(j_lower, j_upper, at=j)
                sed.add_point(self.filters["J"], j, j_error)

                # H
                h = unitconversion.ab_to_jansky(h_mag)
                h_lower = unitconversion.ab_to_jansky(h_mag_upper)
                h_upper = unitconversion.ab_to_jansky(h_mag_lower)
                h_error = ErrorBar(h_lower, h_upper, at=h)
                sed.add_point(self.filters["H"], h, h_error)

                # Ks
                k = unitconversion.ab_to_jansky(k_mag)
                k_lower = unitconversion.ab_to_jansky(k_mag_upper)
                k_upper = unitconversion.ab_to_jansky(k_mag_lower)
                k_error = ErrorBar(k_lower, k_upper, at=k)
                sed.add_point(self.filters["K"], k, k_error)

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

        if table2 is not None:

            # Find galaxy index
            galaxy_index = tables.find_index(table2, self.ngc_name)

            # Check
            if galaxy_index is not None:

                # 3.6 micron
                i1_log = table2[galaxy_index]["logF3.6"]
                i1_log_error = table2[galaxy_index]["e_logF3.6"]
                i1_log_lower = i1_log - i1_log_error
                i1_log_upper = i1_log + i1_log_error

                # 4.5 micron
                i2_log = table2[galaxy_index]["logF4.5"]
                i2_log_error = table2[galaxy_index]["e_logF4.5"]
                i2_log_lower = i2_log - i2_log_error
                i2_log_upper = i2_log + i2_log_error

                # 5.8 micron
                i3_log = table2[galaxy_index]["logF5.8"]
                i3_log_error = table2[galaxy_index]["e_logF5.8"]
                i3_log_lower = i3_log - i3_log_error
                i3_log_upper = i3_log + i3_log_error

                # 8.0 micron
                i4_log = table2[galaxy_index]["logF8.0"]
                i4_log_error = table2[galaxy_index]["e_logF8.0"]
                i4_log_lower = i4_log - i4_log_error
                i4_log_upper = i4_log + i4_log_error

                # print("i4log", i4_log)
                # print("i4_log_error", i4_log_error)
                # print("i4_log_lower", i4_log_lower)
                # print("i4_log_upper", i4_log_upper)

                # 24 micron
                mips24_log = table2[galaxy_index]["logF24"]
                mips24_log_error = table2[galaxy_index]["e_logF24"]
                mips24_log_lower = mips24_log - mips24_log_error
                mips24_log_upper = mips24_log + mips24_log_error

                # 70 micron
                mips70_log = table2[galaxy_index]["logF70"]
                mips70_log_error = table2[galaxy_index]["e_logF70"]
                mips70_log_lower = mips70_log - mips70_log_error
                mips70_log_upper = mips70_log + mips70_log_error

                # 160 micron
                mips160_log = table2[galaxy_index]["logF160"]
                mips160_log_error = table2[galaxy_index]["e_logF160"]
                mips160_log_lower = mips160_log - mips160_log_error
                mips160_log_upper = mips160_log + mips160_log_error

                # Calculate data points and errobars in Janskys, add to the SED

                # 3.6 micron
                i1 = 10. ** i1_log
                i1_lower = 10. ** i1_log_lower
                i1_upper = 10. ** i1_log_upper
                i1_error = ErrorBar(i1_lower, i1_upper, at=i1)
                sed.add_point(self.filters["I1"], i1, i1_error)

                # 4.5 micron
                i2 = 10. ** i2_log
                i2_lower = 10. ** i2_log_lower
                i2_upper = 10. ** i2_log_upper
                i2_error = ErrorBar(i2_lower, i2_upper, at=i2)
                sed.add_point(self.filters["I2"], i2, i2_error)

                # 5.8 micron
                i3 = 10. ** i3_log
                i3_lower = 10. ** i3_log_lower
                i3_upper = 10. ** i3_log_upper
                i3_error = ErrorBar(i3_lower, i3_upper, at=i3)
                sed.add_point(self.filters["I3"], i3, i3_error)

                # 8.0 micron
                i4 = 10. ** i4_log
                i4_lower = 10. ** i4_log_lower
                i4_upper = 10. ** i4_log_upper
                i4_error = ErrorBar(i4_lower, i4_upper, at=i4)
                sed.add_point(self.filters["I4"], i4, i4_error)

                # 24 micron
                mips24 = 10. ** mips24_log
                mips24_lower = 10. ** mips24_log_lower
                mips24_upper = 10. ** mips24_log_upper
                mips24_error = ErrorBar(mips24_lower, mips24_upper, at=mips24)
                sed.add_point(self.filters["MIPS 24"], mips24, mips24_error)

                # 70 micron
                mips70 = 10. ** mips70_log
                mips70_lower = 10. ** mips70_log_lower
                mips70_upper = 10. ** mips70_log_upper
                mips70_error = ErrorBar(mips70_lower, mips70_upper, at=mips70)
                sed.add_point(self.filters["MIPS 70"], mips70, mips70_error)

                # 160 micron
                mips160 = 10. ** mips160_log
                mips160_lower = 10. ** mips160_log_lower
                mips160_upper = 10. ** mips160_log_upper
                mips160_error = ErrorBar(mips160_lower, mips160_upper, at=mips160)
                sed.add_point(self.filters["MIPS 160"], mips160, mips160_error)

        # Check if any points
        if len(sed) == 0:
            log.warning("No photometry found in the SINGS catalog")
            return

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
        sed = ObservedSED(photometry_unit="Jy")

        # "J/MNRAS/445/881": LVL global optical photometry (Cook+, 2014)
        #  - "J/MNRAS/445/881/sample": Galaxies of the Spitzer Local Volume Legacy (LVL): properties (table1) and R25 photometry
        #  - "J/MNRAS/445/881/table3": Photometry within the IR apertures of Dale et al. (2009, Cat. J/ApJ/703/517) (258 rows)
        #  - "J/MNRAS/445/881/table4": Photometry within the UV apertures of Lee et al. (2011, Cat. J/ApJS/192/6) (258 rows)

        #result = self.vizier.query_object(self.config.galaxy_name, catalog="J/MNRAS/445/881/sample")
        # ALL IN AB MAGNITUDE SYSTEM Umag Bmag Vmag Rmag umag gmag rmag imag zmag
        # On SimBad, only sdss bands are used, from /sample ...
        # If nothing is found
        #if len(result) == 0: return

        # Get the table
        table = self.fetch_table("J/MNRAS/445/881/sample", object_name=self.config.galaxy_name)
        if table is None: return

        # Define columns for bands
        relevant_bands = [("U", "U"), ("B", "B"), ("V", "V"), ("R", "R"), ("u", "SDSS u"), ("g", "SDSS g"),
                          ("r", "SDSS r"), ("i", "SDSS i"), ("z", "SDSS z")]

        # Loop over the bands
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = band_prefix_catalog + "mag"
            error_column_name = "e_" + column_name

            # COLUMN DOESN't EXIST?
            if column_name not in table.colnames: continue

            # Skip masked values
            if table[column_name].mask[0]: continue

            # AB magnitude
            magnitude = table[0][column_name]
            magnitude_error = table[0][error_column_name]
            magnitude_lower = magnitude - magnitude_error
            magnitude_upper = magnitude + magnitude_error

            # Convert to Jy
            fluxdensity = unitconversion.ab_to_jansky(magnitude)
            fluxdensity_lower = unitconversion.ab_to_jansky(magnitude_upper)
            fluxdensity_upper = unitconversion.ab_to_jansky(magnitude_lower)
            fluxdensity_error = ErrorBar(fluxdensity_lower, fluxdensity_upper, at=fluxdensity)

            # Add data point to SED
            sed.add_point(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the LVL catalog")
            return

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
        sed = ObservedSED(photometry_unit="Jy")

        #result = self.vizier.query_object(self.config.galaxy_name, catalog="J/ApJ/703/517/sample")
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

        # If no results are found
        #if len(result) == 0: return
        #table = result[0]

        # Get the table
        table = self.fetch_table("J/ApJ/703/517/sample", object_name=self.config.galaxy_name)
        if table is None: return

        # Define column names for bands
        relevant_bands = [("1.25", "J"), ("1.65", "H"), ("2.17", "K"), ("3.6", "I1"), ("4.5", "I2"), ("5.8", "I3"),
                          ("8.0", "I4"), ("24", "MIPS 24"), ("70", "MIPS 70"), ("160", "MIPS 160")]

        # Loop over the bands
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F" + band_prefix_catalog
            error_column_name = "e_" + column_name

            # COLUMN DOESN't EXIST?
            if column_name not in table.colnames: continue

            # Skip masked values
            if table[column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = table[0][column_name]
            fluxdensity_error = ErrorBar(table[0][error_column_name]) if not table[error_column_name].mask[0] else None

            # Add data point to SED
            sed.add_point(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the Spitzer catalog")

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
        sed = ObservedSED(photometry_unit="Jy")

        #result = self.vizier.query_object(self.config.galaxy_name, catalog="J/MNRAS/414/500/catalog")
        # No results found
        #if len(result) == 0: return
        #table = result[0]

        # Get the table
        table = self.fetch_table("J/MNRAS/414/500/catalog", object_name=self.config.galaxy_name)
        if table is None: return

        # Define the column names for the bands
        relevant_bands = [("3.6", "I1"), ("8.0", "I4"), ("24", "MIPS 24"), ("70", "MIPS 70")]

        # Loop over the bands
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F_" + band_prefix_catalog + "_"
            error_column_name = "e_" + column_name

            # COLUMN DOESN't EXIST?
            if column_name not in table.colnames: continue

            # Skip masked values
            if table[column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = table[0][column_name]
            fluxdensity_error = ErrorBar(table[0][error_column_name]) if not table[error_column_name].mask[0] else None

            # Add data point to SED
            sed.add_point(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the Spitzer-IRS catalog")
            return

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
        sed = ObservedSED(photometry_unit="Jy")

        #result = self.vizier.query_object(self.config.galaxy_name, catalog="J/ApJS/178/280/table1")

        #if len(result) != 0:

        # Get table
        table = self.fetch_table("J/ApJS/178/280/table1", object_name=self.config.galaxy_name)
        if table is None: return

        # Define column names for bands
        relevant_bands = [("12", "IRAS 12"), ("25", "IRAS 25"), ("60", "IRAS 60"), ("100", "IRAS 100")]

        # Loop over the bands
        for band_prefix_catalog, filter_name in relevant_bands:

            column_name = "F" + band_prefix_catalog

            # COLUMN DOESN't EXIST?
            if column_name not in table.colnames: continue

            # Skip masked values
            if table[column_name].mask[0]: continue

            # Flux and error already in Jy
            fluxdensity = table[0][column_name]
            fluxdensity_error = ErrorBar(0.0)

            # Add data point to SED
            sed.add_point(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # result = self.vizier.get_catalogs("J/ApJS/178/280/table2")
        # table = result[0]
        #
        # index = tables.find_index(table, self.ngc_name, "Name")
        # # print(index)
        # # index = tables.find_index(table, self.config.galaxy_name, "Name")
        # # print(index)
        # # F170
        # # Jy	(n) 170 micron band flux density
        # # e_F170
        # # Jy	(n) Uncertainty in F170
        # # F158
        # # Jy	(n) 158 micron band flux density
        # # e_F158
        # # Jy	(n) Uncertainty in F158
        # # F145
        # # Jy	(n) 145 micron band flux density
        # # e_F145
        # # Jy	(n) Uncertainty in F145
        # # F122
        # # Jy	(n) 122 micron band flux density
        # # e_F122
        # # Jy	(n) Uncertainty in F122
        # # F88
        # # Jy	(n) 88 micron band flux density
        # # e_F88
        # # Jy	(n) Uncertainty in F88
        # # F63
        # # Jy	(n) 63 micron band flux density
        # # e_F63
        # # Jy	(n) Uncertainty in F63
        # # F57
        # # Jy	(n) 57 micron band flux density
        # # e_F57
        # # Jy	(n) Uncertainty in F57
        # # F52
        # # Jy	(n) 52 micron band flux density
        # # e_F52
        # # Jy	(n) Uncertainty in F52
        # if index is not None: pass
        # if len(sed) == 0: return

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the IRAS catalog")
            return

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
        sed = ObservedSED(photometry_unit="Jy")

        #result = self.vizier.query_object(self.config.galaxy_name, catalog="J/MNRAS/398/109/iifsczv4")
        #if len(result) == 0: return

        # Get table
        table = self.fetch_table("J/MNRAS/398/109/iifsczv4", object_name=self.config.galaxy_name)
        if table is None: return

        # Define column names for bands
        relevant_bands = [("12", "IRAS 12"), ("25", "IRAS 25"), ("60", "IRAS 60"), ("100", "IRAS 100")]

        # Loop over the bands
        for band_prefix_catalog, filter_name in relevant_bands:

            # Define column name
            column_name = "S" + band_prefix_catalog + "um"

            # COLUMN DOESN't EXIST?
            if column_name not in table.colnames: continue

            # Flux and error already in Jy
            fluxdensity = table[0][column_name]
            fluxdensity_error = ErrorBar(0.0)

            # Add data point to SED
            sed.add_point(self.filters[filter_name], fluxdensity, fluxdensity_error)

        # Check the number of points
        if len(sed) == 0:
            log.warning("No photometry found in the IRAS-FSC catalog")
            return

        # Add the SED to the dictionary
        self.seds["IRAS-FSC"] = sed

    # -----------------------------------------------------------------

    def get_s4g(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting fluxes from the S4G catalog ...")

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # Get parameters from S4G catalog
        #result = self.vizier.query_object(self.config.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        #table = result[0]

        # Get table
        table = self.fetch_table("J/PASP/122/1397/s4g", object_name=self.config.galaxy_name)
        if table is None: return

        # I1
        if "__3.6_" in table.colnames:

            i1_mag = table["__3.6_"][0]
            i1_mag_error = table["e__3.6_"][0]
            i1_fluxdensity = unitconversion.ab_to_jansky(i1_mag)
            i1_fluxdensity_lower = unitconversion.ab_to_jansky(i1_mag + i1_mag_error)
            i1_fluxdensity_upper = unitconversion.ab_to_jansky(i1_mag - i1_mag_error)
            i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=i1_fluxdensity)

            # Add data point to SED
            sed.add_point(self.filters["I1"], i1_fluxdensity, i1_error)

        # I2
        if "__4.5_" in table.colnames:

            i2_mag = table["__4.5_"][0]
            i2_mag_error = table["e__4.5_"][0]
            i2_fluxdensity = unitconversion.ab_to_jansky(i2_mag)
            i2_fluxdensity_lower = unitconversion.ab_to_jansky(i2_mag + i2_mag_error)
            i2_fluxdensity_upper = unitconversion.ab_to_jansky(i2_mag - i2_mag_error)
            i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=i2_fluxdensity)

            # Add data point to SED
            sed.add_point(self.filters["I2"], i2_fluxdensity, i2_error)

        # Check number of points
        if len(sed) == 0:
            log.warning("No photometry found in the S4G catalog")
            return

        # Add the SED to the dictionary
        self.seds["S4G"] = sed

    # -----------------------------------------------------------------

    def get_brown(self):

        """
        This function ...
        :return:
        """

        # J/ApJS/212/18/sample
        # AB magnitudes for the sample with neither foreground nor intrinsic dust extinction corrections, and modeled Milky Way foreground dust extinction

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # FUV: [12.5/22.9] GALEX FUV AB band magnitude
        # e_FUV:
        # UVW2:
        # e_UVW2:
        # UVM2:
        # e_UVM2:
        # NUV:
        # e_NUV:
        # UVW1:
        # e_UVW1:
        # Umag: [11.9/15.7] Swift/UVOT U AB band magnitude
        # e_Umag:
        # umag:
        # e_umag:
        # gmag:
        # e_gmag:
        # Vmag:
        # e_Vmag:
        # rmag:
        # e_rmag:
        # imag:
        # e_imag:
        # zmag:
        # e_zmag:
        # Jmag:
        # e_Jmag:
        # Hmag:
        # e_Hmag:
        # Ksmag:
        # e_Ksmag:
        # W1mag:
        # e_W1mag:
        # [3.6]:
        # e_[3.6]:
        # [4.5]:
        # e_[4.5]:
        # W2mag:
        # e_W2mag:
        # [5.8]:
        # e_[5.8]:
        # [8.0]:
        # e_[8.0]:
        # W3mag:
        # e_W3mag:
        # W4mag:
        # e_W4mag:
        # W4'mag: Corrected WISE W4 AB band magnitude
        # e_W4'mag:
        # [24]:
        # e_[24]:

        pass

    # -----------------------------------------------------------------

    def get_planck(self):

        """
        This function ...
        :return:
        """

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # The second release is not yet available ... ??

    # -----------------------------------------------------------------

    def get_emission_lines(self):

        """
        This function ...
        :return:
        """

        # Create an SED
        sed = ObservedSED(photometry_unit="Jy")

        # J/ApJS/190/233/Opt

        # Get result
        #result = self.vizier.get_catalogs("J/ApJS/190/233/Opt")
        #table = result[0]

        # Get result
        table = self.fetch_table("J/ApJS/190/233/Opt")
        if table is None: return
        galaxy_index = tables.find_index(table, self.ngc_name, "Name")

        # FHa: The Hα 6563 Angstrom line flux (aW/m2)
        # e_FHa: Uncertainty in Ha (aW/m2)
        # FNII: The [NII] 6584Å line flux (aW/m2)
        # e_FNII: Uncertainty in NII (aW/m2)

        # Nothing found
        if galaxy_index is None:
            log.warning("No result found for this galaxy in the catalog")
            return

        # H alpha
        ha_flux = table["FHa"][galaxy_index] * u("aW/m2")
        ha_flux_error = table["e_FHa"][galaxy_index] * u("aW/m2")

        # NII
        n2_flux = table["FNII"][galaxy_index] * u("aW/m2")
        n2_flux_error = table["e_FNII"][galaxy_index] * u("aW/m2")

        ha_filter = self.filters["Ha"]
        ha_wavelength = ha_filter.center
        ha_frequency = ha_wavelength.to("Hz", equivalencies=spectral())

        # Calculate flux density
        # print(ha_flux, type(ha_flux))
        # print(ha_frequency, type(ha_frequency))
        #print(ha_flux / ha_frequency)
        ha_fluxdensity = (ha_flux / ha_frequency).to("Jy")
        ha_fluxdensity_error = (ha_flux_error / ha_frequency).to("Jy")
        ha_errorbar = ErrorBar(ha_fluxdensity_error)

        # Add entry
        sed.add_point(ha_filter, ha_fluxdensity, ha_errorbar)

        # Add the SED to the dictionary
        self.seds["Lines"] = sed

    # -----------------------------------------------------------------

    def list(self):

        """
        This function ...
        :return:
        """

        print(fmt.green + fmt.bold + "Found SEDs: (" + str(len(self.seds)) + ")" + fmt.reset)
        print("")

        for label in self.seds:

            print(fmt.underlined + label + fmt.reset + ": " + str(len(self.seds[label])) + " fluxes")
            print("")

            for filter_name in self.seds[label].filter_names():
                print(" * " + filter_name)

            print("")

            # print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the SEDs
        self.write_seds()

    # -----------------------------------------------------------------

    def write_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SEDs ...")

        # Loop over the different SEDs
        for label in self.seds:

            # Debugging
            log.debug("Writing " + label + " SED ...")

            # Determine the path to the new SED file
            sed_path = fs.join(self.config.path, label + ".dat")

            # Save the SED at the specified location
            self.seds[label].saveto(sed_path)

# -----------------------------------------------------------------
