#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.data.properties Contains the PropertyFetcher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.units import dimensionless_angles
from astroquery.ned import Ned

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import DataComponent
from ...magic.basics.coordinate import SkyCoordinate
from ..basics.properties import GalaxyProperties
from ...core.tools import tables
from ...dustpedia.core.database import DustPediaDatabase, get_account
from ...core.units.parsing import parse_unit as u
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class PropertyFetcher(DataComponent):
    
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
        super(PropertyFetcher, self).__init__(*args, **kwargs)

        # The DustPedia database
        self.database = DustPediaDatabase()

        # Info from DustPedia database
        self.info = None

        # The galaxy properties object
        self.properties = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the NGC name
        self.get_ngc_name()

        # Get the basic galaxy information listed on the database
        self.get_dustpedia_info()

        # Get properties from NED
        self.get_ned_properties()

        # Get the basic properties from S4G
        self.get_s4g_properties()

        # Get the inclination of the galaxy
        self.get_inclination()

        # Show
        self.show()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(PropertyFetcher, self).setup()

        # Get username and password for the DustPedia database
        if self.config.database.username is not None:
            username = self.config.database.username
            password = self.config.database.password
        else: username, password = get_account()

        # Login to the DustPedia database
        self.database.login(username, password)

        # Create the galaxy properties object
        self.properties = GalaxyProperties(name=self.galaxy_name)

    # -----------------------------------------------------------------

    def get_ngc_name(self):

        """
        This function ...
        :return:
        """

        # Get the NGC name of the galaxy
        self.properties.ngc_name = self.ngc_name

    # -----------------------------------------------------------------

    def get_dustpedia_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy info from the DustPedia database ...")

        # Get the info
        self.info = self.database.get_galaxy_info_table(self.ngc_name_nospaces)

        # Get the HYPERLEDA (or DustPedia) name
        self.properties.hyperleda_name = self.hyperleda_name

    # -----------------------------------------------------------------

    def get_ned_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the NASA/IPAC Extragalactic Database ...")

        # Search on NED
        ned_result = Ned.query_object(self.galaxy_name)
        ned_entry = ned_result[0]

        # Get a more common name for this galaxy (sometimes, the name obtained from NED is one starting with 2MASX .., use the PGC name in this case)
        if ned_entry["Object Name"].startswith("2MASX "): gal_name = self.ngc_name
        else: gal_name = ned_entry["Object Name"]

        # Get the redshift
        gal_redshift = ned_entry["Redshift"]
        if isinstance(gal_redshift, np.ma.core.MaskedConstant): gal_redshift = None

        # Get the type (G=galaxy, HII ...)
        gal_type = ned_entry["Type"]
        if isinstance(gal_type, np.ma.core.MaskedConstant): gal_type = None

        # Get the distance
        #ned_distance = ned_entry["Distance (arcmin)"]
        #if isinstance(ned_distance, np.ma.core.MaskedConstant): ned_distance = None

        # Set properties
        self.properties.common_name = gal_name
        self.properties.redshift = gal_redshift
        self.properties.galaxy_type = gal_type

    # -----------------------------------------------------------------

    def get_s4g_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the S4G catalog ...")

        # The Vizier querying object
        vizier = Vizier(columns=['Name', 'RAJ2000', 'DEJ2000', 'amaj', 'ell', 'Dmean', "e_Dmean", "PA"])
        vizier.ROW_LIMIT = -1

        # Get parameters from S4G catalog
        result = vizier.query_object(self.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        # Galaxy name for S4G catalog
        self.properties.name = table["Name"][0]

        # Galaxy center from decomposition (?)
        ra_center = table["RAJ2000"][0]
        dec_center = table["DEJ2000"][0]
        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame='fk5')
        self.properties.center = center

        # Center position
        #self.properties.center = SkyCoordinate(ra=self.info["RA"][0], dec=self.info["DEC"][0], unit="deg") # center position from DustPedia

        # Distance
        self.properties.distance = table["Dmean"][0] * u("Mpc")
        self.properties.distance_error = table["e_Dmean"][0] * u("Mpc")

        # Major axis, ellipticity, position angle
        self.properties.major_arcsec = table["amaj"][0] * u("arcsec")
        self.properties.major = (self.properties.distance * self.properties.major_arcsec).to("pc", equivalencies=dimensionless_angles())

        # Ellipticity
        self.properties.ellipticity = table["ell"][0]
        self.properties.position_angle = Angle(table["PA"][0] + 90.0, u("deg"))

        # Magnitudes
        #asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        #asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        #asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        #asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        # I CAN ADD THESE ATTRIBUTES BACK TO THE GALAXYPROPERTIES CLASS IF I WANT
        #self.properties.i1_mag = asymptotic_ab_magnitude_i1
        #self.properties.i1_mag_error = asymptotic_ab_magnitude_i1_error
        #self.properties.i2_mag = asymptotic_ab_magnitude_i2
        #self.properties.i2_mag_error = asymptotic_ab_magnitude_i2_error

        #self.properties.i1_fluxdensity = unitconversion.ab_to_jansky(self.properties.i1_mag) * u("Jy")
        #i1_fluxdensity_lower = unitconversion.ab_to_jansky(
        #    self.properties.i1_mag + self.properties.i1_mag_error) * u("Jy")
        #i1_fluxdensity_upper = unitconversion.ab_to_jansky(
        #    self.properties.i1_mag - self.properties.i1_mag_error) * u("Jy")
        #i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=self.properties.i1_fluxdensity)
        #self.properties.i1_error = i1_error.average

        #self.properties.i2_fluxdensity = unitconversion.ab_to_jansky(self.properties.i2_mag) * u("Jy")
        #i2_fluxdensity_lower = unitconversion.ab_to_jansky(
        #    self.properties.i2_mag + self.properties.i2_mag_error) * u("Jy")
        #i2_fluxdensity_upper = unitconversion.ab_to_jansky(
        #    self.properties.i2_mag - self.properties.i2_mag_error) * u("Jy")
        #i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=self.properties.i2_fluxdensity)
        #self.properties.i2_error = i2_error.average

        # Other ...
        # absolute_magnitude_i1 = table["M3.6"][0]
        # absolute_magnitude_i2 = table["M4.5"][0]
        # stellar_mass = 10.0**table["logM_"][0] * u("Msun")

    # -----------------------------------------------------------------

    # BASED ON FUV/NUV!
    # def get_inclination(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Querying the catalog of radial profiles for 161 face-on spirals ...")
    #
    #     # The Vizier querying object
    #     vizier = Vizier()
    #     vizier.ROW_LIMIT = -1
    #
    #     # Radial profiles for 161 face-on spirals (Munoz-Mateos+, 2007)
    #     radial_profiles_result = vizier.query_object(self.galaxy_name, catalog="J/ApJ/658/1006")
    #
    #     # Catalog doesnt contain data for a lot of galaxies
    #     # If it doesnt, use DustPedia galaxy info as a backup solution
    #     if len(radial_profiles_result) == 0 or len(radial_profiles_result[0]) == 0:
    #
    #         inclination = Angle(self.info["Inclination"][0], "deg")
    #
    #     # We have a table and it is not empty
    #     else:
    #
    #         table = radial_profiles_result[0]
    #         # distance = float(table[0]["Dist"])
    #         inclination = Angle(float(table[0]["i"]), "deg")
    #
    #     # Set the inclination
    #     self.properties.inclination = inclination

    # -----------------------------------------------------------------

    def get_inclination(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the estimated galaxy inlination ...")

        # Take the inclination from the DustPedia info (HYPERLEDA)
        inclination = Angle(self.info["Inclination"][0], "deg")

        # Set
        self.properties.inclination = inclination

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Showing galaxy properties ...")

        print("")
        print("Galaxy properties: ")
        print(" - Name: " + self.properties.name)
        print(" - Distance: " + tostr(self.properties.distance, scientific=True, fancy=True, ndigits=3))
        print(" - Constellation: " + self.properties.center.get_constellation())
        print("")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the DustPedia galaxy info
        self.write_info()

        # Write the galaxy properties
        self.write_properties()

    # -----------------------------------------------------------------

    def write_info(self):

        """
        This function ...
        :return:
        """

        # Infom the user
        log.info("Writing the galaxy info ...")

        # Write the galaxy info table
        tables.write(self.info, self.galaxy_info_path)

    # -----------------------------------------------------------------

    def write_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the galaxy properties ...")

        # Write
        self.properties.saveto(self.galaxy_properties_path)

# -----------------------------------------------------------------
