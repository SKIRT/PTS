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

# Import astronomical modules
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from .component import DataComponent
from ..preparation import unitconversion
from ...core.basics.errorbar import ErrorBar
from ...magic.basics.skygeometry import SkyCoordinate
from ...magic.tools import catalogs
from ..basics.properties import GalaxyProperties
from ...core.tools import tables
from ...dustpedia.core.database import DustPediaDatabase, get_account

# -----------------------------------------------------------------

class PropertyFetcher(DataComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(PropertyFetcher, self).__init__(config)

        # The Vizier querying object
        self.vizier = Vizier()
        self.vizier.ROW_LIMIT = -1

        # The DustPedia database
        self.database = DustPediaDatabase()

        # Info from DustPedia database
        self.info = None

        # The galaxy properties object
        self.properties = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Setup
        self.setup()

        # 2. Get the NGC ID
        self.get_ngc_id()

        # 3. Get the basic galaxy information listed on the database
        self.get_dustpedia_info()

        # 4. Get the basic properties from S4G
        self.get_s4g_properties()

        # 5. Get the inclination of the galaxy
        self.get_inclination()

        # 6. Get spiral properties
        #self.get_spiral_properties()

        # 7. Writing
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
        self.properties = GalaxyProperties()

    # -----------------------------------------------------------------

    def get_ngc_id(self):

        """
        This function ...
        :return:
        """

        # Get the NGC name of the galaxy
        self.properties.ngc_id = catalogs.get_ngc_name(self.galaxy_name)

    # -----------------------------------------------------------------

    def get_dustpedia_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Fetching galaxy info from the DustPedia database ...")

        # Get the info
        self.info = self.database.get_galaxy_info(self.ngc_id_nospaces)

    # -----------------------------------------------------------------

    def get_s4g_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the S4G catalog ...")

        # Get parameters from S4G catalog
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        # Galaxy name for S4G catalog
        self.properties.galaxy_name = table["Name"][0]

        # Galaxy center from decomposition (?)
        ra_center = table["_RAJ2000"][0]
        dec_center = table["_DEJ2000"][0]
        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame='fk5')
        self.properties.center = center

        # Distance
        self.properties.distance = table["Dmean"][0] * Unit("Mpc")
        self.properties.distance_error = table["e_Dmean"][0] * Unit("Mpc")

        # Major axis, ellipticity, position angle
        self.properties.major_arcsec = table["amaj"][0] * Unit("arcsec")
        self.properties.major = (self.properties.distance * self.properties.major_arcsec).to("pc", equivalencies=dimensionless_angles())

        # Ellipticity
        self.properties.ellipticity = table["ell"][0]
        self.properties.position_angle = Angle(table["PA"][0] + 90.0, Unit("deg"))

        # Magnitudes
        asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        self.properties.i1_mag = asymptotic_ab_magnitude_i1
        self.properties.i1_mag_error = asymptotic_ab_magnitude_i1_error
        self.properties.i2_mag = asymptotic_ab_magnitude_i2
        self.properties.i2_mag_error = asymptotic_ab_magnitude_i2_error

        self.properties.i1_fluxdensity = unitconversion.ab_to_jansky(self.properties.i1_mag) * Unit("Jy")
        i1_fluxdensity_lower = unitconversion.ab_to_jansky(
            self.properties.i1_mag + self.properties.i1_mag_error) * Unit("Jy")
        i1_fluxdensity_upper = unitconversion.ab_to_jansky(
            self.properties.i1_mag - self.properties.i1_mag_error) * Unit("Jy")
        i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=self.properties.i1_fluxdensity)
        self.properties.i1_error = i1_error.average

        self.properties.i2_fluxdensity = unitconversion.ab_to_jansky(self.properties.i2_mag) * Unit("Jy")
        i2_fluxdensity_lower = unitconversion.ab_to_jansky(
            self.properties.i2_mag + self.properties.i2_mag_error) * Unit("Jy")
        i2_fluxdensity_upper = unitconversion.ab_to_jansky(
            self.properties.i2_mag - self.properties.i2_mag_error) * Unit("Jy")
        i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=self.properties.i2_fluxdensity)
        self.properties.i2_error = i2_error.average

        # Other ...
        # absolute_magnitude_i1 = table["M3.6"][0]
        # absolute_magnitude_i2 = table["M4.5"][0]
        # stellar_mass = 10.0**table["logM_"][0] * u.Unit("Msun")

    # -----------------------------------------------------------------

    def get_inclination(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the catalog of radial profiles for 161 face-on spirals ...")

        # Radial profiles for 161 face-on spirals (Munoz-Mateos+, 2007)
        radial_profiles_result = self.vizier.query_object(self.galaxy_name, catalog="J/ApJ/658/1006")

        distance = float(radial_profiles_result[0][0]["Dist"])
        inclination = Angle(float(radial_profiles_result[0][0]["i"]), "deg")

        # Set the inclination
        self.properties.inclination = inclination

    # -----------------------------------------------------------------

    def get_spiral_properties(self):

        """
        This function ...
        :return:
        """

        # http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/582/A86

        # J/A+A/582/A86: Catalogue of features in the S4G (Herrera-Endoqui+, 2015)

        # - J/A+A/582/A86/table2: Properties of bars, ring- and lens-structures in the S4G (2387 rows)
        # - J/A+A/582/A86/table3: Properties of spiral arms in the S4G (1854 rows)

        # Get table2
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/A+A/582/A86/table2"])
        table = result[0]

        # Name: Galaxy name
        # Class: Morphological classification
        # Type: Type of feature
        # sma: Semi-major axis [arcsec]
        # PA: Position angle [deg]
        # Ell: Ellipticity
        # smaEll: Semi-major axis from ellipticity [arcsec]
        # dsma: Deprojected semi-major axis [arcsec]
        # dPA: Deprojected position angle [deg]
        # dEll: Deprojected ellipticity
        # dsmaEll: Deprojexted semi-major axis from Ell [arcsec]
        # Qual: [1/3] Quality flag

        # Get table 3
        result = self.vizier.query_object(self.galaxy_name, catalog=["J/A+A/582/A86/table3"])
        table = result[0]

        # Name: Galaxy name
        # Class: Morphological classification
        # Type: Type of arms
        # Segment: Segment
        # Pitchang: Pitch angle [deg]
        # ri: Inner radius [arcsec]
        # ro: Outer radius [arcsec]
        # Qual: [1/2] Quality flag

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
        self.properties.save(self.galaxy_properties_path)

# -----------------------------------------------------------------
