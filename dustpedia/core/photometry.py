#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.dustpedia.core.photometry Contains the DustPediaPhotometry class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.table import Table
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools import tables
from ...core.filter.filter import parse_filter
from ...core.data.sed import ObservedSED
from .sample import DustPediaSample
from ...magic.region.ellipse import SkyEllipseRegion
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.basics.stretch import SkyStretch
from ...core.units.unit import parse_unit as u

# -----------------------------------------------------------------

dustpedia_dat_path = fs.join(introspection.pts_dat_dir("dustpedia"))

# -----------------------------------------------------------------

dustpedia_photometry_path = fs.join(dustpedia_dat_path, "photometry")

# -----------------------------------------------------------------

ap_phot_table_path = fs.join(dustpedia_photometry_path, "DustPedia_Aperture_Photometry.csv")
iras_scanpi_phot_table_path = fs.join(dustpedia_photometry_path, "DustPedia_IRAS_SCANPI.csv")
planck_phot_table_path = fs.join(dustpedia_photometry_path, "DustPedia_Planck_CCS2.csv")

# -----------------------------------------------------------------

# CHRIS:
#
# I'm very pleased to say that the DustPedia photometry is now ready for everyone to enjoy! The files are attached.
#
# There is a short(ish) PDF describing the photometry, along with 3 tables. One main table of aperture-matched
# photometry, one table of supplementary Planck catalogue photometry, and one table of supplementary IRAS photometry.
#
# Also, I have transferred to Manolis the PNG image grids (described in the PDF) showing for each source the
# apertures in all bands for which aperture photometry is performed; these should hopefully be accessible to everyone
# on the database soon.
#
# Please feel free to send me any questions you have. I wouldn't be surprised if there are a few small issues
# still hidden in the photometry, that my checks haven't found. So please report back anything odd you might find!

# -----------------------------------------------------------------

leda_search_object_url = "http://leda.univ-lyon1.fr/ledacat.cgi?"

# -----------------------------------------------------------------

non_flux_columns = ["name", "ra", "dec", "semimaj_arcsec", "axial_ratio", "pos_angle", "global_flag"]

# -----------------------------------------------------------------

class DustPediaPhotometry(object):

    """
    This class ...
    """

    def __init__(self):

        """
        The constructor ...
        """

        # Read the tables
        self.aperture = Table.read(ap_phot_table_path)
        self.iras_scanpi = Table.read(iras_scanpi_phot_table_path)
        self.planck = Table.read(planck_phot_table_path)

        # Initialize dictionaries
        self.aperture_filters = dict()
        self.iras_filters = dict()
        self.planck_filters = dict()

        # Set the filters for which there is photometry
        self.set_filters()

        # Create the DustPedia sample object
        self.sample = DustPediaSample()

    # -----------------------------------------------------------------

    def set_filters(self):

        """
        This function ...
        :return:
        """

        # APERTURE
        for colname in self.aperture.colnames:

            if colname in non_flux_columns: continue
            if colname.endswith("_err") or colname.endswith("_flag"): continue
            fltr = parse_filter(colname)
            self.aperture_filters[colname] = fltr

        # IRAS
        for colname in self.iras_scanpi.colnames:
            if colname in non_flux_columns: continue
            if colname.endswith("_err") or colname.endswith("_flag"): continue
            fltr = parse_filter(colname)
            self.iras_filters[colname] = fltr

        # PLANCK
        for colname in self.planck.colnames:
            if colname in non_flux_columns: continue
            if colname.endswith("_err") or colname.endswith("_flag"): continue
            fltr = parse_filter(colname)
            self.planck_filters[colname] = fltr

    # -----------------------------------------------------------------

    def get_aperture(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Inform the user
        log.info("Getting the aperture for galaxy '" + galaxy_name + "' ...")

        # Get DustPedia name of the galaxy
        objname = self.sample.get_name(galaxy_name)

        # Find in table
        index = tables.find_index(self.aperture, objname, "name")

        # Get properties
        ra = self.aperture["ra"][index]
        dec = self.aperture["dec"][index]
        semimaj_arcsec = self.aperture["semimaj_arcsec"][index]
        axial_ratio = self.aperture["axial_ratio"][index]
        pos_angle = self.aperture["pos_angle"][index]

        # Create center
        center = SkyCoordinate(ra=ra, dec=dec, unit="deg")

        semimaj = semimaj_arcsec * u("arcsec")

        # Create radius
        radius = SkyStretch(semimaj, semimaj/axial_ratio)

        # Create the angle
        angle = Angle(pos_angle, "deg")

        # Create aperture
        aperture = SkyEllipseRegion(center, radius, angle=angle)

        # Return the aperture
        return aperture

    # -----------------------------------------------------------------

    def get_flags(self, galaxy_name):

        """
        This function ...
        :param galaxy_name:
        :return:
        """

        # Inform the user
        log.info("Getting the photometry flags for galaxy '" + galaxy_name + "' ...")

        # Get DustPedia name of galaxy
        objname = self.sample.get_name(galaxy_name)

        # Find in table
        index = tables.find_index(self.aperture, objname, "name")

        # Initialize dictionary for the flags for each filter
        flags = dict()

        # Loop over the filters for which there is aperture photometry
        for colname in self.aperture_filters:

            flag_colname = colname + "_flag"

            # Masked entry
            if hasattr(self.aperture[flag_colname], "mask") and self.aperture[flag_colname].mask[index]: continue

            # Get the flag
            flag = self.aperture[flag_colname][index]

            # Set the flag for this filter
            flags[colname] = flag

        # Return the flags
        return flags

    # -----------------------------------------------------------------

    def get_sed(self, galaxy_name, add_iras=True, add_planck=True, exact_name=False):
            
        """
        This function ...
        :param galaxy_name:
        :param add_iras:
        :param add_planck:
        :param exact_name:
        """

        # Inform the user
        log.info("Getting the DustPedia SED for galaxy '" + galaxy_name + "' ...")

        # Get DustPedia name of galaxy
        if exact_name: objname = galaxy_name
        else: objname = self.sample.get_name(galaxy_name)

        # Find in table
        index = tables.find_index(self.aperture, objname, "name")

        filters = []
        fluxes = []
        errors = []

        for colname in self.aperture_filters:

            # Masked entry
            if hasattr(self.aperture[colname], "mask") and self.aperture[colname].mask[index]: continue

            flux = self.aperture[colname][index]

            filters.append(self.aperture_filters[colname])
            fluxes.append(flux)

            if colname + "_err" in self.aperture.colnames: errors.append(self.aperture[colname + "_err"][index])
            else: errors.append(None)

        index = tables.find_index(self.iras_scanpi, objname, "name")

        # Add IRAS
        if add_iras:

            for colname in self.iras_filters:

                # Masked entry
                if hasattr(self.iras_scanpi[colname], "mask") and self.iras_scanpi[colname].mask[index]: continue

                flux = self.iras_scanpi[colname][index]

                # Filter and flux
                filters.append(self.iras_filters[colname])
                fluxes.append(flux)

                # Add error
                if colname + "_err" in self.iras_scanpi.colnames: errors.append(self.iras_scanpi[colname + "_err"][index])
                else: errors.append(None)

        index = tables.find_index(self.planck, objname, "name")

        # Add Planck fluxes
        if add_planck:

            for colname in self.planck_filters:

                # Masked entry
                if hasattr(self.planck[colname], "mask") and self.planck[colname].mask[index]: continue

                flux = self.planck[colname][index]

                filters.append(self.planck_filters[colname])
                fluxes.append(flux)

                if colname + "_err" in self.planck.colnames: errors.append(self.planck[colname + "_err"][index])
                else: errors.append(None)

        # Create the observed SED
        sed = ObservedSED(photometry_unit="Jy")

        # Add the points
        for i in range(len(filters)):
            fltr = filters[i]
            flux = fluxes[i]
            error = errors[i]
            if error < 0: error = None
            sed.add_point(fltr, flux, error)

        # Return the SED
        return sed

# -----------------------------------------------------------------
