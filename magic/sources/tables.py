#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.sources.tables Contains table classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.table import SmartTable
from ...core.basics.curve import FilterCurve
from ...core.units.parsing import parse_unit as u
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

class FWHMTable(FilterCurve):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Set properties
        kwargs["y_name"] = "FWHM"
        kwargs["y_description"] = "FWHM of the PSF"
        kwargs["y_unit"] = "arcsec"

        # Call the constructor of the base class
        super(FWHMTable, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def add_fwhm(self, fltr, fwhm):

        """
        This function ...
        :param fltr:
        :param fwhm:
        :return:
        """

        self.add_point(fltr, fwhm)

    # -----------------------------------------------------------------

    def fwhm_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.value_for_filter(fltr)

# -----------------------------------------------------------------

class GalaxyTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check
        if "filters" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy: filters = kwargs.pop("filters")
        else: filters = None

        # Call the constructor of the base class
        super(GalaxyTable, self).__init__(*args, **kwargs)

        # Add column info
        if not from_astropy:

            # Add columns
            self.add_column_info("Index", int, None, "index of the extended source in the catalog")
            self.add_column_info("Name", str, None, "name of the galaxy")

            for fltr in filters:

                column_name = str(fltr) + " flux"
                self.add_column_info(column_name, float, u("Jy"), str(fltr) + " flux density")

    # -----------------------------------------------------------------

    def add_galaxy(self, galaxy):

        """
        This function ...
        :param galaxy:
        :return:
        """

        # Setup if necessary
        if len(self.colnames) == 0: self._setup()

        values = []

        index = galaxy.index
        name = galaxy.name

        # Add index and name
        values.append(index)
        values.append(name)

        # Loop over the filters for which we need a flux
        for name in self.colnames:

            # Skip
            if not name.endswith("flux"): continue

            # Filter
            #fltr = BroadBandFilter(name.split(" flux")[0])
            fltr = parse_filter(name.split(" flux")[0])

            # Get flux
            if galaxy.sed is not None and fltr in galaxy.sed.filters(): flux = galaxy.sed.photometry_for_filter(fltr)
            else: flux = None

            # Add the flux to the values
            values.append(flux)

        # Add a row to the table
        self.add_row(values)

# -----------------------------------------------------------------

class StarTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Check
        if "filters" in kwargs: from_astropy = False
        else: from_astropy = True

        # Get properties
        if not from_astropy: filters = kwargs.pop("filters")
        else: filters = None

        # Call the constructor of the base class
        super(StarTable, self).__init__(*args, **kwargs)

        # Add column info
        if not from_astropy:

            self.add_column_info("Index", int, None, "index of the point source in the catalog")
            self.add_column_info("Catalog", str, None, "original catalog")
            self.add_column_info("ID", str, None, "ID of the point source in the original catalog")

            # Loop over the filters
            for fltr in filters:

                column_name = str(fltr) + " FWHM"
                self.add_column_info(column_name, float, u("arcsec"), str(fltr) + " FWHM")

            # Loop over the filters
            for fltr in filters:

                column_name = str(fltr) + " flux"
                self.add_column_info(column_name, float, u("Jy"), str(fltr) + " flux density")

    # -----------------------------------------------------------------

    def add_star(self, star):

        """
        This function ...
        :param star:
        :return:
        """

        if len(self.colnames) == 0: self._setup()

        values = []

        catalog = star.catalog
        id = star.id

        # Add index, catalog and ID
        values.append(star.index)
        values.append(catalog)
        values.append(id)

        # Loop over the filters for which we need a FWHM
        for name in self.colnames:

            if name == "Index": continue
            if name == "Catalog": continue
            if name == "ID": continue

            # FWHM
            if name.endswith("FWHM"):

                filter_name = name.split(" FWHM")[0]

                # Filter
                fltr = parse_filter(filter_name)
                #filter_name = str(fltr)

                #print(star.fwhms)
                if star.fwhms.has_filter(fltr): fwhm = star.fwhms.fwhm_for_filter(fltr)
                #if filter_name in star.fwhms: fwhm = star.fwhms[filter_name]
                else: fwhm = None

                values.append(fwhm)

            # Flux
            elif name.endswith("flux"):

                # Filter
                #fltr = BroadBandFilter(name.split(" flux")[0])
                fltr = parse_filter(name.split(" flux")[0])

                #print(star.sed)
                #print(fltr)

                # Get flux
                flux = star.sed.photometry_for_filter(fltr)

                # Add the flux to the values
                values.append(flux)

            # Unknown
            else: raise ValueError("Don't know what value to fill in for column '" + name + "'")

        # Add the row
        self.add_row(values)

# -----------------------------------------------------------------
