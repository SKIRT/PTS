#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.calibrations.galametz Contains the GalametzTIRCalibration class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.tools import tables
from ...core.filter.filter import parse_filter, Filter
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# The path to the Galametz data directory
galametz_path = fs.join(introspection.pts_dat_dir("magic"), "Galametz")

# The patsh to the tables containing the single band data
single_band_luminosity_table_path = fs.join(galametz_path, "single_luminosity.dat")
single_band_brightness_table_path = fs.join(galametz_path, "single_brightness.dat")

# The path to the table containing the Galametz calibration parameters
multi_band_luminosity_table_path = fs.join(galametz_path, "multi_luminosity.dat")
multi_band_brightness_table_path = fs.join(galametz_path, "multi_brightness.dat")

# -----------------------------------------------------------------

multi_band_column_names = {"MIPS 24mu": "c24",
                         "Pacs blue": "c70",
                         "Pacs green": "c100",
                         "Pacs red": "c160",
                         "SPIRE PSW": "c250"}

# -----------------------------------------------------------------

class GalametzTIRCalibration(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Load the table with the single band data
        self.single_luminosity = tables.from_file(single_band_luminosity_table_path, format="ascii.commented_header")

        # Load the table with the single band data for surface brightness
        self.single_brightness = tables.from_file(single_band_brightness_table_path, format="ascii.commented_header")

        # Load the table with the multi band data
        self.multi_luminosity = tables.from_file(multi_band_luminosity_table_path, format="ascii.commented_header")

        # Load the table with the multi band data for surface brightness
        self.multi_brightness = tables.from_file(multi_band_brightness_table_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    @lazyproperty
    def single_band_filters(self):

        """
        This function ...
        :return: 
        """

        return [parse_filter(filter_name) for filter_name in self.single_luminosity["Band"]]

    # -----------------------------------------------------------------

    def index_for_filter_single(self, fltr):

        """
        This function ...
        :return: 
        """

        return self.single_band_filters.index(fltr)

    # -----------------------------------------------------------------

    def get_parameters_single_luminosity(self, fltr):

        """
        This function ...
        :parma fltr:
        :return: 
        """

        index = self.index_for_filter_single(fltr)
        return [self.single_luminosity["ai"][index], self.single_luminosity["bi"][index]]

    # -----------------------------------------------------------------

    def get_parameters_single_brightness(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        index = self.index_for_filter_single(fltr)
        return [self.single_brightness["ai"][index], self.single_brightness["bi"][index]]

    # -----------------------------------------------------------------

    def get_scatter_single_luminosity(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        index = self.index_for_filter_single(fltr)
        return self.single_luminosity["scatter"][index]

    # -----------------------------------------------------------------

    def get_scatter_single_brightness(self, fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        index = self.index_for_filter_single(fltr)
        return self.single_brightness["scatter"][index]

    # -----------------------------------------------------------------

    @lazyproperty
    def multi_band_filters(self):

        """
        This function ...
        :return: 
        """

        return [parse_filter(filter_name) for filter_name in multi_band_column_names]

    # -----------------------------------------------------------------

    def has_combination_multi_luminosity(self, *args):

        """
        This function ...
        :param args: 
        :return: 
        """

        # Get column names
        needed_column_names, not_needed_column_names = self.get_column_names_multi(*args)

        # Loop over the entries in the galametz table
        for i in range(len(self.multi_luminosity)):
            if is_appropriate_galametz_entry(self.multi_luminosity, i, needed_column_names, not_needed_column_names): return True

        # No row is found
        return False

    # -----------------------------------------------------------------

    def has_combination_multi_brightness(self, *args):

        """
        This function ...
        :param args: 
        :return: 
        """

        # Get column names
        needed_column_names, not_needed_column_names = self.get_column_names_multi(*args)

        # Loop over the entries in the table
        for i in range(len(self.multi_brightness)):
            if is_appropriate_galametz_entry(self.multi_brightness, i, needed_column_names, not_needed_column_names): return True

        # No row is found
        return False

    # -----------------------------------------------------------------

    def get_column_names_multi(self, *args):

        """
        This function ...
        :param args: 
        :return: 
        """

        # Needed column names
        needed_column_names = [] #[multi_band_column_names[filter_name] for filter_name in args]

        for filter_name in args:
            if isinstance(filter_name, Filter): filter_name = str(filter_name)
            name = multi_band_column_names[filter_name]
            needed_column_names.append(name)

        # List of not needed column names
        colnames = self.multi_luminosity.colnames
        not_needed_column_names = [name for name in colnames if name not in needed_column_names]

        not_needed_column_names.remove("R2")
        not_needed_column_names.remove("CV(RMSE)")

        # Return
        return needed_column_names, not_needed_column_names

    # -----------------------------------------------------------------

    def get_parameters_multi_luminosity(self, *args):

        """
        This function ...
        :param args:
        :return:
        """

        # Column names
        needed_column_names, not_needed_column_names = self.get_column_names_multi(*args)

        # The parameters
        parameters = None

        # Loop over the entries in the galametz table
        for i in range(len(self.multi_luminosity)):

            if is_appropriate_galametz_entry(self.multi_luminosity, i, needed_column_names, not_needed_column_names):

                parameters = []
                for name in needed_column_names: parameters.append(self.multi_luminosity[name][i])
                break

            else: continue

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def get_parameters_multi_brightness(self, *args):

        """
        This function ...
        :param args: 
        :return: 
        """

        # Column names
        needed_column_names, not_needed_column_names = self.get_column_names_multi(*args)

        # The parameters
        parameters = None

        # Loop over the entries in the table
        for i in range(len(self.multi_brightness)):

            if is_appropriate_galametz_entry(self.multi_brightness, i, needed_column_names, not_needed_column_names):

                parameters = []
                for name in needed_column_names: parameters.append(self.multi_brightness[name][i])
                break

            else: continue

        # Return the parameters
        return parameters

# -----------------------------------------------------------------

def is_appropriate_galametz_entry(table, i, needed_cols, not_needed_cols):

    """
    This function ...
    :return:
    """

    # Check if masked columns
    for name in not_needed_cols:

        # Verify that this entry is masked
        if not table[name].mask[i]: return False

    # Check needed cols
    for name in needed_cols:

        # Verify that this entry is not masked
        if table[name].mask[i]: return False

    # No mismatches found, thus appropriate entry
    return True

# -----------------------------------------------------------------
