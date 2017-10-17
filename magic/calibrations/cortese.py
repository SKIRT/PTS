#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.calibrations.cortese Contains the CorteseAttenuationCalibration class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.tools import introspection
from ...core.tools import tables
from ...core.filter.filter import parse_filter
from ...core.basics.range import RealRange
from ..tools import colours
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# The path to the Cortese data directory
cortese_path = fs.join(introspection.pts_dat_dir("magic"), "Cortese")

# The patsh to the tables containing the single band data
cortese_table_path = fs.join(cortese_path, "cortese.dat")

# -----------------------------------------------------------------

colour_combinations = {"FUV-H": ("GALEX FUV", "2MASS H"),
                       "FUV-i": ("GALEX FUV", "SDSS i"),
                       "FUV-r": ("GALEX FUV", "SDSS r"),
                       "FUV-g": ("GALEX FUV", "SDSS g"),
                       "FUV-B": ("GALEX FUV", "B")}

# -----------------------------------------------------------------

class CorteseAttenuationCalibration(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

        # Load the table with the single band data
        self.table = tables.from_file(cortese_table_path, format="ascii.commented_header")

    # -----------------------------------------------------------------

    @lazyproperty
    def ssfr_colours(self):

        """
        This function ...
        :return: 
        """

        colours = []
        for name in self.table.colnames:
            if name in ["Tau", "a1", "a2", "a3", "a4", "a5"]: continue
            colours.append(name)
        return colours

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return: 
        """

        filter_names = set()
        for colour in self.ssfr_colours:
            filter_names.add(colour.split("-")[0])
            filter_names.add(colour.split("-")[1])
        return [parse_filter(name) for name in filter_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def taus(self):

        """
        This function ...
        :return: 
        """

        return list(self.table["Tau"])

    # -----------------------------------------------------------------

    def __len__(self):

        """
        THis function ...
        :return: 
        """

        return len(self.table)

    # -----------------------------------------------------------------

    def get_standard_color_name(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour:
        :return:
        """

        # Convert color to standard format
        ssfr_colour = colours.get_colour_name_for_colour(ssfr_colour, short=True, delimiter="-")
        return ssfr_colour

    # -----------------------------------------------------------------

    def get_range_for_tau(self, ssfr_colour, tau):

        """
        This function ...
        :param ssfr_colour:
        :param tau: 
        :return: 
        """

        #ssfr_colour = self.get_standard_color_name(ssfr_colour)

        index = tables.find_index(self.table, tau)
        return self.get_range_for_index(ssfr_colour, index)

    # -----------------------------------------------------------------

    def get_range_for_index(self, ssfr_colour, index):

        """
        This function ...
        :param ssfr_colour: 
        :param index: 
        :return: 
        """

        ssfr_colour = self.get_standard_color_name(ssfr_colour)

        # Get the upper and lower limit
        upper = self.table[ssfr_colour][index]
        if index == len(self.table) - 1: lower = None
        else: lower = self.table[ssfr_colour][index + 1]

        # Return the range
        return RealRange(lower, upper)

    # -----------------------------------------------------------------

    def get_upper_limit(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour: 
        :return: 
        """

        # The absolute upper limit (so 10.5 for FUV-H, 7.5 for FUV-i, 7.3 for FUV-r, 6.7 for FUV-g, and 6.3 for FUV-B
        #absolute_upper_limit = limits[0][1]

        # Return the maximum of the range for the first tau (the first row), and for the corresponding colour column
        return self.get_range_for_index(ssfr_colour, 0).max

    # -----------------------------------------------------------------

    def get_parameters_for_tau(self, tau):

        """
        This function ...
        :return: 
        """

        index = tables.find_index(self.table, tau)
        return self.get_parameters_for_index(index)

    # -----------------------------------------------------------------

    def get_parameters_for_index(self, index):

        """
        This function ...
        :param index: 
        :return: 
        """

        a1 = self.table["a1"][index]
        a2 = self.table["a2"][index]
        a3 = self.table["a3"][index]
        a4 = self.table["a4"][index]
        a5 = self.table["a5"][index]

        return [a1, a2, a3, a4, a5]

    # -----------------------------------------------------------------

    def get_data(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour:
        :return: 
        """

        taus = []
        ranges = []
        parameters = []

        # Loop over all entries in the Cortese et. al
        for index in range(len(self.table)):

            # Get the tau value
            taus.append(self.taus[index])

            # Get the range
            ranges.append(self.get_range_for_index(ssfr_colour, index))

            # Get the parameters
            parameter_list = self.get_parameters_for_index(index)
            parameters.append(parameter_list)

        # Return the ranges and the parameters
        return taus, ranges, parameters

    # -----------------------------------------------------------------

    def taus_ranges_and_parameters(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour: 
        :return: 
        """

        # Get the data and iterate over the entries
        taus, ranges, parameters = self.get_data(ssfr_colour)
        #print(taus, len(taus))
        #print(ranges, len(ranges))
        #print(parameters, len(parameters))
        #print(len(self))
        for index in range(len(self)):
            yield taus[index], ranges[index], parameters[index]

    # -----------------------------------------------------------------

    def minimum_tau_range_and_parameters(self, ssfr_colour):

        """
        This function ...
        :param ssfr_colour:
        :return:
        """

        # Get the data and iterate over the entries
        taus, ranges, parameters = self.get_data(ssfr_colour)
        # print(taus, len(taus))
        # print(ranges, len(ranges))
        # print(parameters, len(parameters))
        # print(len(self))
        #for index in range(len(self)):
        return taus[0], ranges[0], parameters[0]

# -----------------------------------------------------------------
