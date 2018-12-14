#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.absorption.simple Contains the SimpleAbsorption class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from abc import ABCMeta
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ....core.tools.utils import lazyproperty
from ....core.units.parsing import parse_unit as u
from ....core.units.parsing import parse_quantity as q
from ....core.filter.broad import BroadBandFilter

# -----------------------------------------------------------------

cells_name = "cells"
seds_name = "seds"
seds_alt_name = "seds_alt"
best_name = "best"

# -----------------------------------------------------------------

class AbsorptionBase(object):

    """
    This class ...
    """

    # -----------------------------------------------------------------

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, simulations, absorption_curve_cells=None, emission_curve_cells=None, distance=None):

        """
        This function ...
        :param simulations:
        :param absorption_curve_cells:
        :param emission_curve_cells:
        :param distance:
        """

        # Simulations
        self.simulations = simulations

        # Curves from cells
        self.absorption_curve_cells = absorption_curve_cells
        self.emission_curve_cells = emission_curve_cells

        # Set distance
        if distance is not None: self.distance = distance

    # -----------------------------------------------------------------

    @lazyproperty
    def distance(self):
        return self.simulations.distance

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        return self.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_absorption_curve_cells(self):
        return self.absorption_curve_cells is not None

    # -----------------------------------------------------------------

    @property
    def has_emission_curve_cells(self):
        return self.emission_curve_cells is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def specific_luminosity_unit(self):
        return u("W/micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def frequency_luminosity_unit(self):
        return u("W/Hz")

    # -----------------------------------------------------------------

    @lazyproperty
    def bolometric_luminosity_unit(self):
        return u("W")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed_extrapolate_from(self):
        return q("10 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed_fit_from(self):
        return q("4 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def maximum_wavelength(self):
        return q("2000 micron")

    # -----------------------------------------------------------------

    def correct_observed_stellar_sed(self, sed, extrapolate=True):

        # Extrapolate?
        if extrapolate:
            # Get maximum wavelength below 10 micron for which above zero
            start_wavelength = sed.get_max_positive_wavelength(upper=self.observed_stellar_sed_extrapolate_from)
            sed = sed.extrapolated_from(start_wavelength, regression_from_x=self.observed_stellar_sed_fit_from, xlog=True, ylog=True, replace_nan=0.)
        else: sed = sed.copy()

        # Extend and replace invalid values
        sed.set_negatives_to_zero()
        sed = sed.extended_to_right(self.maximum_wavelength, logscale=True, points=self.wavelengths)
        sed.replace_zeros_by_lowest(0.01) # replace by one hundreth of the minimum value (except zero)
        sed.distance = self.distance

        # Return the correct SED
        return sed

    # -----------------------------------------------------------------

    @lazyproperty
    def min_dust_wavelength(self):
        return q("1 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def minimum_wavelength(self):
        return q("0.01 micron")

    # -----------------------------------------------------------------

    def correct_dust_sed(self, sed, trim=True, make_full=False):
        if trim: sed = sed.stripped_negatives_and_zeroes()
        else:
            sed = sed.copy()
            sed.set_negatives_to_zero()
        if trim: sed = sed.splice_right(self.min_dust_wavelength)
        if make_full: sed = sed.extended_to_left(self.minimum_wavelength, logscale=True, points=self.wavelengths)
        sed.distance = self.distance
        return sed

    # -----------------------------------------------------------------

    @lazyproperty
    def max_absorption_wavelength(self):
        return q("2 micron")

    # -----------------------------------------------------------------

    def correct_absorption_sed(self, sed, trim=True):
        sed = sed.stripped_negatives_and_zeroes()
        if trim: sed = sed.splice_left(self.max_absorption_wavelength)
        sed.distance = self.distance
        return sed

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_filter(self):
        return BroadBandFilter("GALEX FUV")

    # -----------------------------------------------------------------

    @lazyproperty
    def fuv_wavelength(self):
        return self.fuv_filter.wavelength

    # -----------------------------------------------------------------

    @property
    def observed_sed(self):
        return self.simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelengths(self):
        return self.observed_sed.wavelengths(unit="micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_luminosity(self):
        return self.observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def stellar_sed(self):
        return self.simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_luminosity(self):
        return self.stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

# -----------------------------------------------------------------

class SimpleAbsorption(AbsorptionBase):

    """
    This class ...
    """

    @lazyproperty
    def observed_stellar_sed(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_stellar_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity(self):
        return self.observed_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed(self):
        return self.simulations.intrinsic_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity(self):
        return self.intrinsic_stellar_sed.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------
    # ABSORPTION SEDs AND LUMINOSITIES
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_cells(self):
        return self.correct_absorption_sed(self.absorption_curve_cells)

    # -----------------------------------------------------------------

    @property
    def has_absorption_sed_cells(self):
        return self.has_absorption_curve_cells

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity_cells(self):
        return self.absorption_sed_cells.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_cells(self):
        return self.absorption_luminosity_cells.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed(self):
        return self.correct_absorption_sed(self.simulations.observed_sed_absorbed)

    # -----------------------------------------------------------------

    @property
    def has_absorption_sed(self):
        return self.simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity(self):
        return self.absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_absorption_luminosity(self):
        return self.has_absorption_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction(self):
        return self.absorption_luminosity.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_absorption_fraction(self):
        return self.has_absorption_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_cells_complete(self):
        return self.correct_dust_sed(self.emission_curve_cells, trim=False, make_full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_cells(self):
        return self.correct_dust_sed(self.emission_curve_cells)

    # -----------------------------------------------------------------

    @property
    def has_dust_sed_cells(self):
        return self.has_emission_curve_cells

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_complete(self):
        return self.correct_dust_sed(self.simulations.observed_diffuse_dust_sed, trim=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed(self):
        return self.correct_dust_sed(self.simulations.observed_diffuse_dust_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity(self):
        return self.dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction(self):
        return self.dust_luminosity.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def best_observed_stellar_sed(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_sed - self.best_dust_sed_complete)

    # -----------------------------------------------------------------

    @property
    def best_absorption_sed(self):
        if self.has_absorption_sed_cells: return self.absorption_sed_cells # cells
        else: return self.absorption_sed # seds

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_luminosity(self):
        return self.best_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_fraction(self):
        return self.best_absorption_luminosity.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_luminosity(self):
        return self.best_absorption_sed.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_fraction(self):
        return self.best_fuv_absorption_luminosity.value / self.intrinsic_fuv_luminosity.value

    # -----------------------------------------------------------------

    @property
    def best_dust_sed_complete(self):
        if self.has_dust_sed_cells: return self.dust_sed_cells_complete
        else: return self.dust_sed_complete

    # -----------------------------------------------------------------

    @property
    def best_dust_sed(self):
        if self.has_dust_sed_cells: return self.dust_sed_cells
        else: return self.dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity(self):
        return self.best_dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction(self):
        return self.best_dust_luminosity.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_seds(self):
        seds = OrderedDict()
        if self.has_absorption_sed: seds[seds_name] = self.absorption_sed
        if self.has_absorption_sed_cells: seds[cells_name] = self.absorption_sed_cells
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def emission_seds(self):
        seds = OrderedDict()
        seds[seds_name] = self.dust_sed
        if self.has_dust_sed_cells: seds[cells_name] = self.dust_sed_cells
        return seds

# -----------------------------------------------------------------
