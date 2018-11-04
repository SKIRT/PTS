#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.absorption.total Contains the TotalAbsorption class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ....core.tools.utils import lazyproperty
from .simple import AbsorptionBase, seds_name, cells_name, seds_alt_name, best_name

# -----------------------------------------------------------------

class TotalAbsorption(AbsorptionBase):

    """
    This class ...
    """

    def __init__(self, simulations, evolved, unevolved, extra=None, absorption_curve_cells=None, emission_curve_cells=None, distance=None):

        """
        This function ...
        :param simulations:
        :param evolved:
        :param unevolved:
        :param absorption_curve_cells:
        :param distance:
        """

        # Call the constructor of the base class
        super(TotalAbsorption, self).__init__(simulations, absorption_curve_cells=absorption_curve_cells, emission_curve_cells=emission_curve_cells, distance=distance)

        # Evolved and unevolved, and star formation
        self.evolved = evolved
        self.unevolved = unevolved
        self.extra = extra

    # -----------------------------------------------------------------

    @property
    def evolved_simulations(self):
        return self.evolved.simulations

    # -----------------------------------------------------------------

    @property
    def unevolved_simulations(self):
        return self.unevolved.simulations

    # -----------------------------------------------------------------

    @property
    def extra_simulations(self):
        return self.extra.simulations

    # -----------------------------------------------------------------

    @property
    def has_extra(self):
        return self.extra is not None

    # -----------------------------------------------------------------
    # STELLAR SEDs and LUMINOSITIES
    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed_diffuse(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_stellar_sed, extrapolate=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity_diffuse(self):
        return self.observed_stellar_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed_all(self):
        sed = self.evolved.observed_stellar_sed + self.unevolved.observed_stellar_sed_all
        if self.has_extra: sed = sed + self.extra.observed_stellar_sed
        return self.correct_observed_stellar_sed(sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity_all(self):
        return self.observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed_diffuse(self):
        sed = self.evolved.intrinsic_stellar_sed + self.unevolved.intrinsic_stellar_sed_diffuse
        if self.has_extra: sed = sed + self.extra.intrinsic_stellar_sed
        return sed
        # or self.model.total_simulations.intrinsic_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed_all(self):
        sed = self.evolved.intrinsic_stellar_sed + self.unevolved.intrinsic_stellar_sed_all
        if self.has_extra: sed = sed + self.extra.intrinsic_stellar_sed
        return sed

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity(self):
        return self.intrinsic_stellar_sed_all.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------
    # ABSORPTION SEDs AND LUMINOSITIES
    #   Diffuse
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_diffuse_cells(self):
        return self.correct_absorption_sed(self.absorption_curve_cells)

    # -----------------------------------------------------------------

    @property
    def has_absorption_sed_diffuse_cells(self):
        return self.has_absorption_curve_cells

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity_diffuse_cells(self):
        return self.absorption_sed_diffuse_cells.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_diffuse_cells(self):
        return self.absorption_luminosity_diffuse_cells.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_diffuse(self):
        return self.correct_absorption_sed(self.simulations.observed_sed_absorbed)

    # -----------------------------------------------------------------

    @property
    def has_absorption_sed_diffuse(self):
        return self.simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity_diffuse(self):
        return self.absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def has_absorption_luminosity_diffuse(self):
        return self.has_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_diffuse(self):
        return self.absorption_luminosity_diffuse.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_absorption_fraction_diffuse(self):
        return self.has_absorption_luminosity_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_diffuse_cells_complete(self):
        return self.correct_dust_sed(self.emission_curve_cells, trim=False, make_full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_diffuse_cells(self):
        return self.correct_dust_sed(self.emission_curve_cells)

    # -----------------------------------------------------------------

    @property
    def has_dust_sed_diffuse_cells(self):
        return self.has_emission_curve_cells

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_diffuse_complete(self):
        return self.correct_dust_sed(self.simulations.observed_diffuse_dust_sed, trim=False)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_diffuse(self):
        return self.correct_dust_sed(self.simulations.observed_diffuse_dust_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_diffuse(self):
        return self.dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_diffuse(self):
        return self.dust_luminosity_diffuse.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_all(self):
        return self.correct_absorption_sed(self.absorption_sed_diffuse + self.unevolved.absorption_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity_all(self):
        return self.absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_all(self):
        return self.absorption_luminosity_all.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_all(self):
        return self.correct_dust_sed(self.dust_sed_diffuse + self.unevolved.dust_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_all(self):
        return self.dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_all(self):
        return self.dust_luminosity_all.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_all_alt(self):
        return self.correct_dust_sed(self.simulations.observed_dust_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_all_alt(self):
        return self.dust_sed_all_alt.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_all_alt(self):
        return self.dust_luminosity_all_alt.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # DIFFUSE
    # -----------------------------------------------------------------

    @lazyproperty
    def best_observed_stellar_sed_diffuse(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_sed - self.best_dust_sed_diffuse_complete, extrapolate=False)

    # -----------------------------------------------------------------

    @property
    def best_absorption_sed_diffuse(self):
        if self.has_absorption_sed_diffuse_cells: return self.absorption_sed_diffuse_cells
        else: return self.absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_luminosity_diffuse(self):
        return self.best_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_fraction_diffuse(self):
        return self.best_absorption_luminosity_diffuse.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_luminosity_diffuse(self):
        return self.best_absorption_sed_diffuse.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_fraction_diffuse(self):
        return self.best_fuv_absorption_luminosity_diffuse.value / self.intrinsic_fuv_luminosity.value

    # -----------------------------------------------------------------

    @property
    def best_dust_sed_diffuse_complete(self):
        if self.has_dust_sed_diffuse_cells: return self.dust_sed_diffuse_cells_complete
        else: return self.dust_sed_diffuse_complete

    # -----------------------------------------------------------------

    @property
    def best_dust_sed_diffuse(self):
        if self.has_dust_sed_diffuse_cells: return self.dust_sed_diffuse_cells
        else: return self.dust_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_diffuse(self):
        return self.best_dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction_diffuse(self):
        return self.best_dust_luminosity_diffuse.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # ALL
    # -----------------------------------------------------------------

    @lazyproperty
    def best_observed_stellar_sed_all(self):
        sed = self.evolved.best_observed_stellar_sed + self.unevolved.best_observed_stellar_sed_all
        if self.has_extra: sed = sed + self.extra.best_observed_stellar_sed
        return self.correct_observed_stellar_sed(sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_sed_all(self):
        return self.correct_absorption_sed(self.best_absorption_sed_diffuse + self.unevolved.absorption_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_luminosity_all(self):
        return self.best_absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_fraction_all(self):
        return self.best_absorption_luminosity_all.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_luminosity_all(self):
        return self.best_absorption_sed_all.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_fraction_all(self):
        return self.best_fuv_absorption_luminosity_all.value / self.intrinsic_fuv_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_sed_all(self):
        return self.correct_dust_sed(self.best_dust_sed_diffuse + self.unevolved.dust_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_all(self):
        return self.best_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction_all(self):
        return self.best_dust_luminosity_all.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # DIFFUSE
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_seds_diffuse(self):
        seds = OrderedDict()
        if self.has_absorption_sed_diffuse: seds[seds_name] = self.absorption_sed_diffuse
        if self.has_absorption_sed_diffuse_cells: seds[cells_name] = self.absorption_sed_diffuse_cells
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def emission_seds_diffuse(self):
        seds = OrderedDict()
        seds[seds_name] = self.dust_sed_diffuse
        if self.has_dust_sed_diffuse_cells: seds[cells_name] = self.dust_sed_diffuse_cells
        return seds

    # -----------------------------------------------------------------
    # ALL
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_seds_all(self):
        seds = OrderedDict()
        seds[seds_name] = self.absorption_sed_all
        seds[best_name] = self.best_absorption_sed_all
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def emission_seds_all(self):
        seds = OrderedDict()
        seds[seds_name] = self.dust_sed_all
        seds[seds_alt_name] = self.dust_sed_all_alt
        seds[best_name] = self.best_dust_sed_all
        return seds

# -----------------------------------------------------------------
