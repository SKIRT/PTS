#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.absorption.sfr Contains the SFRAbsorption class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ....core.tools.utils import lazyproperty
from .simple import AbsorptionBase, seds_name, cells_name, best_name

# -----------------------------------------------------------------

class UnevolvedAbsorption(AbsorptionBase):

    """
    This class ...
    """

    def __init__(self, simulations, young, sfr, absorption_curve_cells=None, emission_curve_cells=None, distance=None):

        """
        This function ...
        :param simulations:
        :param young:
        :param sfr:
        :param absorption_curve_cells:
        :param emission_curve_cells:
        :param distance:
        """

        # Call the constructor of the base class
        super(UnevolvedAbsorption, self).__init__(simulations, absorption_curve_cells=absorption_curve_cells,
                                                  emission_curve_cells=emission_curve_cells, distance=distance)

        # Young and star formation regions absorption
        self.young = young
        self.sfr = sfr

    # -----------------------------------------------------------------

    @property
    def young_simulations(self):
        return self.young.simulations

    # -----------------------------------------------------------------

    @property
    def sfr_simulations(self):
        return self.sfr.simulations

    # -----------------------------------------------------------------
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
        return self.correct_observed_stellar_sed(self.young.observed_stellar_sed + self.sfr.observed_stellar_sed_all)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity_all(self):
        return self.observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed_diffuse(self):
        return self.young.intrinsic_stellar_sed + self.sfr.intrinsic_stellar_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_stellar_sed_all(self):
        return self.young.intrinsic_stellar_sed + self.sfr.intrinsic_stellar_sed_all

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
        #return self.model.unevolved_simulations.observed_sed_absorbed
        return self.correct_absorption_sed(self.young.absorption_sed + self.sfr.absorption_sed_diffuse)

    # -----------------------------------------------------------------

    @property
    def has_absorption_sed_diffuse(self):
        #return self.model.unevolved_simulations.has_observed_sed_absorbed
        return self.young.has_absorption_sed and self.sfr.has_absorption_sed_diffuse

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
    # Internal
    # -----------------------------------------------------------------

    @property
    def absorption_sed_internal(self):
        return self.sfr.absorption_sed_internal

    # -----------------------------------------------------------------

    @property
    def absorption_luminosity_internal(self):
        return self.sfr.absorption_luminosity_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_internal(self):
        return self.absorption_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def dust_sed_internal(self):
        return self.sfr.dust_sed_internal

    # -----------------------------------------------------------------

    @property
    def dust_luminosity_internal(self):
        return self.sfr.dust_luminosity_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_internal(self):
        return self.dust_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_all(self):
        return self.correct_absorption_sed(self.absorption_sed_diffuse + self.sfr.absorption_sed_internal)

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
        return self.correct_dust_sed(self.dust_sed_diffuse + self.sfr.dust_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_all(self):
        return self.dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_all(self):
        return self.dust_luminosity_all.value / self.stellar_luminosity.value

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
        else: return self.absorption_sed_diffuse  # seds

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
        else: return self.dust_sed_diffuse  # seds

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_diffuse(self):
        return self.best_dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction_diffuse(self):
        return self.best_dust_luminosity_diffuse / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # ALL
    # -----------------------------------------------------------------

    @lazyproperty
    def best_observed_stellar_sed_all(self):
        return self.correct_observed_stellar_sed(self.young.best_observed_stellar_sed + self.sfr.best_observed_stellar_sed_all)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_sed_all(self):
        return self.correct_absorption_sed(self.best_absorption_sed_diffuse + self.sfr.absorption_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_luminosity_all(self):
        return self.best_absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, unit=self.specific_luminosity_unit, distance=self.distance)

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
        return self.correct_dust_sed(self.best_dust_sed_diffuse + self.sfr.dust_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_all(self):
        return self.best_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, unit=self.specific_luminosity_unit, distance=self.distance)

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
        seds[best_name] = self.best_dust_sed_all
        return seds

# -----------------------------------------------------------------
