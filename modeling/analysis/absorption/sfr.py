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
from .simple import AbsorptionBase, cells_name, seds_name, seds_alt_name, best_name

# -----------------------------------------------------------------

class SFRAbsorption(AbsorptionBase):

    """
    This class ...
    """

    def __init__(self, simulations, unattenuated_sed, internal_attenuation_sed, intrinsic_dust_sed,
                 absorption_curve_cells=None, emission_curve_cells=None, distance=None):

        """
        This function ...
        :param simulations:
        :param unattenuated_sed:
        :param internal_attenuation_sed:
        :param intrinsic_dust_sed:
        :param absorption_curve_cells:
        :param emission_curve_cells:
        :param distance:
        """

        # Call the constructor of the base class
        super(SFRAbsorption, self).__init__(simulations, absorption_curve_cells=absorption_curve_cells, emission_curve_cells=emission_curve_cells, distance=distance)

        # Set SED of unattenuated SFR and internal attenuation SED
        self.unattenuated_sed = unattenuated_sed
        self.internal_attenuation_sed = internal_attenuation_sed
        self.intrinsic_dust_sed = intrinsic_dust_sed

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
    def observed_stellar_sed_internal(self):
        return self.correct_observed_stellar_sed(self.stellar_sed - self.dust_sed_internal_complete)  # intrinsic here means unaffected by diffuse dust, but still with internal extinction

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity_internal(self):
        return self.observed_stellar_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_sed_all(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_sed - self.dust_sed_all_complete)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_stellar_luminosity_all(self):
        return self.observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_sed_diffuse(self):
        return self.stellar_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_sed_internal(self):
        return self.unattenuated_sed

    # -----------------------------------------------------------------

    @property
    def intrinsic_stellar_sed_all(self):
        return self.unattenuated_sed

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
    #   Internal
    # -----------------------------------------------------------------

    # INCORRECT DUE TO INCORRECT TRANSPARENT SFR SED
    # @lazyproperty
    # def sfr_absorption_sed_internal(self):
    #
    #     # Get stellar SEDs
    #     #intrinsic_stellar = self.model.get_stellar_sed("sfr", "intrinsic")
    #     #intrinsic_stellar = self.model.intrinsic_sfr_sed # same
    #     intrinsic_stellar = self.model.intrinsic_sfr_stellar_sed
    #     # INCORRECT:
    #     transparent_stellar = self.model.intrinsic_transparent_sfr_stellar_sed # NEW FROM TRANSPARENT MAPPINGS
    #
    #     # INTERNALLY ABSORBED
    #     return self.correct_absorption_sed(transparent_stellar - intrinsic_stellar)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_internal(self):
        return self.correct_absorption_sed(self.internal_attenuation_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_luminosity_internal(self):
        return self.absorption_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_fraction_internal(self):
        return self.absorption_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_internal_complete(self):
        return self.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_internal_complete(self):
        return self.dust_sed_internal_complete.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_internal(self):
        return self.correct_dust_sed(self.intrinsic_dust_sed) # NEW FROM TRANSPARENT MAPPINGS, SUBTRACTED

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_internal(self):
        return self.dust_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_internal(self):
        return self.dust_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_internal_alt(self):
        return self.correct_dust_sed(self.simulations.intrinsic_dust_sed)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_luminosity_internal_alt(self):
        return self.dust_sed_internal_alt.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_fraction_internal_alt(self):
        return self.dust_luminosity_internal_alt.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_sed_all(self):
        return self.correct_absorption_sed(self.absorption_sed_diffuse + self.absorption_sed_internal)

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
    def dust_sed_all_complete(self):
        return self.dust_sed_diffuse_complete + self.dust_sed_internal_complete

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_sed_all(self):
        return self.correct_dust_sed(self.dust_sed_diffuse + self.dust_sed_internal)

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
    # INTERNAL
    # -----------------------------------------------------------------

    @property
    def best_observed_stellar_sed_internal(self):
        return self.observed_stellar_sed_internal

    # -----------------------------------------------------------------

    @property
    def best_absorption_sed_internal(self):
        return self.absorption_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_luminosity_internal(self):
        return self.best_absorption_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_fraction_internal(self):
        return self.best_absorption_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_luminosity_internal(self):
        return self.best_absorption_sed_internal.photometry_at(self.fuv_wavelength, unit=self.specific_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_fuv_absorption_fraction_internal(self):
        return self.best_fuv_absorption_luminosity_internal.value / self.intrinsic_fuv_luminosity.value

    # -----------------------------------------------------------------

    @property
    def best_dust_sed_internal(self):
        return self.dust_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_internal(self):
        return self.best_dust_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction_internal(self):
        return self.best_dust_luminosity_internal.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
    # ALL
    # -----------------------------------------------------------------

    @lazyproperty
    def best_observed_stellar_sed_all(self):
        return self.correct_observed_stellar_sed(self.simulations.observed_sed - self.best_dust_sed_all_complete)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_absorption_sed_all(self):
        return self.correct_absorption_sed(self.best_absorption_sed_diffuse + self.absorption_sed_internal)

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
    def best_dust_sed_all_complete(self):
        return self.correct_dust_sed(self.best_dust_sed_diffuse_complete + self.dust_sed_internal_complete, trim=False, make_full=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_sed_all(self):
        return self.correct_dust_sed(self.best_dust_sed_diffuse + self.dust_sed_internal)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_luminosity_all(self):
        return self.best_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def best_dust_fraction_all(self):
        return self.best_dust_luminosity_all.value / self.stellar_luminosity.value

    # -----------------------------------------------------------------
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
        if self.dust_sed_diffuse_cells: seds[cells_name] = self.dust_sed_diffuse_cells
        return seds

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @lazyproperty
    def absorption_seds_internal(self):
        seds = OrderedDict()
        seds[seds_name] = self.absorption_sed_internal
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def emission_seds_internal(self):
        seds = OrderedDict()
        seds[seds_name] = self.dust_sed_internal
        seds[seds_alt_name] = self.dust_sed_internal_alt
        return seds

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
