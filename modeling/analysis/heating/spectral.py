#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.spectral Contains the SpectralDustHeatingAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty, memoize_method, lazyfileproperty
from ....magic.core.frame import Frame
from ....magic.core.datacube import DataCube
from ....magic.core.list import uniformize
from ....core.units.parsing import parse_quantity
from ....core.basics.curve import WavelengthCurve
from ....magic.tools import plotting
from ...core.data import SpectralData3D

# -----------------------------------------------------------------

max_wavelength_absorption = parse_quantity("5 micron")
min_wavelength_emission = parse_quantity("10 micron")

# -----------------------------------------------------------------

cubes_dirname = "cubes"
curves_dirname = "curves"
cells_dirname = "3D"

# -----------------------------------------------------------------

unevolved_name = "unevolved"
evolved_name = "evolved"

unevolved_short = "unev"
evolved_short = "ev"

# -----------------------------------------------------------------

earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

class SpectralDustHeatingAnalyser(DustHeatingAnalysisComponent):
    
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
        super(SpectralDustHeatingAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Writing
        self.write()

        # Plotting
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SpectralDustHeatingAnalyser, self).setup()

    # -----------------------------------------------------------------

    @property
    def do_earth(self):
        return True

    # -----------------------------------------------------------------

    @property
    def do_faceon(self):
        return True

    # -----------------------------------------------------------------

    @property
    def do_edgeon(self):
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def cubes_path(self):
        return fs.create_directory_in(self.spectral_heating_path, cubes_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_path(self):
        return fs.create_directory_in(self.spectral_heating_path, curves_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def cells_path(self):
        return fs.create_directory_in(self.spectral_heating_path, cells_dirname)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def total_simulations(self):
        return self.model.total_simulations

    # -----------------------------------------------------------------

    @property
    def young_simulations(self):
        return self.model.young_simulations

    # -----------------------------------------------------------------

    @property
    def ionizing_simulations(self):
        return self.model.sfr_simulations

    # -----------------------------------------------------------------
    # EARTH
    # -----------------------------------------------------------------

    @property
    def do_cubes_earth_emission(self):
        return self.has_dust_emission_cubes_earth

    # -----------------------------------------------------------------

    @property
    def do_cubes_earth_absorption(self):
        return self.has_dust_absorption_cubes_earth

    # -----------------------------------------------------------------
    # FACEON
    # -----------------------------------------------------------------

    @property
    def do_cubes_faceon_emission(self):
        return self.has_dust_emission_cubes_faceon

    # -----------------------------------------------------------------

    @property
    def do_cubes_faceon_absorption(self):
        return self.has_dust_absorption_cubes_faceon

    # -----------------------------------------------------------------
    # EDGEON
    # -----------------------------------------------------------------

    @property
    def do_cubes_edgeon_emission(self):
        return self.has_dust_emission_cubes_edgeon

    # -----------------------------------------------------------------

    @property
    def do_cubes_edgeon_absorption(self):
        return self.has_dust_absorption_cubes_edgeon

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_earth(self):
        return self.model.old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_earth(self):
        return self.model.has_old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_earth(self):
        return self.model.young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_earth(self):
        return self.model.has_young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_earth(self):
        return self.model.sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_earth(self):
        return self.model.has_sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_earth(self):
        return self.young_dust_emission_cube_earth + self.ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_earth(self):
        return self.has_young_dust_emission_cube_earth and self.has_ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_earth(self):
        return self.old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_earth(self):
        return self.has_old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_earth(self):
        return self.model.total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_earth(self):
        return self.model.has_total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_earth(self):
        return self.has_unevolved_dust_emission_cube_earth and self.has_total_dust_emission_cube_earth and self.has_evolved_dust_emission_cube_earth

    # -----------------------------------------------------------------
    # CUBES: ABSORPTION
    #   EARTH
    # -----------------------------------------------------------------

    @property
    def cube_earth_absorption_path(self):
        return fs.join(self.cubes_path, "earth_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth_absorption(self):
        return fs.is_file(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_earth_absorption_path", True, write=False)
    def cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_cube_earth + (self.total_dust_absorption_cube_earth - self.evolved_dust_absorption_cube_earth)) / self.total_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def cube_earth_absorption_fixed_path(self):
        return fs.join(self.cubes_path, "earth_absorption_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth_absorption_fixed(self):
        return fs.is_file(self.cube_earth_absorption_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_earth_absorption_fixed_path", True)
    def cube_earth_absorption_fixed(self):

        """
        This function ...
        :return:
        """

        # Interpolated
        return get_fixed_cube_absorption(self.cube_earth_absorption)

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def cube_faceon_absorption_path(self):
        return fs.join(self.cubes_path, "faceon_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon_absorption(self):
        return fs.is_file(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_faceon_absorption_path", True, write=False)
    def cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_cube_faceon + (self.total_dust_absorption_cube_faceon - self.evolved_dust_absorption_cube_faceon)) / self.total_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def cube_faceon_absorption_fixed_path(self):
        return fs.join(self.cubes_path, "faceon_absorption_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon_absorption_fixed(self):
        return fs.is_file(self.cube_faceon_absorption_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_faceon_absorption_fixed_path", True)
    def cube_faceon_absorption_fixed(self):

        """
        This function ...
        :return:
        """

        return get_fixed_cube_absorption(self.cube_faceon_absorption)

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def cube_edgeon_absorption_path(self):
        return fs.join(self.cubes_path, "edgeon_absorption.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon_absorption(self):
        return fs.is_file(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_edgeon_absorption_path", True, write=False)
    def cube_edgeon_absorption(self):

        """
        Thisf unction ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_cube_edgeon + (self.total_dust_absorption_cube_edgeon - self.evolved_dust_absorption_cube_edgeon)) / self.total_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def cube_edgeon_absorption_fixed_path(self):
        return fs.join(self.cubes_path, "edgeon_absorption_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon_absorption_fixed(self):
        return fs.is_file(self.cube_edgeon_absorption_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_edgeon_absorption_fixed_path", True)
    def cube_edgeon_absorption_fixed(self):

        """
        This function ...
        :return:
        """

        return get_fixed_cube_absorption(self.cube_edgeon_absorption)

    # -----------------------------------------------------------------
    # CUBES: EMISSION
    #   EARTH
    # -----------------------------------------------------------------

    @property
    def cube_earth_emission_path(self):
        return fs.join(self.cubes_path, "earth_emission.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth_emission(self):
        return fs.is_file(self.cube_earth_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_earth_emission_path", True, write=False)
    def cube_earth_emission(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_cube_earth + (self.total_dust_emission_cube_earth - self.evolved_dust_emission_cube_earth)) / self.total_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def cube_earth_emission_fixed_path(self):
        return fs.join(self.cubes_path, "earth_emission_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_earth_emission_fixed(self):
        return fs.is_file(self.cube_earth_emission_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_earth_emission_fixed_path", True)
    def cube_earth_emission_fixed(self):

        """
        This function ...
        :return:
        """

        return get_fixed_cube_emission(self.cube_earth_emission)

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def cube_faceon_emission_path(self):
        return fs.join(self.cubes_path, "faceon_emission.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon_emission(self):
        return fs.is_file(self.cube_faceon_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_faceon_emission_path", True, write=False)
    def cube_faceon_emission(self):

        """
        THs function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_cube_faceon + (self.total_dust_emission_cube_faceon - self.evolved_dust_emission_cube_faceon)) / self.total_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def cube_faceon_emission_fixed_path(self):
        return fs.join(self.cubes_path, "faceon_emission_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_faceon_emission_fixed(self):
        return fs.is_file(self.cube_faceon_emission_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_faceon_emission_fixed_path", True)
    def cube_faceon_emission_fixed(self):

        """
        This function ...
        :return:
        """

        return get_fixed_cube_emission(self.cube_faceon_emission)

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def cube_edgeon_emission_path(self):
        return fs.join(self.cubes_path, "edgeon.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon_emission(self):
        return fs.is_file(self.cube_edgeon_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_edgeon_emission_path", True, write=False)
    def cube_edgeon_emission(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_cube_edgeon + (self.total_dust_emission_cube_edgeon - self.evolved_dust_emission_cube_edgeon)) / self.total_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def cube_edgeon_emission_fixed_path(self):
        return fs.join(self.cubes_path, "edgeon_fixed.fits")

    # -----------------------------------------------------------------

    @property
    def has_cube_edgeon_emission_fixed(self):
        return fs.is_file(self.cube_edgeon_emission_fixed_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(DataCube, "cube_edgeon_emission_fixed_path", True)
    def cube_edgeon_emission_fixed(self):

        """
        This function ...
        :return:
        """

        return get_fixed_cube_emission(self.cube_edgeon_emission)

    # -----------------------------------------------------------------
    # CURVES: ABSORPTION
    #   EARTH
    # -----------------------------------------------------------------

    @property
    def wavelength_unit(self):
        return "micron"

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_name(self):
        return "Funev_absorption_earth"

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_description(self):
        return "Fraction of spectral luminosity from unevolved stars absorbed by dust (from earth projection)"

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_path(self):
        return fs.join(self.curves_path, "earth_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption(self):
        return fs.is_file(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_earth_absorption_path", True)
    def curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        return self.cube_earth_absorption_fixed.global_curve(self.curve_earth_absorption_name, measure="mean", description=self.curve_earth_absorption_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_wavelengths(self):
        return self.curve_earth_absorption.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_values(self):
        return self.curve_earth_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_evolved_name(self):
        return self.curve_earth_absorption_name.replace(unevolved_short, evolved_short)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_evolved_description(self):
        return self.curve_earth_absorption_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_earth_absorption_evolved_name,
                                                           self.curve_earth_absorption_wavelengths,
                                                           1. - self.curve_earth_absorption_values,
                                                           description=self.curve_earth_absorption_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_name(self):
        return "Funev_absorption_faceon"

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_description(self):
        return self.curve_earth_absorption_description.replace(earth_name, faceon_name)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_path(self):
        return fs.join(self.curves_path, "faceon_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption(self):
        return fs.is_file(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_faceon_absorption_path", True, write=False)
    def curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        return self.cube_faceon_absorption_fixed.global_curve(self.curve_faceon_absorption_name, measure="mean", description=self.curve_faceon_absorption_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_wavelengths(self):
        return self.curve_faceon_absorption.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_values(self):
        return self.curve_faceon_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_evolved_name(self):
        return self.curve_faceon_absorption_name.replace(unevolved_short, evolved_short)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_evolved_description(self):
        return self.curve_faceon_absorption_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_faceon_absorption_evolved_name,
                                                           self.curve_faceon_absorption_wavelengths,
                                                           1. - self.curve_faceon_absorption_values,
                                                           description=self.curve_faceon_absorption_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_name(self):
        return self.curve_earth_absorption_name.replace(earth_name, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_description(self):
        return self.curve_earth_absorption_description.replace(earth_name, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_path(self):
        return fs.join(self.curves_path, "edgeon_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption(self):
        return fs.is_file(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_edgeon_absorption_path", True, write=False)
    def curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        return self.cube_edgeon_absorption_fixed.global_curve(self.curve_edgeon_absorption_name, measure="mean", description=self.curve_edgeon_absorption_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_wavelengths(self):
        return self.curve_edgeon_absorption.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_values(self):
        return self.curve_edgeon_absorption.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_evolved_name(self):
        return self.curve_edgeon_absorption_name.replace(unevolved_short, evolved_short)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_evolved_description(self):
        return self.curve_edgeon_absorption_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_edgeon_absorption_evolved_name,
                                                           self.curve_edgeon_absorption_wavelengths,
                                                           1. - self.curve_edgeon_absorption_values,
                                                           description=self.curve_edgeon_absorption_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    # CURVES: EMISSION
    #   EARTH
    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_name(self):
        return "Funev_emission_earth"

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_description(self):
        return "Fraction of emitted dust spectral luminosity by unevolved stars (from earth projection)"

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_path(self):
        return fs.join(self.curves_path, "earth_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission(self):
        return fs.is_file(self.curve_earth_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_earth_emission_path", True, write=False)
    def curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        return self.cube_earth_emission_fixed.global_curve(self.curve_earth_emission_name, measure="mean", description=self.curve_earth_emission_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_wavelengths(self):
        return self.curve_earth_emission.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_values(self):
        return self.curve_earth_emission.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_evolved_name(self):
        return "Fev_emission_earth"

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_evolved_description(self):
        return self.curve_earth_emission_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_earth_emission_evolved_name,
                                                           self.curve_earth_emission_wavelengths,
                                                           1. - self.curve_earth_emission_values,
                                                           description=self.curve_earth_emission_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    #   FACEON
    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_name(self):
        return self.curve_earth_emission_name.replace(earth_name, faceon_name)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_description(self):
        return self.curve_earth_emission_description.replace(earth_name, faceon_name)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_path(self):
        return fs.join(self.curves_path, "faceon_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission(self):
        return fs.is_file(self.curve_faceon_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_faceon_emission_path", True, write=False)
    def curve_faceon_emission(self):

        """
        This function ...
        :return:
        """

        return self.cube_faceon_emission_fixed.global_curve(self.curve_faceon_emission_name, measure="mean", description=self.curve_faceon_emission_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_wavelengths(self):
        return self.curve_faceon_emission.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_values(self):
        return self.curve_faceon_emission.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_evolved_name(self):
        return "Fev_emission_faceon"

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_evolved_description(self):
        return self.curve_faceon_emission_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_faceon_emission_evolved_name,
                                                           self.curve_faceon_emission_wavelengths,
                                                           1. - self.curve_faceon_emission_values,
                                                           description=self.curve_faceon_emission_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    #   EDGEON
    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_name(self):
        return self.curve_earth_emission_name.replace(earth_name, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_description(self):
        return self.curve_earth_emission_description.replace(earth_name, edgeon_name)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_path(self):
        return fs.join(self.curves_path, "edgeon_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission(self):
        return fs.is_file(self.curve_edgeon_emission_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_edgeon_emission_path", True, write=False)
    def curve_edgeon_emission(self):

        """
        This function ...
        :return:
        """

        return self.cube_edgeon_emission_fixed.global_curve(self.curve_edgeon_emission_name, measure="mean", description=self.curve_edgeon_emission_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_wavelengths(self):
        return self.curve_edgeon_emission.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_values(self):
        return self.curve_edgeon_emission.values(asarray=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_evolved_name(self):
        return self.curve_edgeon_emission_name.replace(unevolved_short, evolved_short)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_evolved_description(self):
        return self.curve_edgeon_emission_description.replace(unevolved_name, evolved_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_evolved(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_edgeon_emission_evolved_name, self.curve_edgeon_emission_wavelengths,
                                                           1. - self.curve_edgeon_emission_values,
                                                           description=self.curve_edgeon_emission_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    # 3D (CELL) DATA
    #   TOTAL
    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_name(self):
        return "Labs_lambda_total"

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the total model"

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_path(self):
        return fs.join(self.cells_path, "total_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_absorption(self):
        return fs.is_file(self.total_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "total_spectral_absorption_path", True, write=True)
    def total_spectral_absorption_data(self):

        """
        This function ...
        :return:
        """

        return SpectralData3D.from_table_file(self.total_contribution_spectral_absorption_filepath, self.cell_x_coordinates,
                                              self.cell_y_coordinates, self.cell_z_coordinates, length_unit=self.length_unit,
                                              name=self.total_spectral_absorption_name, description=self.total_spectral_absorption_description)

    # -----------------------------------------------------------------
    #   UNEVOLVED
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_name(self):
        return "Labs_lambda_unevolved"

    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the unevolved stellar populations"

    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_path(self):
        return fs.join(self.cells_path, "unevolved_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_absorption(self):
        return fs.is_file(self.unevolved_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "unevolved_spectral_absorption_path", True, write=True)
    def unevolved_spectral_absorption_data(self):

        """
        This function ...
        :return:
        """

        return SpectralData3D.from_table_file(self.unevolved_contribution_spectral_absorption_filepath, self.cell_x_coordinates,
                                              self.cell_y_coordinates, self.cell_z_coordinates, length_unit=self.length_unit,
                                              name=self.unevolved_spectral_absorption_name, description=self.unevolved_spectral_absorption_description)

    # -----------------------------------------------------------------
    #   EVOLVED (OLD)
    # -----------------------------------------------------------------

    @property
    def evolved_spectral_absorption_name(self):
        return "Labs_lambda_evolved"

    # -----------------------------------------------------------------

    @property
    def evolved_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the evolved stellar populations"

    # -----------------------------------------------------------------

    @property
    def evolved_spectral_absorption_path(self):
        return fs.join(self.cells_path, "evolved_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_evolved_spectral_absorption(self):
        return fs.is_file(self.evolved_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "evolved_spectral_absorption_path", True, write=True)
    def evolved_spectral_absorption_data(self):

        """
        Thisn function ...
        :return:
        """

        return SpectralData3D.from_table_file(self.old_contribution_spectral_absorption_filepath, self.cell_x_coordinates,
                                              self.cell_y_coordinates, self.cell_z_coordinates, length_unit=self.length_unit,
                                              name=self.evolved_spectral_absorption_name, description=self.evolved_spectral_absorption_description)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_earth(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_earth(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_earth(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_earth(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_earth(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_earth(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_earth(self):
        #return self.young_dust_absorption_cube_earth + self.ionizing_dust_absorption_cube_earth
        return self.model.unevolved_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_earth(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_earth(self):
        return self.old_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_earth(self):
        return self.has_old_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_earth(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_earth(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_earth(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_earth and self.has_total_dust_absorption_cube_earth and self.has_evolved_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_faceon(self):
        return self.model.old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_faceon(self):
        return self.model.has_old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_faceon(self):
        return self.model.young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_faceon(self):
        return self.model.has_young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_faceon(self):
        return self.model.sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_faceon(self):
        return self.model.has_sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def unevolved_dust_emission_cube_faceon(self):
        return self.young_dust_emission_cube_faceon + self.ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_faceon(self):
        return self.has_young_dust_emission_cube_faceon and self.has_ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_faceon(self):
        return self.old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_faceon(self):
        return self.has_old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_faceon(self):
        return self.model.total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_faceon(self):
        return self.model.has_total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_faceon(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_emission_cube_faceon and self.has_total_dust_emission_cube_faceon and self.has_evolved_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_faceon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_faceon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_faceon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_faceon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_faceon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_faceon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_faceon(self):
        return self.young_dust_absorption_cube_faceon + self.ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_faceon(self):
        return self.has_young_dust_absorption_cube_faceon and self.has_ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_faceon(self):
        return self.old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_faceon(self):
        return self.has_old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_faceon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_faceon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_faceon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_faceon and self.has_total_dust_absorption_cube_faceon and self.has_evolved_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_edgeon(self):
        return self.model.old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_edgeon(self):
        return self.model.has_old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_edgeon(self):
        return self.model.young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_edgeon(self):
        return self.model.has_young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_edgeon(self):
        return self.model.sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_edgeon(self):
        return self.model.has_sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_edgeon(self):
        return self.young_dust_emission_cube_edgeon + self.ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_edgeon(self):
        return self.has_young_dust_emission_cube_edgeon and self.has_ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_edgeon(self):
        return self.old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_edgeon(self):
        return self.has_old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_edgeon(self):
        return self.model.total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_edgeon(self):
        return self.model.has_total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.has_unevolved_dust_emission_cube_edgeon and self.has_total_dust_emission_cube_edgeon and self.has_evolved_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_edgeon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_edgeon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_edgeon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_edgeon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_edgeon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_edgeon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_edgeon(self):
        return self.young_dust_absorption_cube_edgeon + self.ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_edgeon(self):
        return self.has_young_dust_absorption_cube_edgeon and self.has_ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_edgeon(self):
        return self.old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_edgeon(self):
        return self.has_old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_edgeon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_edgeon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        return self.has_unevolved_dust_absorption_cube_edgeon and self.has_total_dust_absorption_cube_edgeon and self.has_evolved_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def do_curve_earth_emission(self):
        return self.do_cubes_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_earth_absorption(self):
        return self.do_cubes_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_curve_faceon_emission(self):
        return self.do_cubes_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_faceon_absorption(self):
        return self.do_cubes_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_curve_edgeon_emission(self):
        return self.do_cubes_edgeon_emission

    # -----------------------------------------------------------------

    @property
    def do_curve_edgeon_absorption(self):
        return self.do_cubes_edgeon_absorption

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write cubes
        self.write_cubes()

        # Write the curves
        self.write_curves()

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth_emission(self):
        return self.do_cubes_earth_emission and not self.has_cube_earth

    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth_absorption(self):
        return self.do_cubes_earth_absorption and not self.has_cube_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon_emission(self):
        return self.do_cubes_faceon_emission and not self.has_cube_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon_absorption(self):
        return self.do_cubes_faceon_absorption and not self.has_cube_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_edgeon_emission(self):
        return self.do_cubes_edgeon_emission and not self.has_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def do_write_cube_edgeon_absorption(self):
        return self.do_cubes_edgeon_absorption and not self.has_cube_edgeon_absorption

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def write_cubes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the heating cubes ...")

        # Absorption
        self.write_cubes_absorption()

        # Emission
        self.write_cubes_emission()

    # -----------------------------------------------------------------

    def write_cubes_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption cubes ...")

        # Earth
        if self.do_write_cube_earth_absorption: self.write_cube_earth_absorption()

        # Face-on
        if self.do_write_cube_faceon_absorption: self.write_cube_faceon_absorption()

        # Edge-on
        if self.do_write_cube_edgeon_absorption: self.write_cube_edgeon_absorption()

    # -----------------------------------------------------------------

    def write_cubes_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the emission cubes ...")

        # Earth
        if self.do_write_cube_earth_emission: self.write_cube_earth_emission()

        # Face-on
        if self.do_write_cube_faceon_emission: self.write_cube_faceon_emission()

        # Edge-on
        if self.do_write_cube_edgeon_emission: self.write_cube_edgeon_emission()

    # -----------------------------------------------------------------

    def remove_cube_earth_emission(self):
        fs.remove_file(self.cube_earth_emission_path)

    # -----------------------------------------------------------------

    def write_cube_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_earth_emission.saveto(self.cube_earth_emission_path)

    # -----------------------------------------------------------------

    def remove_cube_earth_absorption(self):
        fs.remove_file(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_earth_absorption.saveto(self.cube_earth_absorption_path)

    # -----------------------------------------------------------------

    def remove_cube_faceon_emission(self):
        fs.remove_file(self.cube_faceon_emission_path)

    # -----------------------------------------------------------------

    def write_cube_faceon_emission(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_faceon_emission.saveto(self.cube_faceon_emission_path)

    # -----------------------------------------------------------------

    def remove_cube_faceon_absorption(self):
        fs.remove_file(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_faceon_absorption.saveto(self.cube_faceon_absorption_path)

    # -----------------------------------------------------------------

    def remove_cube_edgeon_emission(self):
        fs.remove_file(self.cube_edgeon_emission_path)

    # -----------------------------------------------------------------

    def write_cube_edgeon_emission(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_edgeon_emission.saveto(self.cube_edgeon_emission_path)

    # -----------------------------------------------------------------

    def remove_cube_edgeon_absorption(self):
        fs.remove_file(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def write_cube_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.cube_edgeon_absorption.saveto(self.cube_edgeon_absorption_path)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth(self):
        return self.do_curve_earth_emission and not self.has_curve_earth

    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth_absorption(self):
        return self.do_curve_earth_absorption and not self.has_curve_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon(self):
        return self.do_curve_faceon_emission and not self.has_curve_faceon

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_absorption(self):
        return self.do_curve_faceon_absorption and not self.has_curve_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon(self):
        return self.do_curve_edgeon_emission and not self.has_curve_edgeon

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon_absorption(self):
        return self.do_curve_edgeon_absorption and not self.has_curve_edgeon_absorption

    # -----------------------------------------------------------------

    def write_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the curves ...")

        # Absorption
        self.write_curves_absorption()

        # Emission
        self.write_curves_emission()

    # -----------------------------------------------------------------

    def write_curves_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption curves ...")

        # Earth
        if self.do_write_curve_earth_absorption: self.write_curve_earth_absorption()

        # Face-on
        if self.do_write_curve_faceon_absorption: self.write_curve_faceon_absorption()

        # Edge-on
        if self.do_write_curve_edgeon_absorption: self.write_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    def write_curves_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the emission curves ...")

        # Earth
        if self.do_write_curve_earth_emission: self.write_curve_earth_emission()

        # Face-on
        if self.do_write_curve_faceon_emission: self.write_curve_faceon_emission()

        # Edge-on
        if self.do_write_curve_edgeon_emission: self.write_curve_edgeon_emission()

    # -----------------------------------------------------------------

    def remove_curve_earth_emission(self):
        fs.remove_file(self.curve_earth_emission_path)

    # -----------------------------------------------------------------

    def write_curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_earth_emission.saveto(self.curve_earth_emission_path)

    # -----------------------------------------------------------------

    def remove_curve_earth_absorption(self):
        fs.remove_file(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_earth_absorption.saveto(self.curve_earth_absorption_path)

    # -----------------------------------------------------------------

    def remove_curve_faceon_emission(self):
        fs.remove_file(self.curve_faceon_emission_path)

    # -----------------------------------------------------------------

    def write_curve_faceon_emission(self):

        """
        This functino ...
        :return:
        """

        # Save
        self.curve_faceon_emission.saveto(self.curve_faceon_emission_path)

    # -----------------------------------------------------------------

    def remove_curve_faceon_absorption(self):
        fs.remove_file(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_faceon_absorption.saveto(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    def remove_curve_edgeon_emission(self):
        fs.remove_file(self.curve_edgeon_emission_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_edgeon_emission.saveto(self.curve_edgeon_emission_path)

    # -----------------------------------------------------------------

    def remove_curve_edgeon_absorption(self):
        fs.remove_file(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Save
        self.curve_edgeon_absorption.saveto(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps of spectral heating
        self.plot_spectral_maps()

        # Curves of spectral heating
        self.plot_curves()

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_emission_earth(self):
        return self.do_cubes_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_earth(self):
        return self.do_cubes_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_emission_faceon(self):
        return self.do_cubes_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_faceon(self):
        return self.do_cubes_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_emission_edgeon(self):
        return self.do_cubes_edgeon_emission

    # -----------------------------------------------------------------

    @property
    def do_plot_spectral_maps_absorption_edgeon(self):
        return self.do_cubes_edgeon_absorption

    # -----------------------------------------------------------------

    def plot_spectral_maps(self):

        """
        This function ...
        :return:
        """

        # Earth
        if self.do_earth:
            if self.do_plot_spectral_maps_emission_earth: self.plot_spectral_maps_emission_earth()
            if self.do_plot_spectral_maps_absorption_earth: self.plot_spectral_maps_absorption_earth()

        # Face-on
        if self.do_faceon:
            if self.do_plot_spectral_maps_emission_faceon: self.plot_spectral_maps_emission_faceon()
            if self.do_plot_spectral_maps_absorption_faceon: self.plot_spectral_maps_absorption_faceon()

        # Edge-on
        if self.do_edgeon:
            if self.do_plot_spectral_maps_emission_edgeon: self.plot_spectral_maps_emission_edgeon()
            if self.do_plot_spectral_maps_absorption_edgeon: self.plot_spectral_maps_absorption_edgeon()

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_earth.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_earth_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "earth_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_earth_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_earth(self):

        """
        Thsf unction ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_earth(fltr)

            # Get the path
            path = self.get_spectral_map_emission_earth_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_earth(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_earth_absorption.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_earth_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "earth_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_earth(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_earth_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_earth(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_earth_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_faceon.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_faceon_path(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "faceon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_faceon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the face-on projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_emission_faceon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_faceon_absorption.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_faceon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "faceon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_faceon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_faceon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_faceon(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the face-on projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_faceon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_edgeon.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_edgeon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "edgeon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_emission_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_emission_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust emission from the edge-on projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if the map exists
            if self.has_spectral_map_emission_edgeon(fltr): continue
            
            # Get the map
            frame = self.get_spectral_map_emission_edgeon(fltr)
            
            # Get the path
            path = self.get_spectral_map_emission_edgeon_path(fltr)
            
            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.cube_edgeon_absorption.frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_edgeon_path(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.cubes_path, "edgeon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_edgeon(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.is_file(self.get_spectral_map_absorption_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def plot_spectral_maps_absorption_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting maps of the spectral heating by dust absorption from the edge-on projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if the map exists
            if self.has_spectral_map_absorption_edgeon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_edgeon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_edgeon_path(fltr)

            # Plot
            plotting.plot_map(frame, path=path, interval=(0,1,))

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_earth_emission(self):
        return self.do_curve_earth_emission and not self.has_curve_earth_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_earth_absorption(self):
        return self.do_curve_earth_absorption and not self.has_curve_earth_absorption_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_emission(self):
        return self.do_curve_faceon_emission and not self.has_curve_faceon_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_absorption(self):
        return self.do_curve_faceon_absorption and not self.has_curve_faceon_absorption_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_emission(self):
        return self.do_curve_edgeon_emission and not self.has_curve_edgeon_emission_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_absorption(self):
        return self.do_curve_edgeon_absorption and not self.has_curve_edgeon_absorption_plot

    # -----------------------------------------------------------------

    def plot_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curves ...")

        # Absorption
        self.plot_curves_absorption()

        # Emission
        self.plot_curves_emission()

    # -----------------------------------------------------------------

    def plot_curves_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption curves ...")

        # Earth
        if self.do_plot_curve_earth_absorption: self.plot_curve_earth_absorption()

        # Face-on
        if self.do_plot_curve_faceon_absorption: self.plot_curve_faceon_absorption()

        # Edge-on
        if self.do_plot_curve_edgeon_absorption: self.plot_curve_edgeon_absorption()

    # -----------------------------------------------------------------

    def plot_curves_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the emission curves ...")

        # Earth
        if self.do_plot_curve_earth_emission: self.plot_curve_earth_emission()

        # Face-on
        if self.do_plot_curve_faceon_emission: self.plot_curve_faceon_emission()

        # Edge-on
        if self.do_plot_curve_edgeon_emission: self.plot_curve_edgeon_emission()

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_plot_path(self):
        return fs.join(self.curves_path, "earth.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_plot(self):
        return fs.is_file(self.curve_earth_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_emission(self):
        return {unevolved_name: self.curve_earth_emission, evolved_name: self.curve_earth_emission_evolved}

    # -----------------------------------------------------------------

    def plot_curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth_emission, path=self.curve_earth_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_plot_path(self):
        return fs.join(self.curves_path, "earth_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption_plot(self):
        return fs.is_file(self.curve_earth_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_absorption(self):
        return {unevolved_name: self.curve_earth_absorption, evolved_name: self.curve_earth_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth_absorption, path=self.curve_earth_absorption_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_plot_path(self):
        return fs.join(self.curves_path, "faceon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission_plot(self):
        return fs.is_file(self.curve_faceon_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon_emission(self):
        return {unevolved_name: self.curve_faceon_emission, evolved_name: self.curve_faceon_emission_evolved}

    # -----------------------------------------------------------------

    def plot_curve_faceon_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the face-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_faceon_emission, path=self.curve_faceon_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_plot_path(self):
        return fs.join(self.curves_path, "faceon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption_plot(self):
        return fs.is_file(self.curve_faceon_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon_absorption(self):
        return {unevolved_name: self.curve_faceon_absorption, evolved_name: self.curve_faceon_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_faceon_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the face-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_faceon_absorption, path=self.curve_faceon_absorption_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_plot_path(self):
        return fs.join(self.curves_path, "edgeon.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission_plot(self):
        return fs.is_file(self.curve_edgeon_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon_emission(self):
        return {unevolved_name: self.curve_edgeon_emission, evolved_name: self.curve_edgeon_emission_evolved}

    # -----------------------------------------------------------------

    def plot_curve_edgeon_emission(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the edge-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_edgeon_emission, path=self.curve_edgeon_emission_plot_path, xlog=True)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_plot_path(self):
        return fs.join(self.curves_path, "edgeon_absorption.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption_plot(self):
        return fs.is_file(self.curve_edgeon_absorption_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon_absorption(self):
        return {unevolved_name: self.curve_edgeon_absorption, evolved_name: self.curve_edgeon_absorption_evolved}

    # -----------------------------------------------------------------

    def plot_curve_edgeon_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the edge-on projection ...")

        # Plot
        plotting.plot_curves(self.curves_edgeon_absorption, path=self.curve_edgeon_absorption_plot_path, xlog=True)

# -----------------------------------------------------------------

def get_fixed_cube_emission(cube):

    """
    This function ...
    :param cube:
    :return:
    """

    # Create truncated cube
    fixed = cube.truncated(min_wavelength=min_wavelength_emission, copy=True) # copy the frames

    # Fix
    fixed.replace_negatives_by_nans()
    fixed.replace_infs_by_nans()

    # Replace
    fixed.replace_by_nans_where_greater_than(1.1)
    fixed.cutoff_greater(1.)

    # Interpolate nans
    fixed.interpolate_nans(sigma=3.)

    # Set flag
    fixed.metadata["fixed"] = True

    # Return the new cube
    return fixed

# -----------------------------------------------------------------

def get_fixed_cube_absorption(cube):

    """
    This function ...
    :param cube:
    :return:
    """

    # Create truncated cube
    fixed = cube.truncated(max_wavelength=max_wavelength_absorption, copy=True) # copy the frames

    # Fix
    fixed.replace_negatives_by_nans()
    fixed.replace_infs_by_nans()

    # Replace
    fixed.replace_by_nans_where_greater_than(1.1)
    fixed.cutoff_greater(1.)

    # Interpolate nans
    fixed.interpolate_nans(sigma=3.)

    # Set flag
    fixed.metadata["fixed"] = True

    # Return the new cube
    return fixed

# -----------------------------------------------------------------
