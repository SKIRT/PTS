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
from ....magic.core.image import Image
from ....magic.core.frame import Frame
from ....magic.core.datacube import DataCube
from ....core.units.parsing import parse_quantity
from ....core.basics.curve import WavelengthCurve
from ....magic.tools import plotting
from ...core.data import SpectralData3D
from ...projection.data import project_data

# -----------------------------------------------------------------

max_wavelength_absorption = parse_quantity("5 micron")
min_wavelength_emission = parse_quantity("10 micron")

# -----------------------------------------------------------------

cubes_dirname = "cubes"
curves_dirname = "curves"
maps_dirname = "maps"
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
    def absorption_path(self):
        return self.analysis_run.absorption_path

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
    def maps_path(self):
        return fs.create_directory_in(self.spectral_heating_path, maps_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def cells_path(self):
        return fs.create_directory_in(self.spectral_heating_path, cells_dirname)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_absorption_data(self):
        return self.has_total_contribution_spectral_absorption and self.has_old_contribution_spectral_absorption and self.has_unevolved_contribution_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def do_emission_data(self):
        return self.has_total_contribution_spectral_emission and self.has_old_contribution_spectral_emission and self.has_unevolved_contribution_spectral_emission

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

    @property
    def earth_projection(self):
        return self.analysis_run.earth_projection

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):
        return self.analysis_run.faceon_projection

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):
        return self.analysis_run.edgeon_projection

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
    # CUBES: ABSORPTION FRACTION
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
        return fs.join(self.cubes_path, "edgeon_emission.fits")

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
        return fs.join(self.cubes_path, "edgeon_emission_fixed.fits")

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

    @property
    def curve_earth_absorption_seds_path(self):
        return fs.join(self.curves_path, "earth_absorption_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption_seds(self):
        return fs.is_file(self.curve_earth_absorption_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_earth_absorption_seds_path", True, write=False)
    def curve_earth_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_sed_earth + (self.total_dust_absorption_sed_earth
                - self.evolved_dust_absorption_sed_earth)) / self.total_dust_absorption_sed_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_seds_wavelengths(self):
        return self.curve_earth_absorption_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_seds_values(self):
        return self.curve_earth_absorption_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_absorption_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_earth_absorption_evolved_name,
                                                           self.curve_earth_absorption_seds_wavelengths,
                                                           1. - self.curve_earth_absorption_seds_values,
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

    @property
    def curve_faceon_absorption_seds_path(self):
        return fs.join(self.curves_path, "faceon_absorption_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption_seds(self):
        return fs.is_file(self.curve_faceon_absorption_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_faceon_absorption_seds_path", True, write=False)
    def curve_faceon_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_sed_faceon + (self.total_dust_absorption_sed_faceon
                - self.evolved_dust_absorption_sed_faceon)) / self.total_dust_absorption_sed_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_seds_wavelengths(self):
        return self.curve_faceon_absorption_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_seds_values(self):
        return self.curve_faceon_absorption_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_absorption_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_faceon_absorption_evolved_name,
                                                           self.curve_faceon_absorption_seds_wavelengths,
                                                           1. - self.curve_faceon_absorption_seds_values,
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

    @property
    def curve_edgeon_absorption_seds_path(self):
        return fs.join(self.curves_path, "edgeon_absorption_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption_seds(self):
        return fs.is_file(self.curve_edgeon_absorption_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_edgeon_absorption_seds_path", True, write=False)
    def curve_edgeon_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN ABSORBED STELLAR WAVELENGTHS
        return 0.5 * (self.unevolved_dust_absorption_sed_edgeon + (self.total_dust_absorption_sed_edgeon
                - self.evolved_dust_absorption_sed_edgeon)) / self.total_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_seds_wavelengths(self):
        return self.curve_edgeon_absorption_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_seds_values(self):
        return self.curve_edgeon_absorption_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_absorption_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_edgeon_absorption_evolved_name,
                                                           self.curve_edgeon_absorption_seds_wavelengths,
                                                           1. - self.curve_edgeon_absorption_seds_values,
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

    @property
    def curve_earth_emission_seds_path(self):
        return fs.join(self.curves_path, "earth_emission_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_seds(self):
        return fs.is_file(self.curve_earth_emission_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_earth_emission_seds_path", True, write=False)
    def curve_earth_emission_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_sed_earth + (self.total_dust_emission_sed_earth -
                    self.evolved_dust_emission_sed_earth)) / self.total_dust_emission_sed_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_seds_wavelengths(self):
        return self.curve_earth_emission_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_seds_values(self):
        return self.curve_earth_emission_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_earth_emission_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_earth_emission_evolved_name, self.curve_earth_emission_seds_wavelengths,
                                                           1. - self.curve_earth_emission_seds_values, description=self.curve_earth_emission_evolved_description,
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

    @property
    def curve_faceon_emission_seds_path(self):
        return fs.join(self.curves_path, "faceon_emission_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission_seds(self):
        return fs.is_file(self.curve_faceon_emission_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_faceon_emission_seds_path", True)
    def curve_faceon_emission_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_sed_faceon + (self.total_dust_emission_sed_faceon -
                self.evolved_dust_emission_sed_faceon)) / self.total_dust_emission_sed_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_seds_wavelengths(self):
        return self.curve_faceon_emission_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_seds_values(self):
        return self.curve_faceon_emission_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_faceon_emission_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_faceon_emission_evolved_name,
                                                           self.curve_faceon_emission_seds_wavelengths,
                                                           1. - self.curve_faceon_emission_seds_values,
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

    @property
    def curve_edgeon_emission_seds_path(self):
        return fs.join(self.curves_path, "edgeon_emission_seds.dat")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission_seds(self):
        return fs.is_file(self.curve_edgeon_emission_seds_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "curve_edgeon_emission_seds_path", True, write=False)
    def curve_edgeon_emission_seds(self):

        """
        This function ...
        :return:
        """

        # LOOKING IN DUST EMISSION WAVELENGTHS
        return 0.5 * (self.unevolved_dust_emission_sed_edgeon + (self.total_dust_emission_sed_edgeon -
                self.evolved_dust_emission_sed_edgeon)) / self.total_dust_emission_sed_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_seds_wavelengths(self):
        return self.curve_edgeon_emission_seds.wavelengths(asarray=True, unit=self.wavelength_unit)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_seds_values(self):
        return self.curve_edgeon_emission_seds.values(asarray=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def curve_edgeon_emission_evolved_seds(self):

        """
        This function ...
        :return:
        """

        return WavelengthCurve.from_wavelengths_and_values(self.curve_edgeon_emission_evolved_name,
                                                           self.curve_edgeon_emission_seds_wavelengths,
                                                           1. - self.curve_edgeon_emission_seds_values,
                                                           description=self.curve_edgeon_emission_evolved_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    # 3D (CELL) DATA
    #   TOTAL: ABSORPTION
    # -----------------------------------------------------------------

    @property
    def length_unit(self):
        return "pc"

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "total_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_absorption(self):
        return fs.is_file(self.total_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_spectral_absorption_data(self):
        if not self.has_total_spectral_absorption: raise RuntimeError("Cannot find the total spectral absorption data: run the absorption analysis first")
        return SpectralData3D.from_file(self.total_spectral_absorption_path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_total_absorption_data_for_filter(self, fltr):
        return self.total_spectral_absorption_data.get_data3d_for_wavelength(fltr.wavelength)

    # -----------------------------------------------------------------
    #   TOTAL: EMISSION
    # -----------------------------------------------------------------

    @property
    def total_spectral_emission_path(self):
        return fs.join(self.absorption_path, "total_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_emission(self):
        return fs.is_file(self.total_spectral_emission_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_spectral_emission_data(self):
        if not self.has_total_spectral_emission: raise RuntimeError("Cannot find the total spectral emission data: run the absorption analysis first")
        return SpectralData3D.from_file(self.total_spectral_emission_path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_total_emission_data_for_filter(self, fltr):
        return self.total_spectral_emission_data.get_data3d_for_wavelength(fltr.wavelength)

    # -----------------------------------------------------------------
    #   EVOLVED (OLD): ABSORPTION
    # -----------------------------------------------------------------

    @property
    def evolved_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "old_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_evolved_spectral_absorption(self):
        return fs.is_file(self.evolved_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def evolved_spectral_absorption_data(self):
        if not self.has_evolved_spectral_absorption: raise RuntimeError("Cannot find the old spectral absorption data: run the absorption analysis first")
        return SpectralData3D.from_file(self.evolved_spectral_absorption_path)

    # -----------------------------------------------------------------
    #   EVOLVED (OLD): EMISSION
    # -----------------------------------------------------------------

    @property
    def evolved_spectral_emission_path(self):
        return fs.join(self.absorption_path, "old_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_evolved_spectral_emission(self):
        return fs.is_file(self.evolved_spectral_emission_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def evolved_spectral_emission_data(self):
        if not self.has_evolved_spectral_emission: raise RuntimeError("Cannot find the old spectral emission data: run the absorption analysis first")
        return SpectralData3D.from_file(self.evolved_spectral_emission_path)

    # -----------------------------------------------------------------
    #   UNEVOLVED: ABSORPTION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "unevolved_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_absorption(self):
        return fs.is_file(self.unevolved_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_spectral_absorption_data(self):
        if not self.has_unevolved_spectral_absorption: raise RuntimeError("Cannot find the unevolved spectral absorption data: run the absorption analysis first")
        return SpectralData3D.from_file(self.unevolved_spectral_absorption_path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_absorption_data_for_filter(self, fltr):
        return self.unevolved_spectral_absorption_data.get_data3d_for_wavelength(fltr.wavelength)

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_absorption_fraction_data_for_filter(self, fltr):
        return self.get_unevolved_absorption_data_for_filter(fltr) / self.get_total_absorption_data_for_filter(fltr)

    # -----------------------------------------------------------------
    #   UNEVOLVED: EMISSION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_emission_path(self):
        return fs.join(self.absorption_path, "unevolved_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_emission(self):
        return fs.is_file(self.unevolved_spectral_emission_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_spectral_emission_data(self):
        if not self.has_unevolved_spectral_emission: raise RuntimeError("Cannot find the unevolved spectral emission data: run the absorption analysis first")
        return SpectralData3D.from_file(self.unevolved_spectral_emission_path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_emission_data_for_filter(self, fltr):
        return self.unevolved_spectral_emission_data.get_data3d_for_wavelength(fltr.wavelength)

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_emission_fraction_data_for_filter(self, fltr):
        return self.get_unevolved_emission_data_for_filter(fltr) / self.get_total_emission_data_for_filter(fltr)

    # -----------------------------------------------------------------
    # CURVES FROM SPECTRAL 3D DATA
    #   TOTAL: ABSORPTION
    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "total_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_absorption_curve(self):
        return fs.is_file(self.total_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_spectral_absorption_curve(self):
        if not self.has_total_spectral_absorption_curve: raise RuntimeError("Cannot find the total spectral absorption curve: run the absorption analysis first")
        return WavelengthCurve.from_file(self.total_spectral_absorption_curve_path)

    # -----------------------------------------------------------------
    #   TOTAL: EMISSION
    # -----------------------------------------------------------------

    @property
    def total_spectral_emission_curve_path(self):
        return fs.join(self.absorption_path, "total_curve_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_emission_curve(self):
        return fs.is_file(self.total_spectral_emission_curve_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_spectral_emission_curve(self):
        if not self.has_total_spectral_emission_curve: raise RuntimeError("Cannot find the total spectral emission curve: run the emission analysis first")
        return WavelengthCurve.from_file(self.total_spectral_emission_curve_path)

    # -----------------------------------------------------------------
    #   UNEVOLVED: ABSORPTION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "unevolved_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_absorption_curve(self):
        return fs.is_file(self.unevolved_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_spectral_absorption_curve(self):
        if not self.has_unevolved_spectral_absorption_curve: raise RuntimeError("Cannot find the unevolved spectral absorption curve: run the absorption analysis first")
        return WavelengthCurve.from_file(self.unevolved_spectral_absorption_curve_path)

    # -----------------------------------------------------------------
    #   UNEVOLVED: EMISSION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_emission_curve_path(self):
        return fs.join(self.absorption_path, "unevolved_curve_emission.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_emission_curve(self):
        return fs.is_file(self.unevolved_spectral_emission_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "unevolved_spectral_emission_curve_path", True, write=True)
    def unevolved_spectral_emission_curve(self):
        if not self.has_unevolved_spectral_emission_curve: raise RuntimeError("Cannot find the unevolved spectral emission curve: run the absorption analysis first")
        return WavelengthCurve.from_file(self.unevolved_spectral_emission_curve_path)

    # -----------------------------------------------------------------
    #   UNEVOLVED: ABSORPTION FRACTION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_fraction_curve_path(self):
        return fs.join(self.curves_path, "absorption_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_absorption_fraction_curve(self):
        return fs.is_file(self.unevolved_spectral_absorption_fraction_curve_path)

    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_fraction_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_fraction_description(self):
        return "Fraction of dust absorption attributed by unevolved stellar populations"

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_absorption_curve_wavelengths(self):
        return self.unevolved_spectral_absorption_curve.wavelengths(unit=self.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def spectral_absorption_curve_unit(self):
        return self.unevolved_spectral_absorption_curve.unit

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_curve_values(self):
        #print(self.total_spectral_absorption_curve.unit, self.spectral_absorption_curve_unit)
        return self.total_spectral_absorption_curve.values(unit=self.spectral_absorption_curve_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_curve_values(self):
        return self.unevolved_spectral_absorption_curve.values(unit=self.spectral_absorption_curve_unit, asarray=True)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "unevolved_spectral_absorption_fraction_curve_path", True, write=False)
    def unevolved_spectral_absorption_fraction_curve(self):

        """
        This function ...
        :return:
        """

        # Get wavelengths and values
        wavelengths = self.spectral_absorption_curve_wavelengths
        fractions = self.unevolved_spectral_absorption_curve_values / self.total_spectral_absorption_curve_values

        # Create the curve and return
        return WavelengthCurve.from_wavelengths_and_values(self.unevolved_absorption_fraction_name, wavelengths, fractions,
                                                           wavelength_unit=self.wavelength_unit, description=self.unevolved_absorption_fraction_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_spectral_absorption_fraction_curve_values(self):
        return self.unevolved_spectral_absorption_fraction_curve.values(asarray=True)

    # -----------------------------------------------------------------
    #   UNEVOLVED: EMISSION FRACTION
    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_emission_fraction_curve_path(self):
        return fs.join(self.curves_path, "emission_cells.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_emission_fraction_curve(self):
        return fs.is_file(self.unevolved_spectral_emission_fraction_curve_path)

    # -----------------------------------------------------------------

    @property
    def unevolved_emission_fraction_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @property
    def unevolved_emission_fraction_description(self):
        return "Fraction of dust emission attributed by unevolved stellar populations"

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_emission_curve_wavelengths(self):
        return self.unevolved_spectral_emission_curve.wavelengths(unit=self.wavelength_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def spectral_emission_curve_unit(self):
        return self.unevolved_spectral_emission_curve.unit

    # -----------------------------------------------------------------

    @property
    def total_spectral_emission_curve_values(self):
        return self.total_spectral_emission_curve.values(unit=self.spectral_emission_curve_unit, asarray=True)

    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_emission_curve_values(self):
        return self.unevolved_spectral_emission_curve.values(unit=self.spectral_emission_curve_unit, asarray=True)

    # -----------------------------------------------------------------

    @lazyfileproperty(WavelengthCurve, "unevolved_spectral_emission_fraction_curve_path", True, write=False)
    def unevolved_spectral_emission_fraction_curve(self):

        """
        This function ...
        :return:
        """

        # Get wavelengths and values
        wavelengths = self.spectral_emission_curve_wavelengths
        fractions = self.unevolved_spectral_emission_curve_values / self.total_spectral_emission_curve_values

        # Create the curve and return
        return WavelengthCurve.from_wavelengths_and_values(self.unevolved_emission_fraction_name, wavelengths,
                                                           fractions, wavelength_unit=self.wavelength_unit,
                                                           description=self.unevolved_emission_fraction_description)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_spectral_emission_fraction_curve_values(self):
        return self.unevolved_spectral_emission_fraction_curve.values(asarray=True)

    # -----------------------------------------------------------------
    #   EVOLVED: ABSORPTION
    # -----------------------------------------------------------------

    @property
    def evolved_absorption_fraction_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @property
    def evolved_absorption_fraction_description(self):
        return "Fraction of dust absorption attributed by evolved stellar populations"

    # -----------------------------------------------------------------

    @lazyproperty
    def evolved_spectral_absorption_fraction_curve(self):

        """
        This function ...
        :return:
        """

        # Get wavelengths and values
        wavelengths = self.spectral_absorption_curve_wavelengths
        fractions = 1. - self.unevolved_spectral_absorption_fraction_curve_values

        # Create the curve and return
        return WavelengthCurve.from_wavelengths_and_values(self.evolved_absorption_fraction_name, wavelengths, fractions,
                                                           description=self.evolved_absorption_fraction_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    #   EVOLVED: EMISSION
    # -----------------------------------------------------------------

    @property
    def evolved_emission_fraction_name(self):
        return "Heating fraction"

    # -----------------------------------------------------------------

    @property
    def evolved_emission_fraction_description(self):
        return "Fraction of dust emission attributed by evolved stellar populations"

    # -----------------------------------------------------------------

    @lazyproperty
    def evolved_spectral_emission_fraction_curve(self):

        """
        This function ...
        :return:
        """

        # Get wavelengths and values
        wavelengths = self.spectral_emission_curve_wavelengths
        fractions = 1. - self.unevolved_spectral_emission_fraction_curve_values

        # Create the curve and return
        return WavelengthCurve.from_wavelengths_and_values(self.evolved_emission_fraction_name, wavelengths, fractions,
                                                           description=self.evolved_emission_fraction_description,
                                                           wavelength_unit=self.wavelength_unit)

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def absorption_fractions_name(self):
        return "Funev_absorption"

    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_map_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "absorption_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_absorption_fraction_map_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_absorption_fraction_map_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_absorption_fraction_map_for_filter(self, fltr):

        """
        This function ...
        :return:
        """

        if self.has_unevolved_absorption_fraction_map_for_filter(fltr): return Image.from_file(self.get_unevolved_absorption_fraction_map_path_for_filter(fltr))
        else:
            name = self.absorption_fractions_name + "_" + str(fltr)
            data = self.get_unevolved_absorption_fraction_data_for_filter(fltr)
            image = project_data(name, data, self.faceon_projection, return_stddev=True, return_ncells=True, as_image=True) # cell_based=True # NOW DEFAULT
            #image.saveto(self.get_unevolved_absorption_fraction_map_path_for_filter(fltr))
            return image

    # -----------------------------------------------------------------

    @property
    def emission_fractions_name(self):
        return "Funev_emission"

    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_map_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "emission_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_emission_fraction_map_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_emission_fraction_map_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_emission_fraction_map_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if self.has_unevolved_emission_fraction_map_for_filter(fltr): return Image.from_file(self.get_unevolved_emission_fraction_map_path_for_filter(fltr))
        else:
            name = self.emission_fractions_name + "_" + str(fltr)
            data = self.get_unevolved_emission_fraction_data_for_filter(fltr)
            image = project_data(name, data, self.faceon_projection, return_stddev=True, return_ncells=True, as_image=True) #cell_based=True # NOW DEFAULT
            return image

    # -----------------------------------------------------------------

    def fix_map(self, frame, ncells, replace_nans=True):

        """
        This function ...
        :param frame:
        :param ncells:
        :param replace_nans:
        :return:
        """

        # Copy
        fixed = frame.copy()

        # Replace invalid values by NaN
        fixed.replace_negatives_by_nans()
        fixed.replace_infs_by_nans()
        fixed.replace_by_nans_where_greater_than(1.1)
        fixed.cutoff_greater(1.)

        # Get outside nans
        outside_nans = fixed.nans.largest()
        not_nans = outside_nans.inverse()
        not_nans.disk_dilate(radius=self.config.not_nans_dilation_radius)
        do_nans = not_nans.largest().inverse()

        # Get mask
        where = ncells.where_smaller_than(self.config.min_ncells)
        where = where * do_nans.inverse() # don't interpolate outside (where ncells = 0)

        # Replace NaNs to zero that have to stay NaNs (don't interpolate)
        if replace_nans:
            fixed[do_nans] = 0.0
            do_nans.disk_dilate(radius=self.config.not_nans_dilation_radius)
        #plotting.plot_mask(where, title="where")
        #plotting.plot_mask(do_nans, title="other")

        # Put pixels to NaN
        fixed.replace_by_nans(where)

        # Interpolate nans
        fixed.interpolate_nans(sigma=2., error_on_max=replace_nans)
        fixed.replace_by_nans(do_nans)

        # Return the interpolated frame
        return fixed

    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_map_fixed_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "absorption_" + str(fltr) + "_fixed.fits")

    # -----------------------------------------------------------------

    def has_unevolved_absorption_fraction_map_fixed_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_absorption_fraction_map_fixed_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_absorption_fraction_map_fixed_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if self.has_unevolved_absorption_fraction_map_fixed_for_filter(fltr): return Frame.from_file(self.get_unevolved_absorption_fraction_map_fixed_path_for_filter(fltr))
        else:
            #name = self.absorption_fractions_name + "_" + str(fltr) + "_fixed"
            frame = self.get_unevolved_absorption_fraction_frame_for_filter(fltr)
            ncells = self.get_unevolved_absorption_fraction_ncells_for_filter(fltr)
            return self.fix_map(frame, ncells)

    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_map_fixed_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "emission_" + str(fltr) + "_fixed.fits")

    # -----------------------------------------------------------------

    def has_unevolved_emission_fraction_map_fixed_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_emission_fraction_map_fixed_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @memoize_method
    def get_unevolved_emission_fraction_map_fixed_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        if self.has_unevolved_emission_fraction_map_fixed_for_filter(fltr): return Frame.from_file(self.get_unevolved_emission_fraction_map_fixed_path_for_filter(fltr))
        else:
            frame = self.get_unevolved_emission_fraction_frame_for_filter(fltr)
            ncells = self.get_unevolved_emission_fraction_ncells_for_filter(fltr)
            return self.fix_map(frame, ncells)

    # -----------------------------------------------------------------
    # SIMULATION CUBES & SEDS
    #   EARTH
    #     EMISSION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_earth(self):
        return self.model.old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_sed_earth(self):
        return self.model.observed_old_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_earth(self):
        return self.model.has_old_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_sed_earth(self):
        return self.model.has_observed_old_dust_sed_earth

    # -----------------------------------------------------------------
    # 2. Young
    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_earth(self):
        return self.model.young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_sed_earth(self):
        return self.model.observed_young_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_earth(self):
        return self.model.has_young_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_sed_earth(self):
        return self.model.has_observed_young_dust_sed_earth

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_earth(self):
        return self.model.sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_sed_earth(self):
        return self.model.observed_sfr_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_earth(self):
        return self.model.has_sfr_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_sed_earth(self):
        return self.model.has_observed_sfr_dust_sed_earth

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_cube_earth(self):
        return self.has_young_dust_emission_cube_earth and self.has_ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_sed_earth(self):
        return self.has_young_dust_emission_sed_earth and self.has_ionizing_dust_emission_sed_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_earth(self):
        if self.model.has_unevolved_dust_luminosity_cube_earth: return self.model.unevolved_dust_luminosity_cube_earth
        elif self.has_young_and_ionizing_dust_emission_cube_earth: return self.young_dust_emission_cube_earth + self.ionizing_dust_emission_cube_earth # IS THIS EVEN CORRECT?? DUST EMISSION IS NOT LINEAR!
        else: raise IOError("Cannot obtain a cube of the dust emission from unevolved stars in the earth projection")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_sed_earth(self):
        if self.model.has_observed_unevolved_dust_sed_earth: return self.model.observed_unevolved_dust_sed_earth
        elif self.has_young_and_ionizing_dust_emission_sed_earth: return self.young_dust_emission_sed_earth + self.ionizing_dust_emission_sed_earth # IS THIS EVEN CORRECT?? DUST EMISSION IS NOT LINEAR!
        else: raise IOError("Cannot obtain an SED of the dust emission from unevolved stars in the earth projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_earth(self):
        return self.model.has_unevolved_dust_luminosity_cube_earth or self.has_young_and_ionizing_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_sed_earth(self):
        return self.model.has_observed_unevolved_dust_sed_earth or self.has_young_and_ionizing_dust_emission_sed_earth

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_earth(self):
        return self.old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_sed_earth(self):
        return self.old_dust_emission_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_earth(self):
        return self.has_old_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_sed_earth(self):
        return self.has_old_dust_emission_sed_earth

    # -----------------------------------------------------------------
    # 6. Total
    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_earth(self):
        return self.model.total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_sed_earth(self):
        return self.model.observed_total_dust_sed

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_earth(self):
        return self.model.has_total_dust_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_sed_earth(self):
        return self.model.has_observed_total_dust_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_earth(self):
        return self.has_unevolved_dust_emission_cube_earth and self.has_total_dust_emission_cube_earth and self.has_evolved_dust_emission_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_seds_earth(self):
        return self.has_unevolved_dust_emission_sed_earth and self.has_total_dust_emission_sed_earth and self.has_evolved_dust_emission_sed_earth

    # -----------------------------------------------------------------
    #     ABSORPTION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_earth(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_sed_earth(self):
        return self.model.old_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_earth(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_sed_earth(self):
        return self.model.has_old_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------
    # 2. Young
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
    def young_dust_absorption_sed_earth(self):
        return self.model.young_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_sed_earth(self):
        return self.model.has_young_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_earth(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_earth(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_sed_earth(self):
        return self.model.sfr_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_sed_earth(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_cube_earth(self):
        return self.has_young_dust_absorption_cube_earth and self.has_ionizing_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_earth(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth: return self.model.unevolved_absorbed_diffuse_stellar_luminosity_cube_earth
        elif self.has_young_and_ionizing_dust_absorption_cube_earth: return self.young_dust_absorption_cube_earth + self.ionizing_dust_absorption_cube_earth # ABSORPTION IS LINEAR SO OK
        else: raise IOError("Cannot obtain a cube of the dust absorption from unevolved stars in the earth projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_earth(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_earth or self.has_young_and_ionizing_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_sed_earth(self):
        return self.has_young_dust_absorption_sed_earth and self.has_ionizing_dust_absorption_sed_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_sed_earth(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_sed_earth: return self.model.unevolved_absorbed_diffuse_stellar_sed_earth
        elif self.has_young_and_ionizing_dust_absorption_sed_earth: return self.young_dust_absorption_sed_earth + self.ionizing_dust_absorption_sed_earth
        else: raise IOError("Cannot obtain an SED of the dust absorption from unevolved stars in the earth projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_sed_earth(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_sed_earth or self.has_young_and_ionizing_dust_absorption_sed_earth

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
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
    def evolved_dust_absorption_sed_earth(self):
        return self.old_dust_absorption_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_sed_earth(self):
        return self.has_old_dust_absorption_sed_earth

    # -----------------------------------------------------------------
    # 6. Total
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
        return self.has_unevolved_dust_absorption_cube_earth and self.has_total_dust_absorption_cube_earth and self.has_evolved_dust_absorption_cube_earth

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_sed_earth(self):
        return self.model.total_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_sed_earth(self):
        return self.model.has_total_absorbed_diffuse_stellar_sed_earth

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_seds_earth(self):
        return self.has_unevolved_dust_absorption_sed_earth and self.has_total_dust_absorption_sed_earth and self.has_evolved_dust_absorption_cube_earth

    # -----------------------------------------------------------------
    #   FACEON
    #     EMISSION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_faceon(self):
        return self.model.old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_sed_faceon(self):
        return self.model.observed_old_dust_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_faceon(self):
        return self.model.has_old_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_sed_faceon(self):
        return self.model.has_observed_old_dust_sed_faceon

    # -----------------------------------------------------------------
    # 2. Young
    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_faceon(self):
        return self.model.young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_sed_faceon(self):
        return self.model.observed_young_dust_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_faceon(self):
        return self.model.has_young_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_sed_faceon(self):
        return self.model.has_observed_young_dust_sed_faceon

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_faceon(self):
        return self.model.sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_sed_faceon(self):
        return self.model.observed_sfr_dust_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_faceon(self):
        return self.model.has_sfr_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_sed_faceon(self):
        return self.model.has_observed_sfr_dust_sed_faceon

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_cube_faceon(self):
        return self.has_young_dust_emission_cube_faceon and self.has_ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_sed_faceon(self):
        return self.has_young_dust_emission_sed_faceon and self.has_ionizing_dust_emission_sed_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_faceon(self):
        if self.model.has_unevolved_dust_luminosity_cube_faceon: return self.model.unevolved_dust_luminosity_cube_faceon
        elif self.has_young_and_ionizing_dust_emission_cube_faceon: return self.young_dust_emission_cube_faceon + self.ionizing_dust_emission_cube_faceon # IS THIS CORRECT?? NOT LINEAR
        else: raise IOError("Cannot obtain a cube of the dust emission from unevolved stars in the faceon projection")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_sed_faceon(self):
        if self.model.has_observed_unevolved_dust_sed_faceon: return self.model.unevolved_dust_sed_faceon
        elif self.has_young_and_ionizing_dust_emission_sed_faceon: return self.young_dust_emission_sed_faceon + self.ionizing_dust_emission_sed_faceon
        else: raise IOError("Cannot obtain an SED of the dust emission from unevolved stars in the faceon projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_faceon(self):
        return self.has_young_dust_emission_cube_faceon and self.has_ionizing_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_sed_faceon(self):
        return self.has_young_dust_emission_sed_earth and self.has_ionizing_dust_emission_sed_faceon

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_faceon(self):
        return self.old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_sed_faceon(self):
        return self.old_dust_emission_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_faceon(self):
        return self.has_old_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_sed_faceon(self):
        return self.has_old_dust_emission_sed_faceon

    # -----------------------------------------------------------------
    # 6. Total
    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_faceon(self):
        return self.model.total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_sed_faceon(self):
        return self.model.observed_total_dust_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_faceon(self):
        return self.model.has_total_dust_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_sed_faceon(self):
        return self.model.has_observed_total_dust_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_faceon(self):
        return self.has_unevolved_dust_emission_cube_faceon and self.has_total_dust_emission_cube_faceon and self.has_evolved_dust_emission_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_seds_faceon(self):
        return self.has_unevolved_dust_emission_sed_faceon and self.has_total_dust_emission_sed_faceon and self.has_evolved_dust_emission_sed_faceon

    # -----------------------------------------------------------------
    #     ABSORPTION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_faceon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_sed_faceon(self):
        return self.model.old_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_faceon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_sed_faceon(self):
        return self.model.has_old_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------
    # 2. Young
    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_faceon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_sed_faceon(self):
        return self.model.young_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_faceon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_sed_faceon(self):
        return self.model.has_young_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_faceon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_sed_faceon(self):
        return self.model.sfr_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_faceon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_sed_faceon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_cube_faceon(self):
        return self.has_young_dust_absorption_cube_faceon and self.has_ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_sed_faceon(self):
        return self.has_young_dust_absorption_sed_faceon and self.has_ionizing_dust_absorption_sed_faceon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_faceon(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon: return self.model.unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon
        elif self.has_young_and_ionizing_dust_absorption_cube_faceon: return self.young_dust_absorption_cube_faceon + self.ionizing_dust_absorption_cube_faceon
        else: raise IOError("Cannot obtain a cube of the dust absorption from unevolved stars in the faceon projection")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_sed_faceon(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_sed_faceon: return self.model.unevolved_absorbed_diffuse_stellar_sed_faceon
        elif self.has_young_and_ionizing_dust_absorption_sed_faceon: return self.young_dust_absorption_sed_faceon + self.ionizing_dust_absorption_sed_faceon # LINEAR
        else: raise IOError("Cannot obtain an SED of the dust absorption from unevolved stars in the faceon projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_faceon(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_faceon or self.has_young_and_ionizing_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_sed_faceon(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_sed_faceon or self.has_young_and_ionizing_dust_absorption_sed_faceon

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_faceon(self):
        return self.old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_sed_faceon(self):
        return self.old_dust_absorption_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_faceon(self):
        return self.has_old_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_sed_faceon(self):
        return self.has_old_dust_absorption_sed_faceon

    # -----------------------------------------------------------------
    # 6. Total
    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_faceon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_sed_faceon(self):
        return self.model.total_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_faceon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_sed_faceon(self):
        return self.model.has_total_absorbed_diffuse_stellar_sed_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_faceon(self):
        return self.has_unevolved_dust_absorption_cube_faceon and self.has_total_dust_absorption_cube_faceon and self.has_evolved_dust_absorption_cube_faceon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_seds_faceon(self):
        return self.has_unevolved_dust_absorption_sed_faceon and self.has_total_dust_absorption_sed_faceon and self.has_evolved_dust_absorption_sed_faceon

    # -----------------------------------------------------------------
    #   EDGEON
    #     EMISSION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_emission_cube_edgeon(self):
        return self.model.old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_dust_emission_sed_edgeon(self):
        return self.model.observed_old_dust_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_cube_edgeon(self):
        return self.model.has_old_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_emission_sed_edgeon(self):
        return self.model.has_observed_old_dust_sed_edgeon

    # -----------------------------------------------------------------
    # 2. Young
    # -----------------------------------------------------------------

    @property
    def young_dust_emission_cube_edgeon(self):
        return self.model.young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_emission_sed_edgeon(self):
        return self.model.observed_young_dust_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_cube_edgeon(self):
        return self.model.has_young_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_emission_sed_edgeon(self):
        return self.model.has_observed_young_dust_sed_edgeon

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_cube_edgeon(self):
        return self.model.sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_emission_sed_edgeon(self):
        return self.model.observed_sfr_dust_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_cube_edgeon(self):
        return self.model.has_sfr_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_emission_sed_edgeon(self):
        return self.model.has_observed_sfr_dust_sed_edgeon

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_cube_edgeon(self):
        return self.has_young_dust_emission_cube_edgeon and self.has_ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_emission_sed_edgeon(self):
        return self.has_young_dust_emission_sed_edgeon and self.has_ionizing_dust_emission_sed_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_cube_edgeon(self):
        if self.model.has_unevolved_dust_luminosity_cube_edgeon: return self.model.unevolved_dust_luminosity_cube_edgeon
        elif self.has_young_and_ionizing_dust_emission_cube_edgeon: return self.young_dust_emission_cube_edgeon + self.ionizing_dust_emission_cube_edgeon
        else: raise IOError("Cannot obtain a cube of the dust emission from unevolved stars in the edgeon projection")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_emission_sed_edgeon(self):
        if self.model.has_observed_unevolved_dust_sed_edgeon: return self.model.observed_unevolved_dust_sed_edgeon
        elif self.has_young_and_ionizing_dust_emission_sed_edgeon: return self.young_dust_emission_sed_edgeon + self.ionizing_dust_emission_sed_edgeon
        else: raise IOError("Cannot obtain an SED of the dust emission from unevolved stars in the edgeon projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_cube_edgeon(self):
        return self.has_young_dust_emission_cube_edgeon and self.has_ionizing_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_emission_sed_edgeon(self):
        return self.has_young_dust_emission_sed_edgeon and self.has_ionizing_dust_emission_sed_edgeon

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_cube_edgeon(self):
        return self.old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_emission_sed_edgeon(self):
        return self.old_dust_emission_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_cube_edgeon(self):
        return self.has_old_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_emission_sed_edgeon(self):
        return self.has_old_dust_emission_sed_edgeon

    # -----------------------------------------------------------------
    # 6. Total
    # -----------------------------------------------------------------

    @property
    def total_dust_emission_cube_edgeon(self):
        return self.model.total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_emission_sed_edgeon(self):
        return self.model.observed_total_dust_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_cube_edgeon(self):
        return self.model.has_total_dust_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_emission_sed_edgeon(self):
        return self.model.has_observed_total_dust_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_cubes_edgeon(self):
        return self.has_unevolved_dust_emission_cube_edgeon and self.has_total_dust_emission_cube_edgeon and self.has_evolved_dust_emission_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_emission_seds_edgeon(self):
        return self.has_unevolved_dust_emission_sed_edgeon and self.has_total_dust_emission_sed_edgeon and self.has_evolved_dust_emission_sed_edgeon

    # -----------------------------------------------------------------
    #     ABSORPTION
    # -----------------------------------------------------------------
    # 1. Old
    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_cube_edgeon(self):
        return self.model.old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def old_dust_absorption_sed_edgeon(self):
        return self.model.old_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_cube_edgeon(self):
        return self.model.has_old_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_old_dust_absorption_sed_edgeon(self):
        return self.model.has_old_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------
    # 2. Young
    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_cube_edgeon(self):
        return self.model.young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def young_dust_absorption_sed_edgeon(self):
        return self.model.young_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_cube_edgeon(self):
        return self.model.has_young_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_dust_absorption_sed_edgeon(self):
        return self.model.has_young_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------
    # 3. Ionizing
    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_cube_edgeon(self):
        return self.model.sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def ionizing_dust_absorption_sed_edgeon(self):
        return self.model.sfr_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_cube_edgeon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_ionizing_dust_absorption_sed_edgeon(self):
        return self.model.has_sfr_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------
    # 4. Unevolved
    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_cube_edgeon(self):
        return self.has_young_dust_absorption_cube_edgeon and self.has_ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_young_and_ionizing_dust_absorption_sed_edgeon(self):
        return self.has_young_dust_absorption_sed_edgeon and self.has_ionizing_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_cube_edgeon(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon: return self.model.unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon
        elif self.has_young_and_ionizing_dust_absorption_cube_edgeon: return self.young_dust_absorption_cube_edgeon + self.ionizing_dust_absorption_cube_edgeon
        else: raise IOError("Cannot obtain a cube of the dust absorption from unevolved stars in the edgeon projection")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_absorption_sed_edgeon(self):
        if self.model.has_unevolved_absorbed_diffuse_stellar_sed_edgeon: return self.model.unevolved_absorbed_diffuse_stellar_sed_edgeon
        elif self.has_young_and_ionizing_dust_absorption_sed_edgeon: return self.young_dust_absorption_sed_edgeon + self.ionizing_dust_absorption_sed_edgeon
        else: raise IOError("Cannot obtain an SED of the dust absorption from unevolved stars in the edgeon projection")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_cube_edgeon(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_luminosity_cube_edgeon or self.has_young_and_ionizing_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_unevolved_dust_absorption_sed_edgeon(self):
        return self.model.has_unevolved_absorbed_diffuse_stellar_sed_edgeon or self.has_young_and_ionizing_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------
    # 5. Evolved ( == old)
    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_cube_edgeon(self):
        return self.old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def evolved_dust_absorption_sed_edgeon(self):
        return self.old_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_cube_edgeon(self):
        return self.has_old_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_evolved_dust_absorption_sed_edgeon(self):
        return self.has_old_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------
    # 6. Total
    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_cube_edgeon(self):
        return self.model.total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def total_dust_absorption_sed_edgeon(self):
        return self.model.total_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_cube_edgeon(self):
        return self.model.has_total_absorbed_diffuse_stellar_luminosity_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_total_dust_absorption_sed_edgeon(self):
        return self.model.has_total_absorbed_diffuse_stellar_sed_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_cubes_edgeon(self):
        return self.has_unevolved_dust_absorption_cube_edgeon and self.has_total_dust_absorption_cube_edgeon and self.has_evolved_dust_absorption_cube_edgeon

    # -----------------------------------------------------------------

    @property
    def has_dust_absorption_seds_edgeon(self):
        return self.has_unevolved_dust_absorption_sed_edgeon and self.has_total_dust_absorption_sed_edgeon and self.has_evolved_dust_absorption_sed_edgeon

    # -----------------------------------------------------------------
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

        # Write the maps
        self.write_maps()

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth_emission(self):
        return self.do_cubes_earth_emission and not self.has_cube_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_write_cube_earth_absorption(self):
        return self.do_cubes_earth_absorption and not self.has_cube_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon_emission(self):
        return self.do_cubes_faceon_emission and not self.has_cube_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_write_cube_faceon_absorption(self):
        return self.do_cubes_faceon_absorption and not self.has_cube_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_cube_edgeon_emission(self):
        return self.do_cubes_edgeon_emission and not self.has_cube_edgeon_emission

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
    def do_write_curve_earth_emission(self):
        return self.do_curve_earth_emission and not self.has_curve_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth_absorption(self):
        return self.do_curve_earth_absorption and not self.has_curve_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_emission(self):
        return self.do_curve_faceon_emission and not self.has_curve_faceon_emission

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_absorption(self):
        return self.do_curve_faceon_absorption and not self.has_curve_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon_emission(self):
        return self.do_curve_edgeon_emission and not self.has_curve_edgeon_emission

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

        # Absorption from SEDs
        self.write_curves_absorption_seds()

        # Emission from SEDs (to be expected equal)
        self.write_curves_emission_seds()

        # Absorption data
        if self.do_absorption_data: self.write_curves_absorption_data()

        # Emission data
        if self.do_emission_data: self.write_curves_emission_data()

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

        # Debugging
        log.debug("Writing the curve of dust emission from the earth projection ...")

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

        # Debugging
        log.debug("Writing the curve of dust absorption from the earth projection ...")

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

        # Debugging
        log.debug("Writing the curve of dust emission from the faceon projection ...")

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

        # Debugging
        log.debug("Writing the curve of dust absorption from the faceon projection ...")

        # Save
        self.curve_faceon_absorption.saveto(self.curve_faceon_absorption_path)

    # -----------------------------------------------------------------

    def remove_curve_edgeon_emission(self):
        fs.remove_file(self.curve_edgeon_emission_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon_emission(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Writing the curve of dust emission from the edgeon projection ...")

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

        # Debugging
        log.debug("Writing the curve of dust absorption from the edgeon projection ...")

        # Save
        self.curve_edgeon_absorption.saveto(self.curve_edgeon_absorption_path)

    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth_absorption_seds(self):
        return not self.has_curve_earth_absorption_seds

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_absorption_seds(self):
        return not self.has_curve_faceon_absorption_seds

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon_absorption_seds(self):
        return not self.has_curve_edgeon_absorption_seds

    # -----------------------------------------------------------------

    def write_curves_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption curves from SEDs ...")

        # Earth
        if self.do_write_curve_earth_absorption_seds: self.write_curve_earth_absorption_seds()

        # Face-on
        if self.do_write_curve_faceon_absorption_seds: self.write_curve_faceon_absorption_seds()

        # Edge-on
        if self.do_write_curve_edgeon_absorption_seds: self.write_curve_edgeon_absorption_seds()

    # -----------------------------------------------------------------

    def write_curve_earth_absorption_seds(self):
        self.curve_earth_absorption_seds.saveto(self.curve_earth_absorption_seds_path)

    # -----------------------------------------------------------------

    def write_curve_faceon_absorption_seds(self):
        self.curve_faceon_absorption_seds.saveto(self.curve_faceon_absorption_seds_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon_absorption_seds(self):
        self.curve_edgeon_absorption_seds.saveto(self.curve_edgeon_absorption_seds_path)

    # -----------------------------------------------------------------

    @property
    def do_write_curve_earth_emission_seds(self):
        return not self.has_curve_earth_emission_seds

    # -----------------------------------------------------------------

    @property
    def do_write_curve_faceon_emission_seds(self):
        return not self.has_curve_faceon_emission_seds

    # -----------------------------------------------------------------

    @property
    def do_write_curve_edgeon_emission_seds(self):
        return not self.has_curve_edgeon_emission_seds

    # -----------------------------------------------------------------

    def write_curves_emission_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the emission curves from SEDs ...")

        # Earth
        if self.do_write_curve_earth_emission_seds: self.write_curve_earth_emission_seds()

        # Face-on
        if self.do_write_curve_faceon_emission_seds: self.write_curve_faceon_emission_seds()

        # Edge-on
        if self.do_write_curve_edgeon_emission_seds: self.write_curve_edgeon_emission_seds()

    # -----------------------------------------------------------------

    def write_curve_earth_emission_seds(self):
        self.curve_earth_emission_seds.saveto(self.curve_earth_emission_seds_path)

    # -----------------------------------------------------------------

    def write_curve_faceon_emission_seds(self):
        self.curve_faceon_emission_seds.saveto(self.curve_faceon_emission_seds_path)

    # -----------------------------------------------------------------

    def write_curve_edgeon_emission_seds(self):
        self.curve_edgeon_emission_seds.saveto(self.curve_edgeon_emission_seds_path)

    # -----------------------------------------------------------------

    def write_curves_emission_data(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Writing the curve of dust emission from the cell data ...")

        # Save
        #self.unevolved_spectral_emission_curve.saveto(self.unevolved_spectral_emission_curve_path)
        self.unevolved_spectral_emission_fraction_curve.saveto(self.unevolved_spectral_emission_fraction_curve_path)

    # -----------------------------------------------------------------

    def write_curves_absorption_data(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Writing the curve of dust absorption from the cell data ...")

        # Save
        #self.unevolved_spectral_absorption_curve.saveto(self.unevolved_spectral_absorption_curve_path)
        self.unevolved_spectral_absorption_fraction_curve.saveto(self.unevolved_spectral_absorption_fraction_curve_path)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of the spectral heating ...")

        # Absorption
        self.write_maps_absorption()

        # Emission
        self.write_maps_emission()

        # Absorption data
        if self.do_absorption_data: self.write_maps_absorption_data()

        # Emission data
        if self.do_emission_data: self.write_maps_emission_data()

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_absorption_earth(self):
        return self.do_cubes_earth_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_absorption_faceon(self):
        return self.do_cubes_faceon_absorption

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_absorption_edgeon(self):
        return self.do_cubes_edgeon_absorption

    # -----------------------------------------------------------------

    def write_maps_absorption(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption maps for different filters ...")

        # Earth
        if self.do_write_spectral_maps_absorption_earth: self.write_spectral_maps_absorption_earth()

        # Face-on
        if self.do_write_spectral_maps_absorption_faceon: self.write_spectral_maps_absorption_faceon()

        # Edge-on
        if self.do_write_spectral_maps_absorption_edgeon: self.write_spectral_maps_absorption_edgeon()

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_earth(self, fltr):
        frame = self.cube_earth_absorption_fixed.frame_for_filter(fltr, convolve=self.config.spectral_convolution)
        frame.set_meta("SPECCON", self.config.spectral_convolution)
        return frame

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_earth_path(self, fltr):
        return fs.join(self.maps_path, "earth_absorption_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_earth(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_earth_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_absorption_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust absorption from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_earth(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_earth_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_faceon(self, fltr):
        frame = self.cube_faceon_absorption_fixed.frame_for_filter(fltr, convolve=self.config.spectral_convolution)
        frame.set_meta("SPECCON", self.config.spectral_convolution)
        return frame

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_faceon_path(self, fltr):
        return fs.join(self.maps_path, "faceon_absorption_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_faceon(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_faceon_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_absorption_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust absorption from the faceon projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_faceon_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_absorption_edgeon(self, fltr):
        frame = self.cube_edgeon_absorption_fixed.frame_for_filter(fltr, convolve=self.config.spectral_convolution)
        frame.set_meta("SPECCON", self.config.spectral_convolution)
        return frame

    # -----------------------------------------------------------------

    def get_spectral_map_absorption_edgeon_path(self, fltr):
        return fs.join(self.maps_path, "edgeon_absorption_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_edgeon(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_absorption_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust absorption from the edgeon projection ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_spectral_map_absorption_edgeon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_edgeon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_edgeon_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_emission_earth(self):
        return self.do_cubes_earth_emission

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_emission_edgeon(self):
        return self.do_cubes_edgeon_emission

    # -----------------------------------------------------------------

    @property
    def do_write_spectral_maps_emission_faceon(self):
        return self.do_cubes_faceon_emission

    # -----------------------------------------------------------------

    def write_maps_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the emission maps for different filters ...")

        # Earth
        if self.do_write_spectral_maps_emission_earth: self.write_spectral_maps_emission_earth()

        # Face-on
        if self.do_write_spectral_maps_emission_faceon: self.write_spectral_maps_emission_faceon()

        # Edge-on
        if self.do_write_spectral_maps_emission_edgeon: self.write_spectral_maps_emission_edgeon()

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_earth(self, fltr):
        return self.cube_earth_emission_fixed.frame_for_filter(fltr, convolve=self.config.spectral_convolution)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_earth_path(self, fltr):
        return fs.join(self.maps_path, "earth_emission_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_earth(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_earth_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_emission_earth(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust emission from the earth projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_earth(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_earth(fltr)

            # Get the path
            path = self.get_spectral_map_emission_earth_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_faceon(self, fltr):
        return self.cube_faceon_emission_fixed.frame_for_filter(fltr, convolve=False)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_faceon_path(self, fltr):
        return fs.join(self.maps_path, "faceon_emission_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_faceon(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_faceon_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_emission_faceon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust emission from the faceon projection ..")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_faceon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_emission_faceon_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    @memoize_method
    def get_spectral_map_emission_edgeon(self, fltr):
        return self.cube_edgeon_emission_fixed.frame_for_filter(fltr, convolve=self.config.spectral_convolution)

    # -----------------------------------------------------------------

    def get_spectral_map_emission_edgeon_path(self, fltr):
        return fs.join(self.maps_path, "edgeon_emission_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_edgeon(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_edgeon_path(fltr))

    # -----------------------------------------------------------------

    def write_spectral_maps_emission_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing maps of the spectral heating by dust emission from the edgeon projection ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_spectral_map_emission_edgeon(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_edgeon(fltr)

            # Get the path
            path = self.get_spectral_map_emission_edgeon_path(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------
    # ABSORPTION FROM DATA
    # -----------------------------------------------------------------

    # for plotting
    def get_unevolved_absorption_fraction_frame_for_filter(self, fltr):
        return self.get_unevolved_absorption_fraction_map_for_filter(fltr).primary

    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_ncells_for_filter(self, fltr):
        return self.get_unevolved_absorption_fraction_map_for_filter(fltr).frames["ncells"]

    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_stddev_for_filter(self, fltr):
        return self.get_unevolved_absorption_fraction_map_for_filter(fltr).frames["stddev"]

    # -----------------------------------------------------------------

    def write_maps_absorption_data(self):

        """
        This function ...
        :return:
        """

        # Absorption from 3D data
        self.write_maps_absorption_data_fractions()

        # Absorption from 3D data, fixed
        self.write_maps_absorption_data_fixed()

        # Absorption differences
        self.write_maps_absorption_differences()

    # -----------------------------------------------------------------

    def write_maps_absorption_data_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of dust absorption heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_unevolved_absorption_fraction_map_for_filter(fltr): continue

            # Get the map
            image = self.get_unevolved_absorption_fraction_map_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_absorption_fraction_map_path_for_filter(fltr)

            # Write
            image.saveto(path)

    # -----------------------------------------------------------------

    def write_maps_absorption_data_fixed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fixed maps of dust absorption heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if fixed map exists
            if self.has_unevolved_absorption_fraction_map_fixed_for_filter(fltr): continue

            # Get the map
            frame = self.get_unevolved_absorption_fraction_map_fixed_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_absorption_fraction_map_fixed_path_for_filter(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_differences_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "differences_absorption_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_absorption_fraction_differences_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_absorption_fraction_differences_path_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_maps_absorption_differences(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the difference between the dust absorption heating fraction maps from the 3D data and the projected face-on view ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if difference map exists
            if self.has_unevolved_absorption_fraction_differences_for_filter(fltr): continue

            # Get the projected map
            projected_map = self.get_spectral_map_absorption_faceon(fltr)

            # Get the map from the data
            data_map = self.get_unevolved_absorption_fraction_frame_for_filter(fltr)

            # Calculate ratio
            ratio = data_map / projected_map

            # Get the path
            path = self.get_unevolved_absorption_fraction_differences_path_for_filter(fltr)

            # Write ratio
            ratio.saveto(path)

    # -----------------------------------------------------------------
    # EMISSION FROM DATA
    # -----------------------------------------------------------------

    # for plotting
    def get_unevolved_emission_fraction_frame_for_filter(self, fltr):
        return self.get_unevolved_emission_fraction_map_for_filter(fltr).primary

    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_ncells_for_filter(self, fltr):
        return self.get_unevolved_emission_fraction_map_for_filter(fltr).frames["ncells"]

    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_stddev_for_filter(self, fltr):
        return self.get_unevolved_emission_fraction_map_for_filter(fltr).frames["stddev"]

    # -----------------------------------------------------------------

    def write_maps_emission_data(self):

        """
        This function ...
        :return:
        """

        # Emission from 3D data
        self.write_maps_emission_data_fractions()

        # Emission from 3D data, fixed
        self.write_maps_emission_data_fixed()

        # Emission differences
        self.write_maps_emission_differences()

    # -----------------------------------------------------------------

    def write_maps_emission_data_fractions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps of dust emission heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_unevolved_emission_fraction_map_for_filter(fltr): continue

            # Get the map
            image = self.get_unevolved_emission_fraction_map_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_emission_fraction_map_path_for_filter(fltr)

            # Write
            image.saveto(path)

    # -----------------------------------------------------------------

    def write_maps_emission_data_fixed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the fixed maps of dust emission heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if fixed map exists
            if self.has_unevolved_emission_fraction_map_fixed_for_filter(fltr): continue

            # Get the map
            frame = self.get_unevolved_emission_fraction_map_fixed_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_emission_fraction_map_fixed_path_for_filter(fltr)

            # Write
            frame.saveto(path)

    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_differences_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "differences_emission_" + str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def has_unevolved_emission_fraction_differences_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_emission_fraction_differences_path_for_filter(fltr))

    # -----------------------------------------------------------------

    def write_maps_emission_differences(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the difference between the dust emission heating fraction maps from the 3D data and the projected face-on view ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if difference map exists
            if self.has_unevolved_emission_fraction_differences_for_filter(fltr): continue

            # Get the projected map
            projected_map = self.get_spectral_map_emission_faceon(fltr)

            # Get the map from the data
            data_map = self.get_unevolved_emission_fraction_frame_for_filter(fltr)

            # Calculate ratio
            ratio = data_map / projected_map

            # Get the path
            path = self.get_unevolved_emission_fraction_differences_path_for_filter(fltr)

            # Write ratio
            ratio.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps of spectral heating
        self.plot_maps()

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

    @property
    def fraction_limits(self):
        return (0,1,)

    # -----------------------------------------------------------------

    def plot_absorption_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        plotting.plot_map(frame, path=path, interval=self.fraction_limits, background_color="black")

    # -----------------------------------------------------------------

    def plot_emission_map(self, frame, path):

        """
        This function ...
        :param frame:
        :param path:
        :return:
        """

        plotting.plot_map(frame, path=path, interval=self.fraction_limits, background_color="black")

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of the spectral heating ...")

        # Absorption
        self.plot_maps_absorption()

        # Emission
        self.plot_maps_emission()

        # Absorption data
        if self.do_absorption_data: self.plot_maps_absorption_data()

        # Emission data
        if self.do_emission_data: self.plot_maps_emission_data()

    # -----------------------------------------------------------------

    def plot_maps_absorption(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption maps ...")

        # Earth
        if self.do_plot_spectral_maps_absorption_earth: self.plot_spectral_maps_absorption_earth()

        # Face-on
        if self.do_plot_spectral_maps_absorption_faceon: self.plot_spectral_maps_absorption_faceon()

        # Edge-on
        if self.do_plot_spectral_maps_absorption_edgeon: self.plot_spectral_maps_absorption_edgeon()

    # -----------------------------------------------------------------

    def plot_maps_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the emission maps ...")

        # Earth
        if self.do_plot_spectral_maps_emission_earth: self.plot_spectral_maps_emission_earth()

        # Face-on
        if self.do_plot_spectral_maps_emission_faceon: self.plot_spectral_maps_emission_faceon()

        # Edge-on
        if self.do_plot_spectral_maps_emission_edgeon: self.plot_spectral_maps_emission_edgeon()

    # -----------------------------------------------------------------
    # EARTH
    #   EMISSION
    # -----------------------------------------------------------------

    def get_spectral_map_emission_earth_plot_path(self, fltr):
        return fs.join(self.maps_path, "earth_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_earth_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_earth_plot_path(fltr))

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
            if self.has_spectral_map_emission_earth_plot(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_earth(fltr)

            # Get the path
            path = self.get_spectral_map_emission_earth_plot_path(fltr)

            # Plot
            self.plot_emission_map(frame, path)

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    def get_spectral_map_absorption_earth_plot_path(self, fltr):
        return fs.join(self.maps_path, "earth_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_earth_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_earth_plot_path(fltr))

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
            if self.has_spectral_map_absorption_earth_plot(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_earth(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_earth_plot_path(fltr)

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------
    # FACEON
    #   EMISSION
    # -----------------------------------------------------------------

    def get_spectral_map_emission_faceon_plot_path(self, fltr):
        return fs.join(self.maps_path, "faceon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_faceon_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_faceon_plot_path(fltr))

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
            if self.has_spectral_map_emission_faceon_plot(fltr): continue

            # Get the map
            frame = self.get_spectral_map_emission_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_emission_faceon_plot_path(fltr)

            # Plot
            self.plot_emission_map(frame, path)

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    def get_spectral_map_absorption_faceon_plot_path(self, fltr):
        return fs.join(self.maps_path, "faceon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_faceon_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_faceon_plot_path(fltr))

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
            if self.has_spectral_map_absorption_faceon_plot(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_faceon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_faceon_plot_path(fltr)

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------
    # EDGEON
    #   EMISSION
    # -----------------------------------------------------------------

    def get_spectral_map_emission_edgeon_plot_path(self, fltr):
        return fs.join(self.maps_path, "edgeon_emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_emission_edgeon_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_emission_edgeon_plot_path(fltr))

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
            if self.has_spectral_map_emission_edgeon_plot(fltr): continue
            
            # Get the map
            frame = self.get_spectral_map_emission_edgeon(fltr)
            
            # Get the path
            path = self.get_spectral_map_emission_edgeon_plot_path(fltr)
            
            # Plot
            self.plot_emission_map(frame, path)

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    def get_spectral_map_absorption_edgeon_plot_path(self, fltr):
        return fs.join(self.maps_path, "edgeon_absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_spectral_map_absorption_edgeon_plot(self, fltr):
        return fs.is_file(self.get_spectral_map_absorption_edgeon_plot_path(fltr))

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
            if self.has_spectral_map_absorption_edgeon_plot(fltr): continue

            # Get the map
            frame = self.get_spectral_map_absorption_edgeon(fltr)

            # Get the path
            path = self.get_spectral_map_absorption_edgeon_plot_path(fltr)

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------
    # ABSORPTION DATA MAPS
    # -----------------------------------------------------------------

    def get_unevolved_absorption_fraction_map_plot_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "absorption_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_unevolved_absorption_fraction_map_plot_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_absorption_fraction_map_plot_path_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_maps_absorption_data(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of dust absorption heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.absorption_filters:

            # Check if map exists
            if self.has_unevolved_absorption_fraction_map_plot_for_filter(fltr): continue

            # Get the frame
            #frame = self.get_unevolved_absorption_fraction_frame_for_filter(fltr)
            frame = self.get_unevolved_absorption_fraction_map_fixed_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_absorption_fraction_map_plot_path_for_filter(fltr)

            # Plot
            self.plot_absorption_map(frame, path)

    # -----------------------------------------------------------------
    # EMISSION DATA MAPS
    # -----------------------------------------------------------------

    def get_unevolved_emission_fraction_map_plot_path_for_filter(self, fltr):
        return fs.join(self.cells_path, "emission_" + str(fltr) + ".pdf")

    # -----------------------------------------------------------------

    def has_unevolved_emission_fraction_map_plot_for_filter(self, fltr):
        return fs.is_file(self.get_unevolved_emission_fraction_map_plot_path_for_filter(fltr))

    # -----------------------------------------------------------------

    def plot_maps_emission_data(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps of dust emission heating fraction from the 3D data for different filters ...")

        # Loop over the filters
        for fltr in self.config.emission_filters:

            # Check if map exists
            if self.has_unevolved_emission_fraction_map_plot_for_filter(fltr): continue

            # Get the frame
            #frame = self.get_unevolved_absorption_fraction_frame_for_filter(fltr)
            frame = self.get_unevolved_emission_fraction_map_fixed_for_filter(fltr)

            # Get the path
            path = self.get_unevolved_emission_fraction_map_plot_path_for_filter(fltr)

            # Plot
            self.plot_emission_map(frame, path)

    # -----------------------------------------------------------------
    # CURVES
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

        # Absorption from SEDs
        self.plot_curves_absorption_seds()

        # Emission from SEDs
        self.plot_curves_emission_seds()

        # Absorption data
        if self.do_absorption_data: self.plot_curves_absorption_data()

        # Emission data
        if self.do_emission_data: self.plot_curves_emission_data()

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
        return fs.join(self.curves_path, "earth_emission.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_plot(self):
        return fs.is_file(self.curve_earth_emission_plot_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_emission(self):
        return {unevolved_name: self.curve_earth_emission, evolved_name: self.curve_earth_emission_evolved}

    # -----------------------------------------------------------------

    @property
    def curves_emission_y_label(self):
        return "Fraction of emitted dust luminosity"

    # -----------------------------------------------------------------

    def plot_curve_earth_emission(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth_emission, path=self.curve_earth_emission_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

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

    @property
    def curves_absorption_y_label(self):
        return "Fraction of absorbed luminosity"

    # -----------------------------------------------------------------

    @property
    def curve_ylimits(self):
        return (0,1,)

    # -----------------------------------------------------------------

    def plot_curve_earth_absorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the earth projection ...")

        # Plot
        plotting.plot_curves(self.curves_earth_absorption, path=self.curve_earth_absorption_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_plot_path(self):
        return fs.join(self.curves_path, "faceon_emission.pdf")

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
        plotting.plot_curves(self.curves_faceon_emission, path=self.curve_faceon_emission_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

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
        plotting.plot_curves(self.curves_faceon_absorption, path=self.curve_faceon_absorption_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_plot_path(self):
        return fs.join(self.curves_path, "edgeon_emission.pdf")

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
        plotting.plot_curves(self.curves_edgeon_emission, path=self.curve_edgeon_emission_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

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
        plotting.plot_curves(self.curves_edgeon_absorption, path=self.curve_edgeon_absorption_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_earth_absorption_seds(self):
        return not self.has_curve_earth_absorption_seds_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_absorption_seds(self):
        return not self.has_curve_faceon_absorption_seds_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_absorption_seds(self):
        return not self.has_curve_edgeon_absorption_seds_plot

    # -----------------------------------------------------------------

    def plot_curves_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption curves from SEDs ...")

        # Earth
        if self.do_plot_curve_earth_absorption_seds: self.plot_curve_earth_absorption_seds()

        # Face-on
        if self.do_plot_curve_faceon_absorption_seds: self.plot_curve_faceon_absorption_seds()

        # Edge-on
        if self.do_plot_curve_edgeon_absorption_seds: self.plot_curve_edgeon_absorption_seds()

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_absorption_seds(self):
        return {unevolved_name: self.curve_earth_absorption_seds, evolved_name: self.curve_earth_absorption_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_earth_absorption_seds_plot_path(self):
        return fs.join(self.curves_path, "earth_absorption_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_absorption_seds_plot(self):
        return fs.is_file(self.curve_earth_absorption_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_earth_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_earth_absorption_seds, path=self.curve_earth_absorption_seds_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon_absorption_seds(self):
        return {unevolved_name: self.curve_faceon_absorption_seds, evolved_name: self.curve_faceon_absorption_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_faceon_absorption_seds_plot_path(self):
        return fs.join(self.curves_path, "faceon_absorption_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_absorption_seds_plot(self):
        return fs.is_file(self.curve_faceon_absorption_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_faceon_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_faceon_absorption_seds, path=self.curve_faceon_absorption_seds_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon_absorption_seds(self):
        return {unevolved_name: self.curve_edgeon_absorption_seds, evolved_name: self.curve_edgeon_absorption_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_absorption_seds_plot_path(self):
        return fs.join(self.curves_path, "edgeon_absorption_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_absorption_seds_plot(self):
        return fs.is_file(self.curve_edgeon_absorption_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_edgeon_absorption_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_edgeon_absorption_seds, path=self.curve_edgeon_absorption_seds_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_earth_emission_seds(self):
        return not self.has_curve_earth_emission_seds_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_faceon_emission_seds(self):
        return not self.has_curve_faceon_emission_seds_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_curve_edgeon_emission_seds(self):
        return not self.has_curve_edgeon_emission_seds_plot

    # -----------------------------------------------------------------

    def plot_curves_emission_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the emission curves from SEDs ...")

        # Earth
        if self.do_plot_curve_earth_emission_seds: self.plot_curve_earth_emission_seds()

        # Face-on
        if self.do_plot_curve_faceon_emission_seds: self.plot_curve_faceon_emission_seds()

        # Edge-on
        if self.do_plot_curve_edgeon_emission_seds: self.plot_curve_edgeon_emission_seds()

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_earth_emission_seds(self):
        return {unevolved_name: self.curve_earth_emission_seds, evolved_name: self.curve_earth_emission_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_earth_emission_seds_plot_path(self):
        return fs.join(self.curves_path, "earth_emission_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_earth_emission_seds_plot(self):
        return fs.is_file(self.curve_earth_emission_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_earth_emission_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_earth_emission_seds, path=self.curve_earth_emission_seds_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_faceon_emission_seds(self):
        return {unevolved_name: self.curve_faceon_emission_seds, evolved_name: self.curve_faceon_emission_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_faceon_emission_seds_plot_path(self):
        return fs.join(self.curves_path, "faceon_emission_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_faceon_emission_seds_plot(self):
        return fs.is_file(self.curve_faceon_emission_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_faceon_emission_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_faceon_emission_seds, path=self.curve_faceon_emission_seds_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_edgeon_emission_seds(self):
        return {unevolved_name: self.curve_edgeon_emission_seds, evolved_name: self.curve_edgeon_emission_evolved_seds}

    # -----------------------------------------------------------------

    @property
    def curve_edgeon_emission_seds_plot_path(self):
        return fs.join(self.curves_path, "edgeon_emission_seds.pdf")

    # -----------------------------------------------------------------

    @property
    def has_curve_edgeon_emission_seds_plot(self):
        return fs.is_file(self.curve_edgeon_emission_seds_plot_path)

    # -----------------------------------------------------------------

    def plot_curve_edgeon_emission_seds(self):

        """
        This function ...
        :return:
        """

        # Plot
        plotting.plot_curves(self.curves_edgeon_emission_seds, path=self.curve_edgeon_emission_seds_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_absorption_data(self):
        return {unevolved_name: self.unevolved_spectral_absorption_fraction_curve, evolved_name: self.evolved_spectral_absorption_fraction_curve}

    # -----------------------------------------------------------------

    @property
    def curve_absorption_data_plot_path(self):
        return fs.join(self.curves_path, "absorption_cells.pdf")

    # -----------------------------------------------------------------

    def plot_curves_absorption_data(self):
        
        """
        This function ...
        :return: 
        """
        
        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust absorption from the cell data ...")

        # Plot
        plotting.plot_curves(self.curves_absorption_data, path=self.curve_absorption_data_plot_path, xlog=True,
                             y_label=self.curves_absorption_y_label, ylimits=self.curve_ylimits)

    # -----------------------------------------------------------------

    @lazyproperty
    def curves_emission_data(self):
        return {unevolved_name: self.unevolved_spectral_emission_fraction_curve, evolved_name: self.evolved_spectral_emission_fraction_curve}

    # -----------------------------------------------------------------

    @property
    def curve_emission_data_plot_path(self):
        return fs.join(self.curves_path, "emission_cells.pdf")

    # -----------------------------------------------------------------

    def plot_curves_emission_data(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the curve of the spectral heating by dust emission from the cell data ...")

        # Plot
        plotting.plot_curves(self.curves_emission_data, path=self.curve_emission_data_plot_path, xlog=True,
                             y_label=self.curves_emission_y_label, ylimits=self.curve_ylimits)

# -----------------------------------------------------------------

def get_fixed_cube_emission(cube, return_masks=False):

    """
    This function ...
    :param cube:
    :param return_masks:
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

    # EVERYTHING INVALID IS NOW NAN
    # NOW WE WANT TO INTERPOLATE THESE NANS, EXCEPT FOR THE BACKGROUND NANS (LARGEST)

    # Interpolate nans
    #fixed.interpolate_nans(sigma=3.)
    #fixed.interpolate_not_largest_nans(sigma=3., replace_nans=0.)
    masks_image = None
    if return_masks: masks_image = fixed.interpolate_nans_special(sigma=3., replace_nans=0., return_masks=True)
    else: fixed.interpolate_nans_special(sigma=3., replace_nans=0.)

    # Set flag
    fixed.metadata["fixed"] = True

    # Return the new cube
    if return_masks: return fixed, masks_image
    else: return fixed

# -----------------------------------------------------------------

def get_fixed_cube_absorption(cube, return_masks=False):

    """
    This function ...
    :param cube:
    :param return_masks:
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

    # EVERYTHING INVALID IS NOW NAN
    # NOW WE WANT TO INTERPOLATE THESE NANS, EXCEPT FOR THE BACKGROUND NANS (LARGEST)

    # Interpolate nans
    #fixed.interpolate_nans(sigma=3.)
    #fixed.interpolate_not_largest_nans(sigma=3., replace_nans=0.)
    masks_image = None
    if return_masks: masks_image = fixed.interpolate_nans_special(sigma=3., replace_nans=0., return_masks=True)
    else: fixed.interpolate_nans_special(sigma=3., replace_nans=0.)

    # Set flag
    fixed.metadata["fixed"] = True

    # Return the new cube
    if return_masks: return fixed, masks_image
    else: return fixed

# -----------------------------------------------------------------
