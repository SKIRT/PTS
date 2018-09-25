#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.sfr Contains the SFRAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent, AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ..core.data import Data3D, SpectralData3D
from ...core.units.parsing import parse_unit as u
from ...core.data.sed import SED
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

class AbsorptionAnalyser(AnalysisRunComponent):
    
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
        super(AbsorptionAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Writing
        self.write()

        # Show
        self.show()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AbsorptionAnalyser, self).setup()

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

    @property
    def length_unit(self):
        return "pc"

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

    @property
    def absorption_path(self):
        return self.analysis_run.absorption_path

    # -----------------------------------------------------------------
    # SIMULATION DATA
    #   TOTAL
    # -----------------------------------------------------------------

    @property
    def total_contribution_spectral_absorption_filepath(self):
        return self.total_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_total_contribution_spectral_absorption(self):
        return self.total_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------
    #   BULGE
    # -----------------------------------------------------------------

    @property
    def bulge_contribution_spectral_absorption_filepath(self):
        return self.bulge_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_bulge_contribution_spectral_absorption(self):
        return fs.is_file(self.bulge_contribution_spectral_absorption_filepath)

    # -----------------------------------------------------------------
    #   DISK
    # -----------------------------------------------------------------

    @property
    def disk_contribution_spectral_absorption_filepath(self):
        return self.disk_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_disk_contribution_spectral_absorption(self):
        return fs.is_file(self.disk_contribution_spectral_absorption_filepath)

    # -----------------------------------------------------------------
    #   OLD
    # -----------------------------------------------------------------

    @property
    def old_contribution_spectral_absorption_filepath(self):
        return self.old_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_old_contribution_spectral_absorption(self):
        return self.old_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------
    #   YOUNG
    # -----------------------------------------------------------------

    @property
    def young_contribution_spectral_absorption_filepath(self):
        return self.young_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_young_contribution_spectral_absorption(self):
        return self.young_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------
    #   SFR
    # -----------------------------------------------------------------

    @property
    def sfr_contribution_spectral_absorption_filepath(self):
        return self.ionizing_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_sfr_contribution_spectral_absorption(self):
        return self.ionizing_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------
    #   UNEVOLVED
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_spectral_absorption_filepath(self):
        return self.unevolved_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_unevolved_contribution_spectral_absorption(self):
        return self.unevolved_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------
    # 3D CELL ABSORPTION DATA
    #   TOTAL
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
        return fs.join(self.absorption_path, "total_absorption.dat")

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

        # With external xyz
        return SpectralData3D.from_table_file(self.total_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit,
                                              name=self.total_spectral_absorption_name,
                                              description=self.total_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    #   BULGE
    # -----------------------------------------------------------------

    @property
    def bulge_spectral_absorption_name(self):
        return "Labs_lambda_bulge"

    # -----------------------------------------------------------------

    @property
    def bulge_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the bulge model"

    # -----------------------------------------------------------------

    @property
    def bulge_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "bulge_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_bulge_spectral_absorption(self):
        return fs.is_file(self.bulge_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "bulge_spectral_absorption_path", True, write=True)
    def bulge_spectral_absorption_data(self):

        """
        This function ...
        :return:
        """

        # With external xyz
        return SpectralData3D.from_table_file(self.bulge_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit,
                                              name=self.bulge_spectral_absorption_name,
                                              description=self.bulge_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    #   DISK
    # -----------------------------------------------------------------

    @property
    def disk_spectral_absorption_name(self):
        return "Labs_lambda_disk"

    # -----------------------------------------------------------------

    @property
    def disk_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the disk model"

    # -----------------------------------------------------------------

    @property
    def disk_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "disk_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_disk_spectral_absorption(self):
        return fs.is_file(self.disk_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "disk_spectral_absorption_path", True, write=True)
    def disk_spectral_absorption_data(self):

        """
        This function ...
        :return:
        """

        # With external xyz
        return SpectralData3D.from_table_file(self.disk_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit,
                                              name=self.disk_spectral_absorption_name,
                                              description=self.disk_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    #   OLD
    # -----------------------------------------------------------------

    @property
    def old_spectral_absorption_name(self):
        return "Labs_lambda_old"

    # -----------------------------------------------------------------

    @property
    def old_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the old stellar populations"

    # -----------------------------------------------------------------

    @property
    def old_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "old_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_spectral_absorption(self):
        return fs.is_file(self.old_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "old_spectral_absorption_path", True, write=True)
    def old_spectral_absorption_data(self):

        """
        Thisn function ...
        :return:
        """

        # With external xyz
        return SpectralData3D.from_table_file(self.old_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit, name=self.old_spectral_absorption_name,
                                              description=self.old_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    #   YOUNG
    # -----------------------------------------------------------------

    @property
    def young_spectral_absorption_name(self):
        return "Labs_lambda_young"

    # -----------------------------------------------------------------

    @property
    def young_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the young stellar populations"

    # -----------------------------------------------------------------

    @property
    def young_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "young_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_spectral_absorption(self):
        return fs.is_file(self.young_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "young_spectral_absorption_path", True, write=True)
    def young_spectral_absorption_data(self):

        """
        Thisn function ...
        :return:
        """

        # With external xyz
        return SpectralData3D.from_table_file(self.young_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit, name=self.young_spectral_absorption_name,
                                              description=self.young_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    #   SFR
    # -----------------------------------------------------------------

    @property
    def sfr_spectral_absorption_name(self):
        return "Labs_lambda_sfr"

    # -----------------------------------------------------------------

    @property
    def sfr_spectral_absorption_description(self):
        return "Absorbed spectral luminosities in each dust cell for the star formation regions"

    # -----------------------------------------------------------------

    @property
    def sfr_spectral_absorption_path(self):
        return fs.join(self.absorption_path, "sfr_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_sfr_spectral_absorption(self):
        return fs.is_file(self.sfr_spectral_absorption_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SpectralData3D, "sfr_spectral_absorption_path", True, write=True)
    def sfr_spectral_absorption_data(self):

        """
        Thisn function ...
        :return:
        """

        # With external xyz
        return SpectralData3D.from_table_file(self.sfr_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit, name=self.sfr_spectral_absorption_name,
                                              description=self.sfr_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

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
        return fs.join(self.absorption_path, "unevolved_absorption.dat")

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

        # With external xyz
        return SpectralData3D.from_table_file(self.unevolved_contribution_spectral_absorption_filepath,
                                              self.cell_x_coordinates_colname, self.cell_y_coordinates_colname, self.cell_z_coordinates_colname,
                                              length_unit=self.length_unit,
                                              name=self.unevolved_spectral_absorption_name,
                                              description=self.unevolved_spectral_absorption_description,
                                              xyz_filepath=self.cell_coordinates_filepath)

    # -----------------------------------------------------------------
    # CURVES FROM SPECTRAL 3D ABSORPTION DATA
    #   TOTAL
    # -----------------------------------------------------------------

    @property
    def total_absorption_luminosity_name(self):
        return "Absorption luminosity (total)"

    # -----------------------------------------------------------------

    @property
    def total_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the total simulation"

    # -----------------------------------------------------------------

    @property
    def total_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "total_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_total_spectral_absorption_curve(self):
        return fs.is_file(self.total_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "total_spectral_absorption_curve_path", True, write=False)
    def total_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # total_absorption_luminosity_name
        return self.total_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # BULGE
    # -----------------------------------------------------------------

    @property
    def bulge_absorption_luminosity_name(self):
        return "Absorption luminosity (Bulge)"

    # -----------------------------------------------------------------

    @property
    def bulge_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the bulge simulation"

    # -----------------------------------------------------------------

    @property
    def bulge_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "bulge_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_bulge_spectral_absorption_curve(self):
        return fs.is_file(self.bulge_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "bulge_spectral_absorption_curve_path", True, write=False)
    def bulge_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # bulge_absorption_luminosity_name
        return self.bulge_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # DISK
    # -----------------------------------------------------------------

    @property
    def disk_absorption_luminosity_name(self):
        return "Absorption luminosity (Disk)"

    # -----------------------------------------------------------------

    @property
    def disk_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the disk simulation"

    # -----------------------------------------------------------------

    @property
    def disk_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "disk_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_disk_spectral_absorption_curve(self):
        return fs.is_file(self.disk_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "disk_spectral_absorption_curve_path", True, write=False)
    def disk_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # disk_absorption_luminosity_name
        return self.disk_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # OLD
    # -----------------------------------------------------------------

    @property
    def old_absorption_luminosity_name(self):
        return "Absorption luminosity (Old)"

    # -----------------------------------------------------------------

    @property
    def old_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the old simulation"

    # -----------------------------------------------------------------

    @property
    def old_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "old_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_old_spectral_absorption_curve(self):
        return fs.is_file(self.disk_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "old_spectral_absorption_curve_path", True, write=False)
    def old_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # old_absorption_luminosity_name
        return self.old_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # YOUNG
    # -----------------------------------------------------------------

    @property
    def young_absorption_luminosity_name(self):
        return "Absorption luminosity (Young)"

    # -----------------------------------------------------------------

    @property
    def young_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the young simulation"

    # -----------------------------------------------------------------

    @property
    def young_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "young_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_young_spectral_absorption_curve(self):
        return fs.is_file(self.young_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "young_spectral_absorption_curve_path", True, write=False)
    def young_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # young_absorption_luminosity_name
        return self.young_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # SFR
    # -----------------------------------------------------------------

    @property
    def sfr_absorption_luminosity_name(self):
        return "Absorption luminosity (SFR)"

    # -----------------------------------------------------------------

    @property
    def sfr_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the SFR simulation"

    # -----------------------------------------------------------------

    @property
    def sfr_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "sfr_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_sfr_spectral_absorption_curve(self):
        return fs.is_file(self.sfr_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "sfr_spectral_absorption_curve_path", True, write=False)
    def sfr_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # sfr_absorption_luminosity_name
        return self.sfr_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    #   UNEVOLVED
    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_luminosity_name(self):
        return "Absorption luminosity (unevolved)"

    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_luminosity_description(self):
        return "Absorption luminosity in dust cells for the unevolved simulation"

    # -----------------------------------------------------------------

    @property
    def unevolved_spectral_absorption_curve_path(self):
        return fs.join(self.absorption_path, "unevolved_curve_absorption.dat")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_spectral_absorption_curve(self):
        return fs.is_file(self.unevolved_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    @lazyfileproperty(SED, "unevolved_spectral_absorption_curve_path", True, write=True)
    def unevolved_spectral_absorption_curve(self):

        """
        This function ...
        :return:
        """

        # unevolved_absorption_luminosity_name
        return self.unevolved_spectral_absorption_data.get_global_sed()

    # -----------------------------------------------------------------
    # TOTAL STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def total_observed_sed(self):
        return self.model.total_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_luminosity(self):
        return self.total_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def total_stellar_sed(self):
        return self.model.total_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_stellar_luminosity(self):
        return self.total_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # BULGE STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def bulge_observed_sed(self):
        return self.model.bulge_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_observed_luminosity(self):
        return self.bulge_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def bulge_stellar_sed(self):
        return self.model.bulge_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_stellar_luminosity(self):
        return self.bulge_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # DISK STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def disk_observed_sed(self):
        return self.model.disk_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_observed_luminosity(self):
        return self.disk_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def disk_stellar_sed(self):
        return self.model.disk_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_stellar_luminosity(self):
        return self.disk_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # OLD STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def old_observed_sed(self):
        return self.model.old_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_observed_luminosity(self):
        return self.old_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def old_stellar_sed(self):
        return self.model.old_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stellar_luminosity(self):
        return self.old_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # YOUNG STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def young_observed_sed(self):
        return self.model.young_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_observed_luminosity(self):
        return self.young_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def young_stellar_sed(self):
        return self.model.young_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stellar_luminosity(self):
        return self.young_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # SFR STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def sfr_observed_sed(self):
        return self.model.sfr_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_observed_luminosity(self):
        return self.sfr_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def sfr_stellar_sed(self):
        return self.model.sfr_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_stellar_luminosity(self):
        return self.sfr_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # UNEVOLVED STELLAR EMISSION
    # -----------------------------------------------------------------

    @property
    def unevolved_observed_sed(self):
        return self.model.unevolved_simulations.observed_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_observed_luminosity(self):
        return self.unevolved_observed_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def unevolved_stellar_sed(self):
        return self.model.unevolved_simulations.intrinsic_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_stellar_luminosity(self):
        return self.unevolved_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # TOTAL ABSORPTION
    #   Diffuse
    # -----------------------------------------------------------------

    @property
    def total_absorption_sed_diffuse(self):
        return self.model.total_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_luminosity_diffuse(self):
        return self.total_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_fraction_diffuse(self):
        return self.total_absorption_luminosity_diffuse.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def total_dust_sed_diffuse(self):
        return self.model.total_simulations.observed_diffuse_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_diffuse(self):
        return self.total_dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_fraction_diffuse(self):
        return self.total_dust_luminosity_diffuse.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_sed_all(self):
        return self.total_absorption_sed_diffuse + self.sfr_absorption_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_luminosity_all(self):
        return self.total_absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_fraction_all(self):
        return self.total_absorption_luminosity_all.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_sed_all(self):
        return self.total_dust_sed_diffuse + self.sfr_dust_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_all(self):
        return self.total_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_fraction_all(self):
        return self.total_dust_luminosity_all.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def total_dust_sed_all_alt(self):
        return self.model.total_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_luminosity_all_alt(self):
        return self.total_dust_sed_all_alt.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_dust_fraction_all_alt(self):
        return self.total_dust_luminosity_all_alt.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------
    # BULGE ABSORPTION
    # -----------------------------------------------------------------

    @property
    def bulge_absorption_sed(self):
        return self.model.bulge_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_absorption_luminosity(self):
        return self.bulge_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_absorption_fraction(self):
        return self.bulge_absorption_luminosity.value / self.bulge_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def bulge_dust_sed(self):
        return self.model.bulge_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_dust_luminosity(self):
        return self.bulge_dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_dust_fraction(self):
        return self.bulge_dust_luminosity.value / self.bulge_stellar_luminosity.value

    # -----------------------------------------------------------------
    # DISK ABSORPTION
    # -----------------------------------------------------------------

    @property
    def disk_absorption_sed(self):
        return self.model.disk_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_absorption_luminosity(self):
        return self.disk_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_absorption_fraction(self):
        return self.disk_absorption_luminosity.value / self.disk_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def disk_dust_sed(self):
        return self.model.disk_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_dust_luminosity(self):
        return self.disk_dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_dust_fraction(self):
        return self.disk_dust_luminosity.value / self.disk_stellar_luminosity.value

    # -----------------------------------------------------------------
    # OLD ABSORPTION
    # -----------------------------------------------------------------

    @property
    def old_absorption_sed(self):
        return self.model.old_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_absorption_luminosity(self):
        return self.old_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_absorption_fraction(self):
        return self.old_absorption_luminosity.value / self.old_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def old_dust_sed(self):
        return self.model.old_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_dust_luminosity(self):
        return self.old_dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_dust_fraction(self):
        return self.old_dust_luminosity.value / self.old_stellar_luminosity.value

    # -----------------------------------------------------------------
    # YOUNG ABSORPTION
    # -----------------------------------------------------------------

    @property
    def young_absorption_sed(self):
        return self.model.young_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_absorption_luminosity(self):
        return self.young_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_absorption_fraction(self):
        return self.young_absorption_luminosity.value / self.young_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def young_dust_sed(self):
        return self.model.young_simulations.observed_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_dust_luminosity(self):
        return self.young_dust_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_dust_fraction(self):
        return self.young_dust_luminosity.value / self.young_stellar_luminosity.value

    # -----------------------------------------------------------------
    # SFR ABSORPTION
    #   Diffuse
    # -----------------------------------------------------------------

    @property
    def sfr_absorption_sed_diffuse(self):
        return self.model.sfr_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_luminosity_diffuse(self):
        return self.sfr_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_fraction_diffuse(self):
        return self.sfr_absorption_luminosity_diffuse.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def sfr_dust_sed_diffuse(self):
        return self.model.sfr_simulations.observed_diffuse_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_diffuse(self):
        return self.sfr_dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_fraction_diffuse(self):
        return self.sfr_dust_luminosity_diffuse.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------
    #   Internal
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_sed_internal(self):

        # Get stellar SEDs
        intrinsic_stellar = self.model.get_stellar_sed("sfr", "intrinsic")
        transparent_stellar = self.model.intrinsic_sfr_stellar_sed # NEW FROM TRANSPARENT MAPPINGS

        # INTERNALLY ABSORBED
        return transparent_stellar - intrinsic_stellar

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_luminosity_internal(self):
        return self.sfr_absorption_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_fraction_internal(self):
        return self.sfr_absorption_luminosity_internal.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_sed_internal(self):
        return self.model.intrinsic_sfr_dust_sed # NEW FROM TRANSPARENT MAPPINGS, SUBTRACTED

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_internal(self):
        return self.sfr_dust_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_fraction_internal(self):
        return self.sfr_dust_luminosity_internal.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_sed_internal_alt(self):
        return self.model.sfr_simulations.intrinsic_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_internal_alt(self):
        return self.sfr_dust_sed_internal_alt.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_fraction_internal_alt(self):
        return self.sfr_dust_luminosity_internal_alt.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_sed_all(self):
        return self.sfr_absorption_sed_diffuse + self.sfr_absorption_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_luminosity_all(self):
        return self.sfr_absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_fraction_all(self):
        return self.sfr_absorption_luminosity_all.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_sed_all(self):
        return self.sfr_dust_sed_diffuse + self.sfr_dust_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_luminosity_all(self):
        return self.sfr_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_dust_fraction_all(self):
        return self.sfr_dust_luminosity_all.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------
    # UNEVOLVED ABSORPTION
    #   Diffuse
    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_sed_diffuse(self):
        return self.model.unevolved_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_luminosity_diffuse(self):
        return self.unevolved_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_fraction_diffuse(self):
        return self.unevolved_absorption_luminosity_diffuse.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_sed_diffuse(self):
        return self.model.unevolved_simulations.observed_diffuse_dust_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_luminosity_diffuse(self):
        return self.unevolved_dust_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_fraction_diffuse(self):
        return self.unevolved_dust_luminosity_diffuse.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------
    #   All
    # -----------------------------------------------------------------

    @property
    def unevolved_absorption_sed_all(self):
        return self.unevolved_absorption_sed_diffuse + self.sfr_absorption_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_luminosity_all(self):
        return self.unevolved_absorption_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_fraction_all(self):
        return self.unevolved_absorption_luminosity_all.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_sed_all(self):
        return self.unevolved_dust_sed_diffuse + self.sfr_dust_sed_internal

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_luminosity_all(self):
        return self.unevolved_dust_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_dust_fraction_all(self):
        return self.unevolved_dust_luminosity_all.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write curves
        self.write_curves()

    # -----------------------------------------------------------------

    @property
    def do_write_total_curve(self):
        return not self.has_total_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_bulge_curve(self):
        return not self.has_bulge_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_disk_curve(self):
        return not self.has_disk_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_old_curve(self):
        return not self.has_old_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_young_curve(self):
        return not self.has_young_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_sfr_curve(self):
        return not self.has_sfr_spectral_absorption_curve

    # -----------------------------------------------------------------

    @property
    def do_write_unevolved_curve(self):
        return not self.has_unevolved_spectral_absorption_curve

    # -----------------------------------------------------------------

    def write_curves(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the absorption curves ...")

        # Total
        if self.do_write_total_curve: self.write_total_curve()

        # Bulge
        if self.do_write_bulge_curve: self.write_bulge_curve()

        # Disk
        if self.do_write_disk_curve: self.write_disk_curve()

        # Old
        if self.do_write_old_curve: self.write_old_curve()

        # Young
        if self.do_write_young_curve: self.write_young_curve()

        # SFR
        if self.do_write_sfr_curve: self.write_sfr_curve()

        # Unevolved
        if self.do_write_unevolved_curve: self.write_unevolved_curve()

    # -----------------------------------------------------------------

    def write_total_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.total_spectral_absorption_curve.saveto(self.total_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_bulge_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.bulge_spectral_absorption_curve.saveto(self.bulge_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_disk_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.disk_spectral_absorption_curve.saveto(self.disk_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_old_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.old_spectral_absorption_curve.saveto(self.old_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_young_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.young_spectral_absorption_curve.saveto(self.young_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_sfr_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.sfr_spectral_absorption_curve.saveto(self.sfr_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def write_unevolved_curve(self):

        """
        This function ...
        :return:
        """

        # Write
        self.unevolved_spectral_absorption_curve.saveto(self.unevolved_spectral_absorption_curve_path)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Total
        self.show_total()

        # Bulge
        self.show_bulge()

        # Disk
        self.show_disk()

        # Old
        self.show_old()

        # Young
        self.show_young()

        # SFR
        self.show_sfr()

        # Unevolved
        self.show_unevolved()

    # -----------------------------------------------------------------

    def show_total(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.total_observed_luminosity, stellar_bol=self.total_stellar_luminosity)

        # Diffuse
        self.show_total_diffuse()

        # All
        self.show_total_all()

    # -----------------------------------------------------------------

    def show_total_diffuse(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.total_absorption_luminosity_diffuse, dust=self.total_dust_luminosity_diffuse)
        show_values(rel_absorbed=self.total_absorption_fraction_diffuse, rel_dust=self.total_dust_fraction_diffuse)

    # -----------------------------------------------------------------

    def show_total_all(self):

        """
        This function ...
        :return:
        """

        # Show
        show_values(absorbed=self.total_absorption_luminosity_all, dust=self.total_dust_luminosity_all, dust_alt=self.total_dust_luminosity_all_alt)
        show_values(rel_absorbed=self.total_absorption_fraction_all, rel_dust=self.total_dust_fraction_all, rel_dust_alt=self.total_dust_fraction_all_alt)

    # -----------------------------------------------------------------

    def show_bulge(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.bulge_observed_luminosity, stellar_bol=self.bulge_stellar_luminosity)

        # Show absorption
        show_values(absorbed=self.bulge_absorption_luminosity, dust=self.bulge_dust_luminosity)
        show_values(rel_absorbed=self.bulge_absorption_fraction, rel_dust=self.bulge_dust_fraction)

    # -----------------------------------------------------------------

    def show_disk(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.disk_observed_luminosity, stellar_bol=self.disk_stellar_luminosity)

        # Show absorption
        show_values(absorbed=self.disk_absorption_luminosity, dust=self.disk_dust_luminosity)
        show_values(rel_absorbed=self.disk_absorption_fraction, rel_dust=self.disk_dust_fraction)

    # -----------------------------------------------------------------

    def show_old(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.old_observed_luminosity, stellar_bol=self.old_stellar_luminosity)

        # Show absorption
        show_values(absorbed=self.old_absorption_luminosity, dust=self.old_dust_luminosity)
        show_values(rel_absorbed=self.old_absorption_fraction, rel_dust=self.old_dust_fraction)

    # -----------------------------------------------------------------

    def show_young(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.young_observed_luminosity, stellar_bol=self.young_stellar_luminosity)

        # Show absorption
        show_values(absorbed=self.young_absorption_luminosity, dust=self.young_dust_luminosity)
        show_values(rel_absorbed=self.young_absorption_fraction, rel_dust=self.young_dust_fraction)

    # -----------------------------------------------------------------

    def show_sfr(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.sfr_observed_luminosity, stellar_bol=self.sfr_stellar_luminosity)

        # Diffuse
        self.show_sfr_diffuse()

        # Internal
        self.show_sfr_internal()

        # All
        self.show_sfr_all()

    # -----------------------------------------------------------------

    def show_sfr_diffuse(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.sfr_absorption_luminosity_diffuse, dust=self.sfr_dust_luminosity_diffuse)
        show_values(rel_absorbed=self.sfr_absorption_fraction_diffuse, rel_dust=self.sfr_dust_fraction_diffuse)

    # -----------------------------------------------------------------

    def show_sfr_internal(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.sfr_absorption_luminosity_internal, dust=self.sfr_dust_luminosity_internal, dust_alt=self.sfr_dust_luminosity_internal_alt)
        show_values(rel_absorbed=self.sfr_absorption_fraction_internal, rel_dust=self.sfr_dust_fraction_internal, rel_dust_alt=self.sfr_dust_fraction_internal_alt)

    # -----------------------------------------------------------------

    def show_sfr_all(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.sfr_absorption_luminosity_all, dust=self.sfr_dust_luminosity_all)
        show_values(rel_absorbed=self.sfr_absorption_fraction_all, rel_dust=self.sfr_dust_fraction_all)

    # -----------------------------------------------------------------

    def show_unevolved(self):

        """
        This function ...
        :return:
        """

        # Show stellar
        show_values(observed_bol=self.unevolved_observed_luminosity, stellar_bol=self.unevolved_stellar_luminosity)

        # Diffuse
        self.show_unevolved_diffuse()

        # All
        self.show_unevolved_all()

    # -----------------------------------------------------------------

    def show_unevolved_diffuse(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.unevolved_absorption_luminosity_diffuse, dust=self.unevolved_dust_luminosity_diffuse)
        show_values(rel_absorbed=self.unevolved_absorption_fraction_diffuse, rel_dust=self.unevolved_dust_fraction_diffuse)

    # -----------------------------------------------------------------

    def show_unevolved_all(self):

        """
        This function ...
        :return:
        """

        # Show absorption
        show_values(absorbed=self.unevolved_absorption_luminosity_all, dust=self.unevolved_dust_luminosity_all)
        show_values(rel_absorbed=self.unevolved_absorption_fraction_all, rel_dust=self.unevolved_dust_fraction_all)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------

def show_values(**kwargs):
    for label in kwargs: print(" - " + fmt.bold + label + fmt.reset_bold + ": " + tostr(kwargs[label]))

# -----------------------------------------------------------------
