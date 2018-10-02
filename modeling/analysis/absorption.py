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

# Import standard modules
from box import Box
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import AnalysisComponent, AnalysisRunComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, lazyfileproperty
from ..core.data import Data3D, SpectralData3D
from ...core.data.sed import SED
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.basics.configuration import open_mapping, save_mapping
from ...core.basics.configuration import open_box, save_box
from ...core.basics.map import Map
from ...core.tools import types
from ...core.plot.sed import plot_seds
from ...core.units.parsing import parse_quantity as q
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

observed_stellar_name = "observed (stellar)"
intrinsic_stellar_name = "intrinsic (stellar)"
absorbed_name = "absorbed"
dust_name = "dust"

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
        if self.config.show: self.show()

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

    @property
    def total_observed_stellar_sed_diffuse(self):
        return self.model.total_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_stellar_luminosity_diffuse(self):
        return self.total_observed_stellar_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def total_observed_stellar_sed_all(self):
        return self.old_observed_stellar_sed + self.unevolved_observed_stellar_sed_all

    # -----------------------------------------------------------------

    @lazyproperty
    def total_observed_stellar_luminosity_all(self):
        return self.total_observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def bulge_observed_stellar_sed(self):
        return self.model.bulge_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_observed_stellar_luminosity(self):
        return self.bulge_observed_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def disk_observed_stellar_sed(self):
        return self.model.disk_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_observed_stellar_luminosity(self):
        return self.disk_observed_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def old_observed_stellar_sed(self):
        return self.model.old_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_observed_stellar_luminosity(self):
        return self.old_observed_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def young_observed_stellar_sed(self):
        return self.model.young_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_observed_stellar_luminosity(self):
        return self.young_observed_stellar_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def sfr_observed_stellar_sed_diffuse(self):
        return self.model.sfr_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_observed_stellar_luminosity_diffuse(self):
        return self.sfr_observed_stellar_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def sfr_observed_stellar_sed_internal(self):
        return self.model.intrinsic_sfr_sed - self.sfr_dust_sed_internal # intrinsic here means unaffected by diffuse dust, but still with internal extinction

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_observed_stellar_luminosity_internal(self):
        return self.sfr_observed_stellar_sed_internal.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def sfr_observed_stellar_sed_all(self):
        return self.model.sfr_simulations.observed_sed - self.sfr_dust_sed_all

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_observed_stellar_luminosity_all(self):
        return self.sfr_observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

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

    @property
    def unevolved_observed_stellar_sed_diffuse(self):
        return self.model.unevolved_simulations.observed_stellar_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_observed_stellar_luminosity_diffuse(self):
        return self.unevolved_observed_stellar_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def unevolved_observed_stellar_sed_all(self):
        return self.young_observed_stellar_sed + self.sfr_observed_stellar_sed_all

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_observed_stellar_luminosity_all(self):
        return self.unevolved_observed_stellar_sed_all.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------
    # TOTAL ABSORPTION
    #   Diffuse
    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_luminosity_diffuse_cells(self):
        return self.total_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_fraction_diffuse_cells(self):
        return self.total_absorption_luminosity_diffuse_cells.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def total_absorption_sed_diffuse(self):
        return self.model.total_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_sed_diffuse(self):
        return self.model.total_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_luminosity_diffuse(self):
        return self.total_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_luminosity_diffuse(self):
        return self.has_total_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def total_absorption_fraction_diffuse(self):
        return self.total_absorption_luminosity_diffuse.value / self.total_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_total_absorption_fraction_diffuse(self):
        return self.has_total_absorption_luminosity_diffuse

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

    @lazyproperty
    def bulge_absorption_luminosity_cells(self):
        return self.bulge_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_absorption_fraction_cells(self):
        return self.bulge_absorption_luminosity_cells.value / self.bulge_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def bulge_absorption_sed(self):
        return self.model.bulge_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_bulge_absorption_sed(self):
        return self.model.bulge_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_absorption_luminosity(self):
        return self.bulge_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_bulge_absorption_luminosity(self):
        return self.has_bulge_absorption_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_absorption_fraction(self):
        return self.bulge_absorption_luminosity.value / self.bulge_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_bulge_absorption_fraction(self):
        return self.has_bulge_absorption_luminosity

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

    @lazyproperty
    def disk_absorption_luminosity_cells(self):
        return self.disk_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_absorption_fraction_cells(self):
        return self.disk_absorption_luminosity_cells.value / self.disk_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def disk_absorption_sed(self):
        return self.model.disk_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_disk_absorption_sed(self):
        return self.model.disk_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_absorption_luminosity(self):
        return self.disk_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_disk_absorption_luminosity(self):
        return self.has_disk_absorption_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_absorption_fraction(self):
        return self.disk_absorption_luminosity.value / self.disk_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_disk_absorption_fraction(self):
        return self.has_disk_absorption_luminosity

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

    @lazyproperty
    def old_absorption_luminosity_cells(self):
        return self.old_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_absorption_fraction_cells(self):
        return self.old_absorption_luminosity_cells.value / self.old_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def old_absorption_sed(self):
        return self.model.old_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_old_absorption_sed(self):
        return self.model.old_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_absorption_luminosity(self):
        return self.old_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_old_absorption_luminosity(self):
        return self.has_old_absorption_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def old_absorption_fraction(self):
        return self.old_absorption_luminosity.value / self.old_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_old_absorption_fraction(self):
        return self.has_old_absorption_luminosity

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

    @lazyproperty
    def young_absorption_luminosity_cells(self):
        return self.young_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_absorption_fraction_cells(self):
        return self.young_absorption_luminosity_cells.value / self.young_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def young_absorption_sed(self):
        return self.model.young_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_young_absorption_sed(self):
        return self.model.young_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_absorption_luminosity(self):
        return self.young_absorption_sed.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_young_absorption_luminosity(self):
        return self.has_young_absorption_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def young_absorption_fraction(self):
        return self.young_absorption_luminosity.value / self.young_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_young_absorption_fraction(self):
        return self.has_young_absorption_luminosity

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

    @lazyproperty
    def sfr_absorption_luminosity_diffuse_cells(self):
        return self.sfr_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_fraction_diffuse_cells(self):
        return self.sfr_absorption_luminosity_diffuse_cells.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def sfr_absorption_sed_diffuse(self):
        return self.model.sfr_simulations.observed_sed_absorbed

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorption_sed_diffuse(self):
        return self.model.sfr_simulations.has_observed_sed_absorbed

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_luminosity_diffuse(self):
        return self.sfr_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorption_luminosity_diffuse(self):
        return self.has_sfr_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def sfr_absorption_fraction_diffuse(self):
        return self.sfr_absorption_luminosity_diffuse.value / self.sfr_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_sfr_absorption_fraction_diffuse(self):
        return self.has_sfr_absorption_luminosity_diffuse

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
        #intrinsic_stellar = self.model.get_stellar_sed("sfr", "intrinsic")
        intrinsic_stellar = self.model.intrinsic_sfr_sed # same
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

    @lazyproperty
    def unevolved_absorption_luminosity_diffuse_cells(self):
        return self.unevolved_spectral_absorption_curve.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_fraction_diffuse_cells(self):
        return self.unevolved_absorption_luminosity_diffuse_cells.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_sed_diffuse(self):
        #return self.model.unevolved_simulations.observed_sed_absorbed
        return self.young_absorption_sed + self.sfr_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorption_sed_diffuse(self):
        #return self.model.unevolved_simulations.has_observed_sed_absorbed
        return self.has_young_absorption_sed and self.has_sfr_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_luminosity_diffuse(self):
        return self.unevolved_absorption_sed_diffuse.integrate().to(self.bolometric_luminosity_unit, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorption_luminosity_diffuse(self):
        return self.has_unevolved_absorption_sed_diffuse

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_absorption_fraction_diffuse(self):
        return self.unevolved_absorption_luminosity_diffuse.value / self.unevolved_stellar_luminosity.value

    # -----------------------------------------------------------------

    @property
    def has_unevolved_absorption_fraction_diffuse(self):
        return self.has_unevolved_absorption_luminosity_diffuse

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

    @property
    def total_properties_path(self):
        return fs.join(self.absorption_path, "total.txt")

    # -----------------------------------------------------------------

    @property
    def has_total_properties(self):
        return fs.is_file(self.total_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "total_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "total_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def total_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.total_observed_luminosity
        props.stellar_bol = self.total_stellar_luminosity
        props.observed_stellar_bol_diffuse = self.total_observed_stellar_luminosity_diffuse
        props.observed_stellar_bol_all = self.total_observed_stellar_luminosity_all

        # Diffuse
        diffuse = Box(ordered_box=True)
        if self.has_total_absorption_luminosity_diffuse: diffuse.absorbed = self.total_absorption_luminosity_diffuse
        diffuse.absorbed_cells = self.total_absorption_luminosity_diffuse_cells # CELLS
        diffuse.dust = self.total_dust_luminosity_diffuse
        if self.has_total_absorption_fraction_diffuse: diffuse.rel_absorbed = self.total_absorption_fraction_diffuse
        diffuse.rel_absorbed_cells = self.total_absorption_fraction_diffuse_cells # CELLS
        diffuse.rel_dust = self.total_dust_fraction_diffuse

        # All
        all = Box(ordered_box=True)
        all.absorbed = self.total_absorption_luminosity_all
        all.dust = self.total_dust_luminosity_all
        all.dust_alt = self.total_dust_luminosity_all_alt
        all.rel_absorbed = self.total_absorption_fraction_all
        all.rel_dust = self.total_dust_fraction_all
        all.rel_dust_alt = self.total_dust_fraction_all_alt

        # Return
        props.diffuse = diffuse
        props.all = all
        return props

    # -----------------------------------------------------------------

    @property
    def bulge_properties_path(self):
        return fs.join(self.absorption_path, "bulge.txt")

    # -----------------------------------------------------------------

    @property
    def has_bulge_properties(self):
        return fs.is_file(self.bulge_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "bulge_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "bulge_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def bulge_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.bulge_observed_luminosity
        props.stellar_bol = self.bulge_stellar_luminosity
        props.observed_stellar_bol = self.bulge_observed_stellar_luminosity

        # Absorption
        if self.has_bulge_absorption_luminosity: props.absorbed = self.bulge_absorption_luminosity
        props.absorbed_cells = self.bulge_absorption_luminosity_cells
        props.dust = self.bulge_dust_luminosity
        if self.has_bulge_absorption_fraction: props.rel_absorbed = self.bulge_absorption_fraction
        props.rel_absorbed_cells = self.bulge_absorption_fraction_cells
        props.rel_dust = self.bulge_dust_fraction

        # Return
        return props

    # -----------------------------------------------------------------

    @property
    def disk_properties_path(self):
        return fs.join(self.absorption_path, "disk.txt")

    # -----------------------------------------------------------------

    @property
    def has_disk_properties(self):
        return fs.is_file(self.disk_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "disk_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "disk_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def disk_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.disk_observed_luminosity
        props.stellar_bol = self.disk_stellar_luminosity
        props.observed_stellar_bol = self.disk_observed_stellar_luminosity

        # Absorption
        if self.has_disk_absorption_luminosity: props.absorbed = self.disk_absorption_luminosity
        props.absorbed_cells = self.disk_absorption_luminosity_cells
        props.dust = self.disk_dust_luminosity
        if self.has_disk_absorption_fraction: props.rel_absorbed = self.disk_absorption_fraction
        props.rel_absorbed_cells = self.disk_absorption_fraction_cells
        props.rel_dust = self.disk_dust_fraction

        # Return
        return props

    # -----------------------------------------------------------------

    @property
    def old_properties_path(self):
        return fs.join(self.absorption_path, "old.txt")

    # -----------------------------------------------------------------

    @property
    def has_old_properties(self):
        return fs.is_file(self.old_properties_path)

    #-----------------------------------------------------------------

    #@lazyfileproperty(Map, "old_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "old_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def old_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.old_observed_luminosity
        props.stellar_bol = self.old_stellar_luminosity
        props.observed_stellar_bol = self.old_observed_stellar_luminosity

        # Absorption
        if self.has_old_absorption_luminosity: props.absorbed = self.old_absorption_luminosity
        props.absorbed_cells = self.old_absorption_luminosity_cells
        props.dust = self.old_dust_luminosity
        if self.has_old_absorption_fraction: props.rel_absorbed = self.old_absorption_fraction
        props.rel_absorbed_cells = self.old_absorption_fraction_cells
        props.rel_dust = self.old_dust_fraction

        # Return
        return props

    # -----------------------------------------------------------------

    @property
    def young_properties_path(self):
        return fs.join(self.absorption_path, "young.txt")

    # -----------------------------------------------------------------

    @property
    def has_young_properties(self):
        return fs.is_file(self.young_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "young_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "young_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def young_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.young_observed_luminosity
        props.stellar_bol = self.young_stellar_luminosity
        props.observed_stellar_bol = self.young_observed_stellar_luminosity

        # Absorption
        if self.has_young_absorption_luminosity: props.absorbed = self.young_absorption_luminosity
        props.absorbed_cells = self.young_absorption_luminosity_cells
        props.dust = self.young_dust_luminosity
        if self.has_young_absorption_fraction: props.rel_absorbed = self.young_absorption_fraction
        props.rel_absorbed_cells = self.young_absorption_fraction_cells
        props.rel_dust = self.young_dust_fraction

        # Return
        return props

    # -----------------------------------------------------------------

    @property
    def sfr_properties_path(self):
        return fs.join(self.absorption_path, "sfr.txt")

    # -----------------------------------------------------------------

    @property
    def has_sfr_properties(self):
        return fs.is_file(self.sfr_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "sfr_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "sfr_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def sfr_properties(self):

        """
        Thisf unction ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.sfr_observed_luminosity
        props.stellar_bol = self.sfr_stellar_luminosity
        props.observed_stellar_bol_diffuse = self.sfr_observed_stellar_luminosity_diffuse
        props.observed_stellar_bol_internal = self.sfr_observed_stellar_luminosity_internal
        props.observed_stellar_bol_all = self.sfr_observed_stellar_luminosity_all

        # Diffuse
        diffuse = Box(ordered_box=True)
        if self.has_sfr_absorption_luminosity_diffuse: diffuse.absorbed = self.sfr_absorption_luminosity_diffuse
        diffuse.absorbed_cells = self.sfr_absorption_luminosity_diffuse_cells
        diffuse.dust = self.sfr_dust_luminosity_diffuse
        if self.has_sfr_absorption_fraction_diffuse: diffuse.rel_absorbed = self.sfr_absorption_fraction_diffuse
        diffuse.rel_absorbed_cells = self.sfr_absorption_fraction_diffuse_cells
        diffuse.rel_dust = self.sfr_dust_fraction_diffuse

        # Internal
        internal = Box(ordered_box=True)
        internal.absorbed = self.sfr_absorption_luminosity_internal
        internal.dust = self.sfr_dust_luminosity_internal
        internal.dust_alt = self.sfr_dust_luminosity_internal_alt
        internal.rel_absorbed = self.sfr_absorption_fraction_internal
        internal.rel_dust = self.sfr_dust_fraction_internal
        internal.rel_dust_alt = self.sfr_dust_fraction_internal_alt

        # All
        all = Box(ordered_box=True)
        all.absorbed = self.sfr_absorption_luminosity_all
        all.dust = self.sfr_dust_luminosity_all
        all.rel_absorbed = self.sfr_absorption_fraction_all
        all.rel_dust = self.sfr_dust_fraction_all

        # Return
        props.diffuse = diffuse
        props.internal = internal
        props.all = all
        return props

    # -----------------------------------------------------------------

    @property
    def unevolved_properties_path(self):
        return fs.join(self.absorption_path, "unevolved.txt")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_properties(self):
        return fs.is_file(self.unevolved_properties_path)

    # -----------------------------------------------------------------

    #@lazyfileproperty(Map, "unevolved_properties_path", True, write=True, fsave=save_mapping, fload=open_mapping)
    @lazyfileproperty(Box, "unevolved_properties_path", True, write=True, fsave=save_box, fload=open_box)
    def unevolved_properties(self):

        """
        This function ...
        :return:
        """

        # Create
        props = Box(ordered_box=True)

        # Stellar
        props.observed_bol = self.unevolved_observed_luminosity
        props.stellar_bol = self.unevolved_stellar_luminosity
        props.observed_stellar_bol_diffuse = self.unevolved_observed_stellar_luminosity_diffuse
        props.observed_stellar_bol_all = self.unevolved_observed_stellar_luminosity_all

        # Diffuse
        diffuse = Box(ordered_box=True)
        if self.has_unevolved_absorption_luminosity_diffuse: diffuse.absorbed = self.unevolved_absorption_luminosity_diffuse
        diffuse.absorbed_cells = self.unevolved_absorption_luminosity_diffuse_cells
        diffuse.dust = self.unevolved_dust_luminosity_diffuse
        if self.has_unevolved_absorption_fraction_diffuse: diffuse.rel_absorbed = self.unevolved_absorption_fraction_diffuse
        diffuse.rel_absorbed_cells = self.unevolved_absorption_fraction_diffuse_cells
        diffuse.rel_dust = self.unevolved_dust_fraction_diffuse

        # All
        all = Box(ordered_box=True)
        all.absorbed = self.unevolved_absorption_luminosity_all
        all.dust = self.unevolved_dust_luminosity_all
        all.rel_absorbed = self.unevolved_absorption_fraction_all
        all.rel_dust = self.unevolved_dust_fraction_all

        # Return
        props.diffuse = diffuse
        props.all = all
        return props

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

        # Show
        print("TOTAL")
        print("")
        show_properties(self.total_properties)
        print("")

    # -----------------------------------------------------------------

    def show_bulge(self):

        """
        This function ...
        :return:
        """

        # Show
        print("BULGE")
        print("")
        show_properties(self.bulge_properties)
        print("")

    # -----------------------------------------------------------------

    def show_disk(self):

        """
        This function ...
        :return:
        """

        # Show
        print("DISK")
        print("")
        show_properties(self.disk_properties)
        print("")

    # -----------------------------------------------------------------

    def show_old(self):

        """
        This function ...
        :return:
        """

        # Show
        print("OLD")
        print("")
        show_properties(self.old_properties)
        print("")

    # -----------------------------------------------------------------

    def show_young(self):

        """
        This function ...
        :return:
        """

        # Show
        print("YOUNG")
        print("")
        show_properties(self.young_properties)
        print("")

    # -----------------------------------------------------------------

    def show_sfr(self):

        """
        This function ...
        :return:
        """

        # Show
        print("SFR")
        print("")
        show_properties(self.sfr_properties)
        print("")

    # -----------------------------------------------------------------

    def show_unevolved(self):

        """
        This function ...
        :return:
        """

        # Show
        print("UNEVOLVED")
        print("")
        show_properties(self.unevolved_properties)
        print("")

    # -----------------------------------------------------------------

    @property
    def do_plot_total(self):
        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_bulge(self):
        return not self.has_bulge_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_disk(self):
        return not self.has_disk_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_old(self):
        return not self.has_old_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_young(self):
        return not self.has_young_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_sfr(self):
        return True

    # -----------------------------------------------------------------

    @property
    def do_plot_unevolved(self):
        return True

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Total
        if self.do_plot_total: self.plot_total()

        # Bulge
        if self.do_plot_bulge: self.plot_bulge()

        # Disk
        if self.do_plot_disk: self.plot_disk()

        # Old
        if self.do_plot_old: self.plot_old()

        # Young
        if self.do_plot_young: self.plot_young()

        # SFR
        if self.do_plot_sfr: self.plot_sfr()

        # Unevolved
        if self.do_plot_unevolved: self.plot_unevolved()

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_unit(self):
        return u("Jy")

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_options(self):
        options = dict()
        options["absorbed"] = {"above": observed_stellar_name, "above_name": intrinsic_stellar_name}
        options["dust"] = {"above": observed_stellar_name}
        return options

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):
        return q("0.01 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):
        return q("1000 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def min_flux(self):
        return q("1e-13.5 W/m2", density=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_flux(self):
        return q("1e-10 W/m2", density=True)

    # -----------------------------------------------------------------

    def plot_seds(self, observed_stellar, absorbed, dust, path, **extra):

        """
        This function ...
        :param observed_stellar:
        :param absorbed:
        :param dust:
        :param path:
        :param extra:
        :return:
        """

        # Set SEDs
        seds = OrderedDict()
        seds[observed_stellar_name] = observed_stellar
        seds[absorbed_name] = absorbed
        seds[dust_name] = dust

        # Add extra SEDs
        seds.update(**extra)

        # Plot
        plot_seds(seds, options=self.plot_options, distance=self.galaxy_distance, min_wavelength=self.min_wavelength,
                  max_wavelength=self.max_wavelength, min_flux=self.min_flux, max_flux=self.max_flux, unit=self.plot_unit,
                  path=path)

    # -----------------------------------------------------------------

    @property
    def do_plot_total_diffuse(self):
        return not self.has_total_diffuse_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_total_all(self):
        return not self.has_total_all_plot

    # -----------------------------------------------------------------

    def plot_total(self):

        """
        This function ...
        :return:
        """

        # Diffuse
        if self.do_plot_total_diffuse: self.plot_total_diffuse()

        # All
        if self.do_plot_total_all: self.plot_total_all()

    # -----------------------------------------------------------------

    @property
    def total_diffuse_plot_path(self):
        return fs.join(self.absorption_path, "total_diffuse.pdf")

    # -----------------------------------------------------------------

    @property
    def has_total_diffuse_plot(self):
        return fs.is_file(self.total_diffuse_plot_path)

    # -----------------------------------------------------------------

    def plot_total_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the total model (diffuse dust) ...")

        # Plot
        if self.has_total_absorption_sed_diffuse: self.plot_seds(self.total_observed_stellar_sed_diffuse, self.total_absorption_sed_diffuse, self.total_dust_sed_diffuse, self.total_diffuse_plot_path, absorbed_cells=self.total_spectral_absorption_curve)
        else: self.plot_seds(self.total_observed_stellar_sed_diffuse, self.total_spectral_absorption_curve, self.total_dust_sed_diffuse, self.total_diffuse_plot_path)

    # -----------------------------------------------------------------

    @property
    def total_all_plot_path(self):
        return fs.join(self.absorption_path, "total_all.pdf")

    # -----------------------------------------------------------------

    @property
    def has_total_all_plot(self):
        return fs.is_file(self.total_all_plot_path)

    # -----------------------------------------------------------------

    def plot_total_all(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the total model (all dust) ...")

        # Plot
        self.plot_seds(self.total_observed_stellar_sed_all, self.total_absorption_sed_all, self.total_dust_sed_all, self.total_all_plot_path)

    # -----------------------------------------------------------------

    @property
    def bulge_plot_path(self):
        return fs.join(self.absorption_path, "bulge.pdf")

    # -----------------------------------------------------------------

    @property
    def has_bulge_plot(self):
        return fs.is_file(self.bulge_plot_path)

    # -----------------------------------------------------------------

    def plot_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the old stellar bulge ...")

        # Plot
        if self.has_bulge_absorption_sed: self.plot_seds(self.bulge_observed_stellar_sed, self.bulge_absorption_sed, self.bulge_dust_sed, self.bulge_plot_path, absorbed_cells=self.bulge_spectral_absorption_curve)
        else: self.plot_seds(self.bulge_observed_stellar_sed, self.bulge_spectral_absorption_curve, self.bulge_dust_sed, self.bulge_plot_path)

    # -----------------------------------------------------------------

    @property
    def disk_plot_path(self):
        return fs.join(self.absorption_path, "disk.pdf")

    # -----------------------------------------------------------------

    @property
    def has_disk_plot(self):
        return fs.is_file(self.disk_plot_path)

    # -----------------------------------------------------------------

    def plot_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the old stellar disk ...")

        # Plot
        if self.has_disk_absorption_sed: self.plot_seds(self.disk_observed_stellar_sed, self.disk_absorption_sed, self.disk_dust_sed, self.disk_plot_path, absorbed_cells=self.disk_spectral_absorption_curve)
        else: self.plot_seds(self.disk_observed_stellar_sed, self.disk_spectral_absorption_curve, self.disk_dust_sed, self.disk_plot_path)

    # -----------------------------------------------------------------

    @property
    def old_plot_path(self):
        return fs.join(self.absorption_path, "old.pdf")

    # -----------------------------------------------------------------

    @property
    def has_old_plot(self):
        return fs.is_file(self.old_plot_path)

    # -----------------------------------------------------------------

    def plot_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the old stars ...")

        # Plot
        if self.has_old_absorption_sed: self.plot_seds(self.old_observed_stellar_sed, self.old_absorption_sed, self.old_dust_sed, self.old_plot_path, absorbed_cells=self.old_spectral_absorption_curve)
        else: self.plot_seds(self.old_observed_stellar_sed, self.old_spectral_absorption_curve, self.old_dust_sed, self.old_plot_path)

    # -----------------------------------------------------------------

    @property
    def young_plot_path(self):
        return fs.join(self.absorption_path, "young.pdf")

    # -----------------------------------------------------------------

    @property
    def has_young_plot(self):
        return fs.is_file(self.young_plot_path)

    # -----------------------------------------------------------------

    def plot_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the young stars ...")

        # Plot
        if self.has_young_absorption_sed: self.plot_seds(self.young_observed_stellar_sed, self.young_absorption_sed, self.young_dust_sed, self.young_plot_path, absorbed_cells=self.young_spectral_absorption_curve)
        else: self.plot_seds(self.young_observed_stellar_sed, self.young_spectral_absorption_curve, self.young_dust_sed, self.young_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_sfr_diffuse(self):
        return not self.has_sfr_diffuse_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_sfr_internal(self):
        return not self.has_sfr_internal_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_sfr_all(self):
        return not self.has_sfr_all_plot

    # -----------------------------------------------------------------

    def plot_sfr(self):

        """
        This function ...
        :return:
        """

        # Diffuse
        if self.do_plot_sfr_diffuse: self.plot_sfr_diffuse()

        # Internal
        if self.do_plot_sfr_internal: self.plot_sfr_internal()

        # All
        if self.do_plot_sfr_all: self.plot_sfr_all()

    # -----------------------------------------------------------------

    @property
    def sfr_diffuse_plot_path(self):
        return fs.join(self.absorption_path, "sfr_diffuse.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sfr_diffuse_plot(self):
        return fs.is_file(self.sfr_diffuse_plot_path)

    # -----------------------------------------------------------------

    def plot_sfr_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the star formation regions (diffuse dust) ...")

        # Plot
        if self.has_sfr_absorption_sed_diffuse: self.plot_seds(self.sfr_observed_stellar_sed_diffuse, self.sfr_absorption_sed_diffuse, self.sfr_dust_sed_diffuse, self.sfr_diffuse_plot_path, absorbed_cells=self.sfr_spectral_absorption_curve)
        else: self.plot_seds(self.sfr_observed_stellar_sed_diffuse, self.sfr_spectral_absorption_curve, self.sfr_dust_sed_diffuse, self.sfr_diffuse_plot_path)

    # -----------------------------------------------------------------

    @property
    def sfr_internal_plot_path(self):
        return fs.join(self.absorption_path, "sfr_internal.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sfr_internal_plot(self):
        return fs.is_file(self.sfr_internal_plot_path)

    # -----------------------------------------------------------------

    def plot_sfr_internal(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the star formation regions (internal dust) ...")

        # Plot
        self.plot_seds(self.sfr_observed_stellar_sed_internal, self.sfr_absorption_sed_internal, self.sfr_dust_sed_internal, self.sfr_internal_plot_path)

    # -----------------------------------------------------------------

    @property
    def sfr_all_plot_path(self):
        return fs.join(self.absorption_path, "sfr_all.pdf")

    # -----------------------------------------------------------------

    @property
    def has_sfr_all_plot(self):
        return fs.is_file(self.sfr_all_plot_path)

    # -----------------------------------------------------------------

    def plot_sfr_all(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the star formation regions (all dust) ...")

        # Plot
        self.plot_seds(self.sfr_observed_stellar_sed_all, self.sfr_absorption_sed_all, self.sfr_dust_sed_all, self.sfr_all_plot_path)

    # -----------------------------------------------------------------

    @property
    def do_plot_unevolved_diffuse(self):
        return not self.has_unevolved_diffuse_plot

    # -----------------------------------------------------------------

    @property
    def do_plot_unevolved_all(self):
        return not self.has_unevolved_all_plot

    # -----------------------------------------------------------------

    def plot_unevolved(self):

        """
        This function ...
        :return:
        """

        # Diffuse
        if self.do_plot_unevolved_diffuse: self.plot_unevolved_diffuse()

        # All
        if self.do_plot_unevolved_all: self.plot_unevolved_all()

    # -----------------------------------------------------------------

    @property
    def unevolved_diffuse_plot_path(self):
        return fs.join(self.absorption_path, "unevolved_diffuse.pdf")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_diffuse_plot(self):
        return fs.is_file(self.unevolved_diffuse_plot_path)

    # -----------------------------------------------------------------

    def plot_unevolved_diffuse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the unevolved stars (diffuse dust) ...")

        # Plot
        if self.has_unevolved_absorption_sed_diffuse: self.plot_seds(self.unevolved_observed_stellar_sed_diffuse, self.unevolved_absorption_sed_diffuse, self.unevolved_dust_sed_diffuse, self.unevolved_diffuse_plot_path, absorbed_cells=self.unevolved_spectral_absorption_curve)
        else: self.plot_seds(self.unevolved_observed_stellar_sed_diffuse, self.unevolved_spectral_absorption_curve, self.unevolved_dust_sed_diffuse, self.unevolved_diffuse_plot_path)

    # -----------------------------------------------------------------

    @property
    def unevolved_all_plot_path(self):
        return fs.join(self.absorption_path, "unevolved_all.pdf")

    # -----------------------------------------------------------------

    @property
    def has_unevolved_all_plot(self):
        return fs.is_file(self.unevolved_all_plot_path)

    # -----------------------------------------------------------------

    def plot_unevolved_all(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the absorption and dust SEDs for the unevolved stars (all dust) ...")

        # Plot
        self.plot_seds(self.unevolved_observed_stellar_sed_all, self.unevolved_absorption_sed_all, self.unevolved_dust_sed_all, self.unevolved_all_plot_path)

# -----------------------------------------------------------------

def show_properties(properties):
    for label in properties:
        if types.is_dictionary(properties[label]):
            print(" - " + fmt.bold + label + fmt.reset_bold + ":")
            for label2 in properties[label]: print("    * " + fmt.bold + label2 + fmt.reset_bold + ": " + tostr(properties[label][label2]))
        else: print(" - " + fmt.bold + label + fmt.reset_bold + ": " + tostr(properties[label]))

# -----------------------------------------------------------------

def show_values(**kwargs):
    for label in kwargs: print(" - " + fmt.bold + label + fmt.reset_bold + ": " + tostr(kwargs[label]))

# -----------------------------------------------------------------
