#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.component Contains the DustHeatingAnalysisComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ..component import AnalysisRunComponent, earth_name, faceon_name
from ....core.tools import filesystem as fs
from ....core.tools.utils import lazyproperty
from ....core.simulation.logfile import LogFile
from ....core.tools import sequences
from ....core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

cell_dirname = "cell"
projected_dirname = "projected"
spectral_dirname = "spectral"

# -----------------------------------------------------------------

class DustHeatingAnalysisComponent(AnalysisRunComponent):

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
        super(DustHeatingAnalysisComponent, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_heating_path(self):
        return fs.create_directory_in(self.analysis_run.heating_path, cell_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def projected_heating_path(self):
        return fs.create_directory_in(self.analysis_run.heating_path, projected_dirname)

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_heating_path(self):
        return fs.create_directory_in(self.analysis_run.heating_path, spectral_dirname)

    # -----------------------------------------------------------------
    # TOTAL SIMULATION
    #   ABSORPTION: SEE ANALYSISRUNCOMPONENT
    # -----------------------------------------------------------------

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def total_contribution_spectral_absorption_filepath(self):
        return self.total_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_total_contribution_spectral_absorption(self):
        return self.total_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def total_contribution_spectral_emission_filepath(self):
        return self.total_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_total_contribution_spectral_emission(self):
        return self.total_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def total_contribution_logfile_path(self):
        return self.total_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_logfile(self):
        return LogFile.from_file(self.total_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def total_contribution_total_datacube(self):
        return self.total_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_datacube(self):
        return self.total_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_sed(self):
        return self.total_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_sed(self):
        return self.total_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------
    # OLD SIMULATION
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def old_contribution_absorption_filepath(self):
        if self.old_contribution_data.has_absorption: return self.old_contribution_data.absorption_path
        elif self.old_contribution_data.has_isrf: return self.old_contribution_data.isrf_path
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_data(self):
        if self.old_contribution_data.has_absorption: return self.old_contribution_data.absorption
        elif self.old_contribution_data.has_isrf: return self.old_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_column_names(self):
        return fs.get_column_names(self.old_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_column_units(self):
        return fs.get_column_units(self.old_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.old_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_column_index(self):
        return self.old_contribution_absorption_column_names.index(self.old_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_unit(self):
        return u(self.old_contribution_absorption_column_units[self.old_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_luminosities(self):
        return fs.get_column(self.old_contribution_absorption_filepath, self.old_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def old_contribution_spectral_absorption_filepath(self):
        return self.old_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_old_contribution_spectral_absorption(self):
        return self.old_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def old_contribution_spectral_emission_filepath(self):
        return self.old_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_old_contribution_spectral_emission(self):
        return self.old_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def old_contribution_logfile_path(self):
        return self.old_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_logfile(self):
        return LogFile.from_file(self.old_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def old_contribution_total_datacube(self):
        return self.old_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_datacube(self):
        return self.old_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_sed(self):
        return self.old_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_sed(self):
        return self.old_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------
    # YOUNG SIMULATION
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def young_contribution_absorption_filepath(self):
        if self.young_contribution_data.has_absorption: return self.young_contribution_data.absorption_path
        elif self.young_contribution_data.has_isrf: return self.young_contribution_data.isrf_path
        else: raise IOError("No absorption path")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_data(self):
        if self.young_contribution_data.has_absorption: return self.young_contribution_data.absorption
        elif self.young_contribution_data.has_isrf: return self.young_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_column_names(self):
        return fs.get_column_names(self.young_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_column_units(self):
        return fs.get_column_units(self.young_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.young_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_column_index(self):
        return self.young_contribution_absorption_column_names.index(self.young_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_unit(self):
        return u(self.young_contribution_absorption_column_units[self.young_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_luminosities(self):
        return fs.get_column(self.young_contribution_absorption_filepath, self.young_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def young_contribution_spectral_absorption_filepath(self):
        return self.young_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_young_contribution_spectral_absorption(self):
        return self.young_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def young_contribution_spectral_emission_filepath(self):
        return self.young_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_young_contribution_spectral_emission(self):
        return self.young_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def young_contribution_logfile_path(self):
        return self.young_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_logfile(self):
        return LogFile.from_file(self.young_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def young_contribution_total_datacube(self):
        return self.young_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_datacube(self):
        return self.young_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_sed(self):
        return self.young_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_sed(self):
        return self.young_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------
    # IONIZING (SFR) SIMULATION
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_absorption_filepath(self):
        if self.ionizing_contribution_data.has_absorption: return self.ionizing_contribution_data.absorption_path
        elif self.ionizing_contribution_data.has_isrf: return self.ionizing_contribution_data.isrf_path
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_data(self):
        if self.ionizing_contribution_data.has_absorption: return self.ionizing_contribution_data.absorption
        elif self.ionizing_contribution_data.has_isrf: return self.ionizing_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_column_names(self):
        return fs.get_column_names(self.ionizing_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_column_units(self):
        return fs.get_column_units(self.ionizing_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.ionizing_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_column_index(self):
        return self.ionizing_contribution_absorption_column_names.index(self.ionizing_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_unit(self):
        return u(self.ionizing_contribution_absorption_column_units[self.ionizing_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_luminosities(self):
        return fs.get_column(self.ionizing_contribution_absorption_filepath, self.ionizing_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_spectral_absorption_filepath(self):
        return self.ionizing_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_ionizing_contribution_spectral_absorption(self):
        return self.ionizing_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_spectral_emission_filepath(self):
        return self.ionizing_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_ionizing_contribution_spectral_emission(self):
        return self.ionizing_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_logfile_path(self):
        return self.ionizing_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_logfile(self):
        return LogFile.from_file(self.ionizing_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_datacube(self):
        return self.ionizing_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_datacube(self):
        return self.ionizing_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_sed(self):
        return self.ionizing_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_sed(self):
        return self.ionizing_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------
    # EXTRA SIMULATION
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def extra_contribution_absorption_filepath(self):
        if self.extra_contribution_data.has_absorption: return self.extra_contribution_data.absorption_path
        elif self.extra_contribution_data.has_isrf: return self.extra_contribution_data.isrf_path
        else: raise IOError("No absorption path")

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_data(self):
        if self.extra_contribution_data.has_absorption: return self.extra_contribution_data.absorption
        elif self.extra_contribution_data.has_isrf: return self.extra_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_column_names(self):
        return fs.get_column_names(self.extra_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_column_units(self):
        return fs.get_column_units(self.extra_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.extra_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_column_index(self):
        return self.extra_contribution_absorption_column_names.index(self.extra_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_unit(self):
        return u(self.extra_contribution_absorption_column_units[self.extra_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_absorption_luminosities(self):
        return fs.get_column(self.extra_contribution_absorption_filepath, self.extra_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def extra_contribution_spectral_absorption_filepath(self):
        return self.extra_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_extra_contribution_spectral_absorption(self):
        return self.extra_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def extra_contribution_spectral_emission_filepath(self):
        return self.extra_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_extra_contribution_spectral_emission(self):
        return self.extra_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def extra_contribution_logfile_path(self):
        return self.extra_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_contribution_logfile(self):
        return LogFile.from_file(self.extra_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def extra_contribution_total_datacube(self):
        return self.extra_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def extra_contribution_total_faceon_datacube(self):
        return self.extra_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def extra_contribution_total_sed(self):
        return self.extra_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def extra_contribution_total_faceon_sed(self):
        return self.extra_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_absorption_filepath(self):
        if self.unevolved_contribution_data.has_absorption: return self.unevolved_contribution_data.absorption_path
        elif self.unevolved_contribution_data.has_isrf: return self.unevolved_contribution_data.isrf_path
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_data(self):
        if self.unevolved_contribution_data.has_absorption: return self.unevolved_contribution_data.absorption
        elif self.unevolved_contribution_data.has_isrf: return self.unevolved_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_column_names(self):
        return fs.get_column_names(self.unevolved_contribution_absorption_filepath, capitalize=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_column_units(self):
        return fs.get_column_units(self.unevolved_contribution_absorption_filepath)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.unevolved_contribution_absorption_column_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_column_index(self):
        return self.unevolved_contribution_absorption_column_names.index(self.unevolved_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_unit(self):
        return u(self.unevolved_contribution_absorption_column_units[self.unevolved_contribution_absorption_column_index])

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_luminosities(self):
        return fs.get_column(self.unevolved_contribution_absorption_filepath, self.unevolved_contribution_absorption_column_index, float, method="pandas")

    # -----------------------------------------------------------------
    #   SPECTRAL ABSORPTION & EMISSION
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_spectral_absorption_filepath(self):
        return self.unevolved_contribution_data.spectral_absorption_path

    # -----------------------------------------------------------------

    @property
    def has_unevolved_contribution_spectral_absorption(self):
        return self.unevolved_contribution_data.has_spectral_absorption

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_spectral_emission_filepath(self):
        return self.unevolved_contribution_data.spectral_emission_path

    # -----------------------------------------------------------------

    @property
    def has_unevolved_contribution_spectral_emission(self):
        return self.unevolved_contribution_data.has_spectral_emission

    # -----------------------------------------------------------------
    #   LOGFILES
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_logfile_path(self):
        return self.unevolved_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_logfile(self):
        return LogFile.from_file(self.unevolved_contribution_logfile_path)

    # -----------------------------------------------------------------
    #   DATACUBES & SEDs
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_datacube(self):
        return self.unevolved_contribution_data.images[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_datacube(self):
        return self.unevolved_contribution_data.images[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_sed(self):
        return self.unevolved_contribution_data.seds[earth_name]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_sed(self):
        return self.unevolved_contribution_data.seds[faceon_name]["total"]

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):
        return self.analysis_run.wavelength_grid

# -----------------------------------------------------------------
