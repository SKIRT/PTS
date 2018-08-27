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
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty
from ....core.simulation.output import SimulationOutput
from ....core.simulation.logfile import LogFile
from ....core.simulation.data import SimulationData
from ....core.tools import sequences

# -----------------------------------------------------------------

total = "total"
old = "old"
young = "young"
ionizing = "ionizing"
unevolved = "unevolved"
contributions = [total, old, young, ionizing, unevolved]

# -----------------------------------------------------------------

cell_dirname = "cell"
projected_dirname = "projected"
spectral_dirname = "spectral"

# -----------------------------------------------------------------

class DustHeatingAnalysisComponent(AnalysisComponent):
    
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

        # The analysis run
        self.analysis_run = None

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
    #   GENERAL
    # -----------------------------------------------------------------

    @property
    def total_contribution_simulation_path(self):
        #return self.analysis_run.heating_simulation_path_for_contribution(total)
        return self.analysis_run.simulation_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_ski_path(self):
        #return self.analysis_run.heating_ski_path_for_contribution(total)
        return self.analysis_run.ski_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_output_path(self):
        #return self.analysis_run.heating_output_path_for_contribution(total)
        return self.analysis_run.output_path_for_contribution(total)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def total_contribution_output(self):
        #return SimulationOutput.from_directory(self.total_contribution_output_path, self.galaxy_name)
        return self.model.total_simulation_output

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def total_contribution_data(self):
        #return SimulationData.from_output(self.total_contribution_output)
        return self.model.total_simulation_data

    # -----------------------------------------------------------------

    @property
    def total_contribution_cell_properties_filepath(self):
        return self.total_contribution_data.cell_properties_path

    # -----------------------------------------------------------------

    @property
    def cell_properties_path(self):
        return self.total_contribution_cell_properties_filepath

    # -----------------------------------------------------------------

    @property
    def cell_properties(self):
        return self.total_contribution_data.cell_properties

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_mass_fractions(self):
        return np.asarray(self.cell_properties["Mass fraction"])

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def total_contribution_absorption_filepath(self):
        return self.total_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_data(self):
        if self.total_contribution_data.has_absorption: return self.total_contribution_data.absorption
        elif self.total_contribution_data.has_isrf: return self.total_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.total_contribution_absorption_data.colnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_unit(self):
        return self.total_contribution_absorption_data.column_unit(self.total_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_luminosities(self):
        return np.asarray(self.total_contribution_absorption_data[self.total_contribution_absorption_column_name])

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
        return self.total_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_datacube(self):
        return self.total_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_sed(self):
        return self.total_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_sed(self):
        return self.total_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------
    # OLD SIMULATION
    #   GENERAL
    # -----------------------------------------------------------------

    @property
    def old_contribution_simulation_path(self):
        #return self.analysis_run.heating_simulation_path_for_contribution(old)
        return self.analysis_run.simulation_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_ski_path(self):
        #return self.analysis_run.heating_ski_path_for_contribution(old)
        return self.analysis_run.ski_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_output_path(self):
        #return self.analysis_run.heating_output_path_for_contribution(old)
        return self.analysis_run.output_path_for_contribution(old)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def old_contribution_output(self):
        #return SimulationOutput.from_directory(self.old_contribution_output_path, self.galaxy_name)
        return self.model.old_simulation_output

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def old_contribution_data(self):
        #return SimulationData.from_output(self.old_contribution_output)
        return self.model.old_simulation_data

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def old_contribution_absorption_filepath(self):
        return self.old_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_data(self):
        if self.old_contribution_data.has_absorption: return self.old_contribution_data.absorption
        elif self.old_contribution_data.has_isrf: return self.old_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.old_contribution_absorption_data.colnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_unit(self):
        return self.old_contribution_absorption_data.column_unit(self.old_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_luminosities(self):
        return np.asarray(self.old_contribution_absorption_data[self.old_contribution_absorption_column_name])

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
        return self.old_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_datacube(self):
        return self.old_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_sed(self):
        return self.old_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_sed(self):
        return self.old_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------
    # YOUNG SIMULATION
    #   GENERAL
    # -----------------------------------------------------------------

    @property
    def young_contribution_simulation_path(self):
        #return self.analysis_run.heating_simulation_path_for_contribution(young)
        return self.analysis_run.simulation_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_ski_path(self):
        #return self.analysis_run.heating_ski_path_for_contribution(young)
        return self.analysis_run.ski_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_output_path(self):
        #return self.analysis_run.heating_output_path_for_contribution(young)
        return self.analysis_run.output_path_for_contribution(young)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def young_contribution_output(self):
        #return SimulationOutput.from_directory(self.young_contribution_output_path, self.galaxy_name)
        return self.model.young_simulation_output

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def young_contribution_data(self):
        #return SimulationData.from_output(self.young_contribution_output)
        return self.model.young_simulation_data

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def young_contribution_absorption_filepath(self):
        return self.young_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_data(self):
        if self.young_contribution_data.has_absorption: return self.young_contribution_data.absorption
        elif self.young_contribution_data.has_isrf: return self.young_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.young_contribution_absorption_data.colnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_unit(self):
        return self.young_contribution_absorption_data.column_unit(self.young_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_luminosities(self):
        return np.asarray(self.young_contribution_absorption_data[self.young_contribution_absorption_column_name])

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
        return self.young_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_datacube(self):
        return self.young_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_sed(self):
        return self.young_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_sed(self):
        return self.young_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------
    # IONIZING (SFR) SIMULATION
    #   GENERAL
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_simulation_path(self):
        #return self.analysis_run.heating_simulation_path_for_contribution(ionizing)
        return self.analysis_run.simulation_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_ski_path(self):
        #return self.analysis_run.heating_ski_path_for_contribution(ionizing)
        return self.analysis_run.ski_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_output_path(self):
        #return self.analysis_run.heating_output_path_for_contribution(ionizing)
        return self.analysis_run.output_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def ionizing_contribution_output(self):
        #return SimulationOutput.from_directory(self.ionizing_contribution_output_path, self.galaxy_name)
        return self.model.sfr_simulation_output

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def ionizing_contribution_data(self):
        #return SimulationData.from_output(self.ionizing_contribution_output)
        return self.model.sfr_simulation_data

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_absorption_filepath(self):
        return self.ionizing_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_data(self):
        if self.ionizing_contribution_data.has_absorption: return self.ionizing_contribution_data.absorption
        elif self.ionizing_contribution_data.has_isrf: return self.ionizing_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.ionizing_contribution_absorption_data.colnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_unit(self):
        return self.ionizing_contribution_absorption_data.column_unit(self.ionizing_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_luminosities(self):
        return np.asarray(self.ionizing_contribution_absorption_data[self.ionizing_contribution_absorption_column_name])

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
        return self.ionizing_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_datacube(self):
        return self.ionizing_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_sed(self):
        return self.ionizing_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_sed(self):
        return self.ionizing_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------
    # UNEVOLVED SIMULATION
    #   GENERAL
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_simulation_path(self):
        #return self.analysis_run.heating_simulation_path_for_contribution(unevolved)
        return self.analysis_run.simulation_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_ski_path(self):
        #return self.analysis_run.heating_ski_path_for_contribution(unevolved)
        return self.analysis_run.ski_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_output_path(self):
        #return self.analysis_run.heating_output_path_for_contribution(unevolved)
        return self.analysis_run.output_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def unevolved_contribution_output(self):
        #return SimulationOutput.from_directory(self.unevolved_contribution_output_path)
        return self.model.unevolved_simulation_output

    # -----------------------------------------------------------------

    #@lazyproperty
    @property
    def unevolved_contribution_data(self):
        #return SimulationData.from_output(self.unevolved_contribution_output)
        return self.model.unevolved_simulation_data

    # -----------------------------------------------------------------
    #   ABSORPTION
    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_absorption_filepath(self):
        return self.unevolved_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_data(self):
        if self.unevolved_contribution_data.has_absorption: return self.unevolved_contribution_data.absorption
        elif self.unevolved_contribution_data.has_isrf: return self.unevolved_contribution_data.isrf
        else: raise IOError("No absorption data")

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_column_name(self):
        abs_colnames = ["Absorbed bolometric luminosity", "Bolometric luminosity absorbed in cell"]
        return sequences.find_single_in_both(abs_colnames, self.unevolved_contribution_absorption_data.colnames)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_unit(self):
        return self.unevolved_contribution_absorption_data.column_unit(self.unevolved_contribution_absorption_column_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_luminosities(self):
        return np.asarray(self.unevolved_contribution_absorption_data[self.unevolved_contribution_absorption_column_name])

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
        return self.unevolved_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_datacube(self):
        return self.unevolved_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_sed(self):
        return self.unevolved_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_sed(self):
        return self.unevolved_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------
    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DustHeatingAnalysisComponent, self).setup(**kwargs)

        # Load the analysis run
        self.load_run()

    # -----------------------------------------------------------------

    def load_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the analysis run " + self.config.run + " ...")

        # Get the run
        self.analysis_run = self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def model(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.model

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    @property
    def cell_coordinates_filepath(self):

        """
        This function ...
        :return:
        """

        # SKIRT7
        if self.total_contribution_data.has_absorption: return self.total_contribution_absorption_filepath

        # SKIRT8
        else: return self.cell_properties_path

    # -----------------------------------------------------------------

    @property
    def cell_x_coordinates_colname(self):
        return "X coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def cell_y_coordinates_colname(self):
        return "Y coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def cell_z_coordinates_colname(self):
        return "Z coordinate of cell center"

    # -----------------------------------------------------------------

    @property
    def cell_x_coordinates(self):

        """
        This function ...
        :return:
        """

        # SKIRT7
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_x_coordinates_colname])

        # SKIRT8
        else: return self.model.cell_x_coordinates

    # -----------------------------------------------------------------

    @property
    def cell_y_coordinates(self):

        """
        This function ...
        :return:
        """

        # SKIRT7
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_y_coordinates_colname])

        # SKIRT8
        else: return self.model.cell_y_coordinates

    # -----------------------------------------------------------------

    @property
    def cell_z_coordinates(self):

        """
        This function ...
        :return:
        """

        # SKIRT7
        if self.total_contribution_data.has_absorption: return np.asarray(self.total_contribution_absorption_data[self.cell_z_coordinates_colname])

        # SKIRT8
        else: return self.model.cell_z_coordinates

# -----------------------------------------------------------------
