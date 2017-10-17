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

# Import the relevant PTS classes and modules
from ..component import AnalysisComponent
from ....core.tools import filesystem as fs
from ....core.basics.log import log
from ....core.tools.utils import lazyproperty
from ....core.simulation.output import SimulationOutput
from ....core.simulation.logfile import LogFile
from ....core.simulation.data import SimulationData

# -----------------------------------------------------------------

total = "total"
old = "old"
young = "young"
ionizing = "ionizing"
unevolved = "unevolved"
contributions = [total, old, young, ionizing, unevolved]

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

    @property
    def heating_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.wavelength_grid_heating

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.heating_wavelength_grid)

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.heating_path, "cell")

    # -----------------------------------------------------------------

    @lazyproperty
    def projected_heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.heating_path, "projected")

    # -----------------------------------------------------------------

    @property
    def total_contribution_simulation_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_simulation_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_ski_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_ski_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_contribution_output_path(self):

        """
        Thisn function ...
        :return:
        """

        return self.analysis_run.heating_output_path_for_contribution(total)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.total_contribution_output_path, self.galaxy_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.total_contribution_output)

    # -----------------------------------------------------------------

    @property
    def total_contribution_cell_properties_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.cell_properties_path

    # -----------------------------------------------------------------

    @lazyproperty
    def cell_properties(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.cell_properties

    # -----------------------------------------------------------------

    @property
    def total_contribution_absorption_filepath(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_absorption_data(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.absorption

    # -----------------------------------------------------------------

    @property
    def total_contribution_logfile_path(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def total_contribution_logfile(self):

        """
        This ufnction ...
        :return:
        """

        return LogFile.from_file(self.total_contribution_logfile_path)

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_datacube(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_datacube(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_sed(self):

        """
        Thisj function ...
        :return:
        """

        return self.total_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def total_contribution_total_faceon_sed(self):

        """
        This function ...
        :return:
        """

        return self.total_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_simulation_path(self):

        """
        Thisn function ...
        :return:
        """

        return self.analysis_run.heating_simulation_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_ski_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_ski_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def old_contribution_output_path(self):

        """
        Thisf unction ...
        :return:
        """

        return self.analysis_run.heating_output_path_for_contribution(old)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.old_contribution_output_path, self.galaxy_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.old_contribution_output)

    # -----------------------------------------------------------------

    @property
    def old_contribution_absorption_filepath(self):

        """
        This function ...
        :return:
        """

        return self.old_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_absorption_data(self):

        """
        This function ...
        :return:
        """

        return self.old_contribution_data.absorption

    # -----------------------------------------------------------------

    @property
    def old_contribution_logfile_path(self):

        """
        This function ...
        :return:
        """

        return self.old_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def old_contribution_logfile(self):

        """
        This function ...
        :return:
        """

        return LogFile.from_file(self.old_contribution_logfile_path)

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_datacube(self):

        """
        This function ...
        :return:
        """

        return self.old_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_datacube(self):

        """
        Thisj function ...
        :return:
        """

        return self.old_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_sed(self):

        """
        Thisfnctin ...
        :return:
        """

        return self.old_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def old_contribution_total_faceon_sed(self):

        """
        This function ...
        :return:
        """

        return self.old_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_simulation_path(self):

        """
        Thisn function ...
        :return:
        """

        return self.analysis_run.heating_simulation_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_ski_path(self):

        """
        Thisnf unction ...
        :return:
        """

        return self.analysis_run.heating_ski_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def young_contribution_output_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_output_path_for_contribution(young)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_output(self):

        """
        This functio n...
        :return:
        """

        return SimulationOutput.from_directory(self.young_contribution_output_path, self.galaxy_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.young_contribution_output)

    # -----------------------------------------------------------------

    @property
    def young_contribution_absorption_filepath(self):

        """
        Thins function ...
        :return:
        """

        return self.young_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_absorption_data(self):

        """
        This function ...
        :return:
        """

        return self.young_contribution_data.absorption

    # -----------------------------------------------------------------

    @property
    def young_contribution_logfile_path(self):

        """
        This function ...
        :return:
        """

        return self.young_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def young_contribution_logfile(self):

        """
        Thins ufnction ...
        :return:
        """

        return LogFile.from_file(self.young_contribution_logfile_path)

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_datacube(self):

        """
        Thisfunction ...
        :return:
        """

        return self.young_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_datacube(self):

        """
        This function ...
        :return:
        """

        return self.young_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_sed(self):

        """
        Tihs function ...
        :return:
        """

        return self.young_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def young_contribution_total_faceon_sed(self):

        """
        This function ...
        :return:
        """

        return self.young_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_simulation_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_simulation_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_ski_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_ski_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_output_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_output_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.ionizing_contribution_output_path, self.galaxy_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_data(self):

        """
        Thisn function ...
        :return:
        """

        return SimulationData.from_output(self.ionizing_contribution_output)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_absorption_filepath(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_absorption_data(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_contribution_data.absorption

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_logfile_path(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_contribution_logfile(self):

        """
        This functino ...
        :return:
        """

        return LogFile.from_file(self.ionizing_contribution_logfile_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_datacube(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.ionizing_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_datacube(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.ionizing_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_sed(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def ionizing_contribution_total_faceon_sed(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.ionizing_contribution_data.seds["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_simulation_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_simulation_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_ski_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_ski_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_output_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_output_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_output(self):

        """
        This function ...
        :return:
        """

        return SimulationOutput.from_directory(self.unevolved_contribution_output_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_data(self):

        """
        This function ...
        :return:
        """

        return SimulationData.from_output(self.unevolved_contribution_output)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_absorption_filepath(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_contribution_data.absorption_path

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_absorption_data(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_contribution_data.absorption

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_logfile_path(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_contribution_output.logfiles[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def unevolved_contribution_logfile(self):

        """
        Thisfunction ...
        :return:
        """

        return LogFile.from_file(self.unevolved_contribution_logfile_path)

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_datacube(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_contribution_data.images["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_datacube(self):

        """
        This function ...
        :return:
        """

        return self.unevolved_contribution_data.images["faceon"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_sed(self):

        """
        Thisfunction ...
        :return:
        """

        return self.unevolved_contribution_data.seds["earth"]["total"]

    # -----------------------------------------------------------------

    @property
    def unevolved_contribution_total_faceon_sed(self):

        """
        Thisfunction ...
        :return:
        """

        return self.unevolved_contribution_data.seds["faceon"]["total"]

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
