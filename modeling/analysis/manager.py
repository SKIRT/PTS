#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.manager Contains the AnalysisManager class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import numbers
from ...core.launch.manager import SimulationManager
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty
from .component import AnalysisComponent
from ...core.launch.tables import SimulationAssignmentTable
from ...core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

_add_instrument_command_name = "add_instrument" # add an instrument (orientation)
#_add_contributions_command_name = "add_contributions"  # add more contributions (make full instrument) for an existing instrument
_upgrade_instrument_command_name = "upgrade_instrument" # better name for command above?

# -----------------------------------------------------------------

class AnalysisManager(SimulationManager, AnalysisComponent):

    """
    This class ...
    """

    _log_section = "ANALYSIS MANAGER"

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        AnalysisComponent.__init__(self, no_config=True)
        SimulationManager.__init__(self, *args, **kwargs)

        # Add command
        self._commands = self._commands.copy()  # so the dictionary in the core/launch/manager module is not adapted
        self._commands[_add_instrument_command_name] = ("add_instrument_command", True, "simulate an additional orientation for one of the analysis simulations", None)

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.analysis_run.path, "manage")

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_current_path(self):

        """
        This function ...
        :return:
        """

        current_indices = fs.directories_in_path(self.manage_path, returns="name", convert=int)
        return fs.create_directory_in(self.manage_path, str(numbers.lowest_missing_integer(current_indices)))

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        AnalysisComponent.setup(self, **kwargs)

        # Input is shared
        self.config.shared_input = True
        self.config.write_moved = True
        self.config.write_relaunched = True

        # Set paths
        self.config.backup_dir_path = self.manage_current_path
        self.config.backup_dirname = "backup"
        self.config.backup_simulations = True
        self.config.backup_assignment = True

        # Set caching options
        # From generation manager
        # if self.config.cache_volume is not None:
        #
        #     # Get volume path
        #     volume_path = fs.get_volume_path(self.config.cache_volume)
        #     cache_path = fs.join(volume_path, "RT Modeling", self.environment.galaxy_name)
        #
        #     # Set path and root
        #     self.config.cache_path = cache_path
        #     self.config.cache_root = self.environment.path  # set modeling path as cache root path
        #
        #     # Auto-caching
        #     self.config.cache_output = self.config.cache_output
        #     self.config.cache_datacubes = self.config.cache_datacubes
        #     self.config.cache_misc = self.config.cache_misc
        #     self.config.cache_images = self.config.cache_images
        #
        #     # Cache after analysis of simulation
        #     self.config.cache_after_analysis = True

        # Set reference SEDs for plotting simulated SEDS
        reference_sed_paths = OrderedDict()
        reference_sed_paths["Observed clipped fluxes"] = self.environment.observed_sed_path
        reference_sed_paths["Observed truncated fluxes"] = self.environment.truncated_sed_path
        self.config.reference_seds = reference_sed_paths

        # TODO: what to do when assignment table is missing?

        # Load the assignment
        if not self.has_assignment: raise RuntimeError("No assignment table found")
        assignment = SimulationAssignmentTable.from_file(self.assignment_path)

        # Set input
        input = dict()
        input["assignment"] = assignment
        input["timing"] = self.timing_table
        input["memory"] = self.memory_table
        #input["status"] = status
        #input["info_tables"] = [self.parameters_table, self.chi_squared_table]
        #input["remotes"] = remotes
        #input["simulations"] = simulations

        # Interactive
        self.config.interactive = True

        # Set the output path
        self.config.output = self.analysis_run_path

        # Setup
        SimulationManager.setup(self, **input)

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.get_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.path

    # -----------------------------------------------------------------

    @lazyproperty
    def assignment_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.analysis_run_path, "assignment.dat")

    # -----------------------------------------------------------------

    @property
    def has_assignment(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.assignment_path)

    # -----------------------------------------------------------------

    @property
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return self.analysis_context.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return self.analysis_context.memory_table

    # -----------------------------------------------------------------

    @lazyproperty
    def add_instrument_definition(self):

        """
        This function ...
        :return:
        """

        definition = ConfigurationDefinition(write_config=False)
        return definition

    # -----------------------------------------------------------------

    def add_instrument_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # To be implemented
        pass

    # -----------------------------------------------------------------

    def __del__(self):

        """
        This function ...
        :return:
        """

        # Remove output directory if nothing was written
        if fs.is_empty(self.manage_current_path, recursive=True):

            # Remove
            fs.remove_directory(self.manage_current_path)

# -----------------------------------------------------------------
