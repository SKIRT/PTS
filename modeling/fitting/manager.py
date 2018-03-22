#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.generation_status View the status of the simulations of a certain generation.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.basics.configuration import prompt_yn
from pts.core.tools import formatting as fmt
from pts.core.tools.stringify import tostr
from pts.core.remote.ensemble import SKIRTRemotesEnsemble
from pts.core.tools import numbers
from pts.core.launch.manager import SimulationManager, extra_columns
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.core.simulation.remote import is_analysed_status
from pts.core.launch.batchlauncher import SimulationStatusTable
from pts.core.tools import sequences
from pts.modeling.fitting.generation import check_simulation_paths, correct_simulation_paths
from ...core.tools.utils import lazyproperty
from .component import FittingComponent

# -----------------------------------------------------------------

# Number of digits for parameter values
parameters_ndigits = 3

# -----------------------------------------------------------------

class GenerationManager(SimulationManager, FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        #super(GenerationManager, self).__init__(*args, **kwargs)
        FittingComponent.__init__(self, no_config=True)
        SimulationManager.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.fitting_run.generations_path, "manage__" + self.generation.name)

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
        #super(GenerationManager, self).setup(**kwargs)
        FittingComponent.setup(self, **kwargs)
        #SimulationManager.setup(self, **kwargs)

        self.config.shared_input = True
        # manager.config.write_status = config.write_status
        self.config.write_moved = True
        self.config.write_relaunched = True
        # manager.config.write_commands = config.write_commands

        # Set paths
        #backup_path = fs.join(manage_current_path, "backup")

        # Set paths
        # manager.config.path = manage_current_path
        self.config.backup_dir_path = self.manage_current_path
        self.config.backup_dirname = "backup"
        self.config.backup_simulations = True
        self.config.backup_assignment = True

        # Set caching options
        if self.config.cache_volume is not None:

            # Get volume path
            volume_path = fs.get_volume_path(self.config.cache_volume)
            cache_path = fs.join(volume_path, "RT Modeling", self.environment.galaxy_name)

            # Set path and root
            self.config.cache_path = cache_path
            self.config.cache_root = self.environment.path  # set modeling path as cache root path

            # Auto-caching
            self.config.cache_output = self.config.cache_output
            self.config.cache_datacubes = self.config.cache_datacubes
            self.config.cache_misc = self.config.cache_misc
            self.config.cache_images = self.config.cache_images

            # Cache after analysis of simulation
            self.config.cache_after_analysis = True

        # Set reference SEDs for plotting simulated SEDS
        reference_sed_paths = OrderedDict()
        reference_sed_paths["Observed clipped fluxes"] = self.environment.observed_sed_path
        reference_sed_paths["Observed truncated fluxes"] = self.environment.truncated_sed_path
        self.config.reference_seds = reference_sed_paths

        # Set screen script paths
        self.config.screen_scripts = fs.files_in_path(self.generation.path, extension="sh")

        # ANALYSE?
        #self.config.analyse = self.config.analyse
        #self.config.analysis = self.config.analysis

        # Create remotes ensemble
        if not self.config.offline and self.generation.has_assignment_table: remotes = SKIRTRemotesEnsemble(self.generation.host_ids)
        else: remotes = None

        # Check
        if not self.generation.has_assignment_table:

            # raise RuntimeError("No assignment for this generation")
            log.warning("Assignment table is not present for this generation, trying to work without ...")

        if len(self.chi_squared_table) > 0:

            # Get the maximum chi squared value
            min_chi_squared = self.chi_squared_table.best_chi_squared
            max_chi_squared = self.chi_squared_table.worst_chi_squared
            max_chi_squared_magnitude = numbers.order_of_magnitude(max_chi_squared)

            # Set number of digits for chi squared values
            chisq_ndecimal = 1
            chisq_ndigits = max_chi_squared_magnitude + 1 + chisq_ndecimal

            print(fmt.bold + "Best chi squared: " + fmt.reset + tostr(min_chi_squared))
            print(fmt.bold + "Worst chi squared: " + fmt.reset + tostr(max_chi_squared))
            print("")

        # Get status
        status = self.get_generation_status(remotes)

        # Determine input: assignment table or simulation objects
        if self.generation.has_assignment_table:
            assignment = self.generation.assignment_table
            simulations = None
        else:
            assignment = None
            simulations = self.generation.simulations_basic  # basic because there will be no information about which host ID or simulation ID so simulation files cannot be located

            # Fix simulation status
            for simulation in simulations:
                simulation_status = status.get_status(simulation.name)
                if simulation_status == "analysed": simulation.analysed = True

        to_fix_status = dict()

        # Loop over the status entries and check whether analysed simulations actually have their chi squared value set
        for simulation_name in status.simulation_names:

            # Get the simulation status
            simulation_status = status.get_status(simulation_name)

            # Check
            analysed = self.generation.is_analysed(simulation_name)
            if is_analysed_status(simulation_status) and not analysed:

                # Warning
                log.warning("Simulation '" + simulation_name + "' is supposed to be analysed, but chi squared value is missing form the chi squared table")

                # Correct the status
                if self.config.correct_status:

                    # Inform
                    log.warning("Correcting the status ...")

                    # Check the analysis output
                    has_misc = self.generation.has_misc_output(simulation_name)
                    has_plotting = self.generation.has_plotting_output(simulation_name)
                    has_extraction = self.generation.has_extraction_output(simulation_name)

                    # Has any analysis output
                    if has_extraction or has_plotting or has_misc:

                        analysed = []
                        if has_extraction: analysed.append("extraction")
                        if has_plotting: analysed.append("plotting")
                        if has_misc: analysed.append("misc")
                        simulation_status = "analysed: " + ", ".join(analysed)

                    # Has simulation output
                    elif self.generation.is_retrieved(simulation_name): simulation_status = "retrieved"

                    # No simulation output
                    else: simulation_status = "unknown"

                    # Show new status
                    log.warning("New status: '" + simulation_status + "'")

                    # Set the status for the simulation in the table
                    # status.set_status(simulation_name, simulation_status) # NO: GET_STATUS METHOD IS MEMOIZED SO WE NEED NEW SIMULATION STATUS TABLE OBJECT
                    to_fix_status[simulation_name] = simulation_status

                    # Fix the simulation properties
                    if not self.generation.has_simulation(simulation_name): continue
                    else:

                        # Load the simulation object
                        simulation = self.generation.get_simulation(simulation_name)

                        # Needs fixing?
                        if simulation.analysed:

                            # Fix
                            log.warning("Fixing simulation 'analysed' flag ...")

                            # Unset analysed flag
                            simulation.analysed = False

                            # Save the simulation object
                            simulation.save()

                        # Needs fixing?
                        if "SEDFitModelAnalyser" in simulation.analysed_extra:

                            # Fix
                            log.warning("Fixing simulation 'analysed_extra' classes ...")

                            # Unset extra analysis
                            simulation.analysed_extra = sequences.removed(simulation.analysed_extra, "SEDFitModelAnalyser")

                            # Save the simulation object
                            simulation.save()

        # Needs fixing
        if len(to_fix_status) > 0:

            # Get number of simulations for which fixing is necessary
            nfix = len(to_fix_status)

            # Warning
            log.warning("Simulation status table needs fixing for " + str(nfix) + " simulations")

            # Create lists
            simulation_names = []
            status_list = []

            # Loop over the simulations
            for simulation_name in status.simulation_names:

                # Get the correct status
                if simulation_name in to_fix_status: correct_status = to_fix_status[simulation_name]
                else: correct_status = status.get_status(simulation_name)

                # Add to columns
                simulation_names.append(simulation_name)
                status_list.append(correct_status)

            # Create new status table (because status table class is full with lazyproperties and memoized methods)
            status = SimulationStatusTable.from_columns(simulation_names, status_list)

        # Check paths
        if self.config.check_paths:

            # Loop over all simulations
            for simulation_name in self.generation.simulation_names:

                # Get simulation object
                if not self.generation.has_simulation(simulation_name): continue  # simulation objects created on the fly are assumed to have correct paths
                simulation = self.generation.get_simulation(simulation_name)

                # Check paths
                try:
                    check_simulation_paths(simulation)
                except RuntimeError as e:
                    if self.config.correct_paths:
                        log.warning(str(e))
                        log.warning("Fixing simulation paths ...")
                        correct_simulation_paths(simulation, confirm=self.config.confirm_correction)
                    else: raise e

        # Check whether simulations with chi squared (analysed simuations) also have the other analysis output
        if self.config.check_analysis:

            # Simulation names to have status reset
            reset_status = []

            # Loop over the simulations
            for simulation_name in self.generation.simulation_names:

                # Has simulation object?
                if not self.generation.has_simulation(simulation_name): continue

                # Check different analysis output
                has_misc = self.generation.has_misc_output(simulation_name)
                has_plotting = self.generation.has_plotting_output(simulation_name)
                has_extraction = self.generation.has_extraction_output(simulation_name)
                has_timing = self.generation.has_timing(simulation_name)
                has_memory = self.generation.has_memory(simulation_name)
                has_chi_squared = self.generation.is_analysed(simulation_name)

                # Load the simulation
                simulation = self.generation.get_simulation(simulation_name)

                # Check extraction
                if not has_extraction and simulation.analysed_any_extraction:

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' appears to have no extraction output but extraction is supposed to be (partly) executed")

                    # Unset analysed flag
                    if simulation.analysed: simulation.set_analysed(False)
                    simulation.unset_analysed_extraction()
                    reset_status.append(simulation_name)

                # Check plotting
                if not has_plotting and simulation.analysed_any_plotting:

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' appears to have no plotting output but plotting is supposed to be (partly) executed")

                    # Unset analysed flag
                    if simulation.analysed: simulation.set_analysed(False)
                    simulation.unset_analysed_plotting()
                    reset_status.append(simulation_name)

                # Check misc
                if not has_misc and simulation.analysed_any_misc:

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' appears to have no miscellaneous output but miscellaneous analysis is supposed to be (partly) executed")

                    # Unset analysed flag
                    if simulation.analysed: simulation.set_analysed(False)
                    simulation.unset_analysed_misc()
                    reset_status.append(simulation_name)

                # Check batch
                if not (has_timing or has_memory) and simulation.analysed_batch:

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' appears to have no timing or memory entry but batch analysis is supposed to be (partly) executed")

                    # Unset analysed flag
                    if simulation.analysed: simulation.set_analysed(False)
                    simulation.unset_analysed_batch()
                    reset_status.append(simulation_name)

                # Check extra
                if not has_chi_squared and simulation.analysed_all_extra:

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' appears to have no chi squared value but extra analysis is supposed to be executed")

                    # Unset analysed flag
                    if simulation.analysed: simulation.set_analysed(False)
                    simulation.unset_analysed_extra()
                    reset_status.append(simulation_name)

                # Check
                # if has_chi_squared and not (has_extraction and has_plotting and has_misc and (has_timing or has_memory)):
                if has_chi_squared and not (has_extraction and has_plotting and has_misc):  # sometimes timing and memory entries are missing but that is not actually a problem

                    # Show warning
                    log.warning("Simulation '" + simulation_name + "' has a chi squared value but other analysis output seems to be missing")

                    # Get the chi squared value
                    chisq = self.chi_squared_table.chi_squared_for(simulation_name)

                    # Remove from chi squared table and re-analyse?
                    if prompt_yn("reanalyse", "remove the chi squared value of " + str(chisq) + " for simulation '" + simulation_name + "' and re-analyse this simulation?"):

                        # Debugging
                        log.debug("Removing simulation '" + simulation_name + "' simulation from the chi squared table ...")

                        # Remove from the chi squared table
                        self.chi_squared_table.remove_simulation(simulation_name)
                        self.chi_squared_table.save()

                        # Debugging
                        log.debug("Adding re-analyse command for simulation '" + simulation_name + "' ...")

                        # Add re-analyse command
                        reanalyse_command = "reanalyse '" + simulation_name + "' --analysis/basic/local --analysis/basic/ignore_missing_data --analysis/ignore_missing_data --analysis/basic/ignore_bad --analysis/basic/not_skip_ignored_bad_convolution"
                        self.config.commands.append(reanalyse_command)

                # Save the simulation
                simulation.save()

            # Reset status
            for simulation_name in reset_status:

                # Debugging
                log.debug("Resetting status for simulation '" + simulation_name + "' ...")

                # Reset
                status.reset_for_simulation(simulation_name, "retrieved")

        input = dict()
        input["assignment"] = assignment
        input["timing"] = self.timing_table
        input["memory"] = self.memory_table
        input["status"] = status
        input["info_tables"] = [self.parameters_table, self.chi_squared_table]
        input["remotes"] = remotes
        input["simulations"] = simulations

        # Setup
        SimulationManager.setup(self, **input)

    # -----------------------------------------------------------------

    def get_generation_status(self, remotes):

        """
        This function ...
        :param remotes:
        :return:
        """

        # Get the status of the simulations
        return self.generation.get_status(remotes, lazy=self.config.lazy, find_simulations=self.config.find_simulations,
                                       find_remotes=self.config.find_remotes, produce_missing=self.config.produce_missing,
                                       retrieve=self.config.retrieve, check_paths=self.config.check_paths,
                                       correct_paths=self.config.correct_paths, confirm_correction=self.config.confirm_correction,
                                       fix_success=self.config.fix_success)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.load_fitting_run(self.config.run)

    # -----------------------------------------------------------------

    @property
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.timing_table

    # -----------------------------------------------------------------

    @property
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.memory_table

    # -----------------------------------------------------------------

    @lazyproperty
    def generation(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.get_generation(self.config.generation)

    # -----------------------------------------------------------------

    @property
    def chi_squared_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.chi_squared_table

    # -----------------------------------------------------------------

    @property
    def parameters_table(self):

        """
        This function ...
        :return:
        """

        return self.generation.parameters_table

    # -----------------------------------------------------------------

    @property
    def parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.free_parameter_labels

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
