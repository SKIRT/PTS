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
from ...core.basics.configuration import prompt_yn, ConfigurationDefinition
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr, yes_or_no
from ...core.remote.ensemble import SKIRTRemotesEnsemble
from ...core.tools import numbers
from ...core.launch.manager import SimulationManager, extra_columns
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.simulation.remote import is_analysed_status
from ...core.launch.tables import SimulationStatusTable
from ...core.tools import sequences
from .generation import check_simulation_paths, correct_simulation_and_analysis_paths
from ...core.tools.utils import lazyproperty
from .component import FittingComponent
from ...core.launch.batch import MissingSimulation
from ...core.launch.options import AnalysisOptions
from ..config.explore_analysis import definition as analysis_definition

# -----------------------------------------------------------------

# Number of digits for parameter values
parameters_ndigits = 3

# -----------------------------------------------------------------

_simulations_command_name = "simulations"
_set_analysis_command_name = "set_analysis"

# -----------------------------------------------------------------

class GenerationManager(SimulationManager, FittingComponent):

    """
    This class ...
    """

    _log_section = "GENERATION MANAGER"

    # -----------------------------------------------------------------

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

        # Add command
        self._commands = self._commands.copy()  # so the dictionary in the core/launch/manager module is not adapted
        self._commands[_simulations_command_name] = ("show_simulations_command", True, "show info on the simulation files/objects", None)
        self._commands[_set_analysis_command_name] = ("set_analysis_command", True, "set the analysis options for one or multiple simulations", "simulations")

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_path(self):
        return fs.create_directory_in(self.fitting_run.generations_path, self.generation.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def manage_path(self):
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

        # Set options
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
        if self.config.cache_volume is not None: self.set_caching()

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

        # Check whether assignment table is present
        if not self.generation.has_assignment_table: log.warning("Assignment table is not present for this generation, trying to work without ...")

        # Show chi squared values?
        if len(self.chi_squared_table) > 0: self.show_chi_squared()

        # Debugging
        log.debug("Getting the simulation assignment ...")

        # Get assignment
        if self.generation.has_assignment_table: assignment = self.generation.assignment_table
        else: assignment = None

        # Debugging
        log.debug("Loading the simulations ...")

        # Get simulations
        simulations = self.generation.simulations_or_basic

        # Get status?
        if self.config.status is not None: status = SimulationStatusTable.from_file(self.config.status)

        # Get status and fix things according to this real-time status
        else:

            # Debugging
            log.debug("Getting the simulation status ...")

            # Get real-time simulation status
            status = self.get_generation_status(remotes)

            # Check?
            if self.config.check_status:

                # Debugging
                log.debug("Checking the simulation status ...")

                # Check the simulations, adapt status if necessary
                status = self.check_simulations(status)

                # Create the simulations
                _add_to_status = []
                _adapted_assignment = False
                for simulation in simulations:

                    #simulation = simulations[simulation_name]
                    #print(simulation)
                    simulation_name = simulation.name
                    #print(simulation.name)

                    # Simulation in status table?
                    if status.has_simulation(simulation_name):
                        simulation_status = status.get_status(simulation_name)
                        if is_analysed_status(simulation_status): simulation.analysed = True
                    else: _add_to_status.append(simulation_name)

                    # Fix assignment scheme
                    if assignment is not None:
                        #print("here")
                        if not assignment.has_simulation(simulation_name):
                            _adapted_assignment = True
                            assignment.add_simulation_object(simulation)
                        if not assignment.has_host_id(simulation_name):
                            _adapted_assignment = True
                            assignment.set_host_id(simulation_name, simulation.host_id)
                        if not assignment.has_id(simulation_name):
                            _adapted_assignment = True
                            assignment.set_id(simulation_name, simulation.id)

                # Add simulations to status table?
                if len(_add_to_status) > 0:
                    status = SimulationStatusTable.from_previous(status, new_simulations=_add_to_status)
                    for simulation_name in _add_to_status:
                        simulation_status, _ = self.generation.get_simulation_status(simulation_name)
                        status.reset_for_simulation(simulation_name, simulation_status)

                # Save the assignment scheme table
                if _adapted_assignment: assignment.save()

        # Set input for base class
        input = dict()
        input["assignment"] = assignment
        input["timing"] = self.timing_table
        input["memory"] = self.memory_table
        input["status"] = status
        input["info_tables"] = [self.parameters_table, self.chi_squared_table]
        input["remotes"] = remotes
        input["simulations"] = simulations

        # Set the output path
        self.config.output = self.generation_path

        # Setup
        SimulationManager.setup(self, **input)

    # -----------------------------------------------------------------

    def check_simulations(self, status):

        """
        This function ...
        :param status:
        :return:
        """

        # Check supposedly analysed simulations (whether this is correct), passing the current status table
        status = self.check_analysed(status)

        # Check paths of simulations
        if self.config.check_paths: self.check_simulation_paths()

        # Check whether simulations with chi squared (analysed simuations) also have the other analysis output,
        # passing the current status table and adapting it if necessary
        if self.config.check_analysis: status = self.check_analysis(status)

        # Return the new status
        return status

    # -----------------------------------------------------------------

    def check_analysed(self, status):

        """
        Thisn function ...
        :param status:
        :return:
        """

        to_fix_status = dict()

        # Loop over the status entries and check whether analysed simulations actually have their chi squared value set
        for simulation_name in status.simulation_names:

            # Get the simulation status
            simulation_status = status.get_status(simulation_name)

            # Check whether chi-squared value is present (full analysis is completed for simulation)
            analysed = self.generation.is_analysed(simulation_name)
            if is_analysed_status(simulation_status) and not analysed:

                # Warning
                log.warning("Simulation '" + simulation_name + "' is supposed to be analysed, but chi squared value is missing form the chi squared table")

                # Correct the status
                if self.config.correct_status:

                    new_simulation_status = self.correct_analysis_status(simulation_name)
                    # Set the status for the simulation in the table
                    # status.set_status(simulation_name, simulation_status) # NO: GET_STATUS METHOD IS MEMOIZED SO WE NEED NEW SIMULATION STATUS TABLE OBJECT
                    to_fix_status[simulation_name] = new_simulation_status

                else: log.warning("Not correcting the status or simulation properties ...")

        # Needs fixing
        if len(to_fix_status) > 0: status = self.create_new_status_table(status, to_fix_status)

        # Return the status (new or not)
        return status

    # -----------------------------------------------------------------

    def correct_analysis_status(self, simulation_name):

        """
        This function ...
        :param simulation_name:
        :return:
        """

        # Inform
        log.warning("Correcting the analysis status for simulation '" + simulation_name + "' ...")

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

        # Fix the simulation properties
        #if not self.generation.has_simulation(simulation_name): continue
        #else:
        if self.generation.has_simulation(simulation_name):

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

        # Return the new simulation satus
        return simulation_status

    # -----------------------------------------------------------------

    def create_new_status_table(self, status, to_fix):

        """
        This function ...
        :param status:
        :param to_fix:
        :return:
        """

        # Get number of simulations for which fixing is necessary
        nfix = len(to_fix)

        # Warning
        log.warning("Simulation status table needs fixing for " + str(nfix) + " simulations")

        # Create lists
        simulation_names = []
        status_list = []

        # Loop over the simulations
        for simulation_name in status.simulation_names:

            # Get the correct status
            if simulation_name in to_fix: correct_status = to_fix[simulation_name]
            else: correct_status = status.get_status(simulation_name)

            # Add to columns
            simulation_names.append(simulation_name)
            status_list.append(correct_status)

        # Create new status table (because status table class is full with lazyproperties and memoized methods)
        new_status = SimulationStatusTable.from_columns(simulation_names, status_list)

        # Return the new status table
        return new_status

    # -----------------------------------------------------------------

    def check_simulation_paths(self):

        """
        This function ...
        :return:
        """

        # Loop over all simulations
        for simulation_name in self.generation.simulation_names:

            # Get simulation object
            if not self.generation.has_simulation(simulation_name): continue  # simulation objects created on the fly are assumed to have correct paths
            simulation = self.generation.get_simulation(simulation_name)

            # Check paths
            try: check_simulation_paths(simulation)
            except RuntimeError as e:
                if self.config.correct_paths:
                    log.warning(str(e))
                    log.warning("Fixing simulation paths ...")
                    correct_simulation_and_analysis_paths(simulation, confirm=self.config.confirm_correction)
                else: raise e

    # -----------------------------------------------------------------

    def check_analysis(self, status):

        """
        Thisn function ...
        :param status:
        :return:
        """

        # Simulation names to have status reset
        reset_status = []

        # Loop over the simulations
        for simulation_name in self.generation.simulation_names:

            # Has simulation object?
            if not self.generation.has_simulation(simulation_name):
                log.warning("Cannot check analysis options for simulation '" + simulation_name + "': simulation file not found")
                continue

            # Load the simulation
            simulation = self.generation.get_simulation(simulation_name)

            # Check the analysis for this simulation
            adapted = self.check_analysis_for_simulation(simulation)
            if adapted: reset_status.append(simulation_name)

        # Reset status
        for simulation_name in reset_status:

            # Debugging
            log.debug("Resetting status for simulation '" + simulation_name + "' ...")

            # Reset
            status.reset_for_simulation(simulation_name, "retrieved")

        # Return the status (although at this point we don't really create a new one, we just adapt it [is this good enough?])
        return status

    # -----------------------------------------------------------------

    def check_analysis_for_simulation(self, simulation, save=True):

        """
        This function ...
        :param simulation:
        :param save:
        :return:
        """

        # Get the simulation name
        simulation_name = simulation.name

        # Debugging
        log.debug("Checking analysis options for simulation '" + simulation_name + "' ...")

        # Initialize flag
        _adapted = False

        # Check different analysis output
        has_misc = self.generation.has_misc_output(simulation_name)
        has_plotting = self.generation.has_plotting_output(simulation_name)
        has_extraction = self.generation.has_extraction_output(simulation_name)
        has_timing = self.generation.has_timing(simulation_name)
        has_memory = self.generation.has_memory(simulation_name)
        has_chi_squared = self.generation.is_analysed(simulation_name)
        #print(has_misc, has_plotting, has_extraction, has_timing, has_memory, has_chi_squared)
        #print(simulation.analysed_any_extraction, simulation.analysed_any_plotting, simulation.analysed_any_misc, simulation.analysed_batch, simulation.analysed_all_extra)

        # Check extraction
        if not has_extraction and simulation.analysed_any_extraction:

            # Show warning
            log.warning("Simulation '" + simulation_name + "' appears to have no extraction output but extraction is supposed to be (partly) executed")

            # Unset analysed flag
            if simulation.analysed: simulation.set_analysed(False)
            simulation.unset_analysed_extraction()
            #reset_status.append(simulation_name)
            _adapted = True

        # Check plotting
        if not has_plotting and simulation.analysed_any_plotting:

            # Show warning
            log.warning("Simulation '" + simulation_name + "' appears to have no plotting output but plotting is supposed to be (partly) executed")

            # Unset analysed flag
            if simulation.analysed: simulation.set_analysed(False)
            simulation.unset_analysed_plotting()
            #reset_status.append(simulation_name)
            _adapted = True

        # Check misc
        if not has_misc and simulation.analysed_any_misc:

            # Show warning
            log.warning("Simulation '" + simulation_name + "' appears to have no miscellaneous output but miscellaneous analysis is supposed to be (partly) executed")

            # Unset analysed flag
            if simulation.analysed: simulation.set_analysed(False)
            simulation.unset_analysed_misc()
            #reset_status.append(simulation_name)
            _adapted = True

        # Check batch
        if not (has_timing or has_memory) and simulation.analysed_batch:

            # Show warning
            log.warning("Simulation '" + simulation_name + "' appears to have no timing or memory entry but batch analysis is supposed to be (partly) executed")

            # Unset analysed flag
            if simulation.analysed: simulation.set_analysed(False)
            simulation.unset_analysed_batch()
            #reset_status.append(simulation_name)
            _adapted = True

        # Check extra
        if not has_chi_squared and simulation.analysed_all_extra:

            # Show warning
            log.warning("Simulation '" + simulation_name + "' appears to have no chi squared value but extra analysis is supposed to be executed")

            # Unset analysed flag
            if simulation.analysed: simulation.set_analysed(False)
            simulation.unset_analysed_extra()
            #reset_status.append(simulation_name)
            _adapted = True

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

        # Save the simulation if necessary
        if save and _adapted: simulation.save()

        # Return the flag
        return _adapted

    # -----------------------------------------------------------------

    def set_caching(self):

        """
        This function ...
        :return:
        """

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

    # -----------------------------------------------------------------

    def show_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Get the maximum chi squared value
        min_chi_squared = self.chi_squared_table.best_chi_squared
        max_chi_squared = self.chi_squared_table.worst_chi_squared
        max_chi_squared_magnitude = numbers.order_of_magnitude(max_chi_squared)

        # Set number of digits for chi squared values
        chisq_ndecimal = 1
        chisq_ndigits = max_chi_squared_magnitude + 1 + chisq_ndecimal

        print(fmt.bold + "Best chi squared: " + fmt.reset + tostr(min_chi_squared, ndigits=chisq_ndigits))
        print(fmt.bold + "Worst chi squared: " + fmt.reset + tostr(max_chi_squared, ndigits=chisq_ndigits))
        print("")

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

    @lazyproperty
    def show_simulations_definition(self):

        """
        This function ...
        :return:
        """

        # Create the definition
        definition = ConfigurationDefinition(write_config=False)

        # Add settings
        definition.add_optional("move", "directory_path", "move the existing simulation files to this directory")
        definition.add_optional("copy", "directory_path", "copy the existing simulation files to this directory")

        # Return
        return definition

    # -----------------------------------------------------------------

    def show_simulations_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the config
        config = self.get_config_from_command(command, self.show_simulations_definition, **kwargs)

        # Show
        self.show_simulations(move_files_to=config.move, copy_files_to=config.copy)

    # -----------------------------------------------------------------

    def show_simulations(self, move_files_to=None, copy_files_to=None):

        """
        This function ...
        :param move_files_to:
        :param copy_files_to:
        :return:
        """

        #print("simulations")

        has_assignment = self.generation.has_assignment_table

        # Keep track of the simulation filepaths
        filepaths = OrderedDict()

        # Print in columns
        with fmt.print_in_columns() as print_row:

            # Show header
            header = [""]
            header.append("Name")
            header.append("Has file")
            header.append("In assignment")
            header.append("Has chi squared")
            header.append("Host ID (assignment)")
            header.append("Host ID (jobscript)")
            header.append("Host ID (logfile)")
            header.append("Cluster name (assignment)")
            header.append("Cluster name (logfile)")
            header.append("Simulation ID (assignment)")

            # Print the header
            print_row(*header, bold=True)

            # Loop over the simulations
            for simulation_name in self.generation.simulation_names:

                # Initialize the row strings
                parts = []
                parts.append("-")
                parts.append(simulation_name)

                # Has simulation file?
                has_simulation = self.generation.has_simulation(simulation_name)
                parts.append(yes_or_no(has_simulation))

                # Is in assignment?
                in_assignment = self.generation.has_assignment_table and self.generation.in_assignment(simulation_name)
                parts.append(yes_or_no(in_assignment))

                # Check whether chi-squared value is present (full analysis is completed for simulation)
                analysed = self.generation.is_analysed(simulation_name)
                parts.append(yes_or_no(analysed))

                # Get the host ID (from assignment table)
                if has_assignment: host_id_assignment = self.generation._get_host_id_from_assignment(simulation_name)
                else: host_id_assignment = None
                #except MissingSimulation: host_id = None
                parts.append(host_id_assignment if host_id_assignment is not None else "--")

                # Get the host ID from job script
                if self.generation.has_jobscript(simulation_name): host_id_job = self.generation._get_host_id_from_job_script(simulation_name)
                else: host_id_job = None
                parts.append(host_id_job if host_id_job is not None else "--")

                # Get the host ID from logfile
                if self.generation.has_logfile(simulation_name): host_id_log = self.generation._get_host_id_from_logfile(simulation_name)
                else: host_id_log = None
                parts.append(host_id_log if host_id_log is not None else "--")

                # Get the cluster name (assignment)
                if has_assignment: cluster_name_assignment = self.generation._get_cluster_name_from_assignment(simulation_name)
                else: cluster_name_assignment = None
                parts.append(cluster_name_assignment if cluster_name_assignment is not None else "--")

                # Get the cluster name from logfile
                if self.generation.has_logfile(simulation_name): cluster_name_log = self.generation._get_cluster_name_from_logfile(simulation_name)
                else: cluster_name_log = None
                parts.append(cluster_name_log if cluster_name_log is not None else "--")

                # Get the simulation ID
                try: simulation_id = self.generation.get_simulation_id(simulation_name)
                except MissingSimulation: simulation_id = None
                parts.append(str(simulation_id) if simulation_id is not None else "--")

                # Determine the color for the simulation
                #color = "green" if has_simulation else "red"
                if has_simulation and in_assignment: color = "green"
                elif has_simulation: color = "yellow" # not in assignment
                elif in_assignment: color = "magenta" # no simulation file
                else: color = "red" # not in assignment nor in assignment

                # Show the row
                print_row(*parts, color=color)

                # Get the simulation filepath
                if has_simulation: filepaths[simulation_name] = self.generation.get_simulation_filepath(simulation_name)

        nsimulations_with_file = len(filepaths)
        nsimulations = self.generation.nsimulations

        # Copy the simulation files?
        if copy_files_to is not None:
            log.debug("Copying the simulation files (" + str(nsimulations_with_file) + " out of " + str(nsimulations) + " simulations) to '" + copy_files_to + "' ...")
            fs.copy_files(filepaths.values(), copy_files_to)

        # Move the simulation files
        if move_files_to is not None:
            log.debug("Moving the simulation files (" + str(nsimulations_with_file) + " out of " + str(nsimulations) + " simulations) to '" + move_files_to + "' ...")
            fs.move_files(filepaths.values(), move_files_to)

    # -----------------------------------------------------------------

    def create_analysis_options(self, config):

        """
        This function ...
        :param config:
        :return:
        """

        # Create the analysis options
        analysis = AnalysisOptions()

        # Add general options
        self.add_general_analysis_options(analysis, config)

        # Add extraction options
        self.add_extraction_analysis_options(analysis, config)

        # Add plotting options
        self.add_plotting_analysis_options(analysis, config)

        # Add misc options
        self.add_misc_analysis_options(analysis, config)

        # Return
        return analysis

    # -----------------------------------------------------------------

    def add_general_analysis_options(self, analysis, config):

        """
        This function ...
        :param analysis:
        :return:
        """

        # Debugging
        log.debug("Setting general analysis options ...")

        # Set the path to the modeling directory
        analysis.modeling_path = self.config.path

        # Set analyser classes
        if self.is_images_modeling: analysis.analyser_paths = ["pts.modeling.fitting.modelanalyser.ImagesFitModelAnalyser"]
        else: analysis.analyser_paths = ["pts.modeling.fitting.modelanalyser.SEDFitModelAnalyser"]

    # -----------------------------------------------------------------

    @property
    def record_timing(self):

        """
        This function ...
        :return:
        """

        return True
        #if self.config.record_timing: return True
        #elif len(self.remote_host_ids) > 0:
        #    log.warning("Record timing will be enabled because remote execution is used")
        #    return True
        #else: return False

    # -----------------------------------------------------------------

    @property
    def record_memory(self):

        """
        This function ...
        :return:
        """

        return True
        #if self.config.record_memory: return True
        #elif len(self.remote_host_ids) > 0:
        #    log.warning("Record memory will be enabled because remote execution is used")
        #    return True
        #else: return False

    # -----------------------------------------------------------------

    @property
    def extract_timeline(self):

        """
        This
        :return:
        """

        return self.record_timing or self.config.extract_timeline

    # -----------------------------------------------------------------

    @property
    def extract_memory(self):

        """
        This function ...
        :return:
        """

        return self.record_memory or self.config.extract_memory

    # -----------------------------------------------------------------

    def add_extraction_analysis_options(self, analysis, config):

        """
        This function ...
        :param analysis:
        :return:
        """

        # Debugging
        log.debug("Setting extraction options ...")

        # Set
        analysis.extraction.progress = config.extract_progress  # extract progress information
        analysis.extraction.timeline = self.extract_timeline  # extract the simulation timeline
        analysis.extraction.memory = self.extract_memory  # extract memory information

    # -----------------------------------------------------------------

    @property
    def observed_sed_paths(self):

        """
        This function ...
        :return:
        """

        if self.is_galaxy_modeling: return [self.observed_sed_path, self.truncated_sed_path]
        elif self.is_sed_modeling: return [self.observed_sed_path]
        else: raise ValueError("Observed SED not defined for modeling types other than 'galaxy' or 'sed'")

    # -----------------------------------------------------------------

    def add_plotting_analysis_options(self, analysis, config):

        """
        This function ...
        :param analysis:
        :return:
        """

        # Debugging
        log.debug("Setting plotting options ...")

        # Set
        analysis.plotting.seds = config.plot_seds  # Plot the output SEDs
        analysis.plotting.reference_seds = self.observed_sed_paths  # the path to the reference SED (for plotting the simulated SED against the reference points)
        analysis.plotting.format = config.plotting_format  # plot format
        analysis.plotting.progress = config.plot_progress
        analysis.plotting.timeline = config.plot_timeline
        analysis.plotting.memory = config.plot_memory

    # -----------------------------------------------------------------

    @property
    def representation(self):
        return self.generation.model_representation

    # -----------------------------------------------------------------

    @property
    def reference_map_path(self):
        return self.representation.reference_map_path

    # -----------------------------------------------------------------

    @property
    def reference_wcs_path(self):
        return self.reference_map_path

    # -----------------------------------------------------------------

    @property
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.fitting_filters

    # -----------------------------------------------------------------

    @property
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.fitting_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def photometry_image_paths_for_fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        # Create new dictionary
        paths = OrderedDict()

        # Loop over the paths
        for filter_name in self.environment.photometry_image_paths_for_filter_names:

            # Skip filters that are not in the fitting filters list
            if filter_name not in self.fitting_filter_names: continue

            # Add the path
            path = self.environment.photometry_image_paths_for_filter_names[filter_name]
            paths[filter_name] = path

        # Return
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_sed_paths(self):

        """
        This function ...
        :return:
        """

        # Initialize dictionary
        paths = OrderedDict()

        # Set the SED paths
        paths["Observed clipped fluxes"] = self.observed_sed_path
        paths["Observed truncated fluxes"] = self.truncated_sed_path

        # Return the dictionary
        return paths

    # -----------------------------------------------------------------

    def add_misc_analysis_options(self, analysis, config):

        """
        This function ...
        :param analysis:
        :param config:
        :return:
        """

        # Debugging
        log.debug("Setting miscellaneous analysis options ...")

        # From images
        if self.generation.use_images:

            # Calculate observed fluxes from images
            analysis.misc.fluxes_from_images = True
            analysis.misc.fluxes = False
            analysis.misc.images = False

            # Set instrument and coordinate system path
            analysis.misc.fluxes_from_images_instrument = self.earth_instrument_name
            analysis.misc.fluxes_from_images_wcs = self.reference_wcs_path

            # Set mask paths
            analysis.misc.fluxes_from_images_masks = self.photometry_image_paths_for_fitting_filter_names
            analysis.misc.fluxes_from_images_mask_from_nans = True

            # Write the fluxes images
            analysis.misc.write_fluxes_images = True

            # Plot fluxes
            analysis.misc.plot_fluxes_from_images = True
            analysis.misc.plot_fluxes_from_images_reference_seds = self.plot_sed_paths  # self.fit_sed_path

            # Set remote for creating images from datacubes
            #analysis.misc.fluxes_from_images_remote = self.other_host_id
            analysis.misc.fluxes_from_images_remote_spectral_convolution = True

            # Make a plot of the images
            analysis.misc.plot_fluxes_images = True

        # From SEDs
        else:

            # Calculate observed fluxes from SEDs
            analysis.misc.fluxes_from_images = False
            analysis.misc.fluxes = True
            analysis.misc.images = False

            # Plot fluxes
            analysis.misc.plot_fluxes = True
            analysis.misc.plot_fluxes_reference_seds = self.plot_sed_paths

        # Observation filters and observation instruments
        analysis.misc.observation_filters = self.fitting_filter_names
        analysis.misc.observation_instruments = [self.earth_instrument_name]

        # Set spectral convolution flag
        #analysis.misc.fluxes_spectral_convolution = self.spectral_convolution
        analysis.misc.fluxes_spectral_convolution = self.generation.spectral_convolution
        #analysis.misc.fluxes_from_images_spectral_convolution = self.spectral_convolution
        analysis.misc.fluxes_from_images_spectral_convolution = self.generation.spectral_convolution

    # -----------------------------------------------------------------

    @lazyproperty
    def set_analysis_definition(self):

        """
        This function ...
        :return:
        """

        # Create definition
        definition = ConfigurationDefinition(write_config=False)

        # Load settings
        definition.import_settings(analysis_definition)
        definition.add_flag("analysed", "set the analysed flag for the simulations", None)

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    def set_analysis_command(self, command, **kwargs):

        """
        This function ...
        :param command:
        :param kwargs:
        :return:
        """

        # Get the simulation names
        simulation_names, config = self.get_simulation_names_and_config_from_command(command, self.set_analysis_definition, **kwargs)

        # Create the analysis options
        analysis = self.create_analysis_options(config)

        # Loop over the simulations
        for simulation_name in simulation_names:

            # Copy the general analysis options
            options = analysis.copy()

            ## Set simulation-specific options

            # Extraction path
            options.extraction.path = self.get_simulation_extraction_path(simulation_name)

            # Plotting path
            options.plotting.path = self.get_simulation_plotting_path(simulation_name)

            # Misc path
            options.misc.path = self.get_simulation_misc_path(simulation_name)

            # Set analysis for this simulation
            self.set_analysis(simulation_name, options, analysed=config.analysed)

    # -----------------------------------------------------------------

    def set_analysis(self, simulation_name, analysis_options, analysed=None):

        """
        This function sets the analysis options for a particular simulation
        :param simulation_name:
        :param analysis_options:
        :param analysed:
        :return:
        """

        # Debugging
        log.debug("Setting analysis options for simulation '" + simulation_name + "' ...")

        # Get the simulation
        simulation = self.get_simulation(simulation_name)

        # Set the analysis
        simulation.set_analysis_options(analysis_options)

        # Check analyse flags
        adapted = self.check_analysis_for_simulation(simulation, save=False)
        if adapted: log.debug("Adapted analysis status of the simulation")

        # Unset the analysed flag
        if analysed is not None: simulation.set_analysed(analysed, all=True)

        # Save the simulation
        simulation.save()

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
