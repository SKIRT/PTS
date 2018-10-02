#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launcher Contains the AnalysisLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.launch.batch import BatchLauncher
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.basics.log import log
from ...core.launch.options import SchedulingOptions
from ..misc.interface import earth_name, faceon_name, edgeon_name
from ...core.tools import formatting as fmt
from ...core.tools.stringify import tostr
from ...core.simulation.output import output_types as ot
from ...core.tools.utils import memoize_method, lazyproperty
from ...core.remote.remote import Remote
from .initialization import wavelengths_filename, dustgridtree_filename
from ...core.launch.options import AnalysisOptions
from .run import old, young, ionizing, unevolved, total, bulge, disk, contributions
from ...core.prep.smile import SKIRTSmileSchema

# -----------------------------------------------------------------

bulge_component_label = "Evolved stellar bulge"
disk_component_label = "Evolved stellar disk"
young_component_label = "Young stars"
ionizing_component_label = "Ionizing stars"

# -----------------------------------------------------------------

# Set stellar components for different contribution simulations
component_names = {bulge: bulge_component_label,
                   disk: disk_component_label,
                   old: [bulge_component_label, disk_component_label],
                   young: young_component_label,
                   ionizing: ionizing_component_label,
                   unevolved: [young_component_label, ionizing_component_label]}

# -----------------------------------------------------------------

class AnalysisLauncher(AnalysisComponent):

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
        AnalysisComponent.__init__(self, *args, **kwargs)

        # -- Attributes --

        # The model definition
        self.definition = None

        # The ski file for the model
        self.ski = None

        # Input file paths
        self.input_paths = None
        self.has_remote_input_files = False
        self.remote_input_path = None

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski files for the different contributions
        self.ski_contributions = dict()

        # The parameter values
        self.parameter_values = None

        # The parallelization scheme (or the number of processes)
        self.parallelization = None
        self.nprocesses = None

        # The SKIRT smile schema creator
        self.smile = SKIRTSmileSchema()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get the model
        self.get_model()

        # Load the ski file
        self.load_ski()

        # Create the ski files for the different contributions
        self.adjust_ski()

        # Set the simulation input paths
        self.set_input_paths()

        # Set the parallelization scheme
        self.set_parallelization()

        # Set analysis options
        self.set_analysis_options()

        # Writing
        self.write()

        # Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.definition.name

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_run(self):

        """
        This function ...
        :return:
        """

        return self.analysis_runs.load(self.config.run)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisLauncher, self).setup(**kwargs)

        # Set options for the batch launcher
        self.set_launcher_options()

    # -----------------------------------------------------------------

    @property
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    @lazyproperty
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.dust_grid

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):

        """
        This fnction ...
        :return:
        """

        return self.analysis_run.dust_grid_tree

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_cells(self):

        """
        This function ...
        :return:
        """

        if self.config.ncells is not None: return self.config.ncells
        else: return self.analysis_run.ncells

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        if not self.uses_remote: return None
        else:
            remote = Remote()
            if not remote.setup(host_id=self.host_id): raise RuntimeError("Could not connect to the remote host '" + self.host_id + "'")
            else: return remote

    # -----------------------------------------------------------------

    @property
    def local_input_paths(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.input_paths

    # -----------------------------------------------------------------

    @lazyproperty
    def retrieve_types(self):

        """
        This function ...
        :return:
        """

        # Initialize list
        types = []

        # Add types
        types.append(ot.logfiles)
        types.append(ot.seds)
        types.append(ot.total_images)
        types.append(ot.count_images)
        types.append(ot.cell_properties)
        types.append(ot.stellar_density)
        types.append(ot.absorption)

        # Add temperature file retrieval
        # if self.config.temperatures: self.launcher.config.retrieve_types.extend([ot.temperature, ot.cell_temperature])

        # Return the types
        return types

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Basic options
        self.launcher.config.shared_input = True                    # The input directories for the different simulations are shared
        self.launcher.config.group_simulations = self.config.group  # group simulations into larger jobs
        self.launcher.config.group_walltime = self.config.group_walltime  # the preferred walltime for jobs of a group of simulations

        # Set remote host
        if self.config.remote is not None:
            self.launcher.config.remotes = [self.config.remote]         # the remote host on which to run the simulations
            if self.config.cluster_name is not None: self.launcher.set_cluster_for_host(self.config.remote, self.config.cluster_name)

        # Logging options
        self.launcher.config.logging.verbose = True           # verbose logging mode

        ## No simulation analysis options ?

        # SET THE MODELING PATH
        self.launcher.config.analysis.modeling_path = self.config.path

        # SET RETRIEVE TYPES (ONLY RELEVANT FOR REMOTE EXECUTION)
        #self.launcher.config.retrieve_types = self.retrieve_types

        # Don't retrieve and analyse after launching
        self.launcher.config.retrieve = False
        self.launcher.config.analyse = False

        # Write the batch launcher output
        self.launcher.config.write = True
        self.launcher.config.path = self.analysis_run_path

        # Keep remote input and output
        self.launcher.config.keep = self.config.keep

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.path

    # -----------------------------------------------------------------

    def get_model(self):

        """
        Thisnfunction ...
        :return:
        """

        # Inform the user
        log.info("Loading the model ...")

        # Load the model definition
        self.definition = self.static_model_suite.get_model_definition(self.analysis_run.model_name)

        # Get the parameter values
        self.parameter_values = self.analysis_run.parameter_values

        # Show the model name
        print("")
        print(fmt.yellow + fmt.underlined + "Model name" + fmt.reset + ": " + self.model_name)
        if self.generation_name is not None: print(fmt.yellow + fmt.underlined + "Generation name" + fmt.reset + ": " + self.generation_name)
        if self.simulation_name is not None: print(fmt.yellow + fmt.underlined + "Simulation name" + fmt.reset + ": " + self.simulation_name)
        if self.chi_squared is not None: print(fmt.yellow + fmt.underlined + "Chi-squared" + fmt.reset + ": " + tostr(self.chi_squared))

        # Show the parameter values
        print("")
        print("All model parameter values:")
        print("")
        for label in self.parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.generation_name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.simulation_name

    # -----------------------------------------------------------------

    @property
    def chi_squared(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.chi_squared

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.fitting_run_name

    # -----------------------------------------------------------------

    @property
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.fitting_run

    # -----------------------------------------------------------------

    @property
    def from_fitting(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.from_fitting

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        Thisn function ...
        :return:
        """

        #return self.launcher.single_host_id
        return self.config.remote

    # -----------------------------------------------------------------

    @property
    def host(self):

        """
        This function ...
        :return:
        """

        if not self.uses_remote: return None
        else: return self.launcher.single_host

    # -----------------------------------------------------------------

    @property
    def cluster_name(self):

        """
        This function ...
        :return:
        """

        if not self.uses_remote: return None
        return self.launcher.get_clustername_for_host(self.host_id)

    # -----------------------------------------------------------------

    @property
    def uses_remote(self):

        """
        Thisf unction ...
        :return:
        """

        return self.host_id is not None

    # -----------------------------------------------------------------

    @property
    def analysis_run_info(self):

        """
        Thisfunction ...
        :return:
        """

        return self.analysis_run.info

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file ...")

        # Load the ski file
        self.ski = self.analysis_run.ski_file

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski files for simulating the different contributions ...")

        # 1. Set basic settings
        self.set_basic()

        # 3. Set writing options
        self.set_writing_options()

        # 4. Create the ski files for the different contributions
        self.create_ski_contributions()

        # Set the instruments
        self.set_instruments()

    # -----------------------------------------------------------------

    def set_basic(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.config.npackages) + " ...")

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Debugging
        log.debug("Enabling dust self-absorption ..." if self.config.selfabsorption else "Disabling dust self-absorption ...")

        # Set dust self-absorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

        # Debugging
        log.debug("Enabling transient heating ..." if self.config.transient_heating else "Disabling transient heating ...")

        # Set transient heating
        if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
        else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Setting the wavelength grid file ...")

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(fs.name(self.analysis_run.wavelength_grid))

    # -----------------------------------------------------------------

    @lazyproperty
    def local_skirt_version(self):
        return introspection.skirt_main_version()

    # -----------------------------------------------------------------

    @property
    def skirt7(self):
        return self.local_skirt_version == 7

    # -----------------------------------------------------------------

    @property
    def skirt8(self):
        return self.local_skirt_version == 8

    # -----------------------------------------------------------------

    def set_writing_options(self):

        """
        Thisfunction ...
        :return:
        """

        # Debugging
        log.debug("Setting writing options ...")

        # Set dust system writing options
        self.ski.set_write_convergence()

        # EXTRA OUTPUT: ONLY for total simulation
        #if self.config.temperatures: self.ski.set_write_temperature()
        #if self.config.emissivities: self.ski.set_write_emissivity()
        #if self.config.isrf: self.ski.set_write_isrf()

        # Write absorption
        if self.skirt7:
            if not self.smile.supports_writing_absorption: raise RuntimeError("Writing absorption luminosities is not supported in your version of SKIRT7")
            self.ski.set_write_absorption()
        elif self.skirt8: self.ski.set_write_isrf()
        else: raise ValueError("Invalid SKIRT version: " + str(self.local_skirt_version))

        # Write spectral absorption
        if self.smile.supports_writing_spectral_absorption: self.ski.set_write_spectral_absorption()
        else: log.warning("Writing absorption spectra is not supported in your version of SKIRT")

        # Write spectral emission
        if self.smile.supports_writing_spectral_emission: self.ski.set_write_spectral_emission()
        else: log.warning("Writing emission spectra is not supported in your version of SKIRT")

    # -----------------------------------------------------------------

    def create_ski_contributions(self):

        """
        Thisn function ...
        :return:
        """

        # Debugging
        log.debug("Adjusting stellar components for simulating the different contributions ...")

        # Loop over the different contributions, create seperate ski file instance
        for contribution in contributions:

            # Debugging
            log.debug("Adjusting ski file for the contribution of the " + contribution + " stellar population ...")

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # WRITE THE GRID FOR THE TOTAL CONTRIBUTION
            # and write the number of cells crossed by photons
            if contribution == total:

                # Grid
                ski.set_write_grid()
                ski.set_write_density()
                ski.set_write_quality()
                ski.set_write_cell_properties()
                ski.set_write_cells_crossed()
                ski.set_write_depth_map()

                # Other
                if self.config.temperatures: ski.set_write_temperature()
                if self.config.emissivities: ski.set_write_emissivity()
                if self.config.isrf: ski.set_write_isrf()

            # Remove other stellar components, except for the contribution of the total stellar population
            if contribution != total: ski.remove_stellar_components_except(component_names[contribution])

            # For the simulation with only the ionizing stellar component, also write out the stellar density
            # NEW: also for bulge, disk, and young (all individual components)
            if contribution == bulge: ski.set_write_stellar_density()
            if contribution == disk: ski.set_write_stellar_density()
            if contribution == young: ski.set_write_stellar_density()
            if contribution == ionizing: ski.set_write_stellar_density()

            # For the bulge and disk simulations, use less photon packages (only SED instruments)
            if contribution == bulge: ski.setpackages(self.config.npackages_seds)
            if contribution == disk: ski.setpackages(self.config.npackages_seds)

            # Add the ski file instance to the dictionary
            self.ski_contributions[contribution] = ski

    # -----------------------------------------------------------------

    @property
    def sed_earth_instrument(self):
        return self.analysis_run.sed_earth_instrument

    # -----------------------------------------------------------------

    @property
    def simple_earth_instrument(self):
        return self.analysis_run.simple_earth_instrument

    # -----------------------------------------------------------------

    @property
    def full_earth_instrument(self):
        return self.analysis_run.full_earth_instrument

    # -----------------------------------------------------------------

    @property
    def simple_faceon_instrument(self):
        return self.analysis_run.simple_faceon_instrument

    # -----------------------------------------------------------------

    @property
    def full_faceon_instrument(self):
        return self.analysis_run.full_faceon_instrument

    # -----------------------------------------------------------------

    @property
    def simple_edgeon_instrument(self):
        return self.analysis_run.simple_edgeon_instrument

    # -----------------------------------------------------------------

    @property
    def full_edgeon_instrument(self):
        return self.analysis_run.full_edgeon_instrument

    # -----------------------------------------------------------------

    def set_instruments(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adding the instruments ...")

        # Loop over the contributions
        for contribution in self.ski_contributions:

            # Get the ski file
            ski = self.ski_contributions[contribution]

            # Remove the existing instruments
            ski.remove_all_instruments()

            # Total, young, ionizing, and evolved contribution: full instruments for projected heating analysis
            if contribution == total or contribution == young or contribution == ionizing or contribution == old:

                # Add instruments
                ski.add_instrument(earth_name, self.full_earth_instrument)
                ski.add_instrument(faceon_name, self.full_faceon_instrument)
                ski.add_instrument(edgeon_name, self.full_edgeon_instrument)

            # Ionizing contribution: full instrument for earth orientation
            #elif contribution == ionizing: ski.add_instrument(earth_name, self.full_earth_instrument)

            # Bulge and disk: #SED instrument
            # simple instrument
            elif (contribution == bulge) or (contribution == disk): ski.add_instrument(earth_name, self.sed_earth_instrument)

            # Other simulations: only simple earth instrument
            else: ski.add_instrument(earth_name, self.simple_earth_instrument)

    # -----------------------------------------------------------------

    @property
    def input_filenames(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.local_input_paths.keys()

    # -----------------------------------------------------------------

    def set_input_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the simulation input paths ...")

        # No remote
        if not self.uses_remote:
            self.set_input_paths_local()

        # Remote execution
        else:
            self.set_input_paths_remote()

    # -----------------------------------------------------------------

    def set_input_paths_local(self):

        """
        Thisf unction ...
        :return:
        """

        # Debugging
        log.debug("Setting input paths for local execution ...")

        # Check config
        if self.config.remote_input is not None or self.config.remote_input_path is not None:
            raise ValueError("Cannot specifiy remote input path(s) if simulation is not launched remotely")

        # Set input paths
        self.input_paths = self.local_input_paths

        # Set things
        self.has_remote_input_files = False
        self.remote_input_path = None

    # -----------------------------------------------------------------

    def set_input_paths_remote(self):

        """
        Thisn function ...
        :return:
        """

        # Debugging
        log.debug("Setting input paths for remote execution ...")

        # No remote files
        if self.config.remote_input is None and self.config.remote_input_path is None:
            self.find_input_paths_remote()

        # Remote files defined in a dictionary
        elif self.config.remote_input is not None:
            self.set_input_paths_remote_from_dictionary()

        # Remote input directory is specified
        elif self.config.remote_input_path is not None:
            self.set_input_paths_remote_from_directory()

        # We shouldn't get here
        else:
            raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    def find_input_paths_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Searching for input files that are already present on the remote host ...")

        # Find remote directory
        #if not self.heating:
        remote_input_path = self.find_remote_directory_with_all_input()
        #else:
        #    remote_input_path = None

        # Set the original paths
        self.input_paths = self.local_input_paths

        # Found
        if remote_input_path is not None:

            # Set things
            self.has_remote_input_files = False
            self.remote_input_path = remote_input_path

        # No complete remote input directory has been found
        else:

            # Find remote files
            remote_input_paths = self.find_remote_input_files()

            # Remote input files have been found
            if len(remote_input_paths) == 0:

                # Replace files that are on the remote
                for filename in remote_input_paths:
                    # Replace by remote path
                    self.input_paths[filename] = remote_input_paths[filename]

                    # Set things
                    self.has_remote_input_files = True
                    self.remote_input_path = None

            # Not found
            else:

                # Set things
                self.has_remote_input_files = False
                self.remote_input_path = None

    # -----------------------------------------------------------------

    def find_remote_directory_with_all_input(self):

        """
        This function ...
        :return:
        """

        # Get local paths
        paths = self.local_input_paths
        filenames = self.input_filenames

        # Get SKIRT input paths from a previous script in this analysis run with the same host
        remote_input_paths = self.analysis_run.get_remote_script_input_paths_for_host(self.host_id)

        # Loop over the each of the remote input directories
        for remote_input_path in remote_input_paths:

            # Check whether the directory still exists
            if not self.remote.is_directory(remote_input_path): continue

            # Check whether all of the local input files are already in this remote directory
            if self.remote.contains_files(remote_input_path, filenames):

                # Check the wavelengths file
                remote_wavelengths_path = fs.join(remote_input_path, wavelengths_filename)
                # local_wavelengths_path = paths[wavelengths_filename]
                remote_nwavelengths = int(self.remote.get_first_line(remote_wavelengths_path))
                local_nwavelengths = self.nwavelengths
                if remote_nwavelengths != local_nwavelengths: continue  # the number of wavelengths is not equal

                # Check passed, return the directory path
                return remote_input_path

        # No complete input directory from this analysis run is found
        return None

    # -----------------------------------------------------------------

    def find_remote_input_files(self):

        """
        This function ...
        :return:
        """

        # Find from this analysis run
        paths = self.find_remote_input_files_this_analysis_run()
        if len(paths) != 0:
            return paths

        # Find from other analysis runs
        else:
            return self.find_remote_input_files_other_analysis_runs()

    # -----------------------------------------------------------------

    def find_remote_input_files_this_analysis_run(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictinoary for the paths
        remote_paths = dict()

        # Get local paths
        # paths = self.local_input_paths
        filenames = self.input_filenames

        # Get SKIRT input paths from a previous script in this analysis run with the same host
        remote_input_paths = self.analysis_run.get_remote_script_input_paths_for_host(self.host_id)

        # Loop over the each of the remote input directories
        for remote_input_path in remote_input_paths:

            # Check whether the directory still exists
            if not self.remote.is_directory(remote_input_path): continue

            # Loop over the filenames
            for filename in filenames:

                # Determine the remote file path
                remote_filepath = fs.join(remote_input_path, filename)

                # If the file (still) exists
                if self.remote.is_file(remote_filepath):

                    # If it is the wavelength file, check the number of wavelengths
                    if filename == wavelengths_filename:
                        #if self.heating:
                        #    continue  # Don't add remote wavelengths file for heating launcher
                        #else:
                        remote_wavelengths_path = remote_filepath
                        remote_nwavelengths = int(self.remote.get_first_line(remote_wavelengths_path))
                        if remote_nwavelengths != self.nwavelengths: continue
                        remote_paths[filename] = remote_filepath

        # Return the dictionary of paths
        return remote_paths

    # -----------------------------------------------------------------

    def find_remote_input_files_other_analysis_runs(self):

        """
        This function ...
        :return:
        """

        # Look in local analysis runs
        paths = self.find_remote_input_files_other_local_analysis_runs()
        if len(paths) != 0:
            return paths

        # Look in cached analysis runs
        else:
            return self.find_remote_input_files_other_cached_analysis_runs()

    # -----------------------------------------------------------------

    def find_remote_input_files_other_local_analysis_runs(self):

        """
        This function ...
        :return:
        """

        remote_paths = dict()

        filenames = self.input_filenames

        # Loop over the analysis runs
        for analysis_run_name in self.analysis_runs.names:

            # Only other analysis runs
            # if analysis_run.name == self.analysis_run.name: continue
            if analysis_run_name == self.analysis_run.name: continue

            # Load the analysis run
            analysis_run = self.analysis_runs.load(analysis_run_name)

            # Get SKIRT input paths from a previous script in this analysis run with the same host
            remote_input_paths = analysis_run.get_remote_script_input_paths_for_host(self.host_id)

            # Loop over the each of the remote input directories
            for remote_input_path in remote_input_paths:

                # Check whether the directory still exists
                if not self.remote.is_directory(remote_input_path): continue

                # Loop over the filenames
                for filename in filenames:

                    # Skip wavelength files (will not be the same from other analysis run)
                    if filename == wavelengths_filename: continue

                    # Skip dust grid tree files (will not be the same from other analysis run)
                    if filename == dustgridtree_filename: continue

                    # Determine the remote file path
                    remote_filepath = fs.join(remote_input_path, filename)

                    # If the file (still) exists, add it
                    if self.remote.is_file(remote_filepath): remote_paths[filename] = remote_filepath

        # Return the dictionary of paths
        return remote_paths

    # -----------------------------------------------------------------

    def find_remote_input_files_other_cached_analysis_runs(self):

        """
        This function ...
        :return:
        """

        remote_paths = dict()

        filenames = self.input_filenames

        # Get the runs
        # runs = self.cached_analysis_runs

        # Loop over the remote hosts
        # for host_id in runs:

        # Only for the host id used for simulations
        # if host_id != self.host_id: continue

        # Check
        # if self.host_id not in runs: return
        # runs_host = runs[self.host_id]

        # None
        if not self.analysis_context.has_cached_for_host(self.host_id): return []

        # Get the runs
        runs = self.analysis_context.get_cached_runs_for_remote(self.host_id)

        # Loop over the analysis runs
        for run_name in runs.names:

            # Load the analysis run
            run = runs.load(run_name)

            # Get SKIRT input paths from a previous script in this analysis run with the same host
            remote_input_paths = run.get_remote_script_input_paths_for_host(self.host_id)

            # Loop over the each of the remote input directories
            for remote_input_path in remote_input_paths:

                # Check whether the directory still exists
                if not self.remote.is_directory(remote_input_path): continue

                # Loop over the filenames
                for filename in filenames:

                    # Skip wavelength files (will not be the same from other analysis run)
                    if filename == wavelengths_filename: continue

                    # Skip dust grid tree files (will not be the same from other analysis run)
                    if filename == dustgridtree_filename: continue

                    # Determine the remote file path
                    remote_filepath = fs.join(remote_input_path, filename)

                    # If the file (still) exists, add it
                    if self.remote.is_file(remote_filepath): remote_paths[filename] = remote_filepath

        # Return the dictionary of paths
        return remote_paths

    # -----------------------------------------------------------------

    def set_input_paths_remote_from_dictionary(self):

        """
        Thisf unction ...
        :return:
        """

        # Debugging
        log.debug("Setting input paths for remote execution from list of remote files ...")

        # Get local paths
        paths = self.local_input_paths

        # Flag
        has_remote_files = False

        # Replace filepaths by remote filepaths if they have been uploaded to the remote already
        for filename in self.config.remote_input:

            # Check whether valid filename
            if filename not in paths: raise ValueError("The filename '" + filename + "' is not one of the input filenames")

            # Set flag
            has_remote_files = True

            # Replace by remote path
            paths[filename] = self.config.remote_input[filename]

        # Set the input paths
        self.input_paths = paths

        # Set things
        self.has_remote_input_files = has_remote_files
        self.remote_input_path = None

    # -----------------------------------------------------------------

    def set_input_paths_remote_from_directory(self):

        """
        Thisf unction ...
        :return:
        """

        # Debugging
        log.debug("Setting input paths for remote execution from a remote directory path ...")

        # Get local paths
        paths = self.local_input_paths

        # Flag
        are_all_remote = True

        # Search for each file in the directory
        for filename in paths:

            # Determine filepath
            filepath = fs.join(self.config.remote_input_path, filename)

            # Check existence on remote
            if self.remote.is_file(filepath):

                # Add the remote path
                paths[filename] = filepath

            # Doesn't exist
            else:
                are_all_remote = False
                log.warning("Remote version of the '" + filename + "' input file is not found: using local file ...")

        # All files are in the remote directory
        if are_all_remote:

            # Set the original paths
            self.input_paths = paths

            # Set things
            self.has_remote_input_files = False
            self.remote_input_path = self.config.remote_input_path

        else:

            # Set paths
            self.input_paths = paths

            # Set things
            self.has_remote_input_files = True
            self.remote_input_path = None

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        # Remote execution
        if self.uses_remote: self.set_parallelization_remote()

        # Local execution
        else: self.set_parallelization_local()

    # -----------------------------------------------------------------

    def set_parallelization_remote(self):

        """
        Thins function ...
        :param self:
        :return:
        """

        # Parallelization is defined
        if self.config.parallelization_remote is not None:

            self.parallelization = self.config.parallelization_remote
            self.nprocesses = None
            self.launcher.config.check_parallelization = False

            # Set the parallelization
            self.launcher.set_parallelization_for_host(self.host_id, self.parallelization)

        # Parallelization is not defined
        else:

            self.parallelization = None
            self.nprocesses = self.config.nprocesses_remote
            self.launcher.config.data_parallel_remote = self.config.data_parallel_remote

            # Set the number of processes
            self.launcher.set_nprocesses_for_host(self.host_id, self.nprocesses)

    # -----------------------------------------------------------------

    def set_parallelization_local(self):

        """
        Thins function ...
        :param self:
        :return:
        """

        # Parallelization is defined
        if self.config.parallelization_local is not None:

            self.parallelization = self.config.parallelization_local
            self.nprocesses = None
            self.launcher.config.check_parallelization = False

        # Parallelization is not defined
        else:

            self.parallelization = None
            self.nprocesses = self.config.nprocesses_local
            self.launcher.config.data_parallel_local = self.config.data_parallel_local

    # -----------------------------------------------------------------

    @lazyproperty
    def analysis_options(self):

        """
        This function ...
        :return:
        """

        return AnalysisOptions()

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the analysis options ...")

        # 1. Set extraction options
        self.set_extraction_options()

        # 2. Set plotting options
        self.set_plotting_options()

        # 3. Set miscellaneous options
        self.set_misc_options()

        # 4. Set other analysis options
        self.set_other_options()

    # -----------------------------------------------------------------

    @lazyproperty
    def extract_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.total_extract_path

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.total_plot_path

    # -----------------------------------------------------------------

    @lazyproperty
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.total_misc_path

    # -----------------------------------------------------------------

    def set_extraction_options(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Setting extraction analysis options ...")

        # Set options for extraction
        self.analysis_options.extraction.path = self.extract_path
        self.analysis_options.extraction.progress = True
        self.analysis_options.extraction.timeline = True
        self.analysis_options.extraction.memory = True

    # -----------------------------------------------------------------

    def set_plotting_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting plotting analysis options ...")

        # Set options for plotting
        self.analysis_options.plotting.path = self.plot_path
        self.analysis_options.plotting.progress = True
        self.analysis_options.plotting.timeline = True
        self.analysis_options.plotting.seds = True
        self.analysis_options.plotting.grids = False # are already plotted for each initialized analysis run
        self.analysis_options.plotting.reference_seds = [self.observed_sed_path, self.truncated_sed_path]

        # Ignore filters for plotting
        self.analysis_options.plotting.ignore_filters = self.ignore_sed_plot_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.reference_wcs

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.reference_map_path

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range(self):

        """
        This function ...
        :return:
        """

        filters = []
        for fltr in self.observed_filters:
            if not self.wavelength_grid.covers(fltr.wavelength):
                log.warning("The '" + str(fltr) + "' filter is not covered by the wavelength range: not making observations for this filter")
                continue
            filters.append(fltr)
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters_in_range_without_iras(self):

        """
        This function ...
        :return:
        """

        return [fltr for fltr in self.observed_filters_in_range if fltr not in self.iras_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_in_range(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters_in_range]

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names_in_range_without_iras(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters_in_range_without_iras]

    # -----------------------------------------------------------------

    def set_misc_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting miscellaneous analysis options ...")

        # Set miscellaneous options
        self.analysis_options.misc.path = self.misc_path

        # Make RGB images
        #self.analysis_options.misc.rgb = True

        # MAKE datacube ANIMATIONS
        #self.analysis_options.misc.animations = True
        #self.analysis_options.misc.write_animation_frames = True

        # Recreate observed fluxes and images
        self.analysis_options.misc.fluxes = True
        self.analysis_options.misc.images = True

        # WRITE INTERMEDIATE RESULTS
        self.analysis_options.misc.write_intermediate_images = True
        self.analysis_options.misc.write_convolution_kernels = True

        # Use spectral convolution
        self.analysis_options.misc.fluxes_spectral_convolution = self.config.spectral_convolution_fluxes
        self.analysis_options.misc.images_spectral_convolution = self.config.spectral_convolution_images

        # NO SPECTRAL CONVOLUTION FOR CREATING THE PLANCK IMAGES
        self.analysis_options.misc.no_images_spectral_convolution_filters = self.planck_filters

        # For these filters and for the earth instrument
        self.analysis_options.misc.observation_filters = self.observed_filter_names_in_range_without_iras  # the filters for which to create the observations (no IRAS)
        self.analysis_options.misc.observation_instruments = [earth_name]

        # DON'T CREATE OBSERVED IMAGES FOR THE PLANCK FILTERS
        self.analysis_options.misc.no_images_filters = self.planck_filters

        # Group the images per instrument (only when more instruments are being converted into images)
        #self.analysis_options.misc.group_images = True

        # Set WCS path for the images
        self.analysis_options.misc.images_wcs = self.reference_wcs_path
        self.analysis_options.misc.wcs_instrument = earth_name

        # Unit for the images
        self.analysis_options.misc.images_unit = self.config.images_unit

        # Make images on remote host
        # No, not necessary anymore
        #self.analysis_options.misc.make_images_remote = self.images_host_id
        #self.analysis_options.misc.rebin_remote_threshold = self.config.rebin_remote_threshold
        #self.analysis_options.misc.convolve_remote_threshold = self.config.convolve_remote_threshold

        # CONVOLUTION
        # Convolution kernels
        #self.analysis_options.misc.images_kernels = kernel_paths
        self.analysis_options.misc.images_psfs_auto = True # automatically determine the PSF for each filter

        # FWHMS
        self.analysis_options.misc.fwhms_dataset = self.dataset_path

        # REBINNING
        #self.analysis_options.misc.rebin_wcs = # dictionary of FITS files per filter?
        self.analysis_options.misc.rebin_dataset = self.dataset_path # much more convenient to specify
        self.analysis_options.misc.rebin_instrument = earth_name

        # Make images remotely
        self.analysis_options.misc.make_images_remote = self.config.images_remote

        # Nprocesses
        #self.analysis_options.misc.images_nprocesses_local = 2
        #self.analysis_options.misc.images_nprocesses_remote = 8
        self.analysis_options.misc.images_nprocesses_local = 1

    # -----------------------------------------------------------------

    def set_other_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting other analysis options ...")

        # Set the paths of the timing and memory table files
        self.analysis_options.timing_table_path = self.timing_table_path
        self.analysis_options.memory_table_path = self.memory_table_path

        # Set the modeling path
        self.analysis_options.modeling_path = self.config.path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the config
        self.write_config()

        # Write the ski files
        self.write_ski_files()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the configuration ...")

        # Write
        self.config.saveto(self.analysis_run.config_path)

    # -----------------------------------------------------------------

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files ...")

        # Loop over the contributions
        for contribution in self.ski_contributions:

            # Determine the path
            path = self.get_ski_path_for_contribution(contribution)

            # Debugging
            log.debug("Writing the ski file for the " + contribution + " stellar population to '" + path + "' ...")

            # Save the ski file
            self.ski_contributions[contribution].saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_host(self):

        """
        This function ...
        :return:
        """

        hosts = self.launcher.hosts
        assert len(hosts) == 1
        return hosts[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_host_id(self):

        """
        This function ...
        :return:
        """

        return self.remote_host.id

    # -----------------------------------------------------------------

    @lazyproperty
    def remote_cluster_name(self):

        """
        This function ...
        :return:
        """

        return self.remote_host.cluster_name

    # -----------------------------------------------------------------

    @lazyproperty
    def uses_scheduler(self):

        """
        This function ...
        :return:
        """

        return self.launcher.uses_schedulers

    # -----------------------------------------------------------------

    @lazyproperty
    def scheduling_options(self):

        """
        This function ...
        :return:
        """

        return SchedulingOptions(**self.config.scheduling)

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        path = self.analysis_run.simulation_path_for_contribution(contribution)
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @memoize_method
    def get_output_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        path = self.analysis_run.output_path_for_contribution(contribution)
        if not fs.is_directory(path): fs.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @memoize_method
    def get_ski_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return self.analysis_run.ski_path_for_contribution(contribution)

    # -----------------------------------------------------------------

    @memoize_method
    def get_simulation_name_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return self.galaxy_name + "__analysis__" + self.analysis_run.name + "__" + contribution

    # -----------------------------------------------------------------

    @memoize_method
    def get_definition_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        # Get the ski path for this simulation
        ski_path = self.get_ski_path_for_contribution(contribution)

        # Get the local output path for the simulation
        output_path = self.get_output_path_for_contribution(contribution)

        # Create the SKIRT simulation definition
        definition = SingleSimulationDefinition(ski_path, output_path, self.input_paths, name=self.get_simulation_name_for_contribution(contribution))

        # Return the definition
        return definition

    # -----------------------------------------------------------------

    @memoize_method
    def get_scheduling_options_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        # Copy the general options
        options = self.scheduling_options.copy()

        # Set job script path
        simulation_path = self.get_simulation_path_for_contribution(contribution)
        options.local_jobscript_path = fs.join(simulation_path, "job.sh")

        # Return the options
        return options

    # -----------------------------------------------------------------

    @lazyproperty
    def parallelization_local(self):

        """
        This function ...
        :return:
        """

        if self.uses_remote: return None
        else: return self.parallelization

    # -----------------------------------------------------------------

    @lazyproperty
    def nprocesses_local(self):

        """
        Thisfunction ...
        :return:
        """

        if self.uses_remote: return None
        else: return self.nprocesses

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Set the path for storing the batch scripts for manual inspection
        self.launcher.set_script_path(self.remote_host_id, self.analysis_run_path)

        # Enable screen output for remotes without a scheduling system for jobs
        if not self.uses_scheduler: self.launcher.enable_screen_output(self.remote_host_id)

        # Loop over the contributions
        for contribution in contributions:

            # Get definition
            definition = self.get_definition_for_contribution(contribution)

            # Debugging
            log.debug("Adding the simulation of the contribution of the " + contribution + " stellar population to the queue ...")

            # Get the simulation name
            simulation_name = self.get_simulation_name_for_contribution(contribution)

            # Get the analysis options for the total simulation
            if contribution == total: analysis_options = self.analysis_options
            else: analysis_options = None

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name, analysis_options=analysis_options)

            # Set scheduling options for this simulation
            if self.uses_scheduler: self.launcher.set_scheduling_options(self.remote_host_id, simulation_name, self.get_scheduling_options_for_contribution(contribution))

        # Memory usage
        memory = None

        # Set remote input path
        if self.remote_input_path is not None: self.launcher.set_remote_input_path_for_host(self.host_id, self.remote_input_path)

        # Run the launcher, schedules the simulations
        self.launcher.run(memory=memory, ncells=self.ndust_cells, nwavelengths=self.nwavelengths,
                          parallelization_local=self.parallelization_local, nprocesses_local=self.nprocesses_local)

# -----------------------------------------------------------------
