#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launch Contains the AnalysisLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import filesystem as fs
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.basics.log import log
from ...core.launch.options import AnalysisOptions
from ...core.launch.options import SchedulingOptions
from ...core.launch.options import LoggingOptions
from ...core.advanced.runtimeestimator import RuntimeEstimator
from ...core.tools.utils import lazyproperty
from ...core.launch.launcher import SKIRTLauncher
from .initialization import wavelengths_filename, dustgridtree_filename
from ...core.simulation.output import output_types as ot
from ..misc.interface import earth_name
from ...core.remote.remote import Remote
from ...core.prep.deploy import Deployer

# -----------------------------------------------------------------

class AnalysisLauncherBase(AnalysisComponent):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        Thisn function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AnalysisLauncherBase, self).__init__(*args, **kwargs)

        # THE ANALYSIS RUN
        self.analysis_run = None

        # The ski file for the model
        self.ski = None

        # Input file paths
        self.input_paths = None
        self.has_remote_input_files = False
        self.remote_input_path = None

    # -----------------------------------------------------------------

    @abstractproperty
    def heating(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def uses_remote(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def local_input_paths(self):

        """
        This function ...
        :return:
        """

        pass

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
        if not self.uses_remote: self.set_input_paths_local()

        # Remote execution
        else: self.set_input_paths_remote()

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
        if self.config.remote_input is None and self.config.remote_input_path is None: self.find_input_paths_remote()

        # Remote files defined in a dictionary
        elif self.config.remote_input is not None: self.set_input_paths_remote_from_dictionary()

        # Remote input directory is specified
        elif self.config.remote_input_path is not None: self.set_input_paths_remote_from_directory()

        # We shouldn't get here
        else: raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    def find_input_paths_remote(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Searching for input files that are already present on the remote host ...")

        # Find remote directory
        if not self.heating: remote_input_path = self.find_remote_directory_with_all_input()
        else: remote_input_path = None

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
                #local_wavelengths_path = paths[wavelengths_filename]
                remote_nwavelengths = int(self.remote.get_first_line(remote_wavelengths_path))
                local_nwavelengths = self.nwavelengths
                if remote_nwavelengths != local_nwavelengths: continue # the number of wavelengths is not equal

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
        if len(paths) != 0: return paths

        # Find from other analysis runs
        else: return self.find_remote_input_files_other_analysis_runs()

    # -----------------------------------------------------------------

    def find_remote_input_files_this_analysis_run(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictinoary for the paths
        remote_paths = dict()

        # Get local paths
        #paths = self.local_input_paths
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
                        if self.heating: continue # Don't add remote wavelengths file for heating launcher
                        else:
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
        if len(paths) != 0: return paths

        # Look in cached analysis runs
        else: return self.find_remote_input_files_other_cached_analysis_runs()

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
            #if analysis_run.name == self.analysis_run.name: continue
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
        runs = self.cached_analysis_runs

        # Loop over the remote hosts
        for host_id in runs:

            # Loop over the analysis runs
            for run_name in runs[host_id].names:

                # Load the analysis run
                run = runs[host_id].load(run_name)

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

    @abstractproperty
    def host_id(self):

        """
        This function ...
        :return:
        """

        pass

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

class AnalysisLauncher(AnalysisLauncherBase):
    
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
        super(AnalysisLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The remote SKIRT environment
        self.launcher = SKIRTLauncher()

        # The parallelization scheme (or the number of processes)
        self.parallelization = None
        self.nprocesses = None

        # The scheduling options
        self.scheduling_options = None

        # The analysis options
        self.analysis_options = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski file
        self.load_ski()

        # 3. Adjust the ski file
        self.adjust_ski()

        # 4. Set the input paths
        self.set_input_paths()

        # 5. Set the parallelization scheme
        self.set_parallelization()

        # 6. Estimate the runtime for the simulation
        if self.uses_scheduler: self.estimate_runtime()

        # 7. Set the analysis options
        self.set_analysis_options()

        # 8. Writing
        self.write()

        # 9. Launch the simulation
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisLauncher, self).setup(**kwargs)

        # NEW: GET THE RUN
        self.analysis_run = self.get_run(self.config.run)

        # Set remote (and cluster name)
        self.launcher.config.remote = self.config.remote
        self.launcher.config.cluster_name = self.config.cluster_name
        self.launcher.config.attached = self.config.attached
        self.launcher.config.keep = self.config.keep_remote_input_and_output
        self.launcher.config.keep_input = self.config.keep_remote_input

        # Clear remotes
        if self.has_any_host_id and self.config.clear_remotes: self.clear_all_hosts()

        # Deploy SKIRT and PTS
        if self.has_any_host_id and self.config.deploy: self.deploy()

    # -----------------------------------------------------------------

    @property
    def heating(self):

        """
        Thins function ...
        :return:
        """

        return False

    # -----------------------------------------------------------------

    @property
    def local_input_paths(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.input_paths

    # -----------------------------------------------------------------

    def clear_all_hosts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the hosts ...")

        # Setup the remote
        if self.uses_remote:

            remote = Remote()
            if not remote.setup(self.host_id): raise RuntimeError("Could not connect to remote host '" + self.host_id + "'")
            else:

                # Debugging
                log.debug("Claring the " + self.host_id + " remote ...")

                # Clear temporary data and close sessions
                remote.clear_temp_and_sessions()

        # Setup the images remote
        if self.uses_images_remote:

            remote = Remote()
            if not remote.setup(self.images_host_id): log.warning("Could not connect")
            else:

                # Debugging
                log.debug("Clearing the " + self.images_host_id + " remote ...")

                # Clear temporary data and close sessions
                remote.clear_temp_and_sessions()

    # -----------------------------------------------------------------

    def deploy(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Deploying SKIRT and PTS where necessary ...")

        # Create the deployer
        deployer = Deployer()

        # Don't do anything locally
        deployer.config.local = False

        # Set the host ids
        deployer.config.host_ids = self.all_host_ids

        # Set the host id on which PTS should be installed (on the host for extra computations and the fitting hosts
        # that have a scheduling system to launch the pts run_queue command)
        if self.uses_images_remote: deployer.config.pts_on = self.images_host_id #self.moderator.all_host_ids

        # Check versions between local and remote
        deployer.config.check = self.config.check_versions

        # Update PTS dependencies
        deployer.config.update_dependencies = self.config.update_dependencies

        # Do clean install
        deployer.config.clean = self.config.deploy_clean

        # Pubkey pass
        deployer.config.pubkey_password = self.config.pubkey_password

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    @property
    def all_host_ids(self):

        """
        This function ...
        :return:
        """

        host_ids = []
        if self.images_host_id is not None: host_ids.append(self.images_host_id)
        if self.host_id is not None: host_ids.append(self.host_id)
        return host_ids

    # -----------------------------------------------------------------

    @property
    def has_any_host_id(self):

        """
        This function ...
        :return:
        """

        return self.images_host_id is not None or self.host_id is not None

    # -----------------------------------------------------------------

    @property
    def uses_any_remote(self):

        """
        This function ...
        :return:
        """

        return self.has_any_host_id

    # -----------------------------------------------------------------

    @property
    def images_host_id(self):

        """
        This function ...
        :return:
        """

        return self.config.images_remote

    # -----------------------------------------------------------------

    @property
    def host_id(self):

        """
        Thisf unction ...
        :return:
        """

        return self.launcher.host_id

    # -----------------------------------------------------------------

    @property
    def host(self):

        """
        This function ...
        :return:
        """

        return self.launcher.host

    # -----------------------------------------------------------------

    @property
    def cluster_name(self):

        """
        This function ...
        :return:
        """

        return self.launcher.cluster_name

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
    def uses_remote(self):

        """
        Thisf unction ...
        :return:
        """

        return self.host_id is not None

    # -----------------------------------------------------------------

    @property
    def uses_images_remote(self):

        """
        Thisf unction ...
        :return:
        """

        return self.images_host_id is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def uses_scheduler(self):

        """
        This function ...
        :return:
        """

        return self.launcher.uses_scheduler

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.name

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.path

    # -----------------------------------------------------------------

    @property
    def run_output_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.output_path

    # -----------------------------------------------------------------

    @property
    def run_instruments_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.instruments_path

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_file_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.ski_file_path

    # -----------------------------------------------------------------

    @lazyproperty
    def extract_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.extract_path

    # -----------------------------------------------------------------

    @lazyproperty
    def plot_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.plot_path

    # -----------------------------------------------------------------

    @lazyproperty
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.misc_path

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file template ...")

        # Load the ski file template
        self.ski = self.analysis_run.ski_file

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Debugging
        #log.debug("Disabling all writing settings ...")

        # Disable all writing settings
        #self.ski.disable_all_writing_options()

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.config.npackages) + " ...")

        # Set the number of photon packages per wavelength
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

        # Debugging
        log.debug("Setting wavelength grid file ...")

        # Set wavelength grid for ski file
        self.ski.set_file_wavelength_grid(wavelengths_filename)

        # Debugging
        log.debug("Setting file tree dust grid ...")

        # Set dust grid tree file
        self.ski.set_filetree_dust_grid(dustgridtree_filename, write_grid=False)

        # Debugging
        log.debug("Enabling specific writing settings ...")

        # Write temperature data
        if self.config.temperatures: self.ski.set_write_temperature()

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme ...")

        # Set options for parallelization and number of processes
        # Remote execution
        if self.uses_remote:

            # Parallelization is defined
            if self.config.parallelization_remote is not None:

                self.parallelization = self.config.parallelization_remote
                self.nprocesses = None
                self.launcher.config.check_parallelization = False

            # Parallelization is not defined
            else:

                self.parallelization = None
                self.nprocesses = self.config.nprocesses_remote
                self.launcher.config.data_parallel_remote = self.config.data_parallel_remote

        # Local execution
        else:

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

    # def set_parallelization(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Determining the parallelization scheme ...")
    #
    #     # If the host uses a scheduling system
    #     if self.remote.scheduler:
    #
    #         # Debugging
    #         log.debug("Remote host (" + self.remote.host_id + ") uses a scheduling system; determining parallelization scheme based on the requested number of nodes (" + str(self.config.nnodes) + ") ...")
    #
    #         # Create the parallelization scheme from the host configuration and the requested number of nodes
    #         self.parallelization = Parallelization.for_host(self.remote.host, self.config.nnodes, self.config.data_parallel)
    #
    #     # If the remote host does not use a scheduling system
    #     else:
    #
    #         # Debugging
    #         log.debug("Remote host (" + self.remote.host_id + ") does not use a scheduling system; determining parallelization scheme based on the current load of the system and the requested number of cores per process (" + str(self.config.cores_per_process) + ") ...")
    #
    #         # Get the amount of (currently) free cores on the remote host
    #         cores = int(self.remote.free_cores)
    #
    #         # Determine the number of thread to be used per core
    #         threads_per_core = self.remote.threads_per_core if self.remote.use_hyperthreading else 1
    #
    #         # Create the parallelization object
    #         self.parallelization = Parallelization.from_free_cores(cores, self.config.cores_per_process, threads_per_core, self.config.data_parallel)
    #
    #     # Debugging
    #     log.debug("Parallelization scheme that will be used: " + str(self.parallelization))

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        Thisf unction ...
        :return:
        """

        return self.analysis_run.wavelength_grid

    # -----------------------------------------------------------------

    def estimate_runtime(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtime for the simulation based on the timing of previous simulations ...")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator(self.timing_table)

        # Estimate the runtime for the configured number of photon packages and the configured remote host
        runtime = estimator.runtime_for(self.ski, self.parallelization, self.host_id, self.cluster_name, self.config.data_parallel, nwavelengths=self.nwavelengths, ncells=self.ndust_cells)

        # Debugging
        log.debug("The estimated runtime for the simulation is " + str(runtime) + " seconds")

        # Create the scheduling options, set the walltime
        self.scheduling_options = SchedulingOptions()
        self.scheduling_options.walltime = runtime

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the analysis options ...")

        # Initialize the analysis options
        self.analysis_options = AnalysisOptions()

        # 1. Set extraction options
        self.set_extraction_options()

        # 2. Set plotting options
        self.set_plotting_options()

        # 3. Set miscellaneous options
        self.set_misc_options()

        # 4. Set other analysis options
        self.set_other_options()

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
        self.analysis_options.plotting.reference_seds = [self.observed_sed_path]

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

        return [str(fltr) for fltr in self.observed_filter_names_in_range_without_iras]

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
        self.analysis_options.misc.rgb = True

        # Wave movies
        self.analysis_options.misc.wave = False

        # Recreate observed fluxes and images
        self.analysis_options.misc.fluxes = True
        self.analysis_options.misc.images = True

        # Use spectral convolution
        self.analysis_options.misc.fluxes_spectral_convolution = self.config.spectral_convolution_fluxes
        self.analysis_options.misc.images_spectral_convolution = self.config.spectral_convolution_images

        # For these filters and for the earth instrument
        self.analysis_options.misc.observation_filters = self.observed_filter_names_in_range_without_iras  # the filters for which to create the observations (no IRAS)
        self.analysis_options.misc.observation_instruments = [earth_name]

        # Group the images per instrument (only when more instruments are being converted into images)
        #self.analysis_options.misc.group_images = True

        # Set WCS path for the images
        self.analysis_options.misc.images_wcs = self.reference_wcs_path
        self.analysis_options.misc.wcs_instrument = earth_name

        # Unit for the images
        self.analysis_options.misc.images_unit = self.config.images_unit

        # Make images on remote host
        self.analysis_options.misc.make_images_remote = self.images_host_id
        self.analysis_options.misc.rebin_remote_threshold = self.config.rebin_remote_threshold
        self.analysis_options.misc.convolve_remote_threshold = self.config.convolve_remote_threshold

        # CONVOLUTION
        # Convolution kernels
        #self.analysis_options.misc.images_kernels = kernel_paths
        self.analysis_options.misc.images_psfs_auto = True # automatically determine the PSF for each filter

        # REBINNING
        #self.analysis_options.misc.rebin_wcs = # dictionary of FITS files per filter?
        self.analysis_options.misc.rebin_dataset = self.dataset_path # much more convenient to specify
        self.analysis_options.misc.rebin_instrument = earth_name

        # Make images remotely
        self.analysis_options.misc.make_images_remote = self.config.images_remote

        # Nprocesses
        self.analysis_options.misc.images_nprocesses_local = 2
        self.analysis_options.misc.images_nprocesses_remote = 8

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

        # Write the ski file
        self.write_ski()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the configuration used to create this analysis run ...")

        # Write
        self.config.saveto(self.analysis_run.launch_config_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.ski_file_path + "...")

        # Save the ski file
        self.ski.saveto(self.ski_file_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation ...")

        # Set options
        self.launcher.config.show_progress = True
        self.launcher.config.debug_output = self.config.debug_output

        # Create the simulation definition
        definition = SingleSimulationDefinition(self.ski_file_path, self.run_output_path, self.input_paths)

        # Create the logging options
        logging = LoggingOptions(verbose=True, memory=True)

        # Debugging: save the screen output in a text file (for remote execution)
        if self.uses_remote:
            # Determine path, relative to the remote SKIRT directory
            #screen_output_path = fs.join("$SKIRT", "run-debug", self.analysis_run_name + ".txt") # specifying file path is (currently) not possible
            screen_output_path = fs.join("$SKIRT", "run-debug", self.analysis_run_name)
            # Debugging message
            log.debug("Remote simulation output will be written to '" + screen_output_path + "'")
        else: screen_output_path = None

        # Determine the path to the launching script file for manual inspection (for remote execution)
        if self.uses_remote:
            # Determine the path
            local_script_path = fs.join(self.analysis_run_path, self.host_id + ".sh")
            # Debugging
            log.debug("The launching script will be saved locally to '" + local_script_path + "'")
        else: local_script_path = None

        # Set retrieve types (only relevant for remote execution)
        self.launcher.config.retrieve_types = [ot.logfiles, ot.seds, ot.total_images]

        # Add more retrieve types
        if self.config.retrieve_contributions:
            self.launcher.config.retrieve_types.extend([ot.direct_images, ot.transparent_images, ot.scattered_images, ot.dust_images, ot.dust_scattered_images])

        # Add temperature file retrieval
        if self.config.temperatures: self.launcher.config.retrieve_types.extend([ot.temperature, ot.cell_temperature])

        # Other settings
        if log.is_debug(): self.launcher.config.show = True

        # Show the number of dust cells
        log.debug("The number of dust cells is " + str(self.ndust_cells) + "")

        # Debugging
        log.debug("Starting the SKIRT launcher ...")

        # Run the simulation
        self.launcher.run(definition=definition, logging_options=logging, analysis_options=self.analysis_options,
                          scheduling_options=self.scheduling_options, parallelization=self.parallelization,
                          nprocesses=self.nprocesses, local_script_path=local_script_path,
                          screen_output_path=screen_output_path, ncells=self.ndust_cells, remote=self.remote,  # pass remote because possibly already used here (don't connect again in the SKIRTLauncher)
                          remote_input_path=self.remote_input_path, has_remote_input_files=self.has_remote_input_files)

# -----------------------------------------------------------------
