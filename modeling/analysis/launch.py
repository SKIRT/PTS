#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.launch Contains the BestModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools import tables, time
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ...core.tools.logging import log
from ..basics.instruments import FullInstrument
from ...core.launch.options import AnalysisOptions
from ...core.launch.options import SchedulingOptions
from ...core.simulation.arguments import SkirtArguments
from ...core.launch.runtime import RuntimeEstimator
from ...core.launch.parallelization import Parallelization
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.simulation.remote import SkirtRemote
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...magic.misc.kernels import AnianoKernels
from ..basics.projection import GalaxyProjection

# -----------------------------------------------------------------

class BestModelLauncher(AnalysisComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(BestModelLauncher, self).__init__(config)

        # -- Attributes --

        # Create the SKIRT remote execution context
        self.remote = SkirtRemote()

        # The path to the directory with the best model parameters
        self.best_path = None

        # The ski file for the best model
        self.ski = None

        # The wavelength grid
        self.wavelength_grid = None

        # The projection systems
        self.projections = dict()

        # Coordinate system
        self.reference_wcs = None

        # The instruments
        self.instruments = dict()

        # The parallelization scheme
        self.parallelization = None

        # The scheduling options
        self.scheduling_options = None

        # The analysis options
        self.analysis_options = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new BestModelLauncher instance
        launcher = cls()

        # Set the modeling path
        launcher.config.path = arguments.path

        # Return the new instance
        return launcher

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the ski file describing the best model
        self.load_ski()

        # 3. Load the projection systems
        self.load_projections()

        # 4. Create the wavelength grid
        self.create_wavelength_grid()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Adjust the ski file
        self.adjust_ski()

        # 7. Set parallelization
        self.set_parallelization()

        # 8. Estimate the runtime for the simulation
        if self.remote.scheduler: self.estimate_runtime()

        # 9. Set the analysis options
        self.set_analysis_options()

        # 10. Writing
        self.write()

        # 11. Launch the simulation
        self.launch()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(BestModelLauncher, self).setup()

        # Setup the remote execution environment
        self.remote.setup(self.config.remote)

        # The path to the directory with the best model parameters
        self.best_path = fs.join(self.fit_path, "best")

        # Reference coordinate system
        self.reference_wcs = CoordinateSystem.from_file(self.reference_path)

    # -----------------------------------------------------------------

    def load_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the projection systems ...")

        # Load the different projection systems
        for name in ["earth", "faceon", "edgeon"]:

            # Determine the path to the projection file
            path = fs.join(self.components_path, name + ".proj")

            # Load the projection
            projection = GalaxyProjection.from_file(path)

            # Add the projection to the dictionary
            self.projections[name] = projection

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski file for the best fitting model ...")

        # Determine the path to the best model ski file
        path = fs.join(self.best_path, self.galaxy_name + ".ski")

        # Load the ski file
        self.ski = SkiFile(path)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Verify the grid parameters
        if self.config.wavelengths.npoints < 2: raise ValueError(
            "the number of points in the low-resolution grid should be at least 2")
        if self.config.wavelengths.npoints_zoom < 2: raise ValueError(
            "the number of points in the high-resolution subgrid should be at least 2")
        if self.config.wavelengths.min <= 0: raise ValueError("the shortest wavelength should be positive")
        if (self.config.wavelengths.min_zoom <= self.config.wavelengths.min
            or self.config.wavelengths.max_zoom <= self.config.wavelengths.min_zoom
            or self.config.wavelengths.max <= self.config.wavelengths.max_zoom):
            raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

        logmin = np.log10(float(self.config.wavelengths.min))
        logmax = np.log10(float(self.config.wavelengths.max))
        logmin_zoom = np.log10(float(self.config.wavelengths.min_zoom))
        logmax_zoom = np.log10(float(self.config.wavelengths.max_zoom))

        # Build the high- and low-resolution grids independently
        base_grid = np.logspace(logmin, logmax, num=self.config.wavelengths.npoints, endpoint=True, base=10)
        zoom_grid = np.logspace(logmin_zoom, logmax_zoom, num=self.config.wavelengths.npoints_zoom, endpoint=True,
                                base=10)

        # Merge the two grids
        total_grid = []

        # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength < self.config.wavelengths.min_zoom: total_grid.append(wavelength)

        # Add the wavelengths of the high-resolution grid
        for wavelength in zoom_grid: total_grid.append(wavelength)

        # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength > self.config.wavelengths.max_zoom: total_grid.append(wavelength)

        # Create the wavelength grid
        self.wavelength_grid = WavelengthGrid.from_wavelengths(total_grid)

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projection systems
        for name in self.projections:

            # Create the instrument from the projection system
            self.instruments[name] = FullInstrument.from_projection(self.projections[name])

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Enable all writing options for analysis
        #self.ski.enable_all_writing_options()

        # Write out the dust grid data
        self.ski.set_write_grid()

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Determining the parallelization scheme ...")

        # Get the remote host
        host = self.remote.host

        # If the host uses a scheduling system
        if host.scheduler:

            # Get the number of cores per node for this host
            cores_per_node = host.clusters[host.cluster_name].cores

            # Determine the number of cores corresponding to 4 full nodes
            cores = cores_per_node * 4

            # Use 1 core for each process (assume there is enough memory)
            processes = cores

            # Determine the number of threads per core
            if host.use_hyperthreading: threads_per_core = host.clusters[host.cluster_name].threads_per_core
            else: threads_per_core = 1

            # Create a Parallelization instance
            self.parallelization = Parallelization(cores, threads_per_core, processes)

        # If the remote host does not use a scheduling system
        else:

            # Use 4 cores per process
            #cores_per_process = 4

            # Use 10 cores per process
            cores_per_process = 10

            # Get the amount of (currently) free cores on the remote host
            cores = int(self.remote.free_cores)

            # Determine the number of thread to be used per core
            threads_per_core = self.remote.threads_per_core if self.remote.use_hyperthreading else 1

            # Create the parallelization object
            self.parallelization = Parallelization.from_free_cores(cores, cores_per_process, threads_per_core)

        # Debugging
        log.debug("Parallelization scheme that will be used: " + str(self.parallelization))

    # -----------------------------------------------------------------

    def estimate_runtime(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtime for the simulation ...")

        # Debugging
        log.debug("Loading the table with the total runtimes of previous simulations ...")

        # Load the timing table
        timing_table = tables.from_file(self.timing_table_path, format="ascii.ecsv")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator(timing_table)

        # Estimate the runtime for the configured number of photon packages and the configured remote host
        runtime = estimator.runtime_for(self.config.remote, self.config.packages, self.parallelization)

        # Create the scheduling options, set the walltime
        self.scheduling_options = SchedulingOptions()
        self.scheduling_options.walltime = runtime

    # -----------------------------------------------------------------

    def set_analysis_options(self):

        """
        This function ...
        :return:
        """

        # Get the names of the filters for which we have photometry
        filter_names = self.get_observed_filter_names()

        # Set the paths to the for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path("PACS_160")
        for filter_name in filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        # Analysis options
        self.analysis_options = AnalysisOptions()

        # Set options for extraction
        self.analysis_options.extraction.path = self.analysis_extr_path
        self.analysis_options.extraction.progress = True
        self.analysis_options.extraction.timeline = True
        self.analysis_options.extraction.memory = True

        # Set options for plotting
        self.analysis_options.plotting.path = self.analysis_plot_path
        self.analysis_options.plotting.progress = True
        self.analysis_options.plotting.timeline = True
        self.analysis_options.plotting.seds = True
        self.analysis_options.plotting.grids = True
        self.analysis_options.plotting.reference_sed = self.observed_sed_path

        # Set miscellaneous options
        self.analysis_options.misc.path = self.analysis_misc_path
        self.analysis_options.misc.rgb = True
        self.analysis_options.misc.wave = True
        self.analysis_options.misc.fluxes = True
        self.analysis_options.misc.images = True
        self.analysis_options.misc.observation_filters = filter_names
        self.analysis_options.misc.observation_instruments = ["earth"]
        self.analysis_options.misc.make_images_remote = "nancy"
        self.analysis_options.misc.images_wcs = self.reference_path
        self.analysis_options.misc.images_kernels = kernel_paths
        self.analysis_options.misc.images_unit = "MJy/sr"

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

        # Write the input
        self.write_input()

        # Write the ski file
        self.write_ski()

    # -----------------------------------------------------------------

    def write_input(self):

        """
        This function ...
        :return:
        """

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Copy the input map
        self.copy_maps()

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        # Write the wavelength grid
        self.wavelength_grid.to_skirt_input(self.analysis_wavelengths_path)

    # -----------------------------------------------------------------

    def copy_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Copying the input maps ...")

        # Determine the paths to the input maps in the fit/in directory
        fit_in_path = fs.join(self.fit_path, "in")
        old_path = fs.join(fit_in_path, "old_stars.fits")
        young_path = fs.join(fit_in_path, "young_stars.fits")
        ionizing_path = fs.join(fit_in_path, "ionizing_stars.fits")
        dust_path = fs.join(fit_in_path, "dust.fits")

        # Copy the files to the analysis/in directory (if necessary)
        if not fs.has_file(self.analysis_in_path, fs.name(old_path)): fs.copy_file(old_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(young_path)): fs.copy_file(young_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(ionizing_path)): fs.copy_file(ionizing_path, self.analysis_in_path)
        if not fs.has_file(self.analysis_in_path, fs.name(dust_path)): fs.copy_file(dust_path, self.analysis_in_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save the ski file
        self.ski.saveto(self.analysis_ski_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Create the SKIRT arguments object
        arguments = SkirtArguments()

        # Set the arguments
        arguments.ski_pattern = self.analysis_ski_path
        arguments.single = True
        arguments.input_path = self.analysis_in_path
        arguments.output_path = self.analysis_out_path
        arguments.logging.verbose = True
        arguments.logging.memory = True

        # Set the parallelization options
        arguments.parallel.processes = self.parallelization.processes
        arguments.parallel.threads = self.parallelization.threads

        # Debugging: save the screen output in a text file
        remote_skirt_dir_path = self.remote.skirt_dir
        remote_skirt_run_debug_path = fs.join(remote_skirt_dir_path, "run-debug")
        if not self.remote.is_directory(remote_skirt_run_debug_path): self.remote.create_directory(remote_skirt_run_debug_path)
        screen_output_path = fs.join(remote_skirt_run_debug_path, time.unique_name("screen") + ".txt")

        # Save the script file for manual inspection
        host_id = self.config.remote
        scripts_host_path = fs.join(self.analysis_scripts_path, host_id)
        if not fs.is_directory(scripts_host_path): fs.create_directory(scripts_host_path)
        local_script_path = fs.join(scripts_host_path, time.unique_name() + ".sh")

        # Run the simulation
        simulation = self.remote.run(arguments, scheduling_options=self.scheduling_options,
                                     analysis_options=self.analysis_options, screen_output_path=screen_output_path,
                                     local_script_path=local_script_path)

        # Set the retrieve types
        simulation.retrieve_types = ["log", "sed", "image-total"]

        # Save the simulation file
        simulation.save()

# -----------------------------------------------------------------
