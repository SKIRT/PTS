#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.bestmodellauncher Contains the BestModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.launch.batchlauncher import BatchLauncher
from ...core.advanced.runtimeestimator import RuntimeEstimator
from ...core.launch.options import SchedulingOptions
from ...magic.convolution.aniano import AnianoKernels
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from ...core.prep.dustgrids import create_one_dust_grid
from ...core.basics.emissionlines import EmissionLines
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.advanced.dustgridtool import DustGridTool
from ...core.advanced.parallelizationtool import ParallelizationTool
from ...core.basics.configuration import prompt_string

# -----------------------------------------------------------------

contributions = ["total", "old", "young", "ionizing"]
component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                   "young": ["Young stars"],
                   "ionizing": ["Ionizing stars"]}

# -----------------------------------------------------------------

class BestModelLauncher(FittingComponent):
    
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
        super(BestModelLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The fitting run
        self.fitting_run = None

        # The generation name
        self.generation_name = None

        # The ski file template
        self.ski = None

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The analysis options for the simulation of the total stellar contribution
        self.analysis_options_total = dict()

        # The ski files for simulating the contributions of the various stellar components
        self.ski_contributions = dict()

        # The paths to the ski files
        self.ski_paths = dict()

        # The parameter values of the best model
        self.parameter_values = None

        # The wavelength grid and dust grid
        self.wavelength_grid = None
        self.dust_grid = None

        # The path to the directory for the simulations
        self.best_generation_path = None
        self.contributions_simulation_paths = dict()
        self.contributions_output_paths = dict()

        # The path to the wavelength grid file and dust grid file
        self.wavelength_grid_path = None
        self.dust_grid_path = None

        # The paths to the simulation input files
        self.input_paths = None

        # The scheduling options for the different simulations (if using a remote host with scheduling system)
        self.scheduling_options = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load the ski template
        self.load_ski()

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Create the dust grid
        self.create_dust_grid()

        # 4. Set the paths to the input files
        self.set_input()

        # 5. Get the best parameter values
        self.get_parameter_values()

        # 6. Adjust the ski template
        self.adjust_ski()

        # 7. Set the parallelization scheme
        if self.uses_scheduler: self.set_parallelization()

        # 8. Estimate the runtimes, create the scheduling options
        if self.uses_scheduler: self.estimate_runtimes()

        # 9. Write first, then launch the simulations
        self.write()

        # 10. Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BestModelLauncher, self).setup(**kwargs)

        # Load the fitting run
        self.fitting_run = self.load_fitting_run(self.config.fitting_run)

        # Set the default option for the generation name
        last_generation_name = self.fitting_run.last_finished_generation
        if last_generation_name is None: last_generation_name = "first_guess"

        # Set the choices for the generationn name
        generation_names = ["first_guess"] + self.fitting_run.finished_generations

        # Get the generation name
        self.generation_name = prompt_string("generation", "name of the (finished) generation for which to relaunch the "
                                                           "best simulation",
                                             default=last_generation_name,
                                             choices=generation_names)

        # Set options for the batch launcher
        self.set_launcher_options()

        # Set the analysis options for the simulation of the total stellar contribution
        self.set_analysis_options_total()

        # Create a directory for the simulation of the best model of the specified generation
        self.best_generation_path = fs.join(self.fitting_run.best_path, self.generation_name)
        if fs.is_directory(self.best_generation_path): raise RuntimeError("The best model has already been launched for this generation (" + self.generation_name + ")")
        fs.create_directory(self.best_generation_path)

        # Set the paths to the wavelength and dust grid files
        self.wavelength_grid_path = fs.join(self.best_generation_path, "wavelengths.txt")
        self.dust_grid_path = fs.join(self.best_generation_path, "dust_grid.dg")

        # Set the paths to the simulation directories for the different contributions
        for contribution in contributions:

            # Create the directory
            simulation_path = fs.create_directory_in(self.best_generation_path, contribution)

            # Add to the dictionary
            self.contributions_simulation_paths[contribution] = simulation_path

            # Create and set the output directory
            self.contributions_output_paths[contribution] = fs.create_directory_in(simulation_path, "out")

            # Set the ski path
            ski_path = fs.join(simulation_path, self.object_name + ".ski")
            self.ski_paths[contribution] = ski_path

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.name

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Basic options
        self.launcher.config.shared_input = True                    # The input directories (or files) for the different simulations are shared
        self.launcher.config.remotes = [self.config.remote]         # the remote host(s) on which to run the simulations
        self.launcher.config.group_simulations = self.config.group  # group multiple simulations into a single job # TODO: IMPLEMENT THIS
        self.launcher.config.group_walltime = self.config.walltime  # the preferred walltime for jobs of a group of simulations
        #self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        #self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file
        self.launcher.config.cores_per_process = self.config.cores_per_process # the number of cores per process, for non-schedulers
        self.launcher.config.dry = self.config.dry

        # Logging options
        self.launcher.config.logging.verbose = True          # verbose logging

        # Simulation analysis options

        ## General
        self.launcher.config.analysis.relative = True

        ## Extraction
        self.launcher.config.analysis.extraction.path = "extr"     # name for the extraction directory
        self.launcher.config.analysis.extraction.progress = True   # extract progress information
        self.launcher.config.analysis.extraction.timeline = True   # extract the simulation timeline
        self.launcher.config.analysis.extraction.memory = True     # extract memory usage information

        ## Plotting
        self.launcher.config.analysis.plotting.path = "plot"
        self.launcher.config.analysis.plotting.progress = True    # plot progress information
        self.launcher.config.analysis.plotting.timeline = True    # plot timeline
        self.launcher.config.analysis.plotting.memory = True      # plot memory usage
        self.launcher.config.analysis.plotting.seds = True        # plot the simulated SEDs
        self.launcher.config.analysis.plotting.reference_seds = [self.observed_sed_path]  # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.plotting.format = "pdf"     # plot format

        ## Miscellaneous
        self.launcher.config.analysis.misc.path = "misc"          # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        #self.launcher.config.analysis.misc.images = True          # create observed images
        self.launcher.config.analysis.misc.fluxes = True          # calculate observed fluxes
        self.launcher.config.analysis.misc.observation_filters = self.observed_filter_names  # The filters for which to create the observations
        #self.launcher.config.analysis.misc.make_images_remote = self.config.images_remote
        #self.launcher.config.analysis.misc.images_wcs = self.reference_wcs_path
        #self.launcher.config.analysis.misc.images_unit = "MJy/sr"
        #self.launcher.config.analysis.misc.images_kernels = kernel_paths

    # -----------------------------------------------------------------

    def set_analysis_options_total(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the analysis options for the simulation of the total stellar contribution ...")

        # Set the paths to the kernel for each image (except for the SPIRE images)
        kernel_paths = dict()
        aniano = AnianoKernels()
        pacs_red_psf_path = aniano.get_psf_path(self.pacs_red_filter)
        for filter_name in self.observed_filter_names:
            if "SPIRE" in filter_name: continue
            kernel_paths[filter_name] = pacs_red_psf_path

        # Set the analysis options that are different from those used for the simulations of the other stellar contributions
        self.analysis_options_total["misc"] = dict()
        self.analysis_options_total["misc"]["images"] = True
        self.analysis_options_total["misc"]["make_images_remote"] = self.config.images_remote
        self.analysis_options_total["misc"]["images_wcs"] = self.reference_wcs_path
        self.analysis_options_total["misc"]["images_unit"] = "MJy/sr"
        self.analysis_options_total["misc"]["images_kernels"] = kernel_paths

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Create the emission lines instance
        emission_lines = EmissionLines()

        # Fixed wavelengths in the grid
        fixed = [self.i1_filter.pivot, self.fuv_filter.pivot]

        # Create the grid
        # grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added
        grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added = create_one_subgrid_wavelength_grid(self.config.nwavelengths, emission_lines, fixed)

        # Set the grid
        self.wavelength_grid = grid

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        semimajor_angular = self.truncation_ellipse.semimajor  # major axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        x_radius = radius_physical
        y_radius = radius_physical
        z_radius = 3. * Unit("kpc")

        x_min = - x_radius
        x_max = x_radius
        y_min = - y_radius
        y_max = y_radius
        z_min = - z_radius
        z_max = z_radius

        x_extent = x_max - x_min
        y_extent = y_max - y_min
        z_extent = z_max - z_min

        # Set the scale
        scale = self.config.dg.rel_scale * pixelscale

        # Create the grid
        grid = create_one_dust_grid(self.config.dg.grid_type, scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, self.config.dg.min_level, self.config.dg.max_mass_fraction)

        # Set the grid
        self.dust_grid = grid

    # -----------------------------------------------------------------

    def set_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the input paths ...")

        # Set the paths to the input maps
        self.input_paths = self.input_map_paths

        # Determine and set the path to the wavelength grid file
        self.input_paths.append(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    def get_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the parameter values ...")

        # If the first guess model should be used
        if self.generation_name == "first_guess":

            # Inform the user
            log.info("Using the parameter values from the initial guess model ...")

            # Get the parameter values
            self.parameter_values = self.fitting_run.first_guess_parameter_values

        # If the best simulation of a generation has to be used
        else:

            # Inform the user
            log.info("Using the parameter values from the best model in the '" + self.generation_name + "' generation ...")

            # Get the parameter values
            self.parameter_values = self.fitting_run.best_parameter_values_for_generation(self.generation_name)

        # Debugging
        log.debug("The parameter values used for the simulations are: ")
        for label in self.parameter_values: log.debug(" - " + label + ": " + str(self.parameter_values[label]))

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski template ...")

        # 1. Set the best parameter values
        self.set_parameter_values()

        # 2. Set basic properties
        self.set_properties()

        # 3. Adjust the ski files for simulating the contributions of the various stellar components
        self.adjust_ski_contributions()

    # -----------------------------------------------------------------

    def set_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the best parameter values ...")

        # Assert that the parameter values are quantities
        for label in self.parameter_values: assert hasattr(self.parameter_values[label], "unit")

        # Set the parameter values in the ski file template
        self.ski.set_labeled_values(self.parameter_values)

    # -----------------------------------------------------------------

    def set_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting ski file properties ...")

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
        log.debug("Setting the wavelength grid file ...")

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(fs.name(self.wavelength_grid_path))

        # Debugging
        log.debug("Setting the dust grid ...")

        # Set the dust grid
        self.ski.set_dust_grid(self.dust_grid)

    # -----------------------------------------------------------------

    def adjust_ski_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for simulating the contribution of the various stellar components ...")

        # Loop over the different contributions, create seperate ski file instance
        for contribution in contributions:

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # Adjust the ski file for generating simulated images
            if contribution == "total":

                # Debugging
                log.debug("Adjusting ski file for generating simulated images ...")

                # Remove all instruments
                ski.remove_all_instruments()

                # Add the simple instrument
                ski.add_instrument("earth", self.simple_instrument)

            # Remove other stellar components
            else:

                # Debugging
                log.debug("Adjusting ski file for simulating the contribution of the " + contribution + " stellar population ...")
                log.debug("Removing all stellar components other than: " + ", ".join(component_names[contribution]) + " ...")

                # Remove the other components
                ski.remove_stellar_components_except(component_names[contribution])

            # Add the ski file to the dictionary
            self.ski_contributions[contribution] = ski

    # -----------------------------------------------------------------

    def set_parallelization(self):

        """
        This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
        a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
        based on the current load of the corresponding system).
        :return:
        """

        # Inform the user
        log.info("Setting the parallelization scheme for the remote host (" + self.config.remote + ") ...")

        # Create the parallelization tool
        tool = ParallelizationTool()

        # Set configuration options
        tool.config.ski = self.ski
        tool.config.input = self.input_paths

        # Set host properties
        tool.config.nnodes = self.config.nnodes
        tool.config.nsockets = self.remote_host.cluster.sockets_per_node
        tool.config.ncores = self.remote_host.cluster.cores_per_sockets
        tool.config.memory = self.remote_host.cluster.memory

        # MPI available and used
        tool.config.mpi = True
        tool.config.hyperthreading = False
        tool.config.threads_per_core = None

        # Number of dust cells
        tool.config.ncells = None  # number of dust cells (relevant if ski file uses a tree dust grid)

        # Run the parallelization tool
        tool.run()

        # Get the parallelization scheme
        parallelization = tool.parallelization

        # Get the parallelization scheme for this host
        #parallelization = Parallelization.for_host(self.remote_host, self.config.nnodes, self.config.data_parallel)

        # Debugging
        log.debug("Parallelization scheme for host " + self.remote_host_id + ": " + str(parallelization))

        # Set the parallelization for this host
        self.launcher.set_parallelization_for_host(self.remote_host_id, parallelization)

    # -----------------------------------------------------------------

    def estimate_runtimes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the runtimes based on the results of previously finished simulations ...")

        # Create a RuntimeEstimator instance
        estimator = RuntimeEstimator.from_file(self.timing_table_path)

        # Debugging
        log.debug("Estimating the runtime for host '" + self.remote_host_id + "' ...")

        # Get the parallelization scheme that we have defined for this remote host
        parallelization = self.launcher.parallelization_for_host(self.remote_host_id)

        # Dictionary of estimated walltimes for the different simulations
        walltimes = dict()

        # Determine the number of dust cells by building the tree locally
        ncells = self.estimate_ncells()

        # Loop over the different ski files`
        for contribution in self.ski_contributions:

            # Get the ski file
            ski = self.ski_contributions[contribution]

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = estimator.runtime_for(ski, parallelization, self.remote_host_id, self.remote_cluster_name, self.config.data_parallel, nwavelengths=len(self.wavelength_grid), ncells=ncells)

            # Debugging
            log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

            # Set the estimated walltime
            walltimes[contribution] = runtime * 5. # * 5 is temporary

        # Create and set scheduling options for each host that uses a scheduling system
        for contribution in walltimes: self.scheduling_options[contribution] = SchedulingOptions.from_dict({"walltime": walltimes[contribution]})

    # -----------------------------------------------------------------

    def estimate_ncells(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Estimating the number of dust cells ...")

        # Create simulation directory and output directory
        simulation_path = fs.create_directory_in(self.best_generation_path, "temp")

        # Initialize dust grid tool
        tool = DustGridTool()

        # Get the dust grid statistics
        statistics = tool.get_statistics(self.ski_contributions["total"], simulation_path, self.maps_path, self.galaxy_name)
        return statistics.ncells

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the dust grid
        self.write_dust_grid()

        # Write the ski files for simulating the contributions of the various stellar components
        self.write_ski_files()

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        # Write the wavelength grid for SKIRT input
        self.wavelength_grid.to_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid ...")

        # Write the dust grid
        self.dust_grid.saveto(self.dust_grid_path)

    # -----------------------------------------------------------------

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files for simulating the contribution of the various stellar components ...")

        # Loop over the ski files
        for contribution in self.ski_contributions: self.ski_contributions[contribution].saveto(self.ski_paths[contribution])

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

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Set the path to the directory to contain the launch script (job script);
        # just use the directory created for the generation
        self.launcher.set_script_path(self.remote_host_id, self.best_generation_path)

        # Enable screen output logging if using a remote host without a scheduling system for jobs
        if not self.uses_scheduler: self.launcher.enable_screen_output(self.remote_host_id)

        # Loop over the ski paths for the different contributions
        for contribution in self.ski_paths:

            # Debugging
            log.debug("Initiating queuing of simulation for the contribution of the " + contribution + " stellar component ...")

            # Get ski path and output path
            ski_path = self.ski_paths[contribution]
            output_path = self.contributions_output_paths[contribution]

            # Set the simulation name
            simulation_name = self.generation_name.replace(" ", "") + "_" + contribution

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, output_path, self.input_paths)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug(" - name: " + simulation_name)
            log.debug(" - ski path: " + definition.ski_path)
            log.debug(" - output path: " + definition.output_path)

            # Set special analysis options for the simulation of the total stellar contribution
            analysis_options = self.analysis_options_total if contribution == "total" else None

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name, analysis_options=analysis_options)

            # Set scheduling options for this simulation
            if self.uses_scheduler: self.launcher.set_scheduling_options(self.remote_host_id, simulation_name, self.scheduling_options[contribution])

        # Run the launcher, schedules the simulations
        self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in self.launcher.launched_simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Set the path to the fit model analyser
            #analyser_path = "pts.modeling.fitting. ..."

            # Add the analyser class path
            #simulation.add_analyser(analyser_path)

            # Save the simulation object
            simulation.save()

# -----------------------------------------------------------------
