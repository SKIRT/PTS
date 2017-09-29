#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.heating.launch Contains the DustHeatingContributionLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import DustHeatingAnalysisComponent, contributions, total, ionizing
from ....core.tools import filesystem as fs
from ....core.launch.batchlauncher import BatchLauncher
from ....core.simulation.definition import SingleSimulationDefinition
from ....core.basics.log import log
from ....core.advanced.runtimeestimator import RuntimeEstimator
from ....core.launch.options import SchedulingOptions
from ....core.advanced.dustgridtool import DustGridTool
from ....core.tools.utils import lazyproperty
from ...misc.interface import ModelSimulationInterface
from ...basics.instruments import SimpleInstrument
from ....core.remote.remote import Remote
from ..launcher import AnalysisLauncherBase

# -----------------------------------------------------------------

component_names = {"old": ["Evolved stellar bulge", "Evolved stellar disk"],
                    "young": "Young stars",
                    "ionizing": "Ionizing stars",
                    "unevolved": ["Young stars", "Ionizing stars"]}

# -----------------------------------------------------------------

class DustHeatingContributionLauncher(DustHeatingAnalysisComponent, ModelSimulationInterface, AnalysisLauncherBase):
    
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
        #super(DustHeatingContributionLauncher, self).__init__(*args, **kwargs)
        DustHeatingAnalysisComponent.__init__(self, no_config=True)
        ModelSimulationInterface.__init__(self, no_config=True)
        AnalysisLauncherBase.__init__(self, *args, **kwargs)

        # -- Attributes --

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski files for the different contributions
        self.ski_contributions = dict()

        # The paths to the simulation input files
        #self.input_paths = None

        # The wavelength grid
        self.wavelength_grid = None

        # The instruments
        self.instruments = dict()

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

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Load the deprojections
        self.load_deprojections()

        # 4. Create the projections
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # 7. Load the ski file
        self.load_ski()

        # 6. Create the ski files for the different contributions
        self.adjust_ski()

        # 6. Set the simulation input paths
        self.set_input_paths()

        # 7. Set the parallelization scheme
        #if self.uses_scheduler: self.set_parallelization()

        # 8. Estimate the runtimes, create the scheduling options
        if self.uses_scheduler: self.estimate_runtimes()

        # 9. Writing
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
        super(DustHeatingContributionLauncher, self).setup(**kwargs)

        # Set options for the batch launcher
        self.set_launcher_options()

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
        self.launcher.config.group_walltime = self.config.walltime  # the preferred walltime for jobs of a group of simulations

        # Set remote host
        if self.config.remote is not None:
            self.launcher.config.remotes = [self.config.remote]         # the remote host on which to run the simulations
            if self.config.cluster_name is not None: self.launcher.set_cluster_for_host(self.config.remote, self.config.cluster_name)

        # Logging options
        self.launcher.config.logging.verbose = True           # verbose logging mode

        # No simulation analysis options

    # -----------------------------------------------------------------

    # def create_wavelength_grid(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Creating the wavelength grid ...")
    #
    #     # Create the emission lines instance
    #     emission_lines = EmissionLines()
    #
    #     # Fixed wavelengths in the grid
    #     fixed = [self.i1_filter.pivot, self.fuv_filter.pivot] # in micron, for normalization of stellar components
    #
    #     # Range in micron
    #     micron_range = RealRange(self.config.wg.range.min.to("micron").value, self.config.wg.range.max.to("micron").value)
    #
    #     # Create and set the grid
    #     self.wavelength_grid = create_one_logarithmic_wavelength_grid(micron_range, self.config.wg.npoints, emission_lines, fixed)

    # -----------------------------------------------------------------

    # def create_instruments(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Creating the instruments ...")
    #
    #     # Debugging
    #     log.debug("Creating a simple instrument for the earth projection ...")
    #
    #     # Create the instrument and add it to the dictionary
    #     self.instruments["earth"] = self.create_instrument("simple", "earth")
    #
    #     # Debugging
    #     log.debug("Creating a frame instrument for the faceon projection ...")
    #
    #     # Create the instrument and add it to the dictionary
    #     self.instruments["faceon"] = self.create_instrument("frame", "faceon")

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

        #return self.launcher.host

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
    def analysis_run_info(self):

        """
        Thisfunction ...
        :return:
        """

        return self.analysis_run.info

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Create projections
        deprojection_name = self.create_projection_systems(make_faceon=True, make_edgeon=False)

        # Set the deprojection name in the analysis info
        self.analysis_run_info.reference_deprojection = deprojection_name

    # -----------------------------------------------------------------

    @property
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        return SimpleInstrument

    # -----------------------------------------------------------------

    # def set_input(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Setting the input paths ...")
    #
    #     # Set the paths to the input maps
    #     self.input_paths = self.input_map_paths
    #
    #     # Set the path to the wavelength grid file
    #     self.input_paths.append(self.analysis_run.heating_wavelength_grid_path)

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
        self.input_paths = self.analysis_run.input_paths

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

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.config.npackages) + " ...")

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Debugging
        log.debug("Setting the wavelength grid file ...")

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid(fs.name(self.analysis_run.heating_wavelength_grid_path))

        # Debugging
        log.debug("Adding the instruments ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(self.instruments[name])

        # Debugging
        log.debug("Setting writing options ...")

        # Set dust system writing options
        self.ski.set_write_convergence()
        self.ski.set_write_density()
        #self.ski.set_write_depth_map()
        #self.ski.set_write_quality()
        self.ski.set_write_cell_properties()
        #self.ski.set_write_cells_crossed()
        #self.ski.set_write_emissivity()
        #self.ski.set_write_temperature()
        #self.ski.set_write_isrf()
        self.ski.set_write_absorption()
        self.ski.set_write_grid()

        # Debugging
        log.debug("Adjusting stellar components for simulating the different contributions ...")

        # Loop over the different contributions, create seperate ski file instance
        for contribution in contributions:

            # Debugging
            log.debug("Adjusting ski file for the contribution of the " + contribution + " stellar population ...")

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # Remove other stellar components, except for the contribution of the total stellar population
            if contribution != total: ski.remove_stellar_components_except(component_names[contribution])

            # For the simulation with only the ionizing stellar component, also write out the stellar density
            if contribution == ionizing: ski.set_write_stellar_density()

            # Add the ski file instance to the dictionary
            self.ski_contributions[contribution] = ski

    # -----------------------------------------------------------------

    # def set_parallelization(self):
    #
    #     """
    #     This function sets the parallelization scheme for those remote hosts used by the batch launcher that use
    #     a scheduling system (the parallelization for the other hosts is left up to the batch launcher and will be
    #     based on the current load of the corresponding system).
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Setting the parallelization scheme for the remote host (" + self.config.remote + ") ...")
    #
    #     # Get the parallelization scheme for this host
    #     parallelization = Parallelization.for_host(self.remote_host, self.config.nnodes, self.config.data_parallel)
    #
    #     # Debugging
    #     log.debug("Parallelization scheme for host " + self.remote_host_id + ": " + str(parallelization))
    #
    #     # Set the parallelization for this host
    #     self.launcher.set_parallelization_for_host(self.remote_host_id, parallelization)

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
        #ncells = self.estimate_ncells()

        # Loop over the different ski files`
        for contribution in self.ski_contributions:

            # Get the ski file
            ski = self.ski_contributions[contribution]

            # Estimate the runtime for the current number of photon packages and the current remote host
            runtime = estimator.runtime_for(ski, parallelization, self.remote_host_id, self.remote_cluster_name, self.config.data_parallel, nwavelengths=len(self.wavelength_grid), ncells=self.ndust_cells)

            # Debugging
            log.debug("The estimated runtime for this host is " + str(runtime) + " seconds")

            # Set the estimated walltime
            walltimes[contribution] = runtime

        # Create and set scheduling options for each host that uses a scheduling system
        for contribution in walltimes: self.scheduling_options[contribution] = SchedulingOptions.from_dict({"walltime": walltimes[contribution]})

    # -----------------------------------------------------------------

    # def estimate_ncells(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Estimating the number of dust cells ...")
    #
    #     # Create simulation directory and output directory
    #     simulation_path = fs.create_directory_in(self.analysis_run.path, "temp_heating")
    #
    #     # Initialize dust grid tool
    #     tool = DustGridTool()
    #
    #     # Get the dust grid statistics
    #     statistics = tool.get_statistics(self.ski_contributions["total"], simulation_path, self.maps_path, self.galaxy_name)
    #     return statistics.ncells

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

        # Write the instruments
        self.write_instruments()

        # Write the ski files
        self.write_ski_files()
        
    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        # Save the wavelength grid
        self.wavelength_grid.saveto(self.analysis_run.heating_wavelength_grid_path)

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the instruments ...")

        # Loop over the instruments
        for name in self.instruments:

            # Determine path
            path = fs.join(self.analysis_run.heating_instruments_path, name + ".instr")

            # Save
            self.instruments[name].saveto(path)

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
            path = self.analysis_run.heating_ski_path_for_contribution(contribution)

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

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulations ...")

        # Set the path for storing the batch scripts for manual inspection
        self.launcher.set_script_path(self.remote_host_id, self.analysis_run.heating_path)

        # Enable screen output for remotes without a scheduling system for jobs
        if not self.uses_scheduler: self.launcher.enable_screen_output(self.remote_host_id)

        # Loop over the contributions
        for contribution in contributions:

            # Determine a name for this simulation
            simulation_name = self.galaxy_name + "_heating_" + self.analysis_run.name + contribution

            # Create the simulation path
            fs.create_directory(self.analysis_run.heating_simulation_path_for_contribution(contribution))

            # Get the ski path for this simulation
            ski_path = self.analysis_run.heating_ski_path_for_contribution(contribution)

            # Get the local output path for the simulation
            output_path = self.analysis_run.heating_output_path_for_contribution(contribution)
            fs.create_directory(output_path)

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, output_path, self.input_paths)

            # Debugging
            log.debug("Adding the simulation of the contribution of the " + contribution + " stellar population to the queue ...")

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name)

            # Set scheduling options for this simulation
            if self.uses_scheduler: self.launcher.set_scheduling_options(self.remote_host_id, simulation_name, self.scheduling_options[contribution])

        # Run the launcher, schedules the simulations
        self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in self.launcher.launched_simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Set the path to an analyser
            # analyser_path = "pts.modeling.analysis. ..."

            # Add the analyser class path
            # simulation.add_analyser(analyser_path)

            # Save the simulation
            simulation.save()

# -----------------------------------------------------------------
