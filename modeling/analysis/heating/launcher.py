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
from ....core.tools.utils import lazyproperty
from ...misc.interface import ModelSimulationInterface
from ...basics.instruments import SimpleInstrument
from ..launcher import AnalysisLauncherBase
from ....core.tools import formatting as fmt
from ....core.tools.stringify import tostr
from ....core.simulation.output import output_types as ot

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

        # The wavelength grid
        self.wavelength_grid = None

        # The instruments
        self.instruments = dict()

        # The scheduling options for the different simulations (if using a remote host with scheduling system)
        self.scheduling_options = dict()

        # The parameter values
        self.parameter_values = None

        # The parallelization scheme (or the number of processes)
        self.parallelization = None
        self.nprocesses = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the model
        self.get_model()

        # 3. Create the wavelength grid
        self.create_wavelength_grid(check_filters=False)

        # 4. Load the deprojections
        self.load_deprojections()

        # 5. Create the projections
        self.create_projections()

        # 6. Create the instruments
        self.create_instruments()

        # 7. Load the ski file
        self.load_ski()

        # 8. Create the ski files for the different contributions
        self.adjust_ski()

        # 9. Set the simulation input paths
        self.set_input_paths()

        # 10. Set the parallelization scheme
        self.set_parallelization()

        # 11. Estimate the runtimes, create the scheduling options
        if self.uses_scheduler: self.estimate_runtimes()

        # 12. Writing
        self.write()

        # 13. Launch the simulations
        self.launch()

    # -----------------------------------------------------------------

    @property
    def heating(self):

        """
        This function ...
        :return:
        """

        return True

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

    @property
    def local_input_paths(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run.heating_input_paths

    # -----------------------------------------------------------------

    @property
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
        self.launcher.config.group_walltime = self.config.walltime  # the preferred walltime for jobs of a group of simulations

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
        self.launcher.config.retrieve_types = self.retrieve_types

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

        # Show the parameter values
        print("")
        print("All model parameter values:")
        print("")
        for label in self.parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.parameter_values[label]))
        print("")

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

    def create_projections(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Create projections
        deprojection_name = self.create_projection_systems(make_faceon=True, make_edgeon=False, use_grid=self.config.use_grid_resolution, reference_name=self.config.resolution_reference)

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

    @property
    def earth_instrument_properties(self):

        """
        This function ...
        :return:
        """

        return dict()  # no specific properties

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

        # 2. Set the instruments
        self.set_instruments()

        # 3. Set writing options
        self.set_writing_options()

        # 4. Create the ski files for the different contributions
        self.create_ski_contributions()

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
        self.ski.set_file_wavelength_grid(fs.name(self.analysis_run.heating_wavelength_grid_path))

    # -----------------------------------------------------------------

    def set_instruments(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Adding the instruments ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

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
        self.ski.set_write_density()
        # self.ski.set_write_depth_map()
        # self.ski.set_write_quality()
        self.ski.set_write_cell_properties()
        # self.ski.set_write_cells_crossed()

        # EXTRA OUTPUT
        if self.config.temperatures: self.ski.set_write_temperature()
        if self.config.emissivities: self.ski.set_write_emissivity()
        if self.config.isrf: self.ski.set_write_isrf()

        # Write absorption
        self.ski.set_write_absorption()

        #self.ski.set_write_grid()

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

            # Remove other stellar components, except for the contribution of the total stellar population
            if contribution != total: ski.remove_stellar_components_except(component_names[contribution])

            # For the simulation with only the ionizing stellar component, also write out the stellar density
            if contribution == ionizing: ski.set_write_stellar_density()

            # Add the ski file instance to the dictionary
            self.ski_contributions[contribution] = ski

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
        #parallelization = self.launcher.parallelization_for_host(self.remote_host_id)
        parallelization = self.parallelization

        # Dictionary of estimated walltimes for the different simulations
        walltimes = dict()

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the config
        self.write_config()

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the instruments
        self.write_instruments()

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

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid ...")

        # Save the wavelength grid
        self.wavelength_grid.to_skirt_input(self.analysis_run.heating_wavelength_grid_path)

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
            simulation_name = self.galaxy_name + "__heating__" + self.analysis_run.name + "__" + contribution

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

        # Set parallelization and nprocesses
        if self.uses_remote:
            parallelization_local = None
            nprocesses_local = None
        else:
            parallelization_local = self.parallelization
            nprocesses_local = self.nprocesses

        # Memory usage
        memory = None

        # Run the launcher, schedules the simulations
        self.launcher.run(memory=memory, ncells=self.ndust_cells, nwavelengths=self.nwavelengths, parallelization_local=parallelization_local, nprocesses_local=nprocesses_local)

# -----------------------------------------------------------------
