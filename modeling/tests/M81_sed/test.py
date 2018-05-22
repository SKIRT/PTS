#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.modeling.basics.instruments import SEDInstrument
from pts.modeling.tests.base import M81TestBase, fitting_filter_names
from pts.core.launch.options import AnalysisOptions
from pts.do.commandline import Command
from pts.core.simulation.skifile import SkiFile
from pts.modeling.tests.base import instrument_name, free_parameters_absolute_paths
from pts.modeling.tests.base import free_parameters_relative_stellar_component_paths, free_parameters_relative_dust_component_paths
from pts.modeling.tests.base import free_parameters_relative_instruments_paths, free_parameter_types, parameter_ndigits
from pts.core.data.sed import ObservedSED
from pts.core.tools.stringify import tostr
from pts.core.units.parsing import parse_angle
from pts.modeling.tests.base import seds_path, dustpedia_sed_path
from pts.core.tools import sequences
from pts.core.basics.map import Map
from pts.modeling.tests.base import free_parameter_descriptions, free_parameter_units
from pts.modeling.core.environment import SEDModelingEnvironment
from pts.core.tools import introspection
from pts.core.simulation.discover import find_one_simulation_in_path

# -----------------------------------------------------------------

description = "fitting the parameters of a model of M81 based on only a mock observed SED"

# -----------------------------------------------------------------

class M81SEDTest(M81TestBase):

    """
    This class runs the test on M81, but by only adjusting the normalizations (not by creating a model),
    and fitting to a mock observed SED
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor
        :param kwargs:
        """

        # Call the constructor of the base class
        super(M81SEDTest, self).__init__(*args, **kwargs)

        # The galaxy properties
        self.properties = None

        # Bulge and disk model
        self.bulge = None
        self.disk = None

        # The input maps
        self.maps = dict()

        # The instrument
        self.instrument = None

        # The ski template
        self.ski_template = None

        # The observed SED
        self.observed_sed = None

        # The initial parameter values for the fitting
        self.initial_parameter_values = dict()

    # -----------------------------------------------------------------

    @property
    def from_existing_reference(self):

        """
        This function ...
        :return: 
        """

        return self.config.reference_path is not None or self.config.reference_test is not None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the properties
        self.load_properties()

        # 3. Reference simulation
        self.set_reference()

        # 4. Get the real parameter values
        self.get_real_parameter_values()

        # 5. Generate the initial parameter values
        self.generate_initial_parameter_values()

        # 6. Create the ski file template
        self.create_template()

        # 7. Load the observed SED
        self.load_observed_sed()

        # 8. Setup the modeling
        self.setup_modelling()

        # 9. Model
        self.model()

        # 10. Get best parameter values
        self.get_best_parameter_values()

        # 11. Test
        self.test()

        # 12. Plot
        self.plot()

    # -----------------------------------------------------------------

    @property
    def required_input_files(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.simulation_input_path, exact_not_name="wavelengths")

    # -----------------------------------------------------------------

    @property
    def all_input_files(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.simulation_input_path)

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the SED instrument ...")

        # Create the SED instrument
        azimuth = parse_angle("0. deg")
        self.instrument = SEDInstrument(distance=self.galaxy_distance, inclination=self.galaxy_inclination, azimuth=azimuth, position_angle=self.galaxy_position_angle)

    # -----------------------------------------------------------------

    def set_reference(self):

        """
        This function ....
        :return: 
        """

        # Inform the user
        log.info("Setting the reference simulation ...")

        # Load or create
        if self.from_existing_reference: self.load_reference()
        else: self.create_reference()

    # -----------------------------------------------------------------

    def load_reference(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Loading the reference simulation ...")

        # Determine simulation directory
        if self.config.reference_path is not None: simulation_path = self.config.reference_path
        elif self.config.reference_test is not None: simulation_path = fs.join(introspection.pts_tests_dir, self.config.reference_test, "ref")
        else: raise ValueError("Reference path and reference test settings are None")

        # Check whether present
        if not fs.is_directory(simulation_path): raise ValueError("The reference simulation path could not be found")

        # Look for simulation
        prefix, ski_path, in_path, out_path = find_one_simulation_in_path(simulation_path)

        # Load the ski file
        self.ski = SkiFile(ski_path)

        # Other paths
        extr_path = fs.join(simulation_path, "extr")
        plot_path = fs.join(simulation_path, "plot")
        misc_path = fs.join(simulation_path, "misc")
        
        # Check existence
        if not fs.is_directory(extr_path): raise IOError("Extraction directory not found")
        if not fs.is_directory(plot_path): raise IOError("Plotting directory not found")
        if not fs.is_directory(misc_path): raise IOError("Misc directory not found")

        # Copy
        self.copy_reference(ski_path, in_path, out_path, extr_path, plot_path, misc_path)

    # -----------------------------------------------------------------

    def copy_reference(self, ski_path, in_path, out_path, extr_path, plot_path, misc_path):

        """
        This function ...
        :param ski_path:
        :param in_path:
        :param out_path:
        :param extr_path:
        :param plot_path:
        :param misc_path:
        :return: 
        """

        # Copy to new reference path
        fs.copy_file(ski_path, self.reference_ski_path)
        if in_path is not None: fs.copy_from_directory(in_path, self.simulation_input_path)
        fs.copy_from_directory(out_path, self.simulation_output_path)
        fs.copy_from_directory(extr_path, self.simulation_extract_path)
        fs.copy_from_directory(plot_path, self.simulation_plot_path)
        fs.copy_from_directory(misc_path, self.simulation_misc_path)

    # -----------------------------------------------------------------

    def create_reference(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Creating the reference simulation ...")

        # 1. Load the components
        self.load_components()

        # 2. Load the input maps
        self.load_maps()

        # 3. Create instrument
        self.create_instrument()

        # 4. Create deprojection
        self.create_deprojections()

        # 5. Create the wavelength grid
        self.create_wavelength_grid()

        # 6. Create the dust grid
        self.create_dust_grid()

        # 7. Create the ski file
        self.create_ski()

        # 8. Write
        self.write_reference()

        # 9. Plot
        self.plot_reference()

        # 10. Launch the reference simulation
        self.launch_reference()

    # -----------------------------------------------------------------

    def write_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski()

        # Write the input
        self.write_input()

    # -----------------------------------------------------------------

    def plot_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the wavelengths
        self.plot_wavelengths()

        # Plot the filters
        self.plot_filters()

    # -----------------------------------------------------------------

    def launch_reference(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the reference simulation ...")

        # Settings
        settings_launch = dict()
        settings_launch["ski"] = self.reference_ski_path
        settings_launch["input"] = self.simulation_input_path
        settings_launch["output"] = self.simulation_output_path
        settings_launch["create_output"] = True
        settings_launch["remote"] = self.moderator.host_id_for_single("reference")
        settings_launch["attached"] = self.config.attached
        settings_launch["progress_bar"] = True

        # Create the analysis options
        analysis = AnalysisOptions()
        analysis.extraction.path = self.simulation_extract_path
        analysis.plotting.path = self.simulation_plot_path
        analysis.misc.path = self.simulation_misc_path
        analysis.extraction.progress = True
        analysis.extraction.timeline = True
        analysis.extraction.memory = True
        analysis.plotting.progress = True
        analysis.plotting.timeline = True
        analysis.plotting.memory = True
        analysis.plotting.seds = True
        analysis.plotting.grids = True
        analysis.plotting.reference_seds = fs.files_in_path(seds_path)
        analysis.misc.fluxes = True
        analysis.misc.images = False
        analysis.misc.observation_filters = fitting_filter_names
        analysis.misc.observation_instruments = [instrument_name]
        analysis.misc.spectral_convolution = self.config.spectral_convolution

        # Set flux error bars
        dustpedia_sed = ObservedSED.from_file(dustpedia_sed_path)
        filter_names = dustpedia_sed.filter_names()
        errors = dustpedia_sed.errors()
        flux_errors = sequences.zip_into_dict(filter_names, [str(error) for error in errors])
        analysis.misc.flux_errors = flux_errors

        # Input
        input_launch = dict()
        # input_launch["memory"] = MemoryRequirement(serial_memory, parallel_memory)
        input_launch["analysis_options"] = analysis

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        self.launcher = self.run_command(launch)

    # -----------------------------------------------------------------

    def generate_initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating random initial parameter values ...")

        # Get the low and high value
        low_factor = self.config.relative_range_initial.min
        high_factor = self.config.relative_range_initial.max

        # Determine the exponents, to generate random points
        log_low = np.log10(low_factor)
        log_high = np.log10(high_factor)

        # Loop over the real parameter values
        for parameter_name in self.real_parameter_values:

            # Get the parameter value
            value = self.real_parameter_values[parameter_name]

            # Multiply the value with a random number between 1/3 and 3.
            random = np.random.uniform(log_low, log_high)
            random_factor = 10 ** random
            value = value * random_factor # DON'T DO VALUE *= RANDOM_FACTOR HERE: CHANGES THE UNDERLYING QUANTITY OBJECT AS WELL IN SELF.REAL_PARAMETER_VALUES !!

            # Set the value as the initial parameter value
            self.initial_parameter_values[parameter_name] = value

        # Debugging
        log.debug("The initial parameter values are:")
        log.debug("")
        for parameter_name in self.real_parameter_values: log.debug(" - " + parameter_name + ": " + tostr(self.initial_parameter_values[parameter_name], scientific=True, fancy=True, ndigits=parameter_ndigits[parameter_name]))
        log.debug("")

    # -----------------------------------------------------------------

    def create_template(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file template ...")

        # Load as labeled ski file
        self.ski_template = SkiFile(self.reference_ski_path)

        # Add parameter labels
        self.add_labels()

        # Set initial parameter values
        self.set_initial_values()

    # -----------------------------------------------------------------

    def add_labels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding the free parameter labels ...")

        # Add labels for absolute properties
        self.add_labels_absolute()

        # Add labels for stellar component properties
        self.add_labels_stellar_components()

        # Add labels for dust component properties
        self.add_labels_dust_components()

        # Add labels for instruments
        self.add_labels_instruments()

    # -----------------------------------------------------------------

    def add_labels_absolute(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding parameter labels for absolute simulation properties ...")

        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Search in the absolute parameters
            if parameter_name not in free_parameters_absolute_paths: continue

            # Determine path
            path = free_parameters_absolute_paths[parameter_name]

            # label
            self.ski_template.add_label_to_path(path, parameter_name)

    # -----------------------------------------------------------------

    def add_labels_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding parameter labels for stellar component properties ...")

        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Search in the stellar components
            if parameter_name not in free_parameters_relative_stellar_component_paths: continue

            # Determine the relative path to the property and the stellar component name
            path, component_name = free_parameters_relative_stellar_component_paths[parameter_name]

            # Specific component is specified
            if component_name is not None:

                # Get the stellar component
                stellar_component = self.ski_template.get_stellar_component(component_name)

                # label
                self.ski_template.add_label_to_path(path, parameter_name, stellar_component)

            # Non-specific
            else:

                # Loop over the stellar components
                for component_id in self.ski_template.get_stellar_component_ids():

                    # Get the stellar component
                    stellar_component = self.ski_template.get_stellar_component(component_id)

                    # Label
                    self.ski_template.add_label_to_path(path, parameter_name, stellar_component)

    # -----------------------------------------------------------------

    def add_labels_dust_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding labels for dust components ...")

        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Dust components
            if parameter_name not in free_parameters_relative_dust_component_paths: continue

            # Determine the relative path to the property and the dust component name
            path, component_name = free_parameters_relative_dust_component_paths[parameter_name]

            # Specific component is specified
            if component_name is not None:

                # Get the dust component
                dust_component = self.ski_template.get_dust_component(component_name)

                # label
                self.ski_template.add_label_to_path(path, parameter_name, dust_component)

            # Non-specific
            else:

                # Loop over the dust components
                for component_id in self.ski_template.get_dust_component_ids():

                    # Get the dust component
                    dust_component = self.ski_template.get_dust_component(component_id)

                    # Label
                    self.ski_template.add_label_to_path(path, parameter_name, dust_component)

    # -----------------------------------------------------------------

    def add_labels_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding parameter labels for instrument properties ...")

        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Instruments
            if parameter_name not in free_parameters_relative_instruments_paths: continue

            # Determine the relative path to the property and the instrument name
            path, instrument_name = free_parameters_relative_instruments_paths[parameter_name]

            # Specific instrument is specified
            if instrument_name is not None:

                # Get the instrument
                instrument = self.ski.get_instrument(instrument_name)

                # Label
                self.ski_template.add_label_to_path(path, parameter_name, instrument)

            # Non-specific
            else:

                # Loop over the instruments
                for instrument_name in self.ski.get_instrument_names():

                    # Get the instruemnt
                    instrument = self.ski.get_instrument(instrument_name)

                    # Label
                    self.ski_template.add_label_to_path(path, parameter_name, instrument)

    # -----------------------------------------------------------------

    def set_initial_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the initial parameter values for the fitting ...")

        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Get the value
            value = self.initial_parameter_values[parameter_name]

            # Set the value
            self.ski_template.set_labeled_value(parameter_name, value)

    # -----------------------------------------------------------------

    def load_observed_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the mock observed sed ...")

        # Determine path
        paths = fs.files_in_path(self.simulation_misc_path, contains="_fluxes", extension="dat")

        # Check if there is only one SED
        if len(paths) > 1: raise RuntimeError("More than one SED is found")
        path = paths[0]

        # Load the observed SED
        self.observed_sed = ObservedSED.from_file(path)

    # -----------------------------------------------------------------

    def setup_modelling(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Setting up the modelling ...")

        # Settings
        settings_setup = dict()
        settings_setup["type"] = "sed"
        settings_setup["name"] = self.modeling_name
        settings_setup["fitting_host_ids"] = self.moderator.host_ids_for_ensemble("fitting", none_if_none=True)

        # Create input dict for setup
        input_setup = dict()
        input_setup["sed"] = self.observed_sed
        input_setup["ski"] = self.ski_template
        #input_setup["ski_input"] = self.required_input_files
        input_setup["ski_input"] = self.all_input_files # Important so that the fittinginitializer can choose the number
        #  of wavelengths for the fitting based on what the user used as a wavelength grid file for the ski file

        # Create object config
        object_config = dict()
        # object_config["ski"] = ski_path
        input_setup["object_config"] = object_config

        # Construct the command
        stp = Command("setup", "setup the modeling", settings_setup, input_setup, cwd=".")

        # Call the command
        tool = self.run_command(stp)

    # -----------------------------------------------------------------

    def model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the modelling ...")

        # Settings
        settings_model = dict()
        settings_model["ngenerations"] = self.config.ngenerations
        settings_model["nsimulations"] = self.config.nsimulations
        settings_model["fitting_settings"] = {"spectral_convolution": False}

        # For remote execution
        settings_model["attached"] = self.config.attached
        settings_model["fitting_attached"] = self.config.attached

        # Set the random seed
        #settings_model["seed"] = self.config.seed

        # RECURRENCE SETTINGS
        settings_model["check_recurrence"] = self.config.check_recurrence
        settings_model["recurrence_rtol"] = self.config.recurrence_rtol
        settings_model["recurrence_atol"] = self.config.recurrence_atol

        # Input
        input_model = dict()

        # Set galaxy properties
        input_model["properties"] = self.properties

        # Set SEDs
        #input_model["seds"] = dict()

        # Create free parameters config

        # Create descriptions config
        descriptions_config = Map(descriptions=free_parameter_descriptions)
        input_model["descriptions_config"] = descriptions_config

        # Create types config
        types_config = Map(types=free_parameter_types)
        input_model["types_config"] = types_config

        # Create units config
        units_config = Map(units=free_parameter_units)
        input_model["units_config"] = units_config

        # Create the ndigits config
        ndigits_config = Map(ndigits=parameter_ndigits)
        input_model["ndigits_config"] = ndigits_config

        # Create filters config
        filters_config = Map(filters=fitting_filter_names)
        input_model["filters_config"] = filters_config

        # Create genetic config
        input_model["genetic_config"] = Map(genetic=self.config.genetic)

        # Create grid config ??
        #input_model["grid_config"] = Map(grid=self.config.grid)

        # Create ranges config
        ranges_config = Map()
        for parameter_name in self.config.free_parameters: # Define range
            ranges_config[parameter_name + "_range"] = self.config.relative_range_fitting * self.real_parameter_values[parameter_name]
        input_model["ranges_config"] = ranges_config

        # If we're cheating
        if self.config.cheat:
            fixed_parameter_values = defaultdict(list)
            for label in self.config.free_parameters:
                fixed_parameter_values[label].append(self.real_parameter_values[label])
            input_model["fixed_initial_parameters"] = fixed_parameter_values

        # Create initialize config
        initialize_config = Map()
        initialize_config.npackages = self.config.npackages_fitting
        initialize_config.selfabsorption = self.config.selfabsorption
        initialize_config.transient_heating = self.config.transient_heating
        input_model["initialize_config"] = initialize_config

        # Other input
        input_model["fitting_method"] = self.config.fitting_method
        input_model["parameter_grid_scales"] = self.parameter_grid_scales
        input_model["parameter_grid_weights"] = None

        # Construct the command
        command = Command("model", "perform the modelling", settings_model, input_model, cwd=self.modeling_path)

        # Run the command
        self.modeler = self.run_command(command)

    # -----------------------------------------------------------------

    @property
    def parameter_grid_scales(self):

        """
        This function ...
        :return: 
        """

        scales = dict()
        for label in self.config.free_parameters: scales[label] = self.config.scale
        return scales

    # -----------------------------------------------------------------

    @property
    def modeling_environment(self):

        """
        This function ...
        :return: 
        """

        return SEDModelingEnvironment(self.modeling_path)

    # -----------------------------------------------------------------

    def test(self):


        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Testing ...")

        # Check the best value
        #self.check_best() # CANNOT CHECK ANYMORE TEMPORARILY BECAUSE OF RECURRENCE

        # Check the database
        self.check_database()

        # Check the statistics
        #self.check_statistics() # CANNOT CHECK ANYMORE TEMPORARILY BECAUSE OF RECURRENCE

# -----------------------------------------------------------------
