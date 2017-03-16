#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import inspect
import numpy as np
from collections import defaultdict

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools.logging import log
from pts.modeling.basics.instruments import SEDInstrument
from pts.modeling.tests.base import M81TestBase, fitting_filter_names, m81_data_path
from pts.core.launch.options import AnalysisOptions
from pts.do.commandline import Command
from pts.core.simulation.skifile import LabeledSkiFile
from pts.modeling.tests.base import instrument_name, free_parameters_absolute_paths
from pts.modeling.tests.base import free_parameters_relative_stellar_component_paths, free_parameters_relative_dust_component_paths
from pts.modeling.tests.base import free_parameters_relative_instruments_paths, free_parameter_types
from pts.core.data.sed import ObservedSED
from pts.core.tools import parsing, stringify
from pts.core.tools import sequences

# -----------------------------------------------------------------

# Determine path of this directory
this_path = fs.absolute_path(inspect.stack()[0][1])
this_dir_path = fs.directory_of(this_path)

# -----------------------------------------------------------------

description = "fitting the parameters of a model of M81 based on only a mock observed SED"

# -----------------------------------------------------------------

class M81SEDTest(M81TestBase):

    """
    This class runs the test on M81, but by only adjusting the normalizations (not by creating a model),
    and fitting to a mock observed SED
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor
        :param config:
        :param interactive:
        """

        # Call the constructor of the base class
        super(M81SEDTest, self).__init__(config, interactive)

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

        # The real parameter values
        self.real_parameter_values = dict()

        # The initial parameter values for the fitting
        self.initial_parameter_values = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the properties
        self.load_properties()

        # 3. Load the components
        self.load_components()

        # Create instrument
        self.create_instrument()

        # Create deprojection
        self.create_deprojections()

        # Write
        self.write()

        # Plot
        self.plot()

        # Launch the reference simulation
        self.launch_reference()

        # Get the real parameter values
        self.get_real_parameter_values()

        # Generate the initial parameter values
        self.generate_initial_parameter_values()

        # Create the ski file template
        self.create_template()

        # Load the observed SED
        self.load_observed_sed()

        # Setup the modeling
        self.setup_modelling()

        # Model
        self.model()

    # -----------------------------------------------------------------

    def create_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the SED instrument ...")

        azimuth = 0.

        # Create the SED instrument
        self.instrument = SEDInstrument(distance=self.galaxy_distance, inclination=self.galaxy_inclination, azimuth=azimuth, position_angle=self.galaxy_position_angle)

    # -----------------------------------------------------------------

    def write(self):

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

    def plot(self):

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
        #settings_launch["remote"] = self.host_id
        settings_launch["attached"] = True
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
        seds_path = fs.join(m81_data_path, "seds")
        analysis.plotting.reference_seds = fs.files_in_path(seds_path)
        analysis.misc.fluxes = True
        #analysis.misc.images = True
        analysis.misc.observation_filters = fitting_filter_names
        analysis.misc.observation_instruments = [instrument_name]
        #analysis.misc.make_images_remote = self.host_id
        #analysis.misc.images_wcs = self.reference_wcs_path
        #analysis.misc.images_unit = "Jy/pix"
        analysis.misc.spectral_convolution = False

        # Input
        input_launch = dict()
        # input_launch["memory"] = MemoryRequirement(serial_memory, parallel_memory)

        # Launch command
        launch = Command("launch_simulation", "launch the reference simulation", settings_launch, input_launch, cwd=".")
        self.launcher = self.run_command(launch)

    # -----------------------------------------------------------------

    def get_real_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the real parameter values ...")

        # Store the different values encountered in the ski file
        values = defaultdict(list)

        # Add the labels for the free parameters
        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Get the parsing function for this parameter
            parser = getattr(parsing, free_parameter_types[parameter_name])

            # Search in the absolute parameters
            if parameter_name in free_parameters_absolute_paths:

                # Determine the path to the property
                path = free_parameters_absolute_paths

                # Get the current value
                value = parser(self.ski.get_value_for_path(path))

                # Set the value
                values[parameter_name].append(value)

            # Search in the stellar components
            if parameter_name in free_parameters_relative_stellar_component_paths:

                # Determine the relative path to the property and the stellar component name
                path, component_name = free_parameters_relative_stellar_component_paths[parameter_name]

                if component_name is not None:

                    # Get the stellar component
                    stellar_component = self.ski.get_stellar_component(component_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, stellar_component))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the stellar components
                    for component_id in self.ski.get_stellar_component_ids():

                        # Get the stellar component
                        stellar_component = self.ski.get_stellar_component(component_id)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, stellar_component))

                        # Set the value
                        values[parameter_name].append(value)

            # Search in the dust components
            if parameter_name in free_parameters_relative_dust_component_paths:

                # Determine the relative path to the property and the dust component name
                path, component_name = free_parameters_relative_dust_component_paths[parameter_name]

                if component_name is not None:

                    # Get the dust component
                    dust_component = self.ski.get_dust_component(component_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, dust_component))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the dust components
                    for component_id in self.ski.get_dust_component_ids():

                        # Get the dust component
                        dust_component = self.ski.get_dust_component(component_id)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, dust_component))

                        # Set the value
                        values[parameter_name].append(value)

            # Search in instruments
            if parameter_name in free_parameters_relative_instruments_paths:

                # Determine the relative path to the property and the instrument name
                path, instrument_name = free_parameters_relative_instruments_paths[parameter_name]

                if instrument_name is not None:

                    # Get the instrument
                    instrument = self.ski.get_instrument(instrument_name)

                    # Get the current value
                    value = parser(self.ski.get_value_for_path(path, instrument))

                    # Set the value
                    values[parameter_name].append(value)

                else:

                    # Loop over the instruments
                    for instrument_name in self.ski.get_instrument_names():

                        # Get the instruemnt
                        instrument = self.ski.get_instrument(instrument_name)

                        # Get the current value
                        value = parser(self.ski.get_value_for_path(path, instrument))

                        # Set the value
                        values[parameter_name].append(value)

        # Check whether we have only one value for each parameter
        for parameter_name in self.config.free_parameters:

            # Check if any
            if len(values[parameter_name]) == 0: raise ValueError("No parameter values for '" + parameter_name + "' were found in the ski file")

            # Check if all equal
            if not sequences.all_equal(values[parameter_name]): raise ValueError("Parameter values for '" + parameter_name + "' are not equal throughout the ski file")

            # Set the unique real parameter value
            self.real_parameter_values[parameter_name] = values[parameter_name][0]

        # Debugging
        log.debug("The real parameter values are: ")
        log.debug("")
        for parameter_name in self.real_parameter_values: log.debug(" - " + parameter_name + ": " + stringify.stringify(self.real_parameter_values[parameter_name])[1])
        log.debug("")

    # -----------------------------------------------------------------

    def generate_initial_parameter_values(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Generating random initial parameter values ...")

        low_factor = 0.3
        high_factor = 3.

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
            value = value * random_factor

            # Set the value as the initial parameter value
            self.initial_parameter_values[parameter_name] = value

        # Debugging
        log.debug("The initial parameter values are:")
        log.debug("")
        for parameter_name in self.real_parameter_values: log.debug(" - " + parameter_name + ": " + stringify.stringify(self.initial_parameter_values[parameter_name])[1])
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
        self.ski_template = LabeledSkiFile(self.reference_ski_path)

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

        # Add the labels
        # Loop over the free parameters
        for parameter_name in self.config.free_parameters:

            # Search in the absolute parameters
            if parameter_name in free_parameters_absolute_paths:

                # Determine path
                path = free_parameters_absolute_paths[parameter_name]

                # label
                self.ski_template.add_label_to_path(path, parameter_name)

            # Search in the stellar components
            if parameter_name in free_parameters_relative_stellar_component_paths:

                # Determine the relative path to the property and the stellar component name
                path, component_name = free_parameters_relative_stellar_component_paths[parameter_name]

                if component_name is not None:

                    # Get the stellar component
                    stellar_component = self.ski_template.get_stellar_component(component_name)

                    # label
                    self.ski_template.add_label_to_path(path, parameter_name, stellar_component)

                else:

                    # Loop over the stellar components
                    for component_id in self.ski_template.get_stellar_component_ids():

                        # Get the stellar component
                        stellar_component = self.ski_template.get_stellar_component(component_id)

                        # Label
                        self.ski_template.add_label_to_path(path, parameter_name, stellar_component)

            # Dust components
            if parameter_name in free_parameters_relative_dust_component_paths:

                # Determine the relative path to the property and the dust component name
                path, component_name = free_parameters_relative_dust_component_paths[parameter_name]

                if component_name is not None:

                    # Get the dust component
                    dust_component = self.ski_template.get_dust_component(component_name)

                    # label
                    self.ski_template.add_label_to_path(path, parameter_name, dust_component)

                else:

                    # Loop over the dust components
                    for component_id in self.ski_template.get_dust_component_ids():

                        # Get the dust component
                        dust_component = self.ski_template.get_dust_component(component_id)

                        # Label
                        self.ski_template.add_label_to_path(path, parameter_name, dust_component)

            # Instruments
            if parameter_name in free_parameters_relative_instruments_paths:

                # Determine the relative path to the property and the instrument name
                path, instrument_name = free_parameters_relative_instruments_paths[parameter_name]

                if instrument_name is not None:

                    # Get the instrument
                    instrument = self.ski.get_instrument(instrument_name)

                    # Label
                    self.ski_template.add_label_to_path(path, parameter_name, instrument)

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
        settings_setup["name"] = self.galaxy_name
        settings_setup["fitting_host_ids"] = None

        # Create object config
        object_config = dict()
        #object_config["ski"] = ski_path

        # Create input dict for setup
        input_setup = dict()
        input_setup["object_config"] = object_config
        input_setup["sed"] = self.observed_sed
        input_setup["ski"] = self.ski_template

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

# -----------------------------------------------------------------
