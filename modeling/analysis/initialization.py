#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.initialization Contains the AnalysisInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools.logging import log
from ...core.basics.configuration import prompt_string, prompt_string_list, prompt_variable
from ...core.tools import filesystem as fs
from ...core.tools import time
from .run import AnalysisRunInfo
from ..fitting.run import get_best_model_for_generation, get_ski_file_for_simulation
from ..config.parameters import parameter_descriptions, default_units, parsing_types_for_parameter_types
from ..config.parameters import types as parameter_types
from ...core.tools.serialization import write_dict
from ...core.basics.emissionlines import EmissionLines
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from ...core.prep.dustgrids import create_one_dust_grid
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"

# -----------------------------------------------------------------

class AnalysisInitializer(AnalysisComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(AnalysisInitializer, self).__init__(*args, **kwargs)

        # The information about this analysis run
        self.analysis_run_info = None

        # The path to the analysis run info file
        self.run_info_path = None

        # The ski file
        self.ski = None

        # The dictionary of input map paths
        self.input_map_paths = dict()
        self.input_paths = dict()

        # The paths to the input files
        #self.input_paths = None

        # The path to the wavelength grid file
        self.wavelength_grid_path = None

        # The path to the dust grid file
        self.dust_grid_path = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

        # The instruments
        self.instruments = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Load from model
        if self.from_model: self.prompt_model()

        # Prompt for a fitting run
        elif self.from_fitting_run: self.prompt_fitting()

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

        # 3. Create the wavelength grid
        self.create_wavelength_grid()

        # 4. Create the dust grid
        self.create_dust_grid()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Set the input paths
        self.set_input()

        # Adapt ski file
        self.adapt_ski()

        # Write
        self.write()

    # -----------------------------------------------------------------

    @property
    def from_model(self):

        """
        This function ...
        :return:
        """

        return self.config.origin == "model"

    # -----------------------------------------------------------------

    @property
    def from_fitting_run(self):

        """
        This function ...
        :return:
        """

        return self.config.origin == "fitting_run"

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisInitializer, self).setup(**kwargs)

        # Generate a name for this analysis run
        self.analysis_run_name = time.unique_name()

        # Create a directory for this analysis run
        self.analysis_run_path = fs.join(self.analysis_path, self.analysis_run_name)

        # Set the ski file path
        self.ski_file_path = fs.join(self.analysis_run_path, self.galaxy_name + ".ski")

        # Set the path to the wavelength grid file
        # self.wavelength_grid_path = fs.join(self.analysis_run_path, "wavelength_grid.dat")

        # Set the path to the dust grid file
        # self.dust_grid_path = fs.join(self.analysis_run_path, "dust_grid.dg")

        # Set the path to the analysis run info file
        self.run_info_path = fs.join(self.analysis_run_path, "info.dat")

        # Create the info object
        self.analysis_run_info = AnalysisRunInfo()

        # Set the info
        self.analysis_run_info.name = self.analysis_run_name
        self.analysis_run_info.path = self.analysis_run_path

    # -----------------------------------------------------------------

    def prompt_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the initial model to use for analysis ...")

        # Ask for the model name
        model_name = prompt_string("model", "name of the model", choices=self.model_suite.model_names, required=True)

        # Load the model
        #definition = self.model_suite.get_model_definition(model_name)

        # Add the components to the ski file and to the input map paths dictionary
        self.model_suite.add_model_components(model_name, self.ski, self.input_map_paths)

        # PRINT ALL PARAMETERS??

        # Get parameter labels to adapt
        parameter_labels = prompt_string_list("parameters", "names of the parameters to adapt the value", choices=parameter_descriptions)

        # Get the default parameter values
        default_values = dict()
        for label in parameter_labels:

            # Get the value from the ski file
            value = self.ski.get_labeled_value(label)

            # Set the default value
            default_values[label] = value

        # Prompt for parameter values
        parameter_values = dict()
        for label in parameter_labels:

            # Get the parameter type
            parameter_type = parameter_types[label]

            # Get the parsing type
            ptype = parsing_types_for_parameter_types[parameter_type]

            # Get the default unit
            unit = default_units[label]

            # Ask for the value
            value = prompt_variable(label, ptype, "value for the " + parameter_descriptions[label], default=default_values[label], required=False) # can be optional

            # Set the value
            parameter_values[label] = value.to(unit)

        # Set info
        #self.analysis_run_info.fitting_run = None
        #self.analysis_run_info.generation_name = None
        #self.analysis_run_info.simulation_name = None
        self.analysis_run_info.parameter_values = parameter_values
        #self.analysis_run_info.chi_squared = None
        self.analysis_run_info.model_name = model_name

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the origin of the analysis model ...")

        # Prompt for fitting run
        run_id = prompt_string("fitting_run", "name of the fitting run", choices=self.fitting_context.fitting_run_names)

        # Load the fitting run
        fitting_run = self.fitting_context.load_fitting_run(run_id)

        # Check
        if not fitting_run.has_finished_generations: raise ValueError("Fitting run has no finished generations")

        # Prompt for the generation name
        generation_name = prompt_string("generation", "name of the (finished) generation", default=fitting_run.last_finished_generation, choices=fitting_run.finished_generations)

        # Load the best model for the specified generation
        best_model = get_best_model_for_generation(self.config.path, run_id, generation_name)

        # Set info
        self.analysis_run_info.fitting_run = run_id
        self.analysis_run_info.generation_name = generation_name
        self.analysis_run_info.simulation_name = best_model.simulation_name
        self.analysis_run_info.parameter_values = best_model.parameter_values
        self.analysis_run_info.chi_squared = best_model.chi_squared
        self.analysis_run_info.model_name = fitting_run.model_name

        # Load the ski file
        self.ski = get_ski_file_for_simulation(self.config.path, run_id, generation_name, best_model.simulation_name)

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
        pixelscale_angular = self.reference_wcs.average_pixelscale #* u("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        x_radius = radius_physical
        y_radius = radius_physical
        z_radius = 3. * u("kpc")

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

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Loop over the projections
        for projection in ["earth", "faceon", "edgeon"]:

            # Debugging
            log.debug("Creating a full instrument for the " + projection + " projection ...")

            # Create the instrument and add it to the dictionary
            self.instruments[projection] = self.create_instrument("full", projection)

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
        #self.input_paths.append(self.wavelength_grid_path)
        self.input_paths[wavelengths_filename] = self.wavelength_grid_path

    # -----------------------------------------------------------------

    def adapt_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the ski file ...")

        # Set wavelength grid for ski file
        self.ski.set_file_wavelength_grid(wavelengths_filename)

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

        # Write out the dust grid data
        #self.ski.set_write_grid()

        # TODO: GENERATE GRID??!

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the info
        self.write_info()

        # Create the necessary directories
        self.write_directories()

        # Write the ski file
        self.write_ski()

        #self.write_input_map_paths()
        self.write_input_paths()

        self.write_dust_grid()

        self.write_wavelength_grid()

        self.write_instruments()

    # -----------------------------------------------------------------

    def write_info(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the analysis run info to " + self.run_info_path + "...")

        # Write the analysis run info
        self.analysis_run_info.saveto(self.run_info_path)

    # -----------------------------------------------------------------

    def write_directories(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the necessary directories for this analysis run ...")

        # The directory with the instruments
        self.run_instruments_path = fs.create_directory_in(self.analysis_run_path, "instruments")

        # Simulation directories
        self.run_output_path = fs.create_directory_in(self.analysis_run_path, "out")
        self.run_extr_path = fs.create_directory_in(self.analysis_run_path, "extr")
        self.run_plot_path = fs.create_directory_in(self.analysis_run_path, "plot")
        self.run_misc_path = fs.create_directory_in(self.analysis_run_path, "misc")

        # Analysis directories
        self.run_attenuation_path = fs.create_directory_in(self.analysis_run_path, "attenuation")
        self.run_colours_path = fs.create_directory_in(self.analysis_run_path, "colours")
        self.run_residuals_path = fs.create_directory_in(self.analysis_run_path, "residuals")
        self.run_heating_path = fs.create_directory_in(self.analysis_run_path, "heating")

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

    # -----------------------------------------------------------------

    def write_input_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dictionary of input map paths ...")

        # Write
        #write_dict(self.input_map_paths, self.input_maps_file_path)
        write_dict(self.input_paths, self.input_file_path)

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid to " + self.wavelength_grid_path + " ...")

        # Write the wavelength grid
        self.wavelength_grid.to_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid " + self.dust_grid_path + " ...")

        # Write the dust grid
        self.dust_grid.saveto(self.dust_grid_path)

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

            # Determine the path
            path = fs.join(self.run_instruments_path, name + ".instr")

            # Save the instrument
            self.instruments[name].saveto(path)

# -----------------------------------------------------------------
