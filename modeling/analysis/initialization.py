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
from ...core.tools import filesystem as fs
from ...core.tools import time
from .run import AnalysisRunInfo, AnalysisRun
from ...core.tools.serialization import write_dict
from ...core.basics.emissionlines import EmissionLines
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from ...core.prep.dustgrids import create_one_dust_grid
from ...core.units.parsing import parse_unit as u
from ...core.simulation.input import SimulationInput
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.stringify import tostr
from ..misc.select import select_from_model_suite, select_from_fitting_context

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

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

        # The analysis run
        self.analysis_run = None

        # The path to the analysis run info file
        self.run_info_path = None

        # The model definition
        self.definition = None

        # The parameter values
        self.parameter_values = None

        # The ski file
        self.ski = None

        # The dictionary of input map paths
        self.input_paths = None

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

        # 2. Get the model
        self.get_model()

        # 3. Create the analysis run
        self.create_analysis_run()

        # 4. Create the wavelength grid
        self.create_wavelength_grid()

        # 5. Create the dust grid
        self.create_dust_grid()

        # 6. Create the instruments
        self.create_instruments()

        # 7. Adapt ski file
        self.adapt_ski()

        # Set the input
        self.set_input()

        # 8. Write
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

        # Create the info object
        self.analysis_run_info = AnalysisRunInfo()

        # Set the analysis run name and path
        self.analysis_run_info.name = self.analysis_run_name
        self.analysis_run_info.path = self.analysis_run_path

    # -----------------------------------------------------------------

    def get_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the analysis model ...")

        # Load from model
        if self.from_model: self.prompt_model()

        # Prompt for a fitting run
        elif self.from_fitting_run: self.prompt_fitting()

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

        # Show the model parameters
        print("")
        print("Model parameter values:")
        print("")
        for label in self.parameter_values: print(" - " + label + ": " + tostr(self.parameter_values[label]))
        print("")

    # -----------------------------------------------------------------

    def prompt_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the initial model to use for analysis ...")

        # # Ask for the model name
        # model_name = prompt_string("model", "name of the model", choices=self.model_suite.model_names, required=True)
        #
        # # Load the labeled ski template file
        # self.ski = LabeledSkiFile(template_ski_path)
        #
        # # Load the model
        # self.definition = self.model_suite.get_model_definition(model_name)
        #
        # # Add the components to the ski file and to the input map paths dictionary
        # self.model_suite.add_model_components(model_name, self.ski, self.input_paths)
        #
        # # Get the parameter values
        # self.parameter_values = self.prompt_parameters()

        # NEW
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite)
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Set info
        # self.analysis_run_info.fitting_run = None
        # self.analysis_run_info.generation_name = None
        # self.analysis_run_info.simulation_name = None
        self.analysis_run_info.parameter_values = self.parameter_values
        # self.analysis_run_info.chi_squared = None
        self.analysis_run_info.model_name = model_name

    # -----------------------------------------------------------------

    @lazyproperty
    def use_file_tree_dust_grid(self):

        """
        This function ...
        :return:
        """

        smile = SKIRTSmileSchema()
        #return smile.supports_file_tree_grids and self.representation.has_dust_grid_tree
        if not smile.supports_file_tree_grids: raise RuntimeError("A version of SKIRT that supports file tree grids is necessary")
        return True

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the origin of the analysis model ...")

        # # Prompt for fitting run
        # run_id = prompt_string("fitting_run", "name of the fitting run", choices=self.fitting_context.fitting_run_names)
        #
        # # Load the fitting run
        # fitting_run = self.fitting_context.load_fitting_run(run_id)
        #
        # # Prompt for the generation
        # generation_name = self.prompt_generation(fitting_run)
        #
        # # Get the parameter values
        # simulation_name, self.parameter_values, chi_squared = self.get_parameters(fitting_run, generation_name)

        # NEW
        run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values = select_from_fitting_context(self.fitting_context)

        # Set
        self.parameter_values = parameter_values
        self.ski = ski
        self.input_paths = input_paths

        # Set info
        self.analysis_run_info.fitting_run = run_id
        self.analysis_run_info.generation_name = generation_name
        self.analysis_run_info.simulation_name = simulation_name
        self.analysis_run_info.parameter_values = self.parameter_values
        self.analysis_run_info.chi_squared = chi_squared
        self.analysis_run_info.model_name = fitting_run.model_name

        # # Load the ski file
        # if simulation_name is not None: self.ski = get_ski_file_for_simulation(self.config.path, run_id, generation_name, simulation_name)
        # else: self.ski = fitting_run.ski_template
        #
        # # Set parameter values (ALTHOUGH PROBABLY UNNECESSARY)
        # self.ski.set_labeled_values(self.parameter_values)
        #
        # # Load the input paths
        # input_map_paths = fitting_run.input_map_paths
        # self.input_paths = input_map_paths

    # -----------------------------------------------------------------

    def create_analysis_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the analysis run ...")

        # Create the generation object
        self.analysis_run = AnalysisRun(self.galaxy_name, self.analysis_run_info)

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

        # Determine and set the path to the wavelength grid file
        self.input_paths[wavelengths_filename] = self.wavelength_grid_path

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

    def set_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the simulation input ...")

        # Initialize the SimulationInput object
        # self.simulation_input = SimulationInput()

        # Set the paths to the input maps
        # for name in input_map_paths:
        #    path = input_map_paths[name]
        #    #self.simulation_input.add_file(path, name)
        #    self.input_paths[name] = path

        # NEW: DETERMINE AND SET THE PATH TO THE APPROPRIATE DUST GRID TREE FILE
        if self.use_file_tree_dust_grid:

            # self.simulation_input.add_file(self.representation.dust_grid_tree_path)
            self.input_paths[dustgridtree_filename] = self.dust_grid_tree_path

        # Determine and set the path to the appropriate wavelength grid file
        # wavelength_grid_path = self.fitting_run.wavelength_grid_path_for_level(self.generation_info.wavelength_grid_level)
        # self.simulation_input.add_file(wavelength_grid_path)
        self.input_paths[wavelengths_filename] = self.wavelength_grid_path

        # Get the number of wavelengths
        # self.nwavelengths = len(WavelengthGrid.from_skirt_input(wavelength_grid_path))

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

        # Write the input paths
        self.write_input_paths()

        # Write the dust grid
        self.write_dust_grid()

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the instruments
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
