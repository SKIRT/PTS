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

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.tools import time
from .run import AnalysisRunInfo, AnalysisRun
from ...core.tools.serialization import write_dict
from ...core.basics.emissionlines import EmissionLines
from ...core.prep.wavelengthgrids import create_one_subgrid_wavelength_grid
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.stringify import tostr
from ..misc.select import select_from_model_suite, select_from_fitting_context
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection
from ..build.dustgrid import DustGridBuilder
from ..basics.instruments import FullInstrument
from ..build.representation import create_projections_from_dust_grid, create_projections_from_deprojections
from ...core.basics.configuration import prompt_yn

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

        # The ski file
        self.ski = None

        # The model definition
        self.definition = None

        # The parameter values
        self.parameter_values = None

        # The dictionary of input map paths
        self.input_paths = None

        # The wavelength grid
        self.wavelength_grid = None

        # The dust grid
        self.dust_grid = None

        # The deprojections
        self.deprojections = dict()

        # The projections and instruments
        self.projections = dict()
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

        # Load the deprojections
        self.load_deprojections()

        # Create the projections
        self.create_projections()

        # 6. Create the instruments
        self.create_instruments()

        # 7. Adapt ski file
        self.adapt_ski()

        # Set the input
        self.set_input()

        # Build the dust grid
        self.build_dust_grid()

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

        # Select model
        model_name, ski, definition, input_paths, parameter_values = select_from_model_suite(self.model_suite)

        # Set attributes
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
        if not self.analysis_run.has_dust_grid_tree: raise RuntimeError("The dust grid tree is not present at '"  + self.analysis_run.dust_grid_tree_path + "'")
        return True

    # -----------------------------------------------------------------

    def prompt_fitting(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the origin of the analysis model ...")

        # Select model
        run_id, generation_name, simulation_name, fitting_run, chi_squared, ski, input_paths, parameter_values = select_from_fitting_context(self.fitting_context)

        # Load the model definition
        definition = fitting_run.model_definition

        # Set
        self.ski = ski
        self.definition = definition
        self.parameter_values = parameter_values
        self.input_paths = input_paths

        # Set info
        self.analysis_run_info.fitting_run = run_id
        self.analysis_run_info.generation_name = generation_name
        self.analysis_run_info.simulation_name = simulation_name
        self.analysis_run_info.parameter_values = self.parameter_values
        self.analysis_run_info.chi_squared = chi_squared
        self.analysis_run_info.model_name = fitting_run.model_name

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

        # Fixed wavelengths (always in the grid: for normalization)
        fixed = [self.i1_filter.pivot, self.fuv_filter.pivot]

        # Create the emission lines instance
        if self.config.wg.add_emission_lines: emission_lines = EmissionLines()
        else: emission_lines = None

        # Create the grid
        # grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added
        grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added = \
            create_one_subgrid_wavelength_grid(self.config.wg.npoints, emission_lines, fixed, min_wavelength=self.config.wg.range.min, max_wavelength=self.config.wg.range.max)

        # Set the grid
        self.wavelength_grid = grid

        # Debugging
        log.debug("Created a wavelength grid with:")
        log.debug("")
        log.debug(" - number of points: " + str(len(grid)))
        log.debug(" - number of points in subgrids: " + str(subgrid_npoints))
        log.debug(" - number of emission points: " + str(emission_npoints))
        log.debug(" - number of fixed points: " + str(fixed_npoints))
        log.debug(" - filters for which extra sampling was performed: " + str(broad_resampled))
        log.debug(" - narrow band filters for which wavelength was added: " + str(narrow_added))
        log.debug("")
        if log.is_debug():
            print("")
            print(self.wavelength_grid)
            print("")

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """


        # Inform the user
        log.info("Creating the dust grid ...")

        # Load the dust disk deprojection
        deprojection = self.definition.dust_deprojection

        # Set minimum level
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_max_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_max_level
        else: min_level = None

        # Set max ndivisions per pixel
        max_ndivisions_per_pixel = 1. / self.config.dg.scale  # default 1/0.5 = 2 divisions along each direction per pixel

        # Create the dust grid
        # grid_type, deprojection, distance, sky_ellipse, min_level, max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.
        self.dust_grid = create_one_dust_grid_for_galaxy_from_deprojection(self.config.dg.grid_type, deprojection,
                                                                           self.galaxy_distance, self.truncation_ellipse,
                                                                           min_level, self.config.dg.max_mass_fraction,
                                                                           max_ndivisions_per_pixel, self.config.dg.scale_heights)

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.definition.name

    # -----------------------------------------------------------------

    def load_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the deprojections ...")

        # Stellar
        self.load_stellar_deprojections()

        # Dust
        self.load_dust_deprojections()

    # -----------------------------------------------------------------

    def load_stellar_deprojections(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the stellar deprojections ...")

        # Loop over the stellar components
        for name in self.model_suite.get_stellar_component_names(self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.model_suite.load_stellar_component_deprojection(self.model_name, name)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def load_dust_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust deprojections ...")

        # Loop over the dust components
        for name in self.model_suite.get_dust_component_names(self.config.path, self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.model_suite.load_dust_component_deprojection(self.model_name, name)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        azimuth = 0.0

        # Use grid?
        if prompt_yn("grid_resolution", "use the resolution of the dust grid for setting up the instruments?"):

            earth, faceon, edgeon = create_projections_from_dust_grid(self.dust_grid, self.galaxy_distance,
                                                                      self.galaxy_inclination, azimuth,
                                                                      self.disk_position_angle)

        # Use deprojections
        else: earth, faceon, edgeon = create_projections_from_deprojections(self.deprojections, self.galaxy_distance, azimuth)

        # Set the projection systems
        self.projections["earth"] = earth
        self.projections["faceon"] = faceon
        self.projections["edgeon"] = edgeon

    # -----------------------------------------------------------------

    @property
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["earth"]

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["faceon"]

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["edgeon"]

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Create an earth instrument
        self.instruments["earth"] = FullInstrument.from_projection(self.earth_projection)

        # Create a faceon instrument
        self.instruments["faceon"] = FullInstrument.from_projection(self.faceon_projection)

        # Create an edgeon instrument
        self.instruments["edgeon"] = FullInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    @property
    def earth_instrument(self):

        """
        This function ...
        :return:
        """

        return self.instruments["earth"]

    # -----------------------------------------------------------------

    @property
    def faceon_instrument(self):

        """
        This function ...
        :return:
        """

        return self.instruments["faceon"]

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument(self):

        """
        Thsi function ...
        :return:
        """

        return self.instruments["edgeon"]

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

    # -----------------------------------------------------------------

    def build_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the dust grid ...")

        # Create the builder
        builder = DustGridBuilder()

        # Set output path
        builder.config.output = self.analysis_run.dust_grid_build_path

        # Set simulation path
        builder.config.simulation_path = self.analysis_run.dust_grid_simulation_out_path

        # Set whether quality has to be calculated
        builder.config.quality = self.config.check_dust_grid_quality

        # Run the dust grid builder
        builder.run(definition=self.definition, dust_grid=self.dust_grid)

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
            self.input_paths[dustgridtree_filename] = self.analysis_run.dust_grid_tree_path

        # Determine and set the path to the appropriate wavelength grid file
        # wavelength_grid_path = self.fitting_run.wavelength_grid_path_for_level(self.generation_info.wavelength_grid_level)
        # self.simulation_input.add_file(wavelength_grid_path)
        self.input_paths[wavelengths_filename] = self.analysis_run.wavelength_grid_path

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

        # Write the config
        self.write_config()

        # Write the info
        self.write_info()

        # Write the ski file
        self.write_ski()

        # Write the input paths
        self.write_input_paths()

        # Write the dust grid
        self.write_dust_grid()

        # Write the wavelength grid
        self.write_wavelength_grid()

        # Write the projections
        self.write_projections()

        # Write the instruments
        self.write_instruments()

    # -----------------------------------------------------------------

    def write_config(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the configuration used to create this analysis run ...")

        # Write
        self.config.saveto(self.analysis_run.config_path)

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
        write_dict(self.input_paths, self.analysis_run.input_file_path)

    # -----------------------------------------------------------------

    def write_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grid to " + self.analysis_run.wavelength_grid_path + " ...")

        # Write the wavelength grid
        self.wavelength_grid.to_skirt_input(self.analysis_run.wavelength_grid_path)

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid " + self.analysis_run.dust_grid_path + " ...")

        # Write the dust grid
        self.dust_grid.saveto(self.analysis_run.dust_grid_path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections["earth"].saveto(self.analysis_run.earth_projection_path)

        # Write the faceon projection system
        self.projections["faceon"].saveto(self.analysis_run.faceon_projection_path)

        # Write the edgeon projection system
        self.projections["edgeon"].saveto(self.analysis_run.edgeon_projection_path)

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the instruments ...")

        # Write the SED instrument
        self.instruments["earth"].saveto(self.analysis_run.earth_instrument_path)

        # Write the frame instrument
        self.instruments["faceon"].saveto(self.analysis_run.faceon_instrument_path)

        # Write the simple instrument
        self.instruments["edgeon"].saveto(self.analysis_run.edgeon_instrument_path)

# -----------------------------------------------------------------
