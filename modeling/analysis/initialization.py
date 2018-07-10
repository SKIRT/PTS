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
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import time
from .run import AnalysisRunInfo, AnalysisRun
from ...core.tools.serialization import write_dict
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.stringify import tostr
from ..build.dustgrid import DustGridBuilder
from ..basics.instruments import FullInstrument, SimpleInstrument, SEDInstrument, FullSEDInstrument
from ..misc.interface import ModelSimulationInterface, earth_name, edgeon_name, faceon_name
from .run import info_filename
from ...core.tools import formatting as fmt
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..build.definition import ModelDefinition

# -----------------------------------------------------------------

# Input filenames
wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

# -----------------------------------------------------------------

# Define instrument names
full_earth_name = "earth_full"
simple_earth_name = "earth_simple"
sed_earth_name = "earth_sed"
full_sed_earth_name = "earth_full_sed"
simple_faceon_name = "faceon_simple"
full_faceon_name = "faceon_full"
simple_edgeon_name = "edgeon_simple"
full_edgeon_name = "edgeon_full"

# -----------------------------------------------------------------

class AnalysisInitializer(AnalysisComponent, ModelSimulationInterface):

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
        #super(AnalysisInitializer, self).__init__(*args, **kwargs)
        AnalysisComponent.__init__(self, no_config=True)
        ModelSimulationInterface.__init__(self, *args, **kwargs)

        # The information about this analysis run
        self.analysis_run_info = None

        # The analysis run
        self.analysis_run = None

        # The path to the analysis run info file
        self.run_info_path = None

    # -----------------------------------------------------------------

    @property
    def do_build_dust_grid(self):

        """
        This function ...
        :return:
        """

        return not self.grid_from_representation or not self.representation.has_dust_grid_tree

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Get the model
        self.get_model()

        # 2. Create the analysis run
        self.create_analysis_run()

        # 3. Get the wavelength grid
        self.get_wavelength_grid()

        # 4. Get the dust grid
        self.get_dust_grid()

        # 5. Load the deprojections
        self.load_deprojections()

        # 6. Get the projections
        self.get_projections()

        # 7. Get the instruments
        self.create_instruments()

        # 8. Adapt ski file
        self.adapt_ski()

        # 9. Build the dust grid
        if self.do_build_dust_grid: self.build_dust_grid()

        # 10. Set the input
        self.set_input()

        # 11. Write
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
        #super(AnalysisInitializer, self).setup(**kwargs)
        AnalysisComponent.setup(self, **kwargs)
        ModelSimulationInterface.setup(self, **kwargs)

        # Generate a name for this analysis run
        if self.config.name is not None: analysis_run_name = self.config.name
        else: analysis_run_name = time.unique_name()

        # Create a directory for this analysis run
        analysis_run_path = fs.join(self.analysis_path, analysis_run_name)
        if not fs.is_directory(analysis_run_path): fs.create_directory(analysis_run_path)
        elif fs.is_empty(analysis_run_path, recursive=True): fs.clear_directory(analysis_run_path)
        elif self.config.overwrite: fs.clear_directory(analysis_run_path)
        else: raise ValueError("There already exists a directory for this analysis run")

        # Create the info object
        self.analysis_run_info = AnalysisRunInfo()

        # Set the analysis run name and path
        self.analysis_run_info.name = analysis_run_name
        self.analysis_run_info.path = analysis_run_path

        # Set run info path
        self.run_info_path = fs.join(self.analysis_run_path, info_filename)

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):

        """
        Thisf unction ...
        :return:
        """

        return self.analysis_run_info.name

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_run_info.path

    # -----------------------------------------------------------------

    def get_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the analysis model ...")

        # Load from model
        if self.from_model:

            # Get the model
            model_name = self.prompt_model()

            # Update the analysis info
            self.analysis_run_info.model_name = model_name
            self.analysis_run_info.parameter_values = self.parameter_values

        # Prompt for a fitting run
        elif self.from_fitting_run:

            # Get the model
            run_id, generation_name, simulation_name, chi_squared = self.prompt_fitting()

            # Update the analysis info
            self.analysis_run_info.fitting_run = run_id
            self.analysis_run_info.generation_name = generation_name
            self.analysis_run_info.simulation_name = simulation_name
            self.analysis_run_info.chi_squared = chi_squared
            self.analysis_run_info.parameter_values = self.parameter_values

            # Set the name of the corresponding model of the model suite
            self.analysis_run_info.model_name = self.definition.name

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

        # Show the model parameters
        if self.from_fitting_run or self.config.adapt:
            print("")
            print("Adapted model parameter values:")
            print("")
            for label in self.parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.parameter_values[label]))
            print("")

        # Show all model parameters
        else:
            print("")
            print("All model parameter values:")
            print("")
            for label in self.parameter_values: print(" - " + fmt.bold + label + fmt.reset + ": " + tostr(self.parameter_values[label]))
            print("")

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        if not self.from_fitting_run: raise RuntimeError("Not from fitting run")
        return self.fitting_runs.load(self.analysis_run_info.fitting_run)

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

    def create_analysis_run(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the analysis run ...")

        # Create the generation object
        self.analysis_run = AnalysisRun(self.galaxy_name, self.analysis_run_info, hubble_stage=self.hubble_stage)

    # -----------------------------------------------------------------

    def get_projections(self):

        """
        This function ...
        :return:
        """

        # Load from representaiton
        if self.projections_from_representation: self.load_projections()

        # Create new projections
        else: self.create_projections()

    # -----------------------------------------------------------------

    def load_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the projections ...")

        # Set the projection systems
        self.projections[earth_name] = self.representation.earth_projection
        self.projections[faceon_name] = self.representation.faceon_projection
        self.projections[edgeon_name] = self.representation.edgeon_projection

        # Set reference deprojection from representation
        self.analysis_run_info.reference_deprojection = self.representation.reference_deprojection_name

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ..
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Create projections
        deprojection_name = self.create_projection_systems(make_faceon=True, make_edgeon=True, reference_name=self.config.projections_reference)

        # Set the deprojection name in the analysis info
        self.analysis_run_info.reference_deprojection = deprojection_name

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Earth
        self.create_earth_instruments()

        # Faceon
        self.create_faceon_instruments()

        # Edgeon
        self.create_edgeon_instruments()

    # -----------------------------------------------------------------

    def create_earth_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the earth instruments ...")

        # Create full earth instrument
        self.instruments[full_earth_name] = FullInstrument.from_projection(self.earth_projection, **self.full_instrument_properties)

        # Create simple earth instrument
        self.instruments[simple_earth_name] = SimpleInstrument.from_projection(self.earth_projection)

        # Create SED earth instrument
        self.instruments[sed_earth_name] = SEDInstrument.from_projection(self.earth_projection)

        # Create full SED earth instrument
        self.instruments[full_sed_earth_name] = FullSEDInstrument.from_projection(self.earth_projection)

    # -----------------------------------------------------------------

    def create_faceon_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the face-on instruments ...")

        # Create simple faceon instrument
        self.instruments[simple_faceon_name] = SimpleInstrument.from_projection(self.faceon_projection)

        # Full
        self.instruments[full_faceon_name] = FullInstrument.from_projection(self.faceon_projection, **self.full_instrument_properties)

    # -----------------------------------------------------------------

    def create_edgeon_instruments(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on instruments ...")

        # Create simple edgeon instrument
        self.instruments[simple_edgeon_name] = SimpleInstrument.from_projection(self.edgeon_projection)

        # Full
        self.instruments[full_edgeon_name] = FullInstrument.from_projection(self.edgeon_projection, **self.full_instrument_properties)

    # -----------------------------------------------------------------

    @property
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    @property
    def earth_instrument_properties(self):

        """
        This function ...
        :return:
        """

        return None

    # -----------------------------------------------------------------

    @property
    def full_instrument_properties(self):

        """
        This function ...
        :return:
        """

        properties = dict()
        properties["scattering_levels"] = 0 # no scattering levels
        #properties["counts"] = True # record photon counts (for Poisson noise) IS THIS WORKING?
        return properties

    # -----------------------------------------------------------------

    def get_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Load from fitting run
        if self.config.wavelength_grid is not None: self.load_wavelength_grid()

        # Create new
        else: self.create_wavelength_grid()

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the wavelength grid ...")

        # Check whether from fitting run
        if not self.from_fitting_run: raise RuntimeError("Not from fitting run: cannot get wavelength grid")

        # Set filepath
        filepath = fs.join(self.fitting_run.wavelength_grids_path, self.config.wavelength_grid + ".dat")
        if not fs.is_file(filepath): raise ValueError("Wavelength grid '" + self.config.wavelength_grid + "' not found")

        # Load
        self.wavelength_grid = WavelengthGrid.from_skirt_input(filepath)

    # -----------------------------------------------------------------

    def get_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Load from representation
        if self.grid_from_representation: self.load_dust_grid()

        # Create new
        else: self.create_dust_grid()

    # -----------------------------------------------------------------

    @property
    def has_representation(self):
        return self.config.representation is not None

    # -----------------------------------------------------------------

    @property
    def grid_from_representation(self):
        return self.has_representation and self.config.grid_from_representation

    # -----------------------------------------------------------------

    @property
    def projections_from_representation(self):
        return self.has_representation and self.config.projections_from_representation

    # -----------------------------------------------------------------

    @lazyproperty
    def representation(self):

        """
        This function ...
        :return:
        """

        if not self.has_representation: raise RuntimeError("Representation name not specified")
        return self.model_suite.get_representation(self.config.representation)

    # -----------------------------------------------------------------

    def load_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust grid ...")

        # Set the dust grid
        self.dust_grid = self.representation.dust_grid

        # Copy the dust grid building output
        fs.copy_files_from_directory(self.representation.grid_path, self.analysis_run.dust_grid_build_path)

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

        # Set output tree grid file path
        builder.config.tree_path = self.analysis_run.dust_grid_tree_path

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

        # Write the model
        self.write_model()

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

        # Write
        self.ski.saveto(self.analysis_run.ski_file_path)

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
        self.projections[earth_name].saveto(self.analysis_run.earth_projection_path)

        # Write the faceon projection system
        self.projections[faceon_name].saveto(self.analysis_run.faceon_projection_path)

        # Write the edgeon projection system
        self.projections[edgeon_name].saveto(self.analysis_run.edgeon_projection_path)

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the instruments ...")

        # Write the earth instruments
        self.instruments[simple_earth_name].saveto(self.analysis_run.simple_earth_instrument_path)
        self.instruments[full_earth_name].saveto(self.analysis_run.full_earth_instrument_path)
        self.instruments[sed_earth_name].saveto(self.analysis_run.sed_earth_instrument_path)
        self.instruments[full_sed_earth_name].saveto(self.analysis_run.full_sed_earth_instrument_path)

        # Write the faceon instruments
        self.instruments[simple_faceon_name].saveto(self.analysis_run.simple_faceon_instrument_path)
        self.instruments[full_faceon_name].saveto(self.analysis_run.full_faceon_instrument_path)

        # Write the edgeon instruments
        self.instruments[simple_edgeon_name].saveto(self.analysis_run.simple_edgeon_instrument_path)
        self.instruments[full_edgeon_name].saveto(self.analysis_run.full_edgeon_instrument_path)

    # -----------------------------------------------------------------

    @property
    def model_definition_stellar_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.stellar_path

    # -----------------------------------------------------------------

    @property
    def model_definition_dust_path(self):

        """
        This function ...
        :return:
        """

        return self.definition.dust_path

    # -----------------------------------------------------------------

    def write_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model ...")

        # Create
        create_model_definition_in_path(self.model_name, self.analysis_run.model_path, self.model_definition_stellar_path,
                                        self.model_definition_dust_path, parameter_values=self.parameter_values)

# -----------------------------------------------------------------

def create_model_definition_in_path(model_name, model_path, from_stellar_path, from_dust_path, parameter_values=None):

    """
    This function ...
    :param model_name:
    :param model_path:
    :param from_stellar_path:
    :param from_dust_path:
    :param parameter_values:
    :return:
    """

    from ..fitting.configuration import set_definition_values

    # Set the model stellar and dust path
    stellar_path = fs.create_directory_in(model_path, "stellar")
    dust_path = fs.create_directory_in(model_path, "dust")

    # Copy the stellar component directories
    fs.copy_directories_from_directory(from_stellar_path, stellar_path)

    # Copy the dust component directories
    fs.copy_directories_from_directory(from_dust_path, dust_path)

    # Get the paths
    stellar_paths = fs.directories_in_path(stellar_path, returns="dict")
    dust_paths = fs.directories_in_path(dust_path, returns="dict")

    # Load the definition
    definition = ModelDefinition(model_name, model_path, stellar_paths, dust_paths)

    # Adjust parameters
    if parameter_values is not None: set_definition_values(definition, parameter_values)

    # Return the definition
    return definition

# -----------------------------------------------------------------
