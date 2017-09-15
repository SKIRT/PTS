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
from ..basics.instruments import FullInstrument
from ..misc.interface import ModelSimulationInterface, earth_name, edgeon_name, faceon_name

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

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

        # 6. Load the deprojections
        self.load_deprojections()

        # 7. Create the projections
        self.create_projections()

        # 8. Create the instruments
        self.create_instruments()

        # 9. Adapt ski file
        self.adapt_ski()

        # 10. Build the dust grid
        self.build_dust_grid()

        # 11. Set the input
        self.set_input()

        # 12. Write
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
        analysis_run_name = time.unique_name()

        # Create a directory for this analysis run
        analysis_run_path = fs.join(self.analysis_path, analysis_run_name)
        if not fs.is_directory(analysis_run_path): fs.create_directory(analysis_run_path)
        elif fs.is_empty(analysis_run_path, recursive=True): fs.clear_directory(analysis_run_path)
        else: raise ValueError("There already exists a directory for this analysis run")

        # Create the info object
        self.analysis_run_info = AnalysisRunInfo()

        # Set the analysis run name and path
        self.analysis_run_info.name = analysis_run_name
        self.analysis_run_info.path = analysis_run_path

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
        if self.from_model: self.prompt_model()

        # Prompt for a fitting run
        elif self.from_fitting_run: self.prompt_fitting()

        # Invalid
        else: raise ValueError("Invalid value for 'origin'")

        # Show the model parameters
        if self.from_fitting_run or self.config.adapt:
            print("")
            print("Adapted model parameter values:")
            print("")
            for label in self.parameter_values: print(" - " + label + ": " + tostr(self.parameter_values[label]))
            print("")

        # Show all model parameters
        else:
            print("")
            print("All model parameter values:")
            print("")
            for label in self.parameter_values: print(" - " + label + ": " + tostr(self.parameter_values[label]))
            print("")

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
        self.analysis_run = AnalysisRun(self.galaxy_name, self.analysis_run_info)

    # -----------------------------------------------------------------

    @property
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        return FullInstrument

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

        # Write the SED instrument
        self.instruments[earth_name].saveto(self.analysis_run.earth_instrument_path)

        # Write the frame instrument
        self.instruments[faceon_name].saveto(self.analysis_run.faceon_instrument_path)

        # Write the simple instrument
        self.instruments[edgeon_name].saveto(self.analysis_run.edgeon_instrument_path)

# -----------------------------------------------------------------
