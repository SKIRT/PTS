#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.launch Contains the ModelLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from .interface import ModelSimulationInterface
from ...core.tools.stringify import tostr
from ..build.dustgrid import DustGridBuilder
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.utils import lazyproperty
from ...core.tools.serialization import write_dict
from ..basics.instruments import FrameInstrument, SimpleInstrument
from ...core.launch.launcher import SKIRTLauncher
from ...core.simulation.definition import SingleSimulationDefinition

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

# -----------------------------------------------------------------

class ModelLauncher(ModelSimulationInterface):

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
        super(ModelLauncher, self).__init__(*args, **kwargs)

        # Define paths
        self.simulation_path = None
        self.ski_path = None
        self.out_path = None
        self.dust_grid_path = None
        self.wavelength_grid_path = None
        self.dust_grid_build_path = None
        self.dust_grid_simulation_out_path = None
        self.dust_grid_tree_path = None
        self.projections_path = None
        self.instruments_path = None
        self.input_file_path = None

        # The SKIRT launcher
        self.launcher = SKIRTLauncher()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup class
        self.setup(**kwargs)

        # 2. Get the model
        self.get_model()

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

        # Build the dust grid (maybe not necessary since there is only one simulation perfomed?)
        self.build_dust_grid()

        # Set the input
        self.set_input()

        # Write
        self.write()

        # Launch the simulation
        self.launch()

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ModelLauncher, self).setup(**kwargs)

        # Set and create paths
        self.simulation_path = fs.create_directory_in(self.playground_path, self.simulation_name)
        self.ski_path = fs.join(self.simulation_path, self.galaxy_name + ".ski")
        self.out_path = fs.create_directory_in(self.simulation_path, "out")
        self.dust_grid_path = fs.join(self.simulation_path, "dust_grid.dg")
        self.wavelength_grid_path = fs.join(self.simulation_path, "wavelength_grid.dat")
        self.dust_grid_build_path = fs.create_directory_in(self.simulation_path, "dust grid")
        self.dust_grid_simulation_out_path = fs.create_directory_in(self.dust_grid_build_path, "out")
        self.dust_grid_tree_path = fs.join(self.dust_grid_build_path, "tree.dat")
        self.projections_path = fs.create_directory_in(self.simulation_path, "projections")
        self.instruments_path = fs.create_directory_in(self.simulation_path, "instruments")
        self.input_file_path = fs.join(self.simulation_path, "info.dat")

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

    def get_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the model ...")

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

    @property
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        #return FrameInstrument
        return SimpleInstrument

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
        builder.config.output = self.dust_grid_build_path

        # Set simulation path
        builder.config.simulation_path = self.dust_grid_simulation_out_path

        # Set the tree grid file path
        builder.config.tree_path = self.dust_grid_tree_path

        # Set whether quality has to be calculated
        builder.config.quality = self.config.check_dust_grid_quality

        # Run the dust grid builder
        builder.run(definition=self.definition, dust_grid=self.dust_grid)

    # -----------------------------------------------------------------

    @property
    def has_dust_grid_tree(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def use_file_tree_dust_grid(self):

        """
        This function ...
        :return:
        """

        smile = SKIRTSmileSchema()
        if not smile.supports_file_tree_grids: raise RuntimeError("A version of SKIRT that supports file tree grids is necessary")
        if not self.has_dust_grid_tree: raise RuntimeError("The dust grid tree is not present at '"  + self.dust_grid_tree_path + "'")
        return True

    # -----------------------------------------------------------------

    def set_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the simulation input ...")

        # NEW: DETERMINE AND SET THE PATH TO THE APPROPRIATE DUST GRID TREE FILE
        if self.use_file_tree_dust_grid:

            # self.simulation_input.add_file(self.representation.dust_grid_tree_path)
            self.input_paths[dustgridtree_filename] = self.dust_grid_tree_path

        # Determine and set the path to the appropriate wavelength grid file
        # wavelength_grid_path = self.fitting_run.wavelength_grid_path_for_level(self.generation_info.wavelength_grid_level)
        # self.simulation_input.add_file(wavelength_grid_path)
        self.input_paths[wavelengths_filename] = self.wavelength_grid_path

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

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Write
        self.ski.saveto(self.ski_path)

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

    @property
    def earth_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, "earth.proj")

    # -----------------------------------------------------------------

    @property
    def faceon_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, "faceon.proj")

    # -----------------------------------------------------------------

    @property
    def edgeon_projection_path(self):

        """
        This function ..
        :return:
        """

        return fs.join(self.projections_path, "edgeon.proj")

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections["earth"].saveto(self.earth_projection_path)

        # Write the faceon projection system
        self.projections["faceon"].saveto(self.faceon_projection_path)

        # Write the edgeon projection system
        self.projections["edgeon"].saveto(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @property
    def earth_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, "earth.instr")

    # -----------------------------------------------------------------

    @property
    def faceon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, "faceon.instr")

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, "edgeon.instr")

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the instruments ...")

        # Write the SED instrument
        self.instruments["earth"].saveto(self.earth_instrument_path)

        # Write the frame instrument
        self.instruments["faceon"].saveto(self.faceon_instrument_path)

        # Write the simple instrument
        self.instruments["edgeon"].saveto(self.edgeon_instrument_path)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching the simulation ...")

        # Create
        definition = SingleSimulationDefinition(self.ski_path, self.out_path, input_path=self.input_paths, name=self.simulation_name)

        analysis_options = None
        parallelization = None

        # Run the simulation
        self.launcher.run(definition=definition, analysis_options=analysis_options, parallelization=parallelization)

# -----------------------------------------------------------------
