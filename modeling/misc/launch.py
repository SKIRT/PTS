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
from .interface import ModelSimulationInterface, earth_name, faceon_name, edgeon_name
from ...core.tools.stringify import tostr
from ..build.dustgrid import DustGridBuilder
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools.utils import lazyproperty
from ...core.tools.serialization import write_dict
from ..basics.instruments import FrameInstrument, SimpleInstrument, FullInstrument, SEDInstrument
from ...core.launch.launcher import SKIRTLauncher
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools import numbers
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.simulation.grids import load_grid, bintree, octtree, cartesian
from ...core.simulation.grids import CartesianDustGrid, OctTreeDustGrid, BinaryTreeDustGrid
from ..config.launch_model import make_images, make_seds
from ...core.simulation.output import output_types as ot
from ...core.launch.options import AnalysisOptions

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

# -----------------------------------------------------------------

nwavelengths_rel_tolerance = 0.15  # 15 % more or fewer points
nwavelengths_abs_tolerance = 15    # 15 points

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

    @property
    def has_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return self.wavelength_grid is not None

    # -----------------------------------------------------------------

    @property
    def has_dust_grid(self):

        """
        Thisf unction ...
        :return:
        """

        return self.dust_grid is not None

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
        if not self.has_wavelength_grid: self.create_wavelength_grid()

        # 5. Create the dust grid
        if not self.has_dust_grid: self.create_dust_grid()

        # Load the deprojections
        self.load_deprojections()

        # Create the projections
        self.create_projections(make_edgeon=self.config.edgeon, make_faceon=self.config.faceon)

        # 6. Create the instruments
        self.create_instruments()

        # 7. Adapt ski file
        self.adapt_ski()

        # 8. Build the dust grid (to get tree file) (maybe not necessary since there is only one simulation performed?)
        if not self.has_dust_grid_tree: self.build_dust_grid()

        # 9. Set the input
        self.set_input()

        # 10. Write
        self.write()

        # 11. Launch (and analyse) the simulation
        self.launch()

    # -----------------------------------------------------------------

    @property
    def local(self):

        """
        This function ...
        :return:
        """

        return self.config.remote is None

    # -----------------------------------------------------------------

    @property
    def remote(self):

        """
        This function ...
        :return:
        """

        return self.config.remote is not None

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

        # Plotting and misc directories
        self.simulation_plot_path = fs.create_directory_in(self.simulation_path, "plot")
        self.simulation_misc_path = fs.create_directory_in(self.simulation_path, "misc")

        # Load the wavelength grid?
        if self.config.regenerate_wavelength_grid: self.remove_wavelength_grid()
        else: self.load_wavelength_grid()

        # Load dust grid?
        if self.config.regenerate_dust_grid: self.remove_dust_grid()
        else: self.load_dust_grid()

    # -----------------------------------------------------------------

    def remove_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Removing the wavelength grid ...")

        # Remove the wavelength grid file
        if fs.is_file(self.wavelength_grid_path): fs.remove_file(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    def remove_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Removing the dust grid ...")

        # Remove the dust grid file
        if fs.is_file(self.dust_grid_path): fs.remove_file(self.dust_grid_path)

        # Remove the build directory
        if fs.is_directory(self.dust_grid_build_path): fs.remove_directory(self.dust_grid_build_path)

    # -----------------------------------------------------------------

    def load_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the wavelength grid ...")

        # Load the wavelength grid
        self.wavelength_grid = WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

        # Check the number of wavelengths
        if not numbers.is_close(self.wavelength_grid.nwavelengths, self.config.wg.npoints, rtol=nwavelengths_rel_tolerance, atol=nwavelengths_abs_tolerance):
            raise RuntimeError("Wavelength grid has to be regenerated: the number of wavelength points in the existing grid differs too much from the configured value")

    # -----------------------------------------------------------------

    @lazyproperty
    def is_bintree(self):

        """
        This function ...
        :return:
        """

        return isinstance(self.dust_grid, BinaryTreeDustGrid)

    # -----------------------------------------------------------------

    @lazyproperty
    def is_octtree(self):

        """
        This function ...
        :return:
        """

        return isinstance(self.dust_grid, OctTreeDustGrid)

    # -----------------------------------------------------------------

    @lazyproperty
    def is_tree(self):

        """
        This function ...
        :return:
        """

        return self.is_bintree or self.is_octtree

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_min_level(self):

        """
        This function ...
        :return:
        """

        if not self.is_tree: raise ValueError("Not a tree dust grid")
        return self.dust_grid.min_level

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_max_level(self):

        """
        This function ...
        :return:
        """

        if not self.is_tree: raise ValueError("Not a tree dust grid")
        return self.dust_grid.max_level

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_max_mass_fraction(self):

        """
        This function ...
        :return:
        """

        if not self.is_tree: raise ValueError("Not a tree dust grid")
        return self.dust_grid.max_mass_fraction

    # -----------------------------------------------------------------

    def load_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the dust grid ...")

        # Load
        self.dust_grid = load_grid(self.dust_grid_path)

        # Check properties
        if self.config.dg.grid_type == bintree and not self.is_bintree:
            raise RuntimeError("Dust grid has to be regenerated: the existing dust grid is not of type '" + bintree + "'")
        if self.config.dg.grid_type == octtree and not self.is_octtree:
            raise RuntimeError("Dust grid has to be regenerated: the existing dust grid is not of type '" + octtree + "'")
        #if self.config.dg.scale !=
        if self.is_bintree and self.config.dg.bintree_min_level != self.dust_grid.min_level:
            raise RuntimeError("Dust grid has to be regenerated: the existing dust grid doesn't have a minimum level of " + str(self.config.dg.bintree_min_level))
        if self.is_octtree and self.config.dg.octtree_min_level != self.dust_grid.min_level:
            raise RuntimeError("Dust grid has to be regenerated: the existing dust grid doesn't have a minimum level of " + str(self.config.dg.octtree_min_level))
        if self.is_tree and self.config.dg.max_mass_fraction != self.dust_grid.max_mass_fraction:
            raise RuntimeError("Dust grid has to be regenerated: the existing dust grid doesn't have a maximum mass fraction of " + str(self.config.dg.max_mass_fraction))

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

    @property
    def has_parameter_values(self):

        """
        This function ...
        :return:
        """

        return len(self.parameter_values) > 0

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

        # Show the model parameters, if any are chosen
        if self.has_parameter_values:
            print("")
            print("Model parameter values:")
            print("")
            for label in self.parameter_values: print(" - " + label + ": " + tostr(self.parameter_values[label]))
            print("")
        else: log.info("Using the standard parameter values of the model")

    # -----------------------------------------------------------------

    @property
    def make_images(self):

        """
        This function ...
        :return:
        """

        return make_images in self.config.make

    # -----------------------------------------------------------------

    @property
    def make_seds(self):

        """
        This function ...
        :return:
        """

        return make_seds in self.config.make

    # -----------------------------------------------------------------

    @property
    def make_contributions(self):

        """
        This function ...
        :return:
        """

        return self.config.contributions

    # -----------------------------------------------------------------

    @property
    def instrument_class(self):

        """
        This function ...
        :return:
        """

        if self.make_images and self.make_seds:
            if self.make_contributions: return FullInstrument
            else: return SimpleInstrument
        elif self.make_images:
            if self.make_contributions: return FullInstrument
            else: return FrameInstrument
        elif self.make_seds:
            if self.make_contributions: return FullInstrument
            else: return SEDInstrument
        else: raise ValueError("Either SEDs or images should be made")

    # -----------------------------------------------------------------

    def adapt_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the ski file ...")

        # Debugging
        log.debug("Setting the number of photon packages to " + str(self.config.npackages) + " ...")

        # Set the number of photon packages per wavelength
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

            # Set a tree dust grid for the ski file
            self.ski.set_filetree_dust_grid(dustgridtree_filename)

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

        return fs.join(self.projections_path, earth_name + ".proj")

    # -----------------------------------------------------------------

    @property
    def faceon_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, faceon_name + ".proj")

    # -----------------------------------------------------------------

    @property
    def edgeon_projection_path(self):

        """
        This function ..
        :return:
        """

        return fs.join(self.projections_path, edgeon_name + ".proj")

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections[earth_name].saveto(self.earth_projection_path)

        # Write the faceon projection system
        if self.has_faceon_projection: self.projections[faceon_name].saveto(self.faceon_projection_path)

        # Write the edgeon projection system
        if self.has_edgeon_projection: self.projections[edgeon_name].saveto(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @property
    def earth_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, earth_name + ".instr")

    # -----------------------------------------------------------------

    @property
    def faceon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, faceon_name + ".instr")

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, edgeon_name + ".instr")

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the instruments ...")

        # Write the SED instrument
        self.instruments[earth_name].saveto(self.earth_instrument_path)

        # Write the frame instrument
        if self.has_faceon_instrument: self.instruments[faceon_name].saveto(self.faceon_instrument_path)

        # Write the simple instrument
        if self.has_edgeon_instrument: self.instruments[edgeon_name].saveto(self.edgeon_instrument_path)

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

        # Set remote
        self.launcher.config.remote = self.config.remote
        self.launcher.config.attached = self.config.attached

        # Set options
        self.launcher.config.show_progress = True
        self.launcher.config.debug_output = True

        # Set retrieval options
        retrieve_types = []
        if self.make_images:
            retrieve_types.append(ot.total_images)
            if self.make_contributions:
                retrieve_types.append(ot.direct_images)
                retrieve_types.append(ot.transparent_images)
                retrieve_types.append(ot.scattered_images)
                retrieve_types.append(ot.dust_images)
                retrieve_types.append(ot.dust_scattered_images)
        if self.make_seds: retrieve_types.append(ot.seds)
        self.launcher.config.retrieve_types = retrieve_types

        # Set options for parallelization
        if self.local:
            parallelization = None
            nprocesses = self.config.nprocesses_local
            self.launcher.config.data_parallel_local = self.config.data_parallel_local
        elif self.remote:
            parallelization = None
            nprocesses = self.config.nprocesses_remote
            self.launcher.config.data_parallel_remote = self.config.data_parallel_remote
        else: raise RuntimeError("Error")

        # Set analysis options
        analysis_options = AnalysisOptions()
        analysis_options.misc.spectral_convolution = self.config.spectral_convolution
        if self.make_seds: analysis_options.misc.fluxes = True
        if self.make_images: analysis_options.misc.images = True
        analysis_options.misc.observation_filters = self.observed_filter_names
        instruments = [earth_name]
        if self.config.faceon: instruments.append(faceon_name)
        if self.config.edgeon: instruments.append(edgeon_name)
        analysis_options.misc.observation_instruments = instruments

        # Set analysis directories
        analysis_options.plotting.path = self.simulation_plot_path
        analysis_options.misc.path = self.simulation_misc_path

        # Other settings
        if log.is_debug(): self.launcher.config.show = True

        # Run the simulation
        self.launcher.run(definition=definition, analysis_options=analysis_options, parallelization=parallelization, nprocesses=nprocesses)

# -----------------------------------------------------------------
