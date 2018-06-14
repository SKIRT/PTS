#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.misc.contributions Contains the ModelContributionsLauncher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.launch.batch import BatchLauncher
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools.utils import lazyproperty
from ..build.dustgrid import DustGridBuilder
from .interface import ModelSimulationInterface, earth_name, faceon_name, edgeon_name
from ...core.simulation.grids import load_grid, bintree, octtree
from ...core.simulation.grids import OctTreeDustGrid, BinaryTreeDustGrid
from ...core.tools.serialization import write_dict
from ...core.tools import numbers
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ...core.tools.stringify import tostr
from ...core.prep.smile import SKIRTSmileSchema
from ..config.launch_contributions import make_images, make_seds
from ..basics.instruments import FrameInstrument, SimpleInstrument, FullInstrument, SEDInstrument
from ..core.model import bulge_component_name, disk_component_name, young_component_name, ionizing_component_name

# -----------------------------------------------------------------

contributions = ["total", "old", "young", "ionizing"]
component_names = {"old": [bulge_component_name, disk_component_name],
                   "young": [young_component_name],
                   "ionizing": [ionizing_component_name]}

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

# -----------------------------------------------------------------

nwavelengths_rel_tolerance = 0.15  # 15 % more or fewer points
nwavelengths_abs_tolerance = 15    # 15 points

# -----------------------------------------------------------------

class ModelContributionsLauncher(ModelSimulationInterface):
    
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
        super(ModelContributionsLauncher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # Define paths
        self.run_path = None
        self.dust_grid_path = None
        self.wavelength_grid_path = None
        self.dust_grid_build_path = None
        self.dust_grid_simulation_out_path = None
        self.dust_grid_tree_path = None
        self.projections_path = None
        self.instruments_path = None
        self.input_file_path = None

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

        # The ski files for simulating the contributions of the various stellar components
        self.ski_contributions = dict()

        # The paths to the ski files
        self.ski_paths = dict()

        # The simulation and output paths
        self.contributions_simulation_paths = dict()
        self.contributions_output_paths = dict()

        # The scheduling options for the different simulations (if using a remote host with scheduling system)
        self.scheduling_options = dict()

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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the model
        self.get_model()

        # 3. Create the wavelength grid
        if not self.has_wavelength_grid: self.create_wavelength_grid()

        # 4. Create the dust grid
        if not self.has_dust_grid: self.create_dust_grid()

        # 5. Load the deprojections
        self.load_deprojections()

        # 6. Create the projections
        self.create_projections()

        # 7. Create the instruments
        self.create_instruments()

        # 8. Adapt the ski file
        self.adapt_ski()

        # 9. Build the dust grid (to get tree file) (maybe not necessary since there is only one simulation performed?)
        if not self.has_dust_grid_tree: self.build_dust_grid()

        # 10. Set the input
        self.set_input()

        # 11. Write first, then launch the simulations
        self.write()

        # 12. Launch the simulations
        self.launch()

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
    def run_name(self):

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
        super(ModelContributionsLauncher, self).setup(**kwargs)
        #FittingComponent.setup(**kwargs)
        #ModelSimulationInterface.setup(**kwargs)

        # Set and create paths
        self.run_path = fs.create_directory_in(self.playground_path, self.run_name)

        self.simulations_path = fs.create_directory_in(self.run_path, "simulations")

        #self.ski_path = fs.join(self.run_path, self.galaxy_name + ".ski")
        #self.out_path = fs.create_directory_in(self.simulation_path, "out")

        #print(self.run_path)

        self.dust_grid_path = fs.join(self.run_path, "dust_grid.dg")
        self.wavelength_grid_path = fs.join(self.run_path, "wavelength_grid.dat")
        self.dust_grid_build_path = fs.create_directory_in(self.run_path, "dust grid")
        self.dust_grid_simulation_out_path = fs.create_directory_in(self.dust_grid_build_path, "out")
        self.dust_grid_tree_path = fs.join(self.dust_grid_build_path, "tree.dat")
        self.projections_path = fs.create_directory_in(self.run_path, "projections")
        self.instruments_path = fs.create_directory_in(self.run_path, "instruments")
        self.input_file_path = fs.join(self.run_path, "info.dat")

        # Load the wavelength grid?
        if self.config.regenerate_wavelength_grid: self.remove_wavelength_grid()
        else: self.load_wavelength_grid()

        # Load dust grid?
        if self.config.regenerate_dust_grid: self.remove_dust_grid()
        else: self.load_dust_grid()

        # Load the fitting run
        #self.fitting_run = self.load_fitting_run(self.config.fitting_run)

        # # Set the default option for the generation name
        # last_generation_name = self.fitting_run.last_finished_generation
        # if last_generation_name is None: last_generation_name = "first_guess"
        #
        # # Set the choices for the generationn name
        # generation_names = ["first_guess"] + self.fitting_run.finished_generations
        #
        # # Get the generation name
        # self.generation_name = prompt_string("generation", "name of the (finished) generation for which to relaunch the "
        #                                                    "best simulation",
        #                                      default=last_generation_name,
        #                                      choices=generation_names)

        # Set options for the batch launcher
        self.set_launcher_options()

        # Set the analysis options for the simulation of the total stellar contribution
        #self.set_analysis_options_total()

        # Create a directory for the simulation of the best model of the specified generation
        #self.best_generation_path = fs.join(self.fitting_run.best_path, self.generation_name)
        #if fs.is_directory(self.best_generation_path): raise RuntimeError("The best model has already been launched for this generation (" + self.generation_name + ")")
        #fs.create_directory(self.best_generation_path)

        # Set the paths to the wavelength and dust grid files
        #self.wavelength_grid_path = fs.join(self.best_generation_path, "wavelengths.txt")
        #self.dust_grid_path = fs.join(self.best_generation_path, "dust_grid.dg")

        # Set the paths to the simulation directories for the different contributions
        for contribution in contributions:

            # Create the directory
            simulation_path = fs.create_directory_in(self.simulations_path, contribution)

            # Add to the dictionary
            self.contributions_simulation_paths[contribution] = simulation_path

            # Create and set the output directory
            self.contributions_output_paths[contribution] = fs.create_directory_in(simulation_path, "out")

            # Set the ski path
            ski_path = fs.join(simulation_path, self.object_name + ".ski")
            self.ski_paths[contribution] = ski_path

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

    def remove_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Removing the dust grid ...")

        # Remove the dust grid file
        if fs.is_file(self.dust_grid_path): fs.remove_file(self.dust_grid_path)

        # Clear the build directory
        if fs.is_directory(self.dust_grid_build_path): fs.clear_directory(self.dust_grid_build_path, recursive=True)

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

            if self.wavelength_grid.nwavelengths > self.config.wg.npoints and self.config.wg.add_emission_lines:
                log.warning("Number of wavelength points differs much but assuming this is because emission lines were also added")
            else: raise RuntimeError("Wavelength grid has to be regenerated: the number of wavelength points in the existing grid differs too much from the configured value")

    # -----------------------------------------------------------------

    @property
    def remotes(self):

        """
        This function ...
        :return:
        """

        if self.config.remote is not None: return [self.config.remote]
        else: return []

    # -----------------------------------------------------------------

    def set_launcher_options(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting options for the batch simulation launcher ...")

        # Basic options
        self.launcher.config.shared_input = True                    # The input directories (or files) for the different simulations are shared
        self.launcher.config.remotes = self.remotes         # the remote host(s) on which to run the simulations
        self.launcher.config.group_simulations = self.config.group  # group multiple simulations into a single job # TODO: IMPLEMENT THIS
        self.launcher.config.group_walltime = self.config.walltime  # the preferred walltime for jobs of a group of simulations
        #self.launcher.config.timing_table_path = self.timing_table_path  # The path to the timing table file
        #self.launcher.config.memory_table_path = self.memory_table_path  # The path to the memory table file
        self.launcher.config.cores_per_process = self.config.cores_per_process # the number of cores per process, for non-schedulers
        #self.launcher.config.dry = self.config.dry

        # Attached
        #self.launcher.config.remote = self.config.remote
        self.launcher.config.attached = self.config.attached
        self.launcher.config.show_progress = True

        # Logging options
        self.launcher.config.logging.verbose = True          # verbose logging

        # Simulation analysis options

        ## General
        #self.launcher.config.analysis.relative = True

        self.launcher.config.relative_analysis_paths = True

        ## Extraction
        self.launcher.config.analysis.extraction.path = "extr"     # name for the extraction directory
        self.launcher.config.analysis.extraction.progress = True   # extract progress information
        self.launcher.config.analysis.extraction.timeline = True   # extract the simulation timeline
        self.launcher.config.analysis.extraction.memory = True     # extract memory usage information

        ## Plotting
        self.launcher.config.analysis.plotting.path = "plot"
        self.launcher.config.analysis.plotting.progress = True    # plot progress information
        self.launcher.config.analysis.plotting.timeline = True    # plot timeline
        self.launcher.config.analysis.plotting.memory = True      # plot memory usage
        self.launcher.config.analysis.plotting.seds = True        # plot the simulated SEDs
        self.launcher.config.analysis.plotting.reference_seds = [self.observed_sed_path]  # the path to the reference SED (for plotting the simulated SED against the reference points)
        self.launcher.config.analysis.plotting.format = "pdf"     # plot format

        ## Miscellaneous
        self.launcher.config.analysis.misc.path = "misc"          # The base directory where all of the simulations will have a seperate directory with the 'misc' analysis output
        #self.launcher.config.analysis.misc.images = True          # create observed images
        self.launcher.config.analysis.misc.fluxes = True          # calculate observed fluxes
        self.launcher.config.analysis.misc.observation_filters = self.observed_filter_names  # The filters for which to create the observations
        #self.launcher.config.analysis.misc.make_images_remote = self.config.images_remote
        #self.launcher.config.analysis.misc.images_wcs = self.reference_wcs_path
        #self.launcher.config.analysis.misc.images_unit = "MJy/sr"
        #self.launcher.config.analysis.misc.images_kernels = kernel_paths

    # -----------------------------------------------------------------

    # def set_analysis_options_total(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Setting the analysis options for the simulation of the total stellar contribution ...")
    #
    #     # Set the paths to the kernel for each image (except for the SPIRE images)
    #     # HERZIEN:
    #     # kernel_paths = dict()
    #     # aniano = AnianoKernels()
    #     # pacs_red_psf_path = aniano.get_psf_path(self.pacs_red_filter)
    #     # for filter_name in self.observed_filter_names:
    #     #     if "SPIRE" in filter_name: continue
    #     #     kernel_paths[filter_name] = pacs_red_psf_path
    #
    #     # Set the analysis options that are different from those used for the simulations of the other stellar contributions
    #     self.analysis_options_total["misc"] = dict()
    #     self.analysis_options_total["misc"]["images"] = True
    #     self.analysis_options_total["misc"]["make_images_remote"] = self.config.images_remote
    #     self.analysis_options_total["misc"]["images_wcs"] = self.reference_wcs_path
    #     self.analysis_options_total["misc"]["images_unit"] = "MJy/sr"
    #     self.analysis_options_total["misc"]["images_kernels"] = kernel_paths

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

    # def load_ski(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    # def create_wavelength_grid(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Creating the wavelength grid ...")
    #
    #     # Create the emission lines instance
    #     emission_lines = EmissionLines()
    #
    #     # Fixed wavelengths in the grid
    #     fixed = [self.i1_filter.pivot, self.fuv_filter.pivot]
    #
    #     # Create the grid
    #     # grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added
    #     grid, subgrid_npoints, emission_npoints, fixed_npoints, broad_resampled, narrow_added = create_one_subgrid_wavelength_grid(self.config.nwavelengths, emission_lines, fixed)
    #
    #     # Set the grid
    #     self.wavelength_grid = grid

    # -----------------------------------------------------------------

    # def create_dust_grid(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Creating the dust grid ...")
    #
    #     # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
    #     semimajor_angular = self.truncation_ellipse.semimajor  # major axis length of the sky ellipse
    #     radius_physical = (semimajor_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())
    #
    #     # Get the pixelscale in physical units
    #     distance = self.galaxy_properties.distance
    #     pixelscale_angular = self.reference_wcs.average_pixelscale  # in deg
    #     pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())
    #
    #     x_radius = radius_physical
    #     y_radius = radius_physical
    #     z_radius = 3. * Unit("kpc")
    #
    #     x_min = - x_radius
    #     x_max = x_radius
    #     y_min = - y_radius
    #     y_max = y_radius
    #     z_min = - z_radius
    #     z_max = z_radius
    #
    #     x_extent = x_max - x_min
    #     y_extent = y_max - y_min
    #     z_extent = z_max - z_min
    #
    #     # Set the scale
    #     scale = self.config.dg.rel_scale * pixelscale
    #
    #     # Create the grid
    #     grid = create_one_dust_grid(self.config.dg.grid_type, scale, x_extent, y_extent, z_extent, x_min, x_max, y_min, y_max, z_min, z_max, self.config.dg.min_level, self.config.dg.max_mass_fraction)
    #
    #     # Set the grid
    #     self.dust_grid = grid

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

    @property
    def earth_instrument_properties(self):

        """
        This function ...
        :return:
        """

        # Full instrument
        if self.instrument_class == FullInstrument:

            properties = dict()
            properties["scattering_levels"] = self.config.scattering_levels
            properties["counts"] = self.config.counts
            return properties

        # Not full instrument
        else:

            # Check settings
            if self.config.scattering_levels > 0: raise ValueError("Cannot have more scattering levels when instrument is not FullInstrument")
            if self.config.counts: raise ValueError("Cannot record photon counts when instrument is not FullInstrument")

            # Return empty settings
            return dict()

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Create
        self.create_projection_systems(make_edgeon=self.config.edgeon, make_faceon=self.config.faceon)

    # -----------------------------------------------------------------

    # def set_input(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Setting the input paths ...")
    #
    #     # Set the paths to the input maps
    #     self.input_paths = self.input_map_paths
    #
    #     # Determine and set the path to the wavelength grid file
    #     self.input_paths.append(self.wavelength_grid_path)

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

    # FOR FROM FITTING RUN:
    # def get_parameter_values(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Getting the parameter values ...")
    #
    #     # If the first guess model should be used
    #     if self.generation_name == "first_guess":
    #
    #         # Inform the user
    #         log.info("Using the parameter values from the initial guess model ...")
    #
    #         # Get the parameter values
    #         self.parameter_values = self.fitting_run.first_guess_parameter_values
    #
    #     # If the best simulation of a generation has to be used
    #     else:
    #
    #         # Inform the user
    #         log.info("Using the parameter values from the best model in the '" + self.generation_name + "' generation ...")
    #
    #         # Get the parameter values
    #         self.parameter_values = self.fitting_run.best_parameter_values_for_generation(self.generation_name)
    #
    #     # Debugging
    #     log.debug("The parameter values used for the simulations are: ")
    #     for label in self.parameter_values: log.debug(" - " + label + ": " + str(self.parameter_values[label]))

    # -----------------------------------------------------------------

    def adapt_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the ski template ...")

        # 1. Set the best parameter values
        #self.set_parameter_values()

        # 2. Set basic properties
        self.set_properties()

        # 3. Adjust the ski files for simulating the contributions of the various stellar components
        self.adjust_ski_contributions()

    # -----------------------------------------------------------------

    # def set_parameter_values(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Inform the user
    #     log.info("Setting the best parameter values ...")
    #
    #     # Assert that the parameter values are quantities
    #     for label in self.parameter_values: assert hasattr(self.parameter_values[label], "unit")
    #
    #     # Set the parameter values in the ski file template
    #     self.ski.set_labeled_values(self.parameter_values)

    # -----------------------------------------------------------------

    def set_properties(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting ski file properties ...")

        # Debugging
        log.debug("Disabling all writing settings ...")

        # Disable all writing settings
        self.ski.disable_all_writing_options()

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

        # Debugging
        log.debug("Setting the wavelength grid file ...")

        # Set wavelength grid for ski file
        self.ski.set_file_wavelength_grid(wavelengths_filename)

        # # Debugging
        # log.debug("Setting the dust grid ...")
        #
        # # Set the dust grid
        # self.ski.set_dust_grid(self.dust_grid)

        # Set the instruments
        self.set_instruments()

    # -----------------------------------------------------------------

    def set_instruments(self):

        """
        This function ...
        :return:
        """

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: self.ski.add_instrument(name, self.instruments[name])

    # -----------------------------------------------------------------

    def adjust_ski_contributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting ski files for simulating the contribution of the various stellar components ...")

        # Loop over the different contributions, create seperate ski file instance
        for contribution in contributions:

            # Create a copy of the ski file instance
            ski = self.ski.copy()

            # Adjust the ski file for generating simulated images
            if contribution == "total":

                # Debugging
                log.debug("Adjusting ski file for generating simulated images ...")

                # Remove all instruments
                #ski.remove_all_instruments()

                # Add the simple instrument
                #ski.add_instrument("earth", self.simple_instrument)

            # Remove other stellar components
            else:

                # Debugging
                log.debug("Adjusting ski file for simulating the contribution of the " + contribution + " stellar population ...")
                log.debug("Removing all stellar components other than: " + ", ".join(component_names[contribution]) + " ...")

                # Remove the other components
                ski.remove_stellar_components_except(component_names[contribution])

            # Add the ski file to the dictionary
            self.ski_contributions[contribution] = ski

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the ski file
        self.write_ski_files()

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

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files for simulating the contribution of the various stellar components ...")

        # Loop over the ski files
        for contribution in self.ski_contributions: self.ski_contributions[contribution].saveto(self.ski_paths[contribution])

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
        log.info("Writing the wavelength grid ...")

        # Write the wavelength grid for SKIRT input
        self.wavelength_grid.to_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid ...")

        # Write the dust grid
        self.dust_grid.saveto(self.dust_grid_path)

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

        # Set the path to the directory to contain the launch script (job script);
        # just use the directory created for the generation
        if self.remote: self.launcher.set_script_path(self.remote_host_id, self.run_path)

        # Enable screen output logging if using a remote host without a scheduling system for jobs
        if self.remote and not self.uses_scheduler: self.launcher.enable_screen_output(self.remote_host_id)

        #print("INPUT PATHS:", self.input_paths)

        # Loop over the ski paths for the different contributions
        for contribution in self.ski_paths:

            # Debugging
            log.debug("Initiating queuing of simulation for the contribution of the " + contribution + " stellar component ...")

            # Get ski path and output path
            ski_path = self.ski_paths[contribution]
            output_path = self.contributions_output_paths[contribution]

            # Set the simulation name
            #simulation_name = self.generation_name.replace(" ", "") + "_" + contribution
            simulation_name = self.run_name + "_" + contribution

            # Create the SKIRT simulation definition
            definition = SingleSimulationDefinition(ski_path, output_path, self.input_paths)

            # Debugging
            log.debug("Adding a simulation to the queue with:")
            log.debug(" - name: " + simulation_name)
            log.debug(" - ski path: " + definition.ski_path)
            log.debug(" - output path: " + definition.output_path)

            # Set special analysis options for the simulation of the total stellar contribution
            #analysis_options = self.analysis_options_total if contribution == "total" else None
            analysis_options = None

            # Put the parameters in the queue and get the simulation object
            self.launcher.add_to_queue(definition, simulation_name, analysis_options=analysis_options)

            # Set scheduling options for this simulation
            if self.uses_scheduler: self.launcher.set_scheduling_options(self.remote_host_id, simulation_name, self.scheduling_options[contribution])

        # Run the launcher, schedules the simulations
        self.launcher.run()

        # Loop over the scheduled simulations
        for simulation in self.launcher.launched_simulations:

            # Add the path to the modeling directory to the simulation object
            simulation.analysis.modeling_path = self.config.path

            # Set the path to the fit model analyser
            #analyser_path = "pts.modeling.fitting. ..."

            # Add the analyser class path
            #simulation.add_analyser(analyser_path)

            # Save the simulation object
            simulation.save()

# -----------------------------------------------------------------
