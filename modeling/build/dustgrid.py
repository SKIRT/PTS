#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.dustgrid Contains the DustGridBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.basics.configurable import Configurable
from ...core.prep.smile import SKIRTSmileSchema
from ...core.launch.launcher import SKIRTLauncher
from .construct import add_dust_component
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.parallelization import Parallelization
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ...core.simulation.tree import DustGridTree
from ...core.simulation.logfile import LogFile
from ...core.tools import parsing
from ...core.basics.map import Map
from ...core.basics.configuration import save_mapping
from ...core.plot.grids import plotgrids

# -----------------------------------------------------------------

simulation_prefix = "dustgrid"
skifilename = simulation_prefix + ".ski"

# -----------------------------------------------------------------

gridxy_filename = simulation_prefix + "_ds_grhoxy.fits"
geometryxy_filename = simulation_prefix + "_ds_trhoxy.fits"
tree_filename = simulation_prefix + "_ds_tree.dat"
cell_properties_filename = simulation_prefix + "_ds_cellprops.dat"
convergence_filename = simulation_prefix + "_ds_convergence.dat"
quality_filename = simulation_prefix + "_ds_quality.dat"
log_filename = simulation_prefix + "_log.txt"

# -----------------------------------------------------------------

class DustGridBuilder(Configurable):
    
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
        super(DustGridBuilder, self).__init__(*args, **kwargs)

        # Smile
        self.smile = SKIRTSmileSchema()

        # Create the SKIRT launcher
        self.launcher = SKIRTLauncher()

        # The ski file template
        self.ski = None

        # The model
        self.definition = None

        # THe model representation
        #self.representation = None
        # The dust grid
        self.dust_grid = None

        # The dictionary of input map paths
        self.input_map_paths = dict()

        # Paths
        self.ski_path = None
        self.out_path = None

        # ...
        self.ratio = None
        self.mean_ratio = None
        self.median_ratio = None
        self.std = None


        # Quality measures
        self.quality = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the ski file
        self.create_ski()

        # Lauch the simulation
        self.launch()

        # Check
        self.get_quality()

        # 7. Writing
        if self.config.write: self.write()

        # Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DustGridBuilder, self).setup(**kwargs)

        # Determine path
        #self.ski_path = self.output_path_file(skifilename)
        # Determine output path
        #self.out_path = self.output_path #self.output_path_directory("out")

        # Determine paths
        self.ski_path = fs.join(self.config.simulation_path, skifilename)
        self.out_path = self.config.simulation_path

        # Get model definition and representation
        self.definition = kwargs.pop("definition")
        #self.representation = kwargs.pop("representation")
        self.dust_grid = kwargs.pop("dust_grid")

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Create template ski file
        #ski = self.smile.create_panchromatic_template()
        self.ski = self.smile.create_oligochromatic_template()

        # Set components
        self.set_components()

        #self.set_instruments()

        # Set other settings
        self.adjust_ski()

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # 1. Set stellar components
        #self.set_stellar_components()

        # 2. Set dust components
        self.set_dust_components()

    # -----------------------------------------------------------------

    def set_stellar_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar components ...")

    # -----------------------------------------------------------------

    def set_dust_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust components ...")

        # Loop over the dust components
        for name in self.definition.dust_component_names:

            # Load the component
            component = self.definition.get_dust_component(name)

            #print(name, component)

            # Add the dust component
            map_filename = add_dust_component(self.ski, name, component)

            # If map filename is defined, set path in dictionary
            if map_filename is not None: self.input_map_paths[map_filename] = component.map_path

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Remove the stellar system
        self.ski.remove_stellar_system()

        # Add the instrument
        #self.ski.add_instrument("earth", self.representation.sed_instrument)

        # Set the number of photon packages
        self.ski.setpackages(0)

        # Set the name of the wavelength grid file
        #self.ski.set_file_wavelength_grid("wavelengths.txt")
        # NO -> OLIGO

        # Set the dust emissivityex
        #self.set_dust_emissivity()

        # Set the dust grid
        self.ski.set_dust_grid(self.dust_grid)

        # Set all-cells dust library
        #self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        #self.set_selfabsorption()

        # Disable all writing options
        #self.ski.disable_all_writing_options()

        # Enable writing options
        self.ski.enable_all_writing_options()

        # Disable writing stellar density (we don't have a stellar system)
        self.ski.set_write_stellar_density(False)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Write the ski file
        self.write_ski()

        # Create simulation definition
        definition = SingleSimulationDefinition(self.ski_path, self.out_path, self.input_map_paths)

        # Determine parallelization scheme (do singleprocessing-
        ncores = 2
        nthreads_per_core = 2
        nprocesses = 1
        parallelization = Parallelization(ncores, nthreads_per_core, nprocesses)

        # Set settings
        self.launcher.config.progress_bar = True
        self.launcher.config.finish_after = "Writing dust cell properties" # finish after this line has been printed (when the next one comes)
        #self.launcher.config.finish_at = ""

        # Run
        self.launcher.run(definition=definition, parallelization=parallelization)

    # -----------------------------------------------------------------

    @property
    def simulation(self):

        """
        This function ...
        :return:
        """

        return self.launcher.simulation

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save
        self.ski.saveto(self.ski_path, fix=True)

    # -----------------------------------------------------------------

    @property
    def grid_xy_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.out_path, gridxy_filename)

    # -----------------------------------------------------------------

    @property
    def geometry_xy_path(self):
        
        """
        This function ...
        :return: 
        """

        return fs.join(self.out_path, geometryxy_filename)

    # -----------------------------------------------------------------

    def get_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the quality of the dust grid ...")

        # Get the projected quality
        self.get_projected_quality()

        # Get the optical depth quality
        self.get_optical_depth_quality()

        # Density
        self.get_density_quality()

        # Mass
        self.get_dust_mass_quality()

    # -----------------------------------------------------------------

    def get_projected_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the projected quality of the dust grid ...")

        # Load both maps
        gridxy = Frame.from_file(self.grid_xy_path)
        geometryxy = Frame.from_file(self.geometry_xy_path)

        # Determine ratio
        self.ratio = Frame(gridxy / geometryxy)
        self.mean_ratio = Frame.zeros_like(gridxy)
        self.median_ratio = Frame.zeros_like(gridxy)
        self.std = Frame.zeros_like(gridxy)

        # Loop over the unique values and their corresponding patches (masks)
        for value, where in gridxy.unique_values_and_masks:

            # Get original values
            original = geometryxy.data[where]

            # Caculate average
            mean = np.mean(original)
            median = np.median(original)
            std = np.std(original)

            # Print check
            #print(value, mean, median, std)

            self.mean_ratio[where] = mean / value
            self.median_ratio[where] = median / value
            self.std[where] = std

    # -----------------------------------------------------------------

    @property
    def cell_properties_path(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the cell properties file
        cellprops_path = fs.join(self.out_path, cell_properties_filename)
        return cellprops_path

    # -----------------------------------------------------------------

    @property
    def optical_depth_90(self):

        """
        This function ...
        :return:
        """

        # Get the optical depth for which 90% of the cells have a smaller value
        optical_depth = None
        for line in reversed(open(self.cell_properties_path).readlines()):
            if "of the cells have optical depth smaller than" in line:
                optical_depth = float(line.split("than: ")[1])
                break

        # Return the optical depth
        return optical_depth

    # -----------------------------------------------------------------

    def get_optical_depth_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the optical depth quality ...")

        # Debugging
        log.debug("90% of the cells have an optical depth smaller than " + str(self.optical_depth_90))

        # RERUNNING:

        # Adapt the maximal optical depth criterion
        #grid.max_optical_depth = optical_depth

        # Inform the user
        #log.info("Generating the high-resolution grid data ...")

        # Rerun the simulation
        #prefix = generate_grid(grid, out_path)
        #optical_depth = get_optical_depth_criterium(out_path, prefix)

        # Debugging
        #log.debug("90% of the cells have an optical depth smaller than " + str(optical_depth))

        # Create
        self.quality["optical_depth"] = Map()
        self.quality["optical_depth"]["tau90"] = self.optical_depth_90
        self.quality["optical_depth"]["mean"] = self.optical_depth_quality[0]
        self.quality["optical_depth"]["stddev"] = self.optical_depth_quality[1]

    # -----------------------------------------------------------------

    def get_density_quality(self):

        """
        This function ...
        :return:
        """

        # Inform theuser
        log.info("Getting the density quality ...")

        self.quality["density"] = Map()

        self.quality["density"]["mean"] = self.density_quality[0]
        self.quality["density"]["stddev"] = self.density_quality[1]
        self.quality["density"]["surface"] = Map()
        self.quality["density"]["surface"].x = Map()
        self.quality["density"]["surface"].x.expected = self.surface_density_convergence[0]
        self.quality["density"]["surface"].x.actual = self.surface_density_convergence[1]
        self.quality["density"]["surface"].y = Map()
        self.quality["density"]["surface"].y.expected = self.surface_density_convergence[2]
        self.quality["density"]["surface"].y.actual = self.surface_density_convergence[3]
        self.quality["density"]["surface"].z = Map()
        self.quality["density"]["surface"].z.expected = self.surface_density_convergence[4]
        self.quality["density"]["surface"].z.actual = self.surface_density_convergence[5]

    # -----------------------------------------------------------------

    def get_dust_mass_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the dust mass quality ...")

        self.quality["mass"] = Map()
        self.quality["mass"].expected = self.dust_mass_convergence[0]
        self.quality["mass"].actual = self.dust_mass_convergence[1]

    # -----------------------------------------------------------------

    @property
    def convergence_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.out_path, convergence_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def surface_density_convergence(self):

        """
        This function ...
        :return:
        """

        direction = None

        x_expected = x_actual = None
        y_expected = y_actual = None
        z_expected = z_actual = None

        for line in fs.read_lines(self.convergence_path):

            if "X-axis surface density" in line: direction = "x"
            elif "Y-axis surface density" in line: direction = "y"
            elif "Z-axis surface density" in line: direction = "z"
            elif "total dust mass" in line: direction = None # important

            if direction == "x":

                if "expected value" in line: x_expected = parsing.mass_surface_density_quantity(line.split("= ")[1])
                elif "actual value" in line: x_actual = parsing.mass_surface_density_quantity(line.split("= ")[1])

            elif direction == "y":

                if "expected value" in line: y_expected = parsing.mass_surface_density_quantity(line.split("= ")[1])
                elif "actual value" in line: y_actual = parsing.mass_surface_density_quantity(line.split("= ")[1])

            elif direction == "z":

                if "expected value" in line: z_expected = parsing.mass_surface_density_quantity(line.split("= ")[1])
                elif "actual value" in line: z_actual = parsing.mass_surface_density_quantity(line.split("= ")[1])

        # Return
        return x_expected, x_actual, y_expected, y_actual, z_expected, z_actual

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_mass_convergence(self):

        """
        This function ...
        :return:
        """

        triggered = False

        expected = actual = None

        for line in fs.read_lines(self.convergence_path):

            if "X-axis surface density" in line: continue
            elif "Y-axis surface density" in line: continue
            elif "Z-axis surface density" in line: continue
            elif "total dust mass" in line: triggered = True

            if triggered and "expected value" in line:
                expected = parsing.mass_quantity(line.split("= ")[1])

            elif triggered and "actual value" in line:
                actual = parsing.mass_quantity(line.split("= ")[1])

        # Return
        return expected, actual

    # -----------------------------------------------------------------

    @property
    def quality_path(self):

        """
        This fucntion ...
        :return:
        """

        return fs.join(self.out_path, quality_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def optical_depth_quality(self):

        """
        This function ...
        :return:
        """

        mean = None
        stddev = None

        # Search in the lines
        for line in fs.read_lines(self.quality_path):

            if "Mean value of optical depth delta" in line: mean = parsing.real(line.split(": ")[1])
            elif "Standard deviation of optical depth delta" in line: stddev = parsing.real(line.split(": ")[1])

        # Return
        return mean, stddev

    # -----------------------------------------------------------------

    @lazyproperty
    def density_quality(self):

        """
        This function ...
        :return:
        """

        mean = None
        stddev = None

        # Search in the lines
        for line in fs.read_lines(self.quality_path):

            if "Mean value of density delta" in line: mean = parsing.mass_density_quantity(line.split(": ")[1])
            elif "Standard deviation of density delta" in line: stddev = parsing.mass_density_quantity(line.split(": ")[1])

        # Return
        return mean, stddev

    # -----------------------------------------------------------------

    @property
    def log_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.out_path, log_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def log_file(self):

        """
        This function ...
        :return:
        """

        return LogFile(self.log_path)

    # -----------------------------------------------------------------

    @property
    def ncells(self):

        """
        This function ...
        :return:
        """

        return self.log_file.dust_cells_tree

    # -----------------------------------------------------------------

    @property
    def ntree_nodes(self):

        """
        This function ...
        :return:
        """

        return self.log_file.tree_nodes

    # -----------------------------------------------------------------

    @lazyproperty
    def tree_leaf_distribution(self):

        """
        This function ...
        :return:
        """

        return self.log_file.tree_leaf_distribution

    # -----------------------------------------------------------------

    @property
    def tree_levels(self):

        """
        This function ...
        :return:
        """

        return self.log_file.tree_levels

    # -----------------------------------------------------------------

    @property
    def tree_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.out_path, tree_filename)

    # -----------------------------------------------------------------

    @property
    def has_tree(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.tree_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def tree(self):

        """
        This function ...
        :return:
        """

        if not self.has_tree: return None
        return DustGridTree.from_file(self.tree_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write tree
        if self.has_tree: self.write_tree()

        # Write the ratios
        self.write_quality()

        # Write the cell distribution
        self.write_cell_distribution()

    # -----------------------------------------------------------------

    def write_tree(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid tree data ...")

        # Determine path
        path = self.output_path_file("tree.dat")

        # Save the tree
        #self.tree.saveto(path) # THIS CAN TAKE MUCH TOO LONG

        # MUCH QUICKER:
        fs.copy_file(self.tree_path, path)

    # -----------------------------------------------------------------

    def write_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the quality measures ...")

        # Write projected quality
        self.write_projected_quality()

        # Optical depth
        self.write_optical_depth_quality()

        # Density
        self.write_density_quality()

        # Dust mass
        self.write_dust_mass_quality()

    # -----------------------------------------------------------------

    def write_projected_quality(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Writing the quality maps ...")

        # Ratio
        path = self.output_path_file("ratio.fits")
        self.ratio.saveto(path)

        # Ratio of mean in each projected dust cell
        mean_path = self.output_path_file("ratio_mean.fits")
        self.mean_ratio.saveto(mean_path)

        # Ratio of median in each projected dust cell
        median_path = self.output_path_file("ratio_median.fits")
        self.median_ratio.saveto(median_path)

        # Standard deviation of theoretical density in each dust cell
        std_path = self.output_path_file("std.fits")
        self.std.saveto(std_path)

    # -----------------------------------------------------------------

    def write_optical_depth_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the optical depth quality ...")

        # Determine path
        path = self.output_path_file("optical_depth_quality.dat")

        # Save
        save_mapping(path, self.quality["optical_depth"])

    # -----------------------------------------------------------------

    def write_density_quality(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the density quality ...")

        # Determine path
        path = self.output_path_file("density_quality.dat")

        # Save
        save_mapping(path, self.quality["density"])

    # -----------------------------------------------------------------

    def write_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cell tree distribution ...")

        # Determine path
        path = self.output_path_file("tree_distribution.dat")

        # Save the distribution
        self.tree_leaf_distribution.saveto(path)

    # -----------------------------------------------------------------

    def write_dust_mass_quality(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust mass quality ...")

        # Determine path
        path = self.output_path_file("mass_quality.dat")

        # Save
        save_mapping(path, self.quality["mass"])

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the grid
        self.plot_grid()

        # Plot
        self.plot_dust_cell_distribution()

    # -----------------------------------------------------------------

    def plot_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the grid ...")

        # Plot the dust grid for the simulation
        plotgrids(self.simulation, output_path=self.config.output_path(), silent=(not log.is_debug()))

    # -----------------------------------------------------------------

    def plot_dust_cell_distribution(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the dust cell distribution ...")

        # Determine title and path
        title = "Dust cells in each tree level"
        path = self.output_path_file("cells_tree.pdf")

        # Plot
        self.tree_leaf_distribution.plot(title=title, path=path)

# -----------------------------------------------------------------
