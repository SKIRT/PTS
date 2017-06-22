#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.imagesrepresentation Contains the ImagesRepresentationBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

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

# -----------------------------------------------------------------

simulation_prefix = "dustgrid"
skifilename = simulation_prefix + ".ski"

# -----------------------------------------------------------------

gridxy_filename = simulation_prefix + "_ds_grhoxy.fits"
geometryxy_filename = simulation_prefix + "_ds_trhoxy.fits"
tree_filename = simulation_prefix + "_ds_tree.dat"

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
        self.check()

        # 7. Writing
        if self.config.write: self.write()

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

    def check(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the dust grid ...")

        # Load both maps
        gridxy = Frame.from_file(self.grid_xy_path)
        geometryxy = Frame.from_file(self.geometry_xy_path)

        # Determine ratio
        self.ratio = Frame(gridxy / geometryxy)

        self.mean_ratio = Frame.zeros_like(gridxy)
        self.median_ratio = Frame.zeros_like(gridxy)

        self.std = Frame.zeros_like(gridxy)

        # Loop over the unique values in the gridded data
        values = np.unique(gridxy.data)
        for value in values:

            # Check in which pixels this value (get patch each time)
            where = gridxy.where(gridxy) # returns mask

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

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write tree
        if self.has_tree: self.write_tree()

        self.write_ratios()

    # -----------------------------------------------------------------

    def write_tree(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid tree data ...")

        # Copy
        fs.copy_file(self.tree_path, self.output_path_file("tree.dat"))

    # -----------------------------------------------------------------

    def write_ratios(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ratio (quality?) maps ...")

        path = self.output_path_file("ratio.fits")
        self.ratio.saveto(path)

        mean_path = self.output_path_file("ratio_mean.fits")
        self.mean_ratio.saveto(mean_path)

        median_path = self.output_path_file("ratio_median.fits")
        self.median_ratio.saveto(median_path)

        std_path = self.output_path_file("std.fits")
        self.std.saveto(std_path)

# -----------------------------------------------------------------
