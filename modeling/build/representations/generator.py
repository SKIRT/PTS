#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations.generator Contains the RepresentationsGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ..component import BuildComponent
from ...component.galaxy import GalaxyModelingComponent
from ....core.prep.dustgrids import DustGridGenerator
from ....core.basics.log import log
from ....core.basics.range import QuantityRange, RealRange
from .galaxy import GalaxyRepresentationBuilder
from ....core.tools import time
from ....core.tools import tables
from ....core.tools import filesystem as fs
from ....core.prep.templates import get_pan_template
from ....core.advanced.dustgridtool import generate_grid
from ....core.simulation.grids import load_grid

# -----------------------------------------------------------------

class RepresentationGenerator(BuildComponent, GalaxyModelingComponent):
    
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
        #super(RepresentationGenerator, self).__init__(*args, **kwargs)
        BuildComponent.__init__(self, *args, **kwargs)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The model definition
        self.definition = None

        # The dust grid generator
        self.dg_generator = None

        # A name for this representation generation event
        self.event_name = None

        # The representations
        self.representations = OrderedDict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the dust grids
        self.create_dust_grids()

        # 3. Build the representations
        self.build_representations()

        # 5. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(RepresentationGenerator, self).setup(**kwargs)
        BuildComponent.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Create the model definition
        self.definition = self.get_model_definition(self.config.model_name)

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

        # Set the event name
        self.event_name = time.unique_name("generator")

    # -----------------------------------------------------------------

    def create_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        semimajor_angular = self.truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        pixelscale_angular = self.definition.basic_maps_minimum_average_pixelscale.to("deg")
        #pixelscale_angular = self.reference_wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

        # BINTREE: (smallest_cell_pixels, min_level, max_mass_fraction)
        # Low-resolution: 10., 6, 1e-5
        # High-resolution: 0.5, 9, 0.5e-6

        # OCTTREE:
        # Low-resolution: 10., 2, 1e-5
        # High-resolution: 0.5, 3, 0.5e-6

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        min_scale = self.config.dg.scale_range.min * pixelscale
        max_scale = self.config.dg.scale_range.max * pixelscale
        scale_range = QuantityRange(min_scale, max_scale, invert=True)

        # The range of the max mass fraction
        mass_fraction_range = RealRange(self.config.dg.mass_fraction_range.min, self.config.dg.mass_fraction_range.max, invert=True) # must be inverted

        # Set fixed grid properties
        self.dg_generator.grid_type = self.config.dg.grid_type # set grid type
        self.dg_generator.x_radius = radius_physical
        self.dg_generator.y_radius = radius_physical
        self.dg_generator.z_radius = self.definition.dust_scaleheight * self.config.dg.scale_heights

        # Set options
        self.dg_generator.show = False
        self.dg_generator.write = False

        # Set the range of the minimum tree level
        if self.config.dg.grid_type == "bintree": level_range = self.config.dg.bintree_level_range # 6 to 9
        elif self.config.dg.grid_type == "octtree": level_range = self.config.dg.octtree_level_range # 2 to 3
        else: level_range = None

        # Generate the dust grids
        self.dg_generator.run(scale_range=scale_range, level_range=level_range, mass_fraction_range=mass_fraction_range, ngrids=self.config.nrepresentations)

    # -----------------------------------------------------------------

    def build_representations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the representations ...")

        # Loop over the dust grids
        for index, grid in enumerate(self.dg_generator.grids):

            # Create builder
            builder = GalaxyRepresentationBuilder(cwd=self.config.path)

            # Determine a name for this representation
            name = "grid" + str(index)

            # Set name
            builder.config.name = name

            # Set model name
            builder.config.model_name = self.config.model_name

            # Set option to calculate the quality of the dust grid
            builder.config.check_dust_grid_quality = self.config.check_dust_grid_quality

            # Build, passing the dust grid that has been created
            builder.run(dust_grid=grid)

            # Set the path for this representation
            self.representations[name] = builder.representation

    # -----------------------------------------------------------------

    #def generate_dust_grids(self):

        #"""
        #This function ...
        #:return:
        #"""

        ## Inform the user
        #log.info("Generating the dust grids ...")

        ## Loop over the representations
        #for name in self.representations:
            
        #    # Get info
        #    representation = self.representations[name]
        #    dust_grid_path = representation.dust_grid_path
        #    grid_out_path = representation.grid_out_path

        #    # Load the dust grid
        #    grid = load_grid(dust_grid_path)

        #    # Create ski file template
        #    ski = get_pan_template()

        #    # Add the dust grid geometry
        #    deprojection = self.definition.dust_deprojection

        #    # Create list of input paths
        #    input_paths = []

        #    # Add the model input paths
        #    input_paths.extend(self.definition.input_paths)

        #    # Generate the grid
        #    prefix = generate_grid(ski, grid, grid_out_path, input_paths)

        #    # If there is a tree data file, copy it to the main grid directory of the representation
        #    tree_out_path = fs.join(grid_out_path, prefix + "_ds_tree.dat")
        #    if fs.is_file(tree_out_path):

        #        # Debugging
        #        log.debug("Copying dust grid tree data for representation '" + name + "' ...")

        #        # Copy
        #        new_tree_path = fs.copy_file(tree_out_path, representation.grid_path, new_name="tree.dat")

        #    # No tree data
        #    else: log.debug("No tree data for representation '" + name + "'")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the dust grid table
        self.write_dust_grid_table()

    # -----------------------------------------------------------------

    def write_dust_grid_table(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Writing the dust grid table ...")

        # Determine the path
        path = fs.join(self.representations_path, self.event_name + "_dustgrids.dat")

        # Write the dust grids table
        tables.write(self.dg_generator.table, path)

# -----------------------------------------------------------------
