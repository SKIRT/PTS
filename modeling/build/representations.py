#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations Contains the RepresentationsGenerator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ..component.galaxy import GalaxyModelingComponent
from ...core.prep.dustgrids import DustGridGenerator
from ...core.tools.logging import log
from ...core.basics.range import QuantityRange, RealRange
from ...core.units.parsing import parse_unit as u
from .representation import RepresentationBuilder
from ...core.tools import time
from ...core.tools import tables
from ...core.tools import filesystem as fs

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

        # 4. Writing
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
        self.dg_generator.z_radius = 3. * u("kpc")

        # Set options
        self.dg_generator.show = False
        self.dg_generator.write = False

        # Generate the dust grids
        self.dg_generator.run(scale_range=scale_range, level_range=self.config.dg.level_range, mass_fraction_range=mass_fraction_range, ngrids=self.config.nrepresentations)

    # -----------------------------------------------------------------

    def build_representations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the representations ...")

        # Loop over the dust grids

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
