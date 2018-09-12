#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representation Contains the RepresentationBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument, load_instrument
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection, load_projection
from ...core.simulation.grids import load_grid
from ...core.simulation.grids import FileTreeDustGrid
from ...core.simulation.tree import DustGridTree
from ...core.tools.utils import lazyproperty
from ...core.simulation.tree import DustGridTreeDistribution
from ...magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

class Representation(object):

    """
    This class ...
    """

    def __init__(self, name, model_name, path):

        """
        This function ...
        :param name:
        :param model_name:
        :param path:
        """

        # General properties
        self.name = name
        self.model_name = model_name
        self.path = path

        # Directories of the representation
        self.projections_path = fs.create_directory_in(self.path, "projections")
        self.instruments_path = fs.create_directory_in(self.path, "instruments")
        self.grid_path = fs.create_directory_in(self.path, "grid")

        # Individual projection paths
        self.earth_projection_path = fs.join(self.projections_path, "earth.proj")
        self.edgeon_projection_path = fs.join(self.projections_path, "edgeon.proj")
        self.faceon_projection_path = fs.join(self.projections_path, "faceon.proj")

        # Individual instrument paths
        self.sed_instrument_path = fs.join(self.instruments_path, "sed.instr")
        self.frame_instrument_path = fs.join(self.instruments_path, "frame.instr")
        self.simple_instrument_path = fs.join(self.instruments_path, "simple.instr")

        # Dust grid file path
        self.dust_grid_path = fs.join(self.grid_path, "dust_grid.dg")

        # Dust grid SKIRT output path
        self.grid_out_path = fs.create_directory_in(self.grid_path, "out")

        # Dust grid tree path
        self.dust_grid_tree_path = fs.join(self.grid_path, "tree.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):
        return DustGridTree.from_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    def create_file_tree_dust_grid(self, search_method="Neighbor", write=False):

        """
        This function ...
        :param search_method:
        :param write:
        :return: 
        """

        grid = FileTreeDustGrid(filename=self.dust_grid_tree_path, search_method=search_method, write=write)
        return grid

    # -----------------------------------------------------------------

    @lazyproperty
    def has_dust_grid_tree(self):
        return fs.is_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixelscale(self):
        return self.earth_projection.pixelscale

    # -----------------------------------------------------------------

    def get_projection_paths(self):
        return fs.files_in_path(self.projections_path, extension="proj", returns="dict")

    # -----------------------------------------------------------------

    def get_projections(self):

        """
        This function ...
        :return:
        """

        projections = OrderedDict()
        paths = self.get_projection_paths()
        for name in paths: projections[name] = load_projection(paths[name])
        return projections

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_projection(self):
        return GalaxyProjection.from_file(self.earth_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):
        return EdgeOnProjection.from_file(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):
        return FaceOnProjection.from_file(self.faceon_projection_path)

    # -----------------------------------------------------------------

    @property
    def earth_pixelscale(self):
        return self.earth_projection.pixelscale

    # -----------------------------------------------------------------

    @property
    def edgeon_pixelscale(self):
        return self.edgeon_projection.pixelscale

    # -----------------------------------------------------------------

    @property
    def faceon_pixelscale(self):
        return self.faceon_projection.pixelscale

    # -----------------------------------------------------------------

    def get_instrument_paths(self):

        """
        This function ...
        :return:
        """

        return fs.files_in_path(self.instruments_path, extension="instr", returns="dict")

    # -----------------------------------------------------------------

    def get_instruments(self):

        """
        Thisn function ...
        :return:
        """

        instruments = OrderedDict()
        paths = self.get_instrument_paths()
        for name in paths: instruments[name] = load_instrument(paths[name])
        return instruments

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):
        return SEDInstrument.from_file(self.sed_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument(self):
        return FrameInstrument.from_file(self.frame_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_instrument(self):
        return SimpleInstrument.from_file(self.simple_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):
        return load_grid(self.dust_grid_path)

    # -----------------------------------------------------------------

    @property
    def dust_grid_type(self):
        return type(self.dust_grid).__name__

    # -----------------------------------------------------------------

    @property
    def dust_grid_tree_distribution_path(self):
        return fs.join(self.grid_path, "tree_distribution.dat")

    # -----------------------------------------------------------------

    @property
    def has_dust_grid_tree_distribution(self):
        return fs.is_file(self.dust_grid_tree_distribution_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree_distribution(self):
        return DustGridTreeDistribution.from_file(self.dust_grid_tree_distribution_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_min_level(self):
        return self.dust_grid_tree_distribution.min_level

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_max_level(self):
        return self.dust_grid_tree_distribution.max_level

    # -----------------------------------------------------------------

    @lazyproperty
    def ndust_cells(self):
        return self.dust_grid_tree_distribution.ncells

    # -----------------------------------------------------------------

    @lazyproperty
    def properties(self):

        """
        This function ...
        :return:
        """

        properties = dict()
        properties["pixelscale"] = self.pixelscale
        properties["dust_grid_type"] = self.dust_grid_type
        if self.has_dust_grid_tree_distribution:
            properties["dust_grid_min_level"] = self.dust_grid_min_level
            properties["dust_grid_max_level"] = self.dust_grid_max_level
        properties["dust_grid_ncells"] = self.ndust_cells

        return properties

    # -----------------------------------------------------------------

    # BELOW THIS: TOO MUCH OF A HACK, BUT AT THE MOMENT THE REFERENCE DEPROJECTION NAME IS SAVED NOWHERE, AND WE NEED THIS INFO

    @property
    def representations_path(self):
        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def build_path(self):
        return fs.directory_of(self.representations_path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):
        return fs.directory_of(self.build_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_suite(self):
        from .suite import ModelSuite
        return ModelSuite.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojections(self):

        """
        This function ...
        :return:
        """

        from .suite import get_dust_deprojections, get_stellar_deprojections

        # Get deprojections
        deprojections = dict()

        # Load the stellar deprojections
        get_stellar_deprojections(self.model_suite, self.model_name, deprojections=deprojections)

        # Load the dust deprojections
        get_dust_deprojections(self.model_suite, self.model_name, deprojections=deprojections)

        # Return the deprojections
        return deprojections

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):
        return self.earth_projection.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection_name(self):

        """
        This function ...
        :return:
        """

        azimuth = 0.0

        matching_names = []

        # Loop over the deprojections
        for component_name in self.deprojections:

            # Get the deprojection
            deprojection = self.deprojections[component_name]

            # Create the 'earth' projection system
            earth_projection = GalaxyProjection.from_deprojection(deprojection, self.galaxy_distance, azimuth)

            # Check if equivalent
            if self.earth_projection.isclose(earth_projection, rtol=1e-4): matching_names.append(component_name)

        # Check the result
        if len(matching_names) > 1: raise ValueError("Multiple possible deprojections")
        elif len(matching_names) == 0: raise ValueError("No matching deprojection")
        else: return matching_names[0]

    # -----------------------------------------------------------------

    @property
    def reference_deprojection(self):
        return self.deprojections[self.reference_deprojection_name]

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map_path(self):
        return self.reference_deprojection.filepath

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map_wcs(self):
        return CoordinateSystem.from_file(self.reference_map_path)

# -----------------------------------------------------------------
