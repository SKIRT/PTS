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
import math

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools import filesystem as fs
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ...core.tools.logging import log
from ..build.component import get_stellar_component_names, get_dust_component_names
from ..build.component import load_stellar_component_deprojection, load_dust_component_deprojection
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...core.basics.configuration import prompt_string, prompt_yn, prompt_real
from ...core.units.stringify import represent_quantity
from ...core.simulation.grids import load_grid
from ..component.galaxy import GalaxyModelingComponent
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection, smallest_scale_for_dust_grid
from ...core.simulation.grids import FileTreeDustGrid
from ...core.simulation.tree import DustGridTree
from .dustgrid import DustGridBuilder

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

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree_path(self):

        """
        This function ...
        :return: 
        """

        path = fs.join(self.grid_path, "tree.dat")
        return path #fs.is_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):

        """
        This function ...
        :return: 
        """

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

        """
        This function ...
        :return: 
        """

        return fs.is_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.earth_projection.pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return GalaxyProjection.from_file(self.earth_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        return EdgeOnProjection.from_file(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        return FaceOnProjection.from_file(self.faceon_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_instrument(self):

        """
        This function ...
        :return:
        """

        return SEDInstrument.from_file(self.sed_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument(self):

        """
        This function ...
        :return:
        """

        return FrameInstrument.from_file(self.frame_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_instrument(self):

        """
        This function ...
        :return:
        """

        return SimpleInstrument.from_file(self.simple_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):

        """
        This function ...
        :return:
        """

        return load_grid(self.dust_grid_path)

# -----------------------------------------------------------------

class RepresentationBuilderBase(BuildComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RepresentationBuilderBase, self).__init__(*args, **kwargs)

        # The model definition
        self.definition = None

        # The representation
        self.representation = None

        # The projections
        self.projections = dict()

        # The instruments
        self.instruments = dict()

        # The dust grid
        self.dust_grid = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RepresentationBuilderBase, self).setup(**kwargs)

        # Create the model definition
        self.definition = self.get_model_definition(self.config.model_name)

        # Create the representation
        path = fs.create_directory_in(self.representations_path, self.config.name)
        self.representation = Representation(self.config.name, self.config.model_name, path)

    # -----------------------------------------------------------------

    @property
    def representation_name(self):

        """
        This function ...
        :return:
        """

        return self.representation.name

    # -----------------------------------------------------------------

    @property
    def representation_path(self):

        """
        This function ...
        :return:
        """

        return self.representation.path

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.representation.model_name

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
        builder.config.output = self.representation.grid_out_path

        # Run the builder
        builder.run(definition=self.definition, representation=self.representation)

# -----------------------------------------------------------------

class RepresentationBuilder(RepresentationBuilderBase, GalaxyModelingComponent):
    
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
        #super(RepresentationBuilder, self).__init__(*args, **kwargs)
        RepresentationBuilderBase.__init__(self, *args, **kwargs)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The deprojections
        self.deprojections = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the deprojections
        self.load_deprojections()

        # 6. Create the dust grid
        if self.dust_grid is None: self.create_dust_grid()

        # 4. Create the projections
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # Build the dust grid
        self.build_dust_grid()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(RepresentationBuilder, self).setup(**kwargs)
        RepresentationBuilderBase.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Get the dust grid, if passed
        if "dust_grid" in kwargs: self.dust_grid = kwargs.pop("dust_grid")

    # -----------------------------------------------------------------

    def load_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the deprojections used for the model ...")

        # 1. Load stellar deprojections
        self.load_stellar_deprojections()

        # 2. load dust deprojections
        self.load_dust_deprojections()

    # -----------------------------------------------------------------

    def load_stellar_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the stellar deprojections ...")

        # Loop over the stellar components
        for name in get_stellar_component_names(self.config.path, self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = load_stellar_component_deprojection(self.config.path, self.model_name, name)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def load_dust_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the dust deprojections ...")

        # Loop over the dust components
        for name in get_dust_component_names(self.config.path, self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = load_dust_component_deprojection(self.config.path, self.model_name, name)
            if deprojection is not None: self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        azimuth = 0.0

        # Use grid?
        if prompt_yn("grid_resolution", "use the resolution of the dust grid for setting up the instruments?"):

            # Determine smallest scale
            smallest_scale = smallest_scale_for_dust_grid(self.dust_grid)

            # Determine instrument pixelscale
            ratio = prompt_real("pixelscale_to_grid_scale_ratio", "ratio of the instrument pixelscale to the smallest scale of the dust grid (e.g. 10)")
            physical_pixelscale = smallest_scale * ratio

            # Set number of pixels from extent
            extent = self.dust_grid.x_extent

            pixels_x = int(math.ceil(extent/physical_pixelscale))
            pixels_y = pixels_x

            x_center = 0.5 * (pixels_x - 1)
            y_center = 0.5 * (pixels_y - 1)

            # Pixel to physical
            center_x = x_center * physical_pixelscale
            center_y = y_center * physical_pixelscale

            # Create projections
            # distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y
            earth_projection = GalaxyProjection(self.galaxy_distance, self.galaxy_inclination, azimuth, self.disk_position_angle, pixels_x, pixels_y, center_x, center_y, extent, extent)
            faceon_projection = FaceOnProjection.from_projection(earth_projection)
            edgeon_projection = EdgeOnProjection.from_projection(earth_projection)

        # Use deprojections
        else:

            # Get the desired deprojection to base the instruments on
            reference_deprojection = self.prompt_deprojection()

            # Create the 'earth' projection system
            earth_projection = GalaxyProjection.from_deprojection(reference_deprojection, self.galaxy_distance, azimuth)

            # Create the face-on projection system
            faceon_projection = FaceOnProjection.from_deprojection(reference_deprojection, self.galaxy_distance)

            # Create the edge-on projection system
            edgeon_projection = EdgeOnProjection.from_deprojection(reference_deprojection, self.galaxy_distance)

        # Set the projection systems
        self.projections["earth"] = earth_projection
        self.projections["faceon"] = faceon_projection
        self.projections["edgeon"] = edgeon_projection

    # -----------------------------------------------------------------

    def prompt_deprojection(self):

        """
        This function ...
        :return:
        """

        # Dictionary for the options
        options = dict()

        lowest_pixelscale = None
        lowest_pixelscale_name = None
        lowest_pixelscale_title = None

        # Loop over the different deprojection models
        for name, title in self.deprojections:

            # Determine name and description
            option = name
            pixelscale = self.deprojections[(name, title)].pixelscale
            if lowest_pixelscale is None or pixelscale < lowest_pixelscale:
                lowest_pixelscale = pixelscale
                lowest_pixelscale_name = name
                lowest_pixelscale_title = title
            description = "pixelscale of the " + title.lower() + " input map (" + represent_quantity(lowest_pixelscale) + ")"

            # Add the option
            options[option] = description

        # name, description, choices=None, default=None
        answer = prompt_string("reference map", "input map to use as the reference for the spatial resolution (dust grid and instruments) of the model representation", choices=options, default=lowest_pixelscale_name)

        # Set the reference deprojection
        answer_title = None
        for name, title in self.deprojections:
            if name == answer:
                answer_title = title
                break

        # Return the deprojection
        return self.deprojections[(answer, answer_title)]

    # -----------------------------------------------------------------

    @property
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["earth"]

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Create an SED instrument
        self.instruments["SED"] = SEDInstrument.from_projection(self.earth_projection)

        # Create a frame instrument to generate datacube
        self.instruments["frame"] = FrameInstrument.from_projection(self.earth_projection)

        # Create a simple instrument (SED + frame)
        self.instruments["simple"] = SimpleInstrument.from_projection(self.earth_projection)

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Load the dust disk deprojection
        deprojection = self.definition.dust_deprojection

        # Set minimum level
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_max_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_max_level
        else: min_level = None

        # Set max ndivisions per pixel
        max_ndivisions_per_pixel = 1. / self.config.dg.scale  # default 1/0.5 = 2 divisions along each direction per pixel

        # Create the dust grid
        # grid_type, deprojection, distance, sky_ellipse, min_level, max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.
        self.dust_grid = create_one_dust_grid_for_galaxy_from_deprojection(self.config.dg.grid_type, deprojection,
                                                                           self.galaxy_distance, self.truncation_ellipse,
                                                                           min_level, self.config.dg.max_mass_fraction,
                                                                           max_ndivisions_per_pixel, self.config.dg.scale_heights)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 3. Write the dust grids
        self.write_dust_grid()

        # 1. Write the projections
        self.write_projections()

        # 2. Write the instruments
        self.write_instruments()

        # 4. Write the representations table
        self.write_table()

    # -----------------------------------------------------------------

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid ...")

        # Write the dust grid
        self.dust_grid.saveto(self.representation.dust_grid_path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections["earth"].saveto(self.representation.earth_projection_path)

        # Write the edgeon projection system
        self.projections["edgeon"].saveto(self.representation.edgeon_projection_path)

        # Write the faceon projection system
        self.projections["faceon"].saveto(self.representation.faceon_projection_path)

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SED, frame and simple instruments ...")

        # Write the SED instrument
        self.instruments["SED"].saveto(self.representation.sed_instrument_path)

        # Write the frame instrument
        self.instruments["frame"].saveto(self.representation.frame_instrument_path)

        # Write the simple instrument
        self.instruments["simple"].saveto(self.representation.simple_instrument_path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing the representations table ...")

        # Open the table, add an entry, and save the table
        table = self.representations_table
        table.add_entry(self.representation_name, self.model_name)
        table.save()

# -----------------------------------------------------------------
