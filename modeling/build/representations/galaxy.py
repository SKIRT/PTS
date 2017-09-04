#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations.galaxy Contains the GalaxyRepresentationBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from ...basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ....core.basics.log import log
from ...basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ....core.basics.configuration import prompt_string, prompt_yn, prompt_real
from ....core.units.stringify import represent_quantity
from ...component.galaxy import GalaxyModelingComponent
from ....core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection, smallest_scale_for_dust_grid
from .base import RepresentationBuilderBase

# -----------------------------------------------------------------

class GalaxyRepresentationBuilder(RepresentationBuilderBase, GalaxyModelingComponent):
    
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
        RepresentationBuilderBase.__init__(self, no_config=True)
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
        for name in self.suite.get_stellar_component_names(self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.suite.load_stellar_component_deprojection(self.model_name, name)
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
        for name in self.suite.get_dust_component_names(self.config.path, self.model_name):

            # Load the deprojection of the component, if applicable
            title, deprojection = self.suite.load_dust_component_deprojection(self.model_name, name)
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

            # Create the projections
            # dust_grid, galaxy_distance, galaxy_inclination, azimuth, disk_position_angle
            earth, faceon, edgeon = create_projections_from_dust_grid(self.dust_grid, self.galaxy_distance, self.galaxy_inclination, azimuth, self.disk_position_angle)

        # Use deprojections
        # galaxy_distance, azimuth
        else: earth, faceon, edgeon = create_projections_from_deprojections(self.deprojections, self.galaxy_distance, azimuth)

        # Set the projection systems
        self.projections["earth"] = earth
        self.projections["faceon"] = faceon
        self.projections["edgeon"] = edgeon

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
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_min_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_min_level
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

def create_projections_from_dust_grid(dust_grid, galaxy_distance, galaxy_inclination, azimuth, disk_position_angle):

    """
    This function ...
    :param dust_grid:
    :param galaxy_distance:
    :param galaxy_inclination:
    :param azimuth:
    :param disk_position_angle
    :return:
    """

    # Determine smallest scale
    smallest_scale = smallest_scale_for_dust_grid(dust_grid)

    # Determine instrument pixelscale
    ratio = prompt_real("pixelscale_to_grid_scale_ratio", "ratio of the instrument pixelscale to the smallest scale of the dust grid (e.g. 10)")
    physical_pixelscale = smallest_scale * ratio

    # Set number of pixels from extent
    extent = dust_grid.x_extent

    pixels_x = int(math.ceil(extent / physical_pixelscale))
    pixels_y = pixels_x

    x_center = 0.5 * (pixels_x - 1)
    y_center = 0.5 * (pixels_y - 1)

    # Pixel to physical
    center_x = x_center * physical_pixelscale
    center_y = y_center * physical_pixelscale

    # Create projections
    # distance, inclination, azimuth, position_angle, pixels_x, pixels_y, center_x, center_y, field_x, field_y
    earth_projection = GalaxyProjection(galaxy_distance, galaxy_inclination, azimuth,
                                        disk_position_angle, pixels_x, pixels_y, center_x, center_y, extent,
                                        extent)
    faceon_projection = FaceOnProjection.from_projection(earth_projection)
    edgeon_projection = EdgeOnProjection.from_projection(earth_projection)

    # Return the projections
    return earth_projection, faceon_projection, edgeon_projection

# -----------------------------------------------------------------

def create_projections_from_deprojections(deprojections, galaxy_distance, azimuth):

    """
    This function ...
    :param deprojections:
    :param galaxy_distance:
    :param azimuth:
    :return:
    """

    # Get the desired deprojection to base the instruments on
    reference_deprojection = prompt_deprojection(deprojections)

    # Create the 'earth' projection system
    earth_projection = GalaxyProjection.from_deprojection(reference_deprojection, galaxy_distance, azimuth)

    # Create the face-on projection system
    faceon_projection = FaceOnProjection.from_deprojection(reference_deprojection, galaxy_distance)

    # Create the edge-on projection system
    edgeon_projection = EdgeOnProjection.from_deprojection(reference_deprojection, galaxy_distance)

    # Return the projections
    return earth_projection, faceon_projection, edgeon_projection

# -----------------------------------------------------------------

def prompt_deprojection(deprojections):

    """
    This function ...
    :param deprojections:
    :return:
    """

    # Dictionary for the options
    options = dict()

    lowest_pixelscale = None
    lowest_pixelscale_name = None
    lowest_pixelscale_title = None

    # Loop over the different deprojection models
    for name, title in deprojections:

        # Determine name and description
        option = name
        pixelscale = deprojections[(name, title)].pixelscale
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
    for name, title in deprojections:
        if name == answer:
            answer_title = title
            break

    # Return the deprojection
    return deprojections[(answer, answer_title)]

# -----------------------------------------------------------------
