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

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ...basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ....core.basics.log import log
from ...basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ....core.basics.configuration import prompt_string, prompt_yn, prompt_real
from ....core.units.stringify import represent_quantity
from ...component.galaxy import GalaxyModelingComponent
from ....core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection, smallest_scale_for_dust_grid
from .base import RepresentationBuilderBase
from ....magic.basics.vector import PixelShape
from ...basics.projection import get_center, get_physical_center
from ....magic.basics.stretch import PhysicalExtent
from ....core.tools.stringify import tostr
from ....core.tools import numbers

# -----------------------------------------------------------------

default_scale_heights = 15. # number of times to take the scale height as the vertical radius of the model

# -----------------------------------------------------------------

class GalaxyRepresentationBuilder(RepresentationBuilderBase, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

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
        for name in self.suite.get_dust_component_names(self.model_name):

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
            earth, faceon, edgeon = create_projections_from_dust_grid(self.dust_grid, self.galaxy_distance, self.disk_inclination, azimuth, self.disk_position_angle)

        # Use deprojections
        # galaxy_distance, azimuth
        else: earth, faceon, edgeon = create_projections_from_deprojections(self.deprojections, self.galaxy_distance, azimuth, self.config.old_scale_heights, scale_heights_reference="old")

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
    earth_projection = GalaxyProjection(distance=galaxy_distance, inclination=galaxy_inclination, azimuth=azimuth,
                                        position_angle=disk_position_angle, pixels_x=pixels_x, pixels_y=pixels_y,
                                        center_x=center_x, center_y=center_y, field_x=extent, field_y=extent)
    faceon_projection = FaceOnProjection.from_projection(earth_projection)
    edgeon_projection = EdgeOnProjection.from_projection(earth_projection)

    # Return the projections
    return earth_projection, faceon_projection, edgeon_projection

# -----------------------------------------------------------------

def create_projections_from_deprojections(deprojections, galaxy_distance, azimuth, scale_heights=default_scale_heights,
                                          return_deprojection_name=False, scale_heights_reference=None,
                                          reference_deprojection_name=None, radial_factor=1, from_projection=True):

    """
    This function ...
    :param deprojections:
    :param galaxy_distance:
    :param azimuth:
    :param scale_heights:
    :param return_deprojection_name:
    :param scale_heights_reference:
    :param reference_deprojection_name:
    :param radial_factor:
    :param from_projection:
    :return:
    """

    # Reference deprojection name not given
    # Get the desired deprojection to base the instruments on
    if reference_deprojection_name is None: reference_deprojection, deprojection_name = prompt_deprojection(deprojections, return_name=True)

    # Reference deprojection name is given
    else: reference_deprojection, deprojection_name = get_deprojection(deprojections, reference_deprojection_name), reference_deprojection_name

    # CREATE EARTH PROJECTION
    earth_projection = create_projection_from_deprojection(reference_deprojection, galaxy_distance, azimuth)

    # CREATE FACEON PROJECTION
    if from_projection: faceon_projection = create_faceon_projection_from_earth_projection(earth_projection, radial_factor=radial_factor)
    else: faceon_projection = create_faceon_projection_from_deprojection(reference_deprojection, radial_factor=radial_factor)

    # Determine the reference scale height for determining the physical z scale
    if scale_heights_reference is not None:
        deprojection = get_deprojection(deprojections, scale_heights_reference)
        scale_height = deprojection.scale_height
    else: scale_height = reference_deprojection.scale_height

    # Determine the z extent of the model based on the given number of scale heights
    z_extent = 2. * scale_height * scale_heights

    # CREATE EDGEON PROJECTION
    if from_projection: edgeon_projection = create_edgeon_projection_from_earth_projection(earth_projection, z_extent, radial_factor=radial_factor)
    else: edgeon_projection = create_edgeon_projection_from_deprojection(reference_deprojection, z_extent, radial_factor=radial_factor)

    # Return the projections
    if return_deprojection_name: return earth_projection, faceon_projection, edgeon_projection, deprojection_name
    else: return earth_projection, faceon_projection, edgeon_projection

# -----------------------------------------------------------------

def get_physical_pixelscale_from_map(the_map, distance, downsample_factor=1.):

    """
    This function ...
    :param the_map:
    :param distance:
    :param downsample_factor:
    :return:
    """

    pixelscale = the_map.average_pixelscale * downsample_factor
    return (abs(pixelscale) * distance).to("pc", equivalencies=dimensionless_angles())

# -----------------------------------------------------------------

def create_projection_from_deprojection(deprojection, distance=None, azimuth=0.0):

    """
    This function ...
    :param deprojection:
    :param distance:
    :param azimuth:
    :return:
    """

    # Get the distance
    if distance is None:
        if deprojection.distance is None: raise ValueError("Distance is not defined")
        distance = deprojection.distance

    # Create
    return GalaxyProjection.from_deprojection(deprojection, distance, azimuth)

# -----------------------------------------------------------------

def create_faceon_projection_from_deprojection(deprojection, radial_factor=1):

    """
    This function ...
    :param deprojection:
    :param radial_factor:
    :return:
    """

    # Get properties
    distance = deprojection.distance
    physical_pixelscale = deprojection.pixelscale

    # Determine extent in the radial direction
    radial_extent = max(deprojection.x_range.span, deprojection.y_range.span)

    # Determine number of pixels
    npixels = int(round(radial_extent / physical_pixelscale)) * radial_factor

    # Create and return
    return create_faceon_projection(npixels, physical_pixelscale, distance)

# -----------------------------------------------------------------

def create_faceon_projection_from_earth_projection(earth_projection, radial_factor=1):

    """
    This function ...
    :param earth_projection:
    :param radial_factor:
    :return:
    """

    # Determine extent in the radial direction
    radial_extent = max(earth_projection.field_x, earth_projection.field_y)

    # Get pixelscale
    #physical_pixelscale = earth_projection.physical_pixelscale
    physical_pixelscale = earth_projection.physical_pixelscale.average

    # Determine number of pixels
    npixels = int(round(radial_extent / physical_pixelscale)) * radial_factor

    # Create and return
    return create_faceon_projection(npixels, physical_pixelscale, earth_projection.distance)

# -----------------------------------------------------------------

def create_faceon_projection(npixels, physical_pixelscale, distance):

    """
    This function ...
    :param npixels:
    :param physical_pixelscale:
    :param distance:
    :return:
    """

    # Make sure that the number of pixels is odd
    npixels = numbers.nearest_odd_integer(npixels)
    npixels = PixelShape.square(npixels)

    # Get the center pixel
    center = get_center(npixels)

    # Get field of view
    field = PhysicalExtent(physical_pixelscale * npixels.x, physical_pixelscale * npixels.y)

    # Get physical center
    center_physical = get_physical_center(field, npixels, center)

    # Create the face-on projection system
    faceon_projection = FaceOnProjection(distance=distance, pixels_x=npixels.x, pixels_y=npixels.y,
                                         center_x=center_physical.x, center_y=center_physical.y,
                                         field_x=field.x, field_y=field.y)

    # Return the projection
    return faceon_projection

# -----------------------------------------------------------------

def create_edgeon_projection_from_deprojection(deprojection, z_extent, radial_factor=1):

    """
    Thisf unction ...
    :param deprojection:
    :param z_extent:
    :param radial_factor:
    :return:
    """

    # Get properties
    distance = deprojection.distance
    physical_pixelscale = deprojection.pixelscale

    # Determine extent in the radial and in the vertical direction
    radial_extent = max(deprojection.x_range.span, deprojection.y_range.span)

    # Determine number of pixels
    nx = int(round(radial_extent / physical_pixelscale)) * radial_factor
    nz = int(round(z_extent / physical_pixelscale))

    # Create and return
    return create_edgeon_projection(nx, nz, physical_pixelscale, distance)

# -----------------------------------------------------------------

def create_edgeon_projection_from_earth_projection(earth_projection, z_extent, radial_factor=1):

    """
    This function ...
    :param earth_projection:
    :param z_extent:
    :param radial_factor:
    :return:
    """

    # Get pixelscale
    #physical_pixelscale = earth_projection.physical_pixelscale
    physical_pixelscale = earth_projection.physical_pixelscale.average

    # Determine extent in the radial and in the vertical direction
    radial_extent = max(earth_projection.field_x, earth_projection.field_y)

    # Determine the z extent
    #z_extent = 2. * scaleheight * scale_heights

    # Determine number of pixels
    nx = int(round(radial_extent / physical_pixelscale)) * radial_factor
    nz = int(round(z_extent / physical_pixelscale))

    # Create and return
    return create_edgeon_projection(nx, nz, physical_pixelscale, earth_projection.distance)

# -----------------------------------------------------------------

def create_edgeon_projection(nx, nz, physical_pixelscale, distance):

    """
    This function ...
    :param nx:
    :param nz:
    :param physical_pixelscale:
    :param distance:
    :return:
    """

    # Make sure that the number of pixels is odd
    nx = numbers.nearest_odd_integer(nx)
    nz = numbers.nearest_odd_integer(nz)

    # Return the pixel shape
    npixels = PixelShape.from_xy(nx, nz)

    # Get the center pixel
    center = get_center(npixels)

    # Get field of view
    field = PhysicalExtent(physical_pixelscale * npixels.x, physical_pixelscale * npixels.y)

    # Get physical center
    center_physical = get_physical_center(field, npixels, center)

    # Create the projection
    # edgeon_projection = EdgeOnProjection.from_deprojection(reference_deprojection, galaxy_distance)
    edgeon_projection = EdgeOnProjection(distance=distance, pixels_x=npixels.x, pixels_y=npixels.y,
                                         center_x=center_physical.x, center_y=center_physical.y,
                                         field_x=field.x, field_y=field.y)

    # Return the projection
    return edgeon_projection

# -----------------------------------------------------------------

def prompt_deprojection(deprojections, return_name=False):

    """
    This function ...
    :param deprojections:
    :param return_name:
    :return:
    """

    # Dictionary for the options
    options = dict()

    lowest_pixelscale = None
    lowest_pixelscale_name = None
    #lowest_pixelscale_title = None

    # Loop over the different deprojection models
    for name, title in deprojections:

        # Determine name and description
        option = name
        pixelscale = deprojections[(name, title)].pixelscale
        if lowest_pixelscale is None or pixelscale < lowest_pixelscale:
            lowest_pixelscale = pixelscale
            lowest_pixelscale_name = name
            #lowest_pixelscale_title = title
        description = "pixelscale of the " + title.lower() + " input map (" + represent_quantity(lowest_pixelscale) + ")"

        # Add the option
        options[option] = description

    # name, description, choices=None, default=None
    answer = prompt_string("reference_map", "input map to use as the reference for the spatial resolution (dust grid and instruments) of the model representation", choices=options, default=lowest_pixelscale_name)

    # Set the reference deprojection
    answer_title = None
    for name, title in deprojections:
        if name == answer:
            answer_title = title
            break

    # Return the deprojection
    deprojection = deprojections[(answer, answer_title)]
    if return_name: return deprojection, answer
    else: return deprojection

# -----------------------------------------------------------------

def get_deprojection(deprojections, name):

    """
    This function ...
    :param deprojections:
    :param name:
    :return:
    """

    matching_titles = []
    valid_names = set()

    # Loop over the deprojections
    for name_i, title_i in deprojections:
        if name_i == name: matching_titles.append(title_i)
        valid_names.add(name_i)
    valid_names = list(valid_names)

    # Return the deprojection if possible
    if len(matching_titles) == 1: return deprojections[(name, matching_titles[0])]
    elif len(matching_titles) == 0: raise ValueError("Deprojection with the name '" + name + "' not found, valid names are: '" + tostr(valid_names) + "'")
    else: raise ValueError("Ambigious result: matching titles are '" + tostr(matching_titles) + "'")

# -----------------------------------------------------------------
