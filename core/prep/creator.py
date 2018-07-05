#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.prep.creator Contains the ModelCreator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..basics.configurable import Configurable
from ...modeling.basics.models import create_3d_model, DeprojectionModel3D
from ...modeling.build.representations.galaxy import create_projection_from_deprojection, create_faceon_projection_from_deprojection, create_edgeon_projection_from_deprojection
from ...modeling.basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from ...core.tools import filesystem as fs
from ...modeling.basics.models import deprojection
from ...core.tools.utils import lazyproperty
from ...magic.core.frame import Frame
from ..basics.configuration import prompt_variable
from ...magic.basics.coordinate import PixelCoordinate
from ..tools import browser

# -----------------------------------------------------------------

widget_width = 400
widget_height = 500
widget_style = "dark" # minimal, dark or light

# -----------------------------------------------------------------

class ModelCreator(Configurable):

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
        super(ModelCreator, self).__init__(*args, **kwargs)

        # The model
        self.model = None

        # The projections
        self.projections = None

        # The deprojected map
        self.deprojected = None

        # The projected images
        self.projected = None

        # The view
        self.view = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

        # 2. Create the model
        self.create()

        # 3. Prompt for the model parameters
        if not self.from_map: self.prompt()

        # 4. Deproject the deprojection model
        if self.is_deprojection and self.config.deproject: self.deproject()

        # 5. Project the model
        if self.config.project: self.project()

        # 6. Create the view
        if self.config.view: self.create_view()

        # 7. Write
        self.write()

        # Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup fucntion of the base class
        super(ModelCreator, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @property
    def from_map(self):

        """
        Thisnfunction ...
        :return:
        """

        return self.config.map is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def map(self):

        """
        Thisn function ...
        :return:
        """

        if self.from_map: return Frame.from_file(self.config.map)
        elif self.is_deprojection: return self.model.map
        else: return None

    # -----------------------------------------------------------------

    @property
    def map_filepath(self):

        """
        This function ...
        :return:
        """

        if self.from_map: return self.config.map
        elif self.is_deprojection: return self.model.filename
        else: return None

    # -----------------------------------------------------------------

    @property
    def map_wcs(self):

        """
        This function ...
        :return:
        """

        return self.map.wcs

    # -----------------------------------------------------------------

    @property
    def center_pix(self):

        """
        This function ...
        :return:
        """

        if self.is_deprojection: return self.model.center
        else: return None

    # -----------------------------------------------------------------

    @property
    def center_sky(self):

        """
        This function ...
        :return:
        """

        from ...magic.basics.coordinate import SkyCoordinate
        if self.is_deprojection: return SkyCoordinate.from_pixel(self.center_pix, self.map_wcs)
        else: return None

    # -----------------------------------------------------------------

    def create(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the model ...")

        # Create model from map
        if self.from_map:

            # Model must be deprojection
            if self.config.model_type != deprojection: raise ValueError("Cannot specify map for models other than deprojection models")

            # Prompt for the additional properties
            center, distance, position_angle, inclination, scaleheight = prompt_deprojection_parameters(center=self.config.center, distance=self.config.distance, position_angle=self.config.position_angle, inclination=self.config.inclination, scale_height=self.config.scale_height)
            if isinstance(center, PixelCoordinate): center_sky = center.to_sky(self.map_wcs)
            else: center_sky = center

            # Create deprojection
            self.model = DeprojectionModel3D.from_wcs(self.map_wcs, center_sky, distance, position_angle, inclination, self.map_filepath, scaleheight)

        # Initialize the model
        else: self.model = create_3d_model(self.config.model_type)

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the model parameters ...")

        # Prompt the properties
        self.model.prompt_properties()

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojection_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory(self.config.name + "_deprojection", create=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def deprojection_temp_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.deprojection_path, "temp")

    # -----------------------------------------------------------------

    def deproject(self):

        """
        This funtion ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the map ...")

        from ...magic.misc.deproject import deproject

        # Create deprojected map
        self.deprojected = deproject(self.model, method=self.config.deprojection_method, output_path=self.deprojection_path, simulation_path=self.deprojection_temp_path, npackages=self.config.npackages, parallelization=self.config.parallelization)

    # -----------------------------------------------------------------

    @lazyproperty
    def projection_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory(self.config.name + "_projection", create=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def projection_temp_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.projection_path, "temp")

    # -----------------------------------------------------------------

    def project(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the model ...")

        from ...modeling.projection.projector import project

        # Create projections
        self.create_projections()

        # Project
        self.projected = project(self.model, self.projections, output_path=self.projection_path, simulation_path=self.projection_temp_path, npackages=self.config.npackages, parallelization=self.config.parallelization, coordinate_systems={"earth":self.map_wcs})

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Initialize
        self.projections = dict()

        # From deprojection?
        if self.is_deprojection:

            # Create projections
            earth, faceon, edgeon = create_projections_from_deprojection(self.model, scale_heights=self.config.scale_heights)

            # Set projections
            self.projections["earth"] = earth
            self.projections["faceon"] = faceon
            self.projections["edgeon"] = edgeon

        # Prompt
        else: self.prompt_projections()

    # -----------------------------------------------------------------

    def prompt_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for projections ...")

        # Create earth projection
        earth = GalaxyProjection.prompt(center=self.center_pix, distance=self.config.distance, position_angle=self.config.position_angle, inclination=self.config.inclination)

        # Create faceon projection
        faceon = FaceOnProjection.from_projection(earth)

        # Create edgeon projection
        edgeon = EdgeOnProjection.from_projection(earth)

        # Set
        self.projections["earth"] = earth
        self.projections["faceon"] = faceon
        self.projections["edgeon"] = edgeon

    # -----------------------------------------------------------------

    @property
    def has_projections(self):

        """
        Thisfunction ...
        :return:
        """

        return self.projections is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def view_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory(self.config.name + "_view", create=True)

    # -----------------------------------------------------------------

    def create_view(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the 3D view ...")

        from ...modeling.plotting.model import plot_galaxy_components, generate_html

        # Plot settings
        plot_kwargs = dict()
        plot_kwargs["width"] = widget_width
        plot_kwargs["height"] = widget_height
        plot_kwargs["style"] = widget_style

        # Render settings
        render_kwargs = dict()
        #render_kwargs["only_body"] = True

        # Add the component
        components = {"model": self.model}
        title = self.config.name

        box = plot_galaxy_components(components, draw=True, show=False, **plot_kwargs)
        html = generate_html(box, title, self.view_path, **render_kwargs)

        # Set the view
        self.view = html

    # -----------------------------------------------------------------

    @property
    def has_view(self):

        """
        Thisfunction ...
        :return:
        """

        return self.view is not None

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the model
        self.write_model()

        # Write the projections
        if self.has_projections: self.write_projections()

        # Write the view
        if self.has_view: self.write_view()

    # -----------------------------------------------------------------

    @property
    def is_deprojection(self):

        """
        Thisn function ...
        :return:
        """

        return isinstance(self.model, DeprojectionModel3D)

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.model.distance

    # -----------------------------------------------------------------

    @property
    def has_projected(self):

        """
        This function ...
        :return:
        """

        return self.projected is not None

    # -----------------------------------------------------------------

    def write_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model ...")

        # Determine the path
        path = self.output_path_file(self.config.name + ".mod")

        # Save the model
        self.model.saveto(path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projections ...")

        # Create directory
        path = self.output_path_directory(self.config.name + "_projections", create=True)

        # Save
        for name in self.projections:

            # Save
            filepath = fs.join(path, name + ".proj")
            self.projections[name].saveto(filepath)

    # -----------------------------------------------------------------

    @property
    def view_filepath(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.view_path, "model.html")

    # -----------------------------------------------------------------

    def write_view(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the view ...")

        # Write
        import io
        with io.open(self.view_filepath, "w", encoding='utf8') as f:
            f.write(self.view)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        # Show the view
        if self.has_view: self.show_view()

    # -----------------------------------------------------------------

    def show_view(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Showing the view ...")

        # Open the view
        browser.open_path(self.view_filepath)

# -----------------------------------------------------------------

def create_projections_from_deprojection(deprojection, distance=None, azimuth=0.0, scale_heights=15):

    """
    This function ...
    :param deprojection:
    :param distance:
    :param azimuth:
    :param scale_heights:
    :return:
    """

    # Determine the distance
    if deprojection.distance is None:
        if distance is None: raise ValueError("Distance is not set")
    else: distance = deprojection.distance

    # Create the 'earth' projection system
    earth_projection = create_projection_from_deprojection(deprojection, distance, azimuth)

    # Create the face-on projection system
    faceon_projection = create_faceon_projection_from_deprojection(deprojection)

    # Determine the reference scale height for determining the physical z scale
    # Determine the z extent of the model based on the given number of scale heights
    scale_height = deprojection.scale_height
    z_extent = 2. * scale_height * scale_heights

    # Create the edge-on projection system
    edgeon_projection = create_edgeon_projection_from_deprojection(deprojection, z_extent)

    # Return the projections
    return earth_projection, faceon_projection, edgeon_projection

# -----------------------------------------------------------------

def prompt_deprojection_parameters(**properties):

    """
    This function ...
    :param properties:
    :return:
    """

    # Prompt center
    if properties.get("center", None) is None: center = prompt_variable("center", "sky_or_pixel_coordinate", "center position of the galaxy (in sky or image coordinates)")
    else: center = properties.pop("center")

    # Prompt distance
    if properties.get("distance", None) is None: distance = prompt_variable("distance", "length_quantity", "distance of the galaxy")
    else: distance = properties.pop("distance")

    # Prompt position angle
    if properties.get("position_angle", None) is None: position_angle = prompt_variable("position_angle", "angle", "position angle of the plane of the galaxy in the image")
    else: position_angle = properties.pop("position_angle")

    # Prompt inclination angle
    if properties.get("inclination", None) is None: inclination = prompt_variable("inclination", "angle", "inclination angle of the plane of the galaxy in the image")
    else: inclination = properties.pop("inclination")

    # Prompt scale height
    if properties.get("scale_height", None) is None: scale_height = prompt_variable("scale_height", "length_quantity", "vertical scale height")
    else: scale_height = properties.pop("scale_height")

    # Return the parameters
    return center, distance, position_angle, inclination, scale_height

# -----------------------------------------------------------------
