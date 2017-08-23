#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.deprojection.deprojector Contains the Deprojector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.simulation.execute import SkirtExec
from ...core.tools import introspection
from ...core.simulation.skifile import LabeledSkiFile
from ..basics.instruments import FrameInstrument
from ...magic.core.frame import Frame
from ...core.units.parsing import parse_unit as u
from ..component.galaxy import GalaxyModelingComponent
from ..build.suite import create_deprojection_for_map
from ...core.tools.stringify import tostr
from ..plotting.model import xy

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")

# -----------------------------------------------------------------

class Deprojector(GalaxyModelingComponent):
    
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
        super(Deprojector, self).__init__(*args, **kwargs)

        # The SKIRT execution environment
        self.skirt = SkirtExec()

        ##

        # Paths
        self.root_path = None

        # The maps
        self.maps = dict()

        # The map paths
        self.map_paths = None

        # The scale heights
        self.scale_heights = dict()

        # The deprojection models
        self.deprojections = dict()

        # The deprojected maps
        self.deprojected = dict()

        # The edge-on maps
        self.edgeon = dict()

        ##

        # The output paths
        self.output_paths = dict()

        ##

        # The ski files
        self.ski_files = dict()

        # The ski file paths
        self.ski_paths = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the deprojection models
        self.create_models()

        # 4. Write
        self.write()

        # 5. Deproject
        self.deproject()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Deprojector, self).setup(**kwargs)

        # Get paths
        self.root_path = kwargs.pop("root_path", None)

        # Checks
        if "map" in kwargs and "maps" in kwargs: raise ValueError("Cannot specify 'map' and 'maps' simultaneously")

        # Get the single map
        if "map" in kwargs:
            if "name" not in kwargs: raise ValueError("When passing only one map, a name must be specified")
            name = kwargs.pop("name")
            self.maps[name] = kwargs.pop("name")

        # Get the maps
        elif "maps" in kwargs: self.maps = kwargs.pop("maps")

        # Checks
        if "map_path" in kwargs and "map_paths" in kwargs: raise ValueError("Cannot specify 'map_path' and 'map_paths' simultaneously")

        # Path(s) are given
        if "map_path" in kwargs or "map_paths" in kwargs:

            # No maps -> load the maps from file
            if self.no_maps:

                if "map_path" in kwargs:
                    if "name" not in kwargs: raise ValueError("When passing only one map path, a name must be specified")
                    name = kwargs.pop("name")
                    self.maps[name] = Frame.from_file(kwargs.pop("map_path"))

                if "map_paths" in kwargs:
                    paths = kwargs.pop("map_paths")
                    for name in paths:
                        self.maps[name] = Frame.from_file(paths[name])

            else:

                # Check
                if "map_path" in kwargs and not self.has_single_map: raise ValueError("Specified a single map path but not single map")

                # Set the map path
                if "map_path" in kwargs: self.map_paths[self.single_map_name] = kwargs.pop("map_path")

                # Set the map paths
                if "map_paths" in kwargs: self.map_paths = kwargs.pop("map_paths")

        # Check that each map has a path, otherwise set it via the map's path attribute
        if self.map_paths is not None:
            for name in self.maps:
                if name not in self.map_paths or self.map_paths[name] is None:
                    path = self.maps[name].path
                    if path is None: raise ValueError("The path of the '" + name + "' map is undefined")
                    self.map_paths[name] = path

        # Check if each map has a path, then set them
        elif self.all_maps_have_path:
            self.map_paths = dict()
            for name in self.maps: self.map_paths[name] = self.maps[name].path

        # Checks
        if "scale_height" in kwargs and "scale_heights" in kwargs: raise ValueError("Cannot specify 'scale_height' and 'scale_heights' simultaneously")

        # Set the same scale height for each map
        if "scale_height" in kwargs:
            scale_height = kwargs.pop("scale_height")
            for name in self.maps: self.scale_heights[name] = scale_height

        # Set the scale heights for each map
        elif "scale_heights" in kwargs: self.scale_heights = kwargs.pop("scale_heights")

        # Check that each map has a scale height
        for name in self.maps:
            if name not in self.scale_heights: raise ValueError("Scale height for '" + name + "' map is not defined")

        # Make directories
        if self.root_path is not None: self.create_directories()

        # Check leftover arguments
        if len(kwargs) > 0: raise ValueError("Could not resolve all input: " + tostr(kwargs))

    # -----------------------------------------------------------------

    @property
    def all_maps_have_path(self):

        """
        This function ...
        :return:
        """

        for name in self.maps:
            if self.maps[name].path is None: return False
        return True

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the directories ...")

        # Loop over the maps
        for name in self.map_names: self.output_paths[name] = fs.create_directory_in(self.root_path, name)

    # -----------------------------------------------------------------

    @property
    def map_names(self):

        """
        This function ...
        :return:
        """

        return self.maps.keys()

    # -----------------------------------------------------------------

    @property
    def nmaps(self):

        """
        This function ...
        :return:
        """

        return len(self.maps)

    # -----------------------------------------------------------------

    @property
    def no_maps(self):

        """
        This function ...
        :return:
        """

        return self.nmaps == 0

    # -----------------------------------------------------------------

    @property
    def has_single_map(self):

        """
        This function ...
        :return:
        """

        return self.nmaps == 1

    # -----------------------------------------------------------------

    @property
    def single_map_name(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_map: raise ValueError("Not a single map")
        return self.map_names[0]

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return:
        """

        return self.maps[self.single_map_name]

    # -----------------------------------------------------------------

    @property
    def single_deprojected(self):

        """
        This function ...
        :return:
        """

        return self.deprojected[self.single_map_name]

    # -----------------------------------------------------------------

    def create_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection models ...")

        # Loop over the maps
        for name in self.map_names:

            # Debugging
            log.debug("Creating the deprojection model for the '" + name + "' map ...")

            # Get the file path
            if self.map_paths is not None: filename = self.map_paths[name]
            else: filename = None

            # Get the scale height
            scaleheight = self.scale_heights[name]

            # Create deprojection model
            deprojection = create_deprojection_for_map(self.galaxy_properties, self.disk_position_angle, self.maps[name], filename, scaleheight)

            # Set the deprojection
            self.deprojections[name] = deprojection

    # -----------------------------------------------------------------

    @property
    def has_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.map_paths is None

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        if not self.has_map_paths: self.write_maps()

        # Write the deprojection
        if self.config.writing.deprojections: self.write_deprojections()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Loop over the maps
        for name in self.maps:

            # Determine the path
            path = fs.join(self.output_paths[name], "map.fits")

            # Save the map
            self.maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_deprojections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojection models ...")

        # Loop over the maps
        for name in self.deprojections:

            # Determine the path
            path = fs.join(self.output_paths[name], "deprojection.mod")

            # Save the deprojection
            self.deprojections[name].saveto(path)

    # -----------------------------------------------------------------

    @property
    def deproject_with_pts(self):

        """
        This function ...
        :return:
        """

        return self.config.method == "pts"

    # -----------------------------------------------------------------

    @property
    def deproject_with_skirt(self):

        """
        This function ...
        :return:
        """

        return self.config.method == "skirt"

    # -----------------------------------------------------------------

    def deproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the maps ...")

        # Deproject within PTS
        if self.deproject_with_pts: self.deproject_pts()

        # Deproject with SKIRT
        elif self.deproject_with_skirt: self.deproject_skirt()

    # -----------------------------------------------------------------

    def deproject_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the maps with PTS ...")

        unit = "pc"

        shape = 128

        # Loop over the components
        for name in self.deprojections:

            # Debugging
            log.debug("Computing the deprojected surface density of the '" + name + "' map ...")

            # Set the model limits
            x_range_scalar = self.deprojections[name].x_range.to(unit).value.as_tuple()
            y_range_scalar = self.deprojections[name].y_range.to(unit).value.as_tuple()
            #z_range_scalar = self.deprojections[name].z_range.to(unit).value.as_tuple()
            limits = [x_range_scalar, y_range_scalar]

            # Create coordinate data
            x, y = xy(shape=shape, limits=limits)

            # Calculate the surface density
            density = self.deprojections[name].surface_density_function(normalize=True)(x, y)

            # Create deprojected map
            deprojected = Frame(density)

            # Set
            self.deprojected[name] = deprojected

    # -----------------------------------------------------------------

    def deproject_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting by launching SKIRT simulations ...")

        # Creat the ski files
        self.create_ski()

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski files ...")

        # Create instrument
        # self.instrument = FrameInstrument.from_projection(self.earth_projection)
        self.instrument = FrameInstrument.from_projection(self.faceon_projection)

        # Load ski template
        self.ski_template = LabeledSkiFile(template_ski_path)

        # Convert to oligochromatic simulation
        self.ski_template.to_oligochromatic(1. * u("micron"))

        # Remove the dust system
        self.ski_template.remove_dust_system()

        # Set number of packages per wavelength
        self.ski_template.setpackages(1e6)

        # Add one instrument
        self.ski_template.remove_all_instruments()
        self.ski_template.add_instrument("faceon", self.instrument)

        # Old

        # old_ski = self.ski_template.copy()
        # old_ski.remove_stellar_components_except("Evolved stellar disk")
        # # old_ski.remove_dust_system()
        # old_ski.set_stellar_component_geometry("Evolved stellar disk", self.deprojections["old stars"])
        # self.ski_files["old stars"] = old_ski
        #
        # young_ski = self.ski_template.copy()
        # young_ski.remove_stellar_components_except("Young stars")
        # young_ski.set_stellar_component_geometry("Young stars", self.deprojections["young stars"])
        # self.ski_files["young stars"] = young_ski
        #
        # ionizing_ski = self.ski_template.copy()
        # ionizing_ski.remove_stellar_components_except("Ionizing stars")
        # ionizing_ski.set_stellar_component_geometry("Ionizing stars", self.deprojections["ionizing stars"])
        # self.ski_files["ionizing stars"] = ionizing_ski
        #
        # dust_ski = self.ski_template.copy()
        # dust_ski.remove_stellar_components_except("Ionizing stars")
        # dust_ski.set_stellar_component_geometry("Ionizing stars", self.deprojections["dust"])
        # self.ski_files["dust"] = dust_ski

    # -----------------------------------------------------------------

    def launch_skirt(self):

        """
        This function ...
        :return:
        """

        # Loop over the ski files
        for name in self.ski_files:

            #
            ski_path = self.ski_paths[name]
            out_path = self.output_paths[name]

            #prefix = fs.strip_extension(fs.name(ski_path))

            # Perform the SKIRT simulation
            simulation = self.skirt.execute(ski_path, inpath=self.maps_path, outpath=out_path, single=True)

            # Determine path
            frame_path = fs.join(out_path, simulation.prefix() + "_faceon_total.fits")

            # Open the output frame
            frame = Frame.from_file(frame_path)

            # Set wcs
            frame.wcs = self.reference_wcs

            # Save frame
            frame.saveto(frame_path)

            # Set the deprojected map
            self.deprojected[name] = frame

# -----------------------------------------------------------------

# DEPROJECT WITH PYTHON:

# Import standard modules
#import math
#from skimage import transform as tf

# -----------------------------------------------------------------

#tform = tf.SimilarityTransform(scale=1, rotation=math.pi / 4,
#                               translation=(text.shape[0] / 2, -100))
#rotated = tf.warp(text, tform)

# -----------------------------------------------------------------
