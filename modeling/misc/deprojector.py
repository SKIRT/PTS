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
from ...core.launch.launcher import SKIRTLauncher
from ...core.tools import introspection
from ...core.simulation.skifile import LabeledSkiFile
from ..basics.instruments import FrameInstrument
from ...magic.core.frame import Frame
from ...core.units.parsing import parse_unit as u
from ..component.galaxy import GalaxyModelingComponent
from ..build.suite import create_deprojection_for_map
from ...core.tools.stringify import tostr
from ..plotting.model import xy
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.parallelization import Parallelization
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")

# -----------------------------------------------------------------

faceon_name = "faceon"
edgeon_name = "edgeon"

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
        #self.skirt = SkirtExec()

        # Create the SKIRT launcher
        self.launcher = SKIRTLauncher()

        # Smile
        self.smile = SKIRTSmileSchema()

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
        self.deprojections = None

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

        # The dust grids
        self.dust_grids = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the deprojection models if necessary
        if not self.has_models: self.create_models()

        # 4. Write
        self.write()

        # 5. Deproject
        self.deproject()

        # Clear
        self.clear()

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
        if "deprojections" in kwargs and "maps" in kwargs: raise ValueError("Cannot specify 'deprojections' and 'maps' simultaneously")
        if "deprojections" in kwargs and "deprojection" in kwargs: raise ValueError("Cannot specify 'deprojections' and 'deprojection' simultaneously")
        if "deprojections" in kwargs and "map" in kwargs: raise ValueError("Cannot specify 'deprojections' and 'map' simultaneously")

        # Get the deprojection
        if "deprojection" in kwargs:
            if "name" not in kwargs: raise ValueError("When passing only one deprojection, a name must be specified")
            name = kwargs.pop("name")
            self.deprojections = dict()
            self.deprojections[name] = kwargs.pop("deprojection")

        # Get multiple deprojections
        elif "deprojections" in kwargs:
            self.deprojections = kwargs.pop("deprojections")

        # Load maps
        else: self.load_maps(**kwargs)

        # Make directories
        if self.root_path is not None: self.create_directories()

        # Check leftover arguments
        if len(kwargs) > 0: raise ValueError("Could not resolve all input: " + tostr(kwargs))

    # -----------------------------------------------------------------

    def load_maps(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Checks
        if "map" in kwargs and "maps" in kwargs: raise ValueError("Cannot specify 'map' and 'maps' simultaneously")

        # Get the single map
        if "map" in kwargs:
            if "name" not in kwargs: raise ValueError("When passing only one map, a name must be specified")
            name = kwargs.pop("name")
            self.maps[name] = kwargs.pop("map")

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

    # -----------------------------------------------------------------

    @property
    def has_models(self):

        """
        This function ...
        :return:
        """

        return self.deprojections is not None

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

        # Initialize the dictionary
        self.deprojections = dict()

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

            # Set the map path in the deprojection model
            self.deprojections[name].filename = path

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

        # Loop over the components
        for name in self.deprojections:

            # Debugging
            log.debug("Computing the deprojected surface density of the '" + name + "' map ...")

            # Determine the number of pixels
            npixels = int(round(float(max(self.maps[name].xsize, self.maps[name].ysize)) / self.config.downsample_factor))
            shape = (npixels, npixels)

            # Debugging
            log.debug("Using " + str(npixels) + " x " + str(npixels) + " pixels for the deprojected map")

            # Set the model limits
            x_range_scalar = self.deprojections[name].x_range.to(unit).value * 2.
            y_range_scalar = self.deprojections[name].y_range.to(unit).value * 2.
            #z_range_scalar = self.deprojections[name].z_range.to(unit).value.as_tuple()
            limits = [x_range_scalar.as_tuple(), y_range_scalar.as_tuple()]

            # Create coordinate data
            x, y = xy(shape=shape, limits=limits)

            # Calculate the surface density
            density = self.deprojections[name].surface_density_function(normalize=True)(x, y)

            # Create deprojected map
            deprojected = Frame(density)

            # Set
            self.deprojected[name] = deprojected

    # -----------------------------------------------------------------

    # @lazyproperty
    # def stellar_ski_template(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     # Load ski template
    #     ski_template = LabeledSkiFile(template_ski_path)
    #
    #     # Convert to oligochromatic simulation
    #     ski_template.to_oligochromatic(1. * u("micron"))
    #
    #     # Remove the dust system
    #     ski_template.remove_dust_system()
    #
    #     # Set number of packages per wavelength
    #     ski_template.setpackages(0)
    #
    #     # Add one instrument
    #     ski_template.remove_all_instruments()
    #
    #     # Return
    #     return ski_template

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        ski =  self.smile.create_oligochromatic_template()

        # Remove the existing instruments
        ski.remove_all_instruments()

        # Remove the stellar system
        ski.remove_stellar_system()

        # Set the number of photon packages
        ski.setpackages(0)

        # Enable writing options
        #ski.enable_all_writing_options()

        # Disable writing stellar density (we don't have a stellar system)
        # BUT DON'T CALL THE FUNCTION WHEN THE SKIRT VERSION DOES NOT SUPPORT WRITING STELLAR DENSITY
        #if self.smile.supports_writing_stellar_density: ski.set_write_stellar_density(False)

        # Return the ski template
        return ski

    # -----------------------------------------------------------------

    def deproject_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting by launching SKIRT simulations ...")

        # Create the dust grids
        self.create_dust_grids()

        # Create the ski files
        self.create_ski_files()

        # Write the ski files
        self.write_ski_files()

        # Launch SKIRT
        self.launch()

    # -----------------------------------------------------------------

    def create_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grids ...")

        # Loop over the deprojections
        for name in self.deprojections:

            deprojection = self.deprojections[name]

            # Set minimum level
            if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_min_level
            elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_min_level
            else: min_level = None

            # Set max ndivisions per pixel
            max_ndivisions_per_pixel = 1. / self.config.dg.rel_scale  # default 1/0.5 = 2 divisions along each direction per pixel

            # Create the dust grid
            # grid_type, deprojection, distance, sky_ellipse, min_level, max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.
            dust_grid = create_one_dust_grid_for_galaxy_from_deprojection(self.config.dg.grid_type, deprojection,
                                                                               self.galaxy_distance,
                                                                               self.truncation_ellipse,
                                                                               min_level, self.config.dg.max_mass_fraction,
                                                                               max_ndivisions_per_pixel,
                                                                               self.config.dg.scale_heights)

            # Set the dust grid
            self.dust_grids[name] = dust_grid

    # -----------------------------------------------------------------

    def create_ski_files(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski files ...")

        # Loop over the components
        for name in self.deprojections:

            # Debugging
            log.debug("Creating a ski file for the deprojection of the '" + name + "' map ...")

            # Make copy
            ski = self.ski_template.copy()

            # FOR STELLAR GEOMETRIES:
            # Add faceon instrument
            #faceon = FrameInstrument.from_deprojection_faceon(self.deprojections[name])
            #ski.add_instrument(faceon_name, faceon)

            # Add edgeon instrument
            #edgeon = FrameInstrument.from_deprojection_edgeon(self.deprojections[name])
            #ski.add_instrument(edgeon_name, edgeon)

            # Set component
            #ski.remove_stellar_components_except("Evolved stellar disk")
            #ski.remove_dust_system()
            #ski.set_stellar_component_geometry("Evolved stellar disk", self.deprojections[name])

            # Add the dust component
            #map_filename = add_dust_component(self.ski, name, component)

            # If map filename is defined, set path in dictionary
            #if map_filename is not None: self.input_map_paths[map_filename] = component.map_path

            # USE DUST GEOMETRY

            # Set filename for deprojection model
            #map_filename = "map.fits"
            deprojection = self.deprojections[name].copy()
            deprojection.filename = fs.name(deprojection.filename)

            # Get title
            title = name

            # Set the geometry
            #ski.set_dust_component_geometry(title, deprojection)

            dust_mass = 1e7 * u("Msun")
            mix = "themis"
            ski.create_new_dust_component(title, deprojection, normalization_value=dust_mass, mix=mix)

            dust_grid = self.dust_grids[name]

            # Set the dust grid
            ski.set_dust_grid(dust_grid)

            # Enable writing options
            ski.enable_all_writing_options()

            # Add the ski file
            self.ski_files[name] = ski

    # -----------------------------------------------------------------

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski files ...")

        # Loop over the ski files
        for name in self.ski_files:

            # Determine path
            filepath = fs.join(self.output_paths[name], name + ".ski")

            # Save the ski file
            self.ski_files[name].saveto(filepath, fix=True)

            # Set the path
            self.ski_paths[name] = filepath

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Loop over the maps
        for name in self.maps:

            # Write the ski file
            #self.write_ski()

            ski_path = self.ski_paths[name]
            out_path = self.output_paths[name]
            #input_map_paths = [self.deprojections[name].filepath]
            in_path = fs.directory_of(self.deprojections[name].filepath)

            # Debugging
            log.debug("Launching SKIRT for the '" + name + "' map ...")

            # Create simulation definition
            definition = SingleSimulationDefinition(ski_path, out_path, in_path)

            # Determine parallelization scheme (do singleprocessing)
            ncores = 2
            nthreads_per_core = 2
            nprocesses = 1
            parallelization = Parallelization(ncores, nthreads_per_core, nprocesses)

            # Set settings
            self.launcher.config.progress_bar = True
            self.launcher.config.finish_after = "Writing dust cell properties"  # finish after this line has been printed (when the next one comes)
            # self.launcher.config.finish_at = ""

            # Run
            self.launcher.run(definition=definition, parallelization=parallelization)
            simulation = self.launcher.simulation

            # FOR STELLAR:
            # Determine path
            #frame_path = fs.join(out_path, simulation.prefix() + "_faceon_total.fits")

            # Open the output frame
            #frame = Frame.from_file(frame_path)

            ### FACEON

            gridxy_filename = simulation.prefix() + "_ds_grhoxy.fits"
            geometryxy_filename = simulation.prefix() + "_ds_trhoxy.fits"

            #grid_xy_path = fs.join(out_path, gridxy_filename)
            geometry_xy_path = fs.join(out_path, geometryxy_filename)

            # Open the output frame
            frame = Frame.from_file(geometry_xy_path)

            # Set wcs NO DOESN'T MAKE ANY SENSE ON THE FACEON MAP!!
            #frame.wcs = self.maps[name].wcs

            # Save frame
            #frame.saveto(frame_path)

            # Set the deprojected map
            self.deprojected[name] = frame

            ### EDGEON

            gridxz_filename = simulation.prefix() + "_ds_grhoxz.fits"
            geometryxz_filename = simulation.prefix() + "_ds_trhoxz.fits"

            geometry_xz_path = fs.join(out_path, geometryxz_filename)

            # Open the frame
            edgeon = Frame.from_file(geometry_xz_path)

            # Set the edgeon map
            self.edgeon[name] = edgeon

    # -----------------------------------------------------------------

    # def launch(self):
    #
    #     """
    #     This function ...
    #     :return:
    #     """
    #
    #     #
    #
    #     # Loop over the ski files
    #     for name in self.ski_files:
    #
    #         #
    #         ski_path = self.ski_paths[name]
    #         out_path = self.output_paths[name]
    #
    #         #prefix = fs.strip_extension(fs.name(ski_path))
    #
    #         # Perform the SKIRT simulation
    #         simulation = self.skirt.execute(ski_path, inpath=self.maps_path, outpath=out_path, single=True)
    #
    #         # Determine path
    #         frame_path = fs.join(out_path, simulation.prefix() + "_faceon_total.fits")
    #
    #         # Open the output frame
    #         frame = Frame.from_file(frame_path)
    #
    #         # Set wcs
    #         frame.wcs = self.reference_wcs
    #
    #         # Save frame
    #         frame.saveto(frame_path)
    #
    #         # Set the deprojected map
    #         self.deprojected[name] = frame

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing ...")

        # Clear maps if not desired as output
        if not self.config.writing.maps: self.clear_maps()

    # -----------------------------------------------------------------

    def clear_maps(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.maps:

            # Determine the path
            path = fs.join(self.output_paths[name], "map.fits")

            # Remove if existing
            if fs.is_file(path): fs.remove_file(path)

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
