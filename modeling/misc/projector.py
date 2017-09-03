#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.deprojection.projector Contains the Projector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.launch.launcher import SKIRTLauncher
from ...magic.core.frame import Frame
from ...core.units.parsing import parse_unit as u
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.parallelization import Parallelization
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...magic.basics.coordinate import SkyCoordinate
from ...core.tools.stringify import tostr
from ..basics.instruments import FrameInstrument

# -----------------------------------------------------------------

faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

class Projector(GalaxyModelingComponent):
    
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
        super(Projector, self).__init__(*args, **kwargs)

        # Create the SKIRT launcher
        self.launcher = SKIRTLauncher()

        # Smile
        self.smile = SKIRTSmileSchema()

        ##

        # WCS
        self.wcs = None

        # Coordinate systems for creating different projections
        self.coordinate_systems = dict()

        ##

        # The models
        self.models = dict()

        ##

        # The projections
        self.projections = None
        self.faceon_projection = None
        self.edgeon_projection = None

        ##

        # The instruments
        self.instruments = dict()
        self.faceon_instrument = None
        self.edgeon_instrument = None

        ##

        # Paths
        self.root_path = None

        # The face-on maps
        self.faceon = dict()

        # The edge-on maps
        self.edgeon = dict()

        # The deprojected maps
        self.projected = dict()

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

        # Create the projections
        self.set_projections()

        # 4. Write
        self.write()

        # 5. Project
        self.project()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Projector, self).setup(**kwargs)

        # Get paths
        self.root_path = kwargs.pop("root_path", None)

        # Check
        if "name" in kwargs and "coordinate_systems" in kwargs:
            raise ValueError("Cannot specify 'name' and 'coordinate_systems' simultaneously")

        # Check
        if "name" in kwargs and "projections" in kwargs:
            raise ValueError("Cannot specify 'name' and 'projections' simultaneously")

        # Check
        if "projections" in kwargs and "coordinate_systems" in kwargs:
            raise ValueError("Cannot specify 'projections' and 'coordinate_systems' simultaneously")

        # Check
        if "wcs" in kwargs and "coordinate_system" in kwargs: raise ValueError("Cannot specify 'wcs' and 'coordinate_system': these represent the same input")

        # Check whether WCS is passed
        if "wcs" in kwargs or "coordinate_system" in kwargs:

            # Use the WCS to create faceon and edgeon projection
            wcs = kwargs.pop("wcs") if "wcs" in kwargs else kwargs.pop("coordinate_system")
            self.wcs = wcs

            # Check whether name is defined, if so, the WCS is also used to create a seperate projection
            if "name" in kwargs: self.coordinate_systems[kwargs.pop("name")] = wcs

            # Not not specified
            else: log.warning("Name for the coordinate system is not defined. WCS will only be used to create face-on and/or edge-on projections")

        # Get coordinate systems
        if "coordinate_systems" in kwargs: self.coordinate_systems = kwargs.pop("coordinate_systems")

        # Get projections
        if "projections" in kwargs: self.projections = kwargs.pop("projections")

        # Make directories
        if self.root_path is not None: self.create_directories()

        # Check leftover arguments
        if len(kwargs) > 0: raise ValueError("Could not resolve all input: " + tostr(kwargs))

    # -----------------------------------------------------------------

    def create_directories(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the directories ...")

        # Loop over the maps
        for name in self.model_names: self.output_paths[name] = fs.create_directory_in(self.root_path, name)

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):

        """
        This function ...
        :return:
        """

        return self.wcs is not None

    # -----------------------------------------------------------------

    @property
    def model_names(self):

        """
        This function ...
        :return:
        """

        return self.models.keys()

    # -----------------------------------------------------------------

    @property
    def nmodels(self):

        """
        This function ...
        :return:
        """

        return len(self.models)

    # -----------------------------------------------------------------

    @property
    def no_models(self):

        """
        This function ...
        :return:
        """

        return self.nmodels == 0

    # -----------------------------------------------------------------

    @property
    def has_single_model(self):

        """
        This function ...
        :return:
        """

        return self.nmodels == 1

    # -----------------------------------------------------------------

    @property
    def single_model_name(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_model: raise ValueError("Not a single model")
        return self.model_names[0]

    # -----------------------------------------------------------------

    @property
    def projection_names(self):

        """
        This function ...
        :return:
        """

        return self.projections.keys()

    # -----------------------------------------------------------------

    @property
    def nprojections(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.projection_names)

    # -----------------------------------------------------------------

    @property
    def no_projections(self):

        """
        This function ...
        :return:
        """

        return self.nprojections == 0

    # -----------------------------------------------------------------

    @property
    def has_single_projection(self):

        """
        This function ...
        :return:
        """

        return self.nprojections == 1

    # -----------------------------------------------------------------

    @property
    def single_projection_name(self):

        """
        This function ...
        :return:
        """

        if not self.has_single_projection: raise ValueError("Not a single projection")
        return self.projection_names[0]

    # -----------------------------------------------------------------

    @property
    def single_model(self):

        """
        This function ...
        :return:
        """

        return self.models[self.single_model_name]

    # -----------------------------------------------------------------

    @property
    def single_projected(self):

        """
        This function ...
        :return:
        """

        return self.projected[self.single_projection_name][self.single_model_name]

    # -----------------------------------------------------------------

    @property
    def single_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.edgeon[self.single_model_name]

    # -----------------------------------------------------------------

    @property
    def single_faceon(self):

        """
        This function ...
        :return:
        """

        return self.faceon[self.single_model_name]

    # -----------------------------------------------------------------

    @property
    def has_projections(self):

        """
        Thisf unction ...
        :return:
        """

        return self.projections is not None

    # -----------------------------------------------------------------

    def set_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the projections ...")

        # Create faceon projection
        if self.config.faceon: self.set_faceon_projection()

        # Set edgeon projection
        if self.config.edgeon: self.set_edgeon_projection()

        # Create other projections
        if not self.has_projections: self.set_other_projections()

    # -----------------------------------------------------------------

    @property
    def center_sky(self):

        """
        Thisfunction
        :return:
        """

        if not self.has_wcs: raise RuntimeError("Cannot calculate the sky center coordinate: no WCS")
        return SkyCoordinate.from_pixel(self.config.center, self.wcs)

    # -----------------------------------------------------------------

    @property
    def has_distance(self):

        """
        This function ...
        :return:
        """

        return self.config.distance is not None

    # -----------------------------------------------------------------

    @property
    def has_center(self):

        """
        This function ...
        :return:
        """

        return self.config.center is not None

    # -----------------------------------------------------------------

    @property
    def has_npixels(self):

        """
        This function ...
        :return:
        """

        return self.config.npixels is not None

    # -----------------------------------------------------------------

    @property
    def has_field(self):

        """
        This function ...
        :return:
        """

        return self.config.field is not None

    # -----------------------------------------------------------------

    @property
    def has_basic_projection_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_distance and self.has_center and self.has_npixels and self.has_field

    # -----------------------------------------------------------------

    def set_faceon_projection(self):

        """
        This function ....
        :return:
        """

        # Inform the user
        log.info("Setting the face-on projection ...")

        # Set the faceon projection
        if self.has_wcs: self.faceon_projection = FaceOnProjection.from_wcs(self.wcs, self.center_sky, self.config.distance)
        elif self.has_basic_projection_properties: self.faceon_projection = FaceOnProjection(self.config.distance, self.config.npixels.x, self.config.npixels.y, self.config.center.x, self.config.center.y, self.config.field.x, self.config.field.y)
        else: raise ValueError("Cannot set the face-on projection: projection properties not set and coordinate system not specified")

    # -----------------------------------------------------------------

    def set_edgeon_projection(self):

        """
        This function ...
        :return:
        """

        # Infrom the user
        log.info("Setting the edge-on projection ...")

        # Set the edgeon projection
        if self.has_wcs: self.edgeon_projection = EdgeOnProjection.from_wcs(self.wcs, self.center_sky, self.config.distance)
        elif self.has_basic_projection_properties: self.edgeon_projection = EdgeOnProjection(self.config.distance, self.config.npixels.x, self.config.npixels.y, self.config.center.x, self.config.center.y, self.config.field.x, self.config.field.y)
        else: raise ValueError("Cannot set the edge-on projection: projection properties not set and coordinate system not specified")

    # -----------------------------------------------------------------

    @property
    def has_inclination(self):

        """
        This function ...
        :return:
        """

        return self.config.inclination is not None

    # -----------------------------------------------------------------

    @property
    def has_azimuth(self):

        """
        This function ...
        :return:
        """

        return self.config.azimuth is not None

    # -----------------------------------------------------------------

    @property
    def has_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.config.position_angle is not None

    # -----------------------------------------------------------------

    @property
    def has_galaxy_projection_properties(self):

        """
        This function ...
        :return:
        """

        return self.has_center and self.has_distance and self.has_inclination and self.has_azimuth and self.has_position_angle

    # -----------------------------------------------------------------

    def set_other_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the other projections ...")

        # Check
        if not self.has_galaxy_projection_properties: raise RuntimeError("Cannot create projections from coordinate systems: galaxy projection properties are not specified")

        # Initialize
        self.projections = dict()

        # Loop over the coordinate systems
        for name in self.coordinate_systems:

            # Get wcs
            wcs = self.coordinate_systems[name]

            # Create projection
            projection = GalaxyProjection.from_wcs(wcs, self.center_sky, self.config.distance, self.config.inclination, self.config.azimuth, self.config.position_angle)

            # Set the projection
            self.projections[name] = projection

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the projections
        self.write_projections()

    # -----------------------------------------------------------------

    @property
    def has_faceon(self):

        """
        This function ...
        :return:
        """

        return self.faceon_projection is not None

    # -----------------------------------------------------------------

    @property
    def has_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.edgeon_projection is not None

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projections ...")

        # Face-on
        if self.has_faceon: self.write_faceon_projection()

        # Edge-on
        if self.has_edgeon: self.write_edgeon_projection()

        # Write other projections
        self.write_other_projections()

    # -----------------------------------------------------------------

    def write_faceon_projection(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_edgeon_projection(self):

        """
        This funtion ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_other_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the other projections ...")

        # Loop over the projections
        for name in self.projections:

            pass

    # -----------------------------------------------------------------

    def project(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the maps ...")

        # Project with SKIRT
        self.project_skirt()

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

    def project_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting by launching SKIRT simulations ...")

        # Create the instruments
        self.create_instruments()

        # Create the ski files
        self.create_ski_files()

        # Write the ski files
        self.write_ski_files()

        # Launch SKIRT
        self.launch()

    # -----------------------------------------------------------------

    def create_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the instruments ...")

        # Face-on
        if self.has_faceon: self.create_faceon_instrument()

        # Edge-on
        if self.has_edgeon: self.create_edgeon_instrument()

        # Other
        self.create_other_instruments()

    # -----------------------------------------------------------------

    def create_faceon_instrument(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the face-on instrument ...")

        # Create the faceon instrument
        self.faceon_instrument = FrameInstrument.from_projection(self.faceon_projection)

    # -----------------------------------------------------------------

    def create_edgeon_instrument(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating the edge-on instrument ...")

        # Create the edgeon instrument
        self.edgeon_instrument = FrameInstrument.from_projection(self.edgeon_projection)

    # -----------------------------------------------------------------

    def create_other_instruments(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating the other instruments ...")

        # Loop over the projections
        for name in self.projections:

            # Create frame instrument
            instrument = FrameInstrument.from_projection(self.projections[name])

            # Set the instrument
            self.instruments[name] = instrument

    # -----------------------------------------------------------------

    def create_ski_files(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski files ...")

        # Loop over the components
        for name in self.models:

            # Debugging
            log.debug("Creating a ski file for the deprojection of the '" + name + "' map ...")

            # Make copy
            ski = self.ski_template.copy()

            # Add faceon
            if self.has_faceon: ski.add_instrument(faceon_name, self.faceon_instrument)

            # Add edgeon
            if self.has_edgeon: ski.add_instrument(edgeon_name, self.edgeon_instrument)

            # Add other
            for name in self.instruments: ski.add_instrument(name, self.instruments[name])

            # Set component
            ski.remove_stellar_components_except("Evolved stellar disk")
            ski.remove_dust_system()
            ski.set_stellar_component_geometry("Evolved stellar disk", self.deprojections[name])

            # Add the stellar component
            ski.create_new_stellar_component()

            # Add the dust component
            #map_filename = add_dust_component(self.ski, name, component)

            # If map filename is defined, set path in dictionary
            #if map_filename is not None: self.input_map_paths[map_filename] = component.map_path

            # Enable writing options
            #ski.enable_all_writing_options()

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

        # Loop over the models
        for name in self.models:

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
            self.projected[name] = frame

            ### EDGEON

            gridxz_filename = simulation.prefix() + "_ds_grhoxz.fits"
            geometryxz_filename = simulation.prefix() + "_ds_trhoxz.fits"

            geometry_xz_path = fs.join(out_path, geometryxz_filename)

            # Open the frame
            edgeon = Frame.from_file(geometry_xz_path)

            # Set the edgeon map
            self.edgeon[name] = edgeon

# -----------------------------------------------------------------
