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

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.launch.launcher import SKIRTLauncher
from ...magic.core.frame import Frame
from ..component.galaxy import GalaxyModelingComponent
from ...core.tools.utils import lazyproperty
from ...core.prep.smile import SKIRTSmileSchema
from ...core.simulation.definition import SingleSimulationDefinition
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...magic.basics.coordinate import SkyCoordinate
from ...magic.basics.coordinate import PhysicalCoordinate, PixelCoordinate
from ...core.tools.stringify import tostr
from ..basics.instruments import FrameInstrument
from ..basics.models import DeprojectionModel3D
from ..basics.projection import get_physical_center

# -----------------------------------------------------------------

faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

map_filename = "map.fits"

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
        self.coordinate_systems = None

        ##

        # The models
        self.models = None

        ##

        self.map_aliases = None

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

        # The projected maps
        self.projected = defaultdict(dict)

        ##

        # The input paths (per model)
        self.input_paths = dict()

        # The output paths (per model)
        self.output_paths = dict()

        ##

        # The ski files
        self.ski_files = dict()

        # The ski file paths
        self.ski_paths = dict()

        # The dust grids
        self.dust_grids = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Create the projections
        self.set_projections()

        # 3. Write
        self.write()

        # 4. Project
        self.project()

        # 5. Clear
        self.clear()

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

        # Get the models
        self.models = kwargs.pop("models")

        # Map aliases
        self.map_aliases = kwargs.pop("map_aliases")

        # Check names
        if "name" in kwargs:
            if kwargs["name"] == edgeon_name or kwargs["name"] == faceon_name: raise ValueError("Name cannot be '" + faceon_name + "' or '" + edgeon_name + "'")
        if "coordinate_systems" in kwargs:
            for name in kwargs["coordinate_systems"]:
                if name == edgeon_name or name == faceon_name: raise ValueError("No name can be '" + faceon_name + "' or '" + edgeon_name + "'")
        if "projections" in kwargs:
            for name in kwargs["projections"]:
                if name == edgeon_name or name == faceon_name: raise ValueError("No name can be '" + faceon_name + "' or '" + edgeon_name + "'")

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
        else: log.warning("Root path not specified")

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
        for name in self.model_names:

            # Debugging
            log.debug("Creating output path for the '" + name + "' model ...")

            # Set output path
            self.output_paths[name] = fs.create_directory_in(self.root_path, name)

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
    def has_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        return self.coordinate_systems is not None

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
        if not self.has_projections and self.has_coordinate_systems: self.set_other_projections()

    # -----------------------------------------------------------------

    @property
    def center_sky(self):

        """
        Thisfunction
        :return:
        """

        if not self.has_wcs: raise RuntimeError("Cannot calculate the sky center coordinate: no coordinate system")
        return SkyCoordinate.from_pixel(self.config.center, self.wcs)

    # -----------------------------------------------------------------

    @property
    def center_physical(self):

        """
        This function ...
        :return:
        """

        #if not self.has_wcs: raise RuntimeError("Cannot calculate the physical center coordinate: no coordinate system")

        if self.has_wcs: return PhysicalCoordinate.from_pixel(self.config.center, self.wcs, self.galaxy_distance, from_center=True)
        else: return get_physical_center(self.config.field, self.config.npixels, self.config.center)

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
        elif self.has_basic_projection_properties: self.faceon_projection = FaceOnProjection(distance=self.config.distance, pixels_x=self.config.npixels.x,
                                                                                             pixels_y=self.config.npixels.y, center_x=self.center_physical.x,
                                                                                             center_y=self.center_physical.y, field_x=self.config.field.x,
                                                                                             field_y=self.config.field.y)
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
        elif self.has_basic_projection_properties: self.edgeon_projection = EdgeOnProjection(distance=self.config.distance, pixels_x=self.config.npixels.x,
                                                                                             pixels_y=self.config.npixels.y, center_x=self.center_physical.x,
                                                                                             center_y=self.center_physical.y, field_x=self.config.field.x,
                                                                                             field_y=self.config.field.y)
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
        if self.config.writing.projections: self.write_projections()

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
        if self.has_projections: self.write_other_projections()

    # -----------------------------------------------------------------

    def write_faceon_projection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the face-on projection ...")

        # Loop over the model names
        for name in self.models:

            # Determine the path
            path = fs.join(self.output_paths[name], "faceon.proj")

            # Save the face-on projection
            self.faceon_projection.saveto(path)

    # -----------------------------------------------------------------

    def write_edgeon_projection(self):

        """
        This funtion ...
        :return:
        """

        # Inform the user
        log.info("Writing the edge-on projection ...")

        # Loop over the model names
        for name in self.models:

            # Determine the path
            path = fs.join(self.output_paths[name], "edgeon.proj")

            # Save the edge-on projection
            self.edgeon_projection.saveto(path)

    # -----------------------------------------------------------------

    def write_other_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the other projections ...")

        # Loop over the model names
        for name in self.models:

            # Loop over the projections
            for projection_name in self.projections:

                # Determine the path
                path = fs.join(self.output_paths[name], projection_name + ".proj")

                # Save the projection
                self.projections[projection_name].saveto(path)

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

        # Debugging
        log.debug("Creating the ski file template ...")

        # Create
        ski =  self.smile.create_oligochromatic_template()

        # Remove the existing instruments
        ski.remove_all_instruments()

        # Remove the stellar system
        #ski.remove_stellar_system()

        # Remove the dust system
        if ski.has_dust_system: ski.remove_dust_system()

        # Set the number of photon packages
        ski.setpackages(self.config.npackages)

        # Enable writing options
        #ski.enable_all_writing_options()

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

        # Set the input paths per model
        self.set_input_paths()

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
        if self.has_projections: self.create_other_instruments()

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

    def set_input_paths(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the ski input paths ...")

        # Loop over the models
        for name in self.models:

            # Only for deprojection models
            if isinstance(self.models[name], DeprojectionModel3D):

                # Map is present
                if self.models[name].has_map:

                    # .. because file is present?
                    map_path = self.models[name].filepath
                    if fs.is_file(map_path):

                        inpath = fs.directory_of(map_path) # directory of input map
                        filename = fs.name(map_path)

                    # .. or because map is loaded
                    elif self.models[name].map_is_loaded:

                        inpath = self.output_paths[name]
                        map_path = fs.join(inpath, map_filename)
                        filename = map_filename
                        self.models[name].map.saveto(map_path) # save the map

                    # We shouldn't get here
                    else: raise RuntimeError("We shouldn't get here")

                # Map is defined in map aliases
                elif self.map_aliases is not None:

                    # Search through the aliases
                    map_name = self.models[name].filename
                    if map_name in self.map_aliases:

                        the_map = self.map_aliases[map_name]
                        inpath = self.output_paths[name]
                        map_path = fs.join(inpath, map_filename)
                        filename = map_filename
                        the_map.saveto(map_path) # save the map

                    # Not defined in aliases
                    else: raise ValueError("Map alias '" + map_name + "' not defined in map aliases")

                # Cannot find the map
                else: raise ValueError("Map does not exist")

                # Set the model map filename
                self.models[name].filename = filename

            # Other kinds of models
            else: inpath = None

            # Set the input path for this model
            self.input_paths[name] = inpath

    # -----------------------------------------------------------------

    def create_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski files ...")

        # Loop over the components
        for name in self.models:

            # Debugging
            log.debug("Creating a ski file for the projection of the '" + name + "' map ...")

            # Make copy
            ski = self.ski_template.copy()

            # Add faceon
            if self.has_faceon: ski.add_instrument(faceon_name, self.faceon_instrument)

            # Add edgeon
            if self.has_edgeon: ski.add_instrument(edgeon_name, self.edgeon_instrument)

            # Add other
            if self.has_projections:
                for instrument_name in self.instruments: ski.add_instrument(instrument_name, self.instruments[name])

            # Remove the dust system
            if ski.has_dust_system: ski.remove_dust_system()

            # Set title
            title = name

            # Add the stellar component
            ski.create_new_stellar_component(title, geometry=self.models[name], luminosities=[1])

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

        #projected_models = dict()

        # Loop over the models
        for name in self.models:

            # Determine ski path and output path for the simulation
            ski_path = self.ski_paths[name]
            out_path = self.output_paths[name]
            in_path = self.input_paths[name]

            # Debugging
            log.debug("Launching SKIRT for the '" + name + "' map ...")

            # Create simulation definition
            definition = SingleSimulationDefinition(ski_path, out_path, in_path)

            # Set settings
            self.launcher.config.show_progress = True

            # Run
            self.launcher.run(definition=definition, parallelization=self.config.parallelization)
            simulation = self.launcher.simulation

            # Get the faceon map
            if self.has_faceon:
                faceon_path = fs.join(out_path, simulation.prefix() + "_" + faceon_name + "_total.fits")
                faceon = Frame.from_file(faceon_path)
            else: faceon = None

            # Get the edgeon map
            if self.has_edgeon:
                edgeon_path = fs.join(out_path, simulation.prefix() + "_" + edgeon_name + "_total.fits")
                edgeon = Frame.from_file(edgeon_path)
            else: edgeon = None

            # Get the other maps
            #other = dict()
            if self.has_projections:
                for instrument_name in self.instruments:
                    other_path = fs.join(out_path, simulation.prefix() + "_" + instrument_name + "_total.fits")
                    other_map = Frame.from_file(other_path)
                    if self.has_coordinate_systems: other_map.wcs = self.coordinate_systems[instrument_name] # set WCS
                    #other[instrument_name] = other_map
                    # Set the map
                    self.projected[instrument_name][name] = other_map

            #projected_models[name] = other

            # Set faceon and edgeon
            self.faceon[name] = faceon
            self.edgeon[name] = edgeon

        # Set the projected maps for each model per projection
        # ACTUALLY, because self.projected is defaultdict, we can as well do all of this inside the loop above..
        # for projection_name in self.projections:
        #     for model_name in projected_models:
        #         self.projected[projection_name][model_name] = projected_models[model_name][projection_name]

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing ...")

        # Clear maps
        # Loop over the models
        for name in self.models:

            inpath = self.output_paths[name]
            map_path = fs.join(inpath, map_filename)
            if fs.is_file(map_path): fs.remove_file(map_path)

# -----------------------------------------------------------------
