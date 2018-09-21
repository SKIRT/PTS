#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.projection.model Contains the ComponentProjections class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import createsimulations, SkirtSimulation, StaticSkirtSimulation
from ...core.simulation.output import SimulationOutput
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.execute import run_simulation
from ...core.prep.smile import get_oligochromatic_template
from ...core.tools import introspection
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from ..basics.instruments import FrameInstrument, FullInstrument
from ..build.representations.galaxy import create_projection_from_deprojection, create_faceon_projection_from_earth_projection, create_edgeon_projection_from_earth_projection
from ..basics.models import DeprojectionModel3D
from .general import default_scale_heights

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

output_filename = "output.txt"

# -----------------------------------------------------------------

def project_model(name, model, projection=None, description=None, path=None, center=None, wcs=None, distance=None):

    """
    This function ...
    :param name:
    :param model:
    :param projection:
    :param description:
    :param path:
    :param center:
    :param wcs:
    :param distance:
    :return:
    """

    # Create the projections object
    projections = ComponentProjections(name, model, path=path, description=description,
                                projection=projection, center=center,
                                earth_wcs=wcs, distance=distance, earth=True, faceon=False, edgeon=False)

    # Return the projected map
    return projections.earth

# -----------------------------------------------------------------

class ComponentProjections(object):

    """
    This class ...
    """

    def __init__(self, name, model, projection=None, projection_faceon=None, projection_edgeon=None,
                 path=None, earth=True, faceon=True, edgeon=True, npackages=default_npackages,
                 description=None, input_filepaths=None, distance=None, wcs=None, center=None, radial_factor=1,
                 earth_wcs=None, scale_heights=default_scale_heights):

        """
        This function ...
        :param name: name of the component
        :param model: the geometric model of the component
        :param projection: can be created from deprojection model (if model IS), or from wcs if given
        :param path: path for the projections, if none is given a temporary directory is created
        :param earth: perform earth projection
        :param faceon: perform faceon projection
        :param npackages:
        :param description:
        :param input_filepaths:
        :param distance:
        :param wcs:
        :param center: galaxy center, as sky coordinate
        :param radial_factor:
        :param earth_wcs:
        :param
        """

        # Set the name and description
        self.name = name
        self.description = description

        # Set the model
        self.model = model

        # Extra info: set lazy properties
        if distance is not None: self.distance = distance
        if wcs is not None: self.wcs = wcs
        if center is not None: self.center = center

        # Set the earth projection (if necessary)
        if projection is None: projection = self.create_projection_earth()
        self.projection_earth = projection

        # Set the face-on projection (if necessary)
        if projection_faceon is None and faceon: projection_faceon = self.create_projection_faceon(radial_factor=radial_factor)
        self.projection_faceon = projection_faceon

        # Set the edge-on projection (if necessary)
        if projection_edgeon is None and edgeon: projection_edgeon = self.create_projection_edgeon(radial_factor=radial_factor, scale_heights=scale_heights)
        self.projection_edgeon = projection_edgeon

        # Set the path
        if path is None: path = introspection.create_temp_dir("projections__" + name)
        self.path = path

        # Set the number of photon packages
        self.npackages = npackages

        # Set the input filepaths
        self.input_paths = input_filepaths

        # Set the WCS for the earth map
        self.earth_wcs = earth_wcs

        # Get the earth simulation
        if earth: self.simulation_earth = self.get_simulation_earth()
        else: self.simulation_earth = None

        # Get the face-on simulation
        if faceon: self.simulation_faceon = self.get_simulation_faceon()
        else: self.simulation_faceon = None

        # Get the edge-on simulation
        if edgeon: self.simulation_edgeon = self.get_simulation_edgeon()
        else: self.simulation_edgeon = None

    # -----------------------------------------------------------------

    @property
    def has_deprojection(self):
        return isinstance(self.model, DeprojectionModel3D)

    # -----------------------------------------------------------------

    @property
    def deprojection(self):

        """
        This function ...
        :return:
        """

        if self.has_deprojection: return self.model
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_map(self):
        return self.has_deprojection and self.deprojection.has_map

    # -----------------------------------------------------------------

    @lazyproperty
    def distance(self):

        """
        This function ...
        :return:
        """

        if not self.has_deprojection: return None
        else: return self.deprojection.distance

    # -----------------------------------------------------------------

    @property
    def has_distance(self):
        #return self.has_deprojection and self.distance is not None
        return self.distance is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def map(self):

        """
        This function ...
        :return:
        """

        if not self.has_map: return None
        else: return self.deprojection.map

    # -----------------------------------------------------------------

    @lazyproperty
    def wcs(self):

        """
        This function ...
        :return:
        """

        if not self.has_map: return None
        else: return self.map.wcs

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):
        #return self.has_map and self.wcs is not None
        return self.wcs is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def inclination(self):

        """
        This function ...
        :return:
        """

        if not self.has_deprojection: return None
        else: return self.deprojection.inclination

    # -----------------------------------------------------------------

    @property
    def has_inclination(self):
        return self.inclination is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def position_angle(self):

        """
        This function ...
        :return:
        """

        if not self.has_deprojection: return None
        else: return self.deprojection.position_angle

    # -----------------------------------------------------------------

    @property
    def has_position_angle(self):
        return self.position_angle is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def center(self):

        """
        This function ...
        :return:
        """

        if not (self.has_wcs and self.has_deprojection): return None
        else: return self.deprojection.center.to_sky(self.wcs)

    # -----------------------------------------------------------------

    @property
    def has_center(self):
        return self.center is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def scaleheight(self):

        """
        This function ...
        :return:
        """

        if self.has_deprojection: return self.deprojection.scale_height
        elif hasattr(self.model, "axial_scale"): return self.model.axial_scale
        elif hasattr(self.model, "effective_radius"):
            radius = self.model.effective_radius
            if hasattr(self.model, "z_flattening"):
                flattening = self.model.z_flattening
                return radius * flattening
            else: return radius
        else: return None

    # -----------------------------------------------------------------

    @property
    def has_scaleheight(self):
        return self.scaleheight is not None

    # -----------------------------------------------------------------

    def create_projection_earth(self):

        """
        This function ...
        :return:
        """

        # Check whether deprojection is defined
        if self.has_deprojection: return create_projection_from_deprojection(self.deprojection)
        else:

            # Checks
            if not self.has_wcs: raise ValueError("Coordinate system is not defined")
            if not self.has_distance: raise ValueError("Galaxy distance is not defined")
            if not self.has_center: raise ValueError("Galaxy center coordinate is not defined")
            if not self.has_inclination: raise ValueError("Galaxy inclination is not defined")
            if not self.has_position_angle: raise ValueError("Galaxy position angle is not defined")

            azimuth = 0.0
            return GalaxyProjection.from_wcs(self.wcs, self.center, self.distance, self.inclination, azimuth, self.position_angle)

    # -----------------------------------------------------------------

    def create_projection_faceon(self, radial_factor=1):

        """
        This function ...
        :param radial_factor:
        :return:
        """

        return create_faceon_projection_from_earth_projection(self.projection_earth, radial_factor=radial_factor)

    # -----------------------------------------------------------------

    def create_projection_edgeon(self, radial_factor=1, scale_heights=default_scale_heights):

        """
        This function ...
        :param radial_factor:
        :param scale_heights:
        :return:
        """

        z_extent = 2. * self.scaleheight * scale_heights
        return create_edgeon_projection_from_earth_projection(self.projection_earth, z_extent, radial_factor=radial_factor)

    # -----------------------------------------------------------------

    @property
    def has_earth(self):
        return self.simulation_earth is not None

    # -----------------------------------------------------------------

    @property
    def has_faceon(self):
        return self.simulation_faceon is not None

    # -----------------------------------------------------------------

    @property
    def has_edgeon(self):
        return self.simulation_edgeon is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_path(self):
        return fs.create_directory_in(self.path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_out_path(self):
        return fs.create_directory_in(self.earth_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_path(self):
        return fs.create_directory_in(self.path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_out_path(self):
        return fs.create_directory_in(self.faceon_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_path(self):
        return fs.create_directory_in(self.path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_out_path(self):
        return fs.create_directory_in(self.edgeon_path, "out")

    # -----------------------------------------------------------------

    @property
    def has_earth_simulation(self):
        return fs.has_files_in_path(self.earth_out_path, extension="fits")

    # -----------------------------------------------------------------

    @property
    def has_faceon_simulation(self):
        return fs.has_files_in_path(self.faceon_out_path, extension="fits")

    # -----------------------------------------------------------------

    @property
    def has_edgeon_simulation(self):
        return fs.has_files_in_path(self.edgeon_out_path, extension="fits")

    # -----------------------------------------------------------------

    @property
    def earth_output_filepath(self):
        return fs.join(self.earth_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_earth_output_file(self):
        return fs.is_file(self.earth_output_filepath)

    # -----------------------------------------------------------------

    def get_simulation_earth(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        if self.has_earth_simulation:

            simulation = StaticSkirtSimulation(prefix=self.simulation_prefix, inpath=self.input_paths, outpath=self.earth_out_path, ski_path=self.ski_path_earth, name=self.simulation_name)
            if self.has_earth_output_file: simulation.output = SimulationOutput.from_file(self.earth_output_filepath) # static simulation so lazy property
            return simulation

        # Run the simulation
        else: return self.run_earth_simulation()

    # -----------------------------------------------------------------

    @property
    def faceon_output_filepath(self):
        return fs.join(self.faceon_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_faceon_output_file(self):
        return fs.is_file(self.faceon_output_filepath)

    # -----------------------------------------------------------------

    def get_simulation_faceon(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        if self.has_faceon_simulation:

            simulation = StaticSkirtSimulation(prefix=self.simulation_prefix, inpath=self.input_paths, outpath=self.faceon_out_path, ski_path=self.ski_path_faceon, name=self.simulation_name)
            if self.has_faceon_output_file: simulation.output = SimulationOutput.from_file(self.faceon_output_filepath)
            return simulation

        # Run the simulation
        else: return self.run_faceon_simulation()

    # -----------------------------------------------------------------

    @property
    def edgeon_output_filepath(self):
        return fs.join(self.edgeon_path, output_filename)

    # -----------------------------------------------------------------

    @property
    def has_edgeon_output_file(self):
        return fs.is_file(self.edgeon_output_filepath)

    # -----------------------------------------------------------------

    def get_simulation_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        # Simulation already performed?
        if self.has_edgeon_simulation:

            simulation = StaticSkirtSimulation(prefix=self.simulation_prefix, inpath=self.input_paths, outpath=self.edgeon_out_path, ski_path=self.ski_path_edgeon, name=self.simulation_name)
            if self.has_edgeon_output_file: simulation.output = SimulationOutput.from_file(self.edgeon_output_filepath)
            return simulation

        # Run the simulation
        else: return self.run_edgeon_simulation()

    # -----------------------------------------------------------------

    @lazyproperty
    def description_string(self):

        """
        This function ...
        :return:
        """

        if self.description is None: return ""
        else: return self.description + " "

    # -----------------------------------------------------------------

    def run_earth_simulation(self):

        """
        This function ...
        :return:
        """

        # Show message
        log.info("Running SKIRT for the " + self.description_string + "earth projection ...")

        # Run simulation
        return run_simulation(self.earth_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    def run_faceon_simulation(self):

        """
        This function ...
        :return:
        """

        # Show message
        log.info("Running SKIRT for the " + self.description_string + "face-on projection ...")

        # Run simulation
        return run_simulation(self.faceon_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    def run_edgeon_simulation(self):

        """
        This function ...
        :return:
        """

        # Show message
        log.info("Running SKIRT for the " + self.description_string + "edge-on projection ...")

        # Run simulation
        return run_simulation(self.edgeon_definition, show_progress=True, debug_output=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_skifile_earth: self.create_skifile_earth()

        # Create the definition and return
        return SingleSimulationDefinition(self.ski_path_earth, self.earth_out_path, input_path=self.input_paths, name=self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_skifile_faceon: self.create_skifile_faceon()

        # Create the definition and return
        return SingleSimulationDefinition(self.ski_path_faceon, self.faceon_out_path, input_path=self.input_paths, name=self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_definition(self):

        """
        This function ...
        :return:
        """

        # Create the skifile if necessary
        if not self.has_skifile_edgeon: self.create_skifile_edgeon()

        # Create the definition and return
        return SingleSimulationDefinition(self.ski_path_edgeon, self.edgeon_out_path, input_path=self.input_paths, name=self.name)

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):
        return self.name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):
        return self.name

    # -----------------------------------------------------------------

    @property
    def ski_path_earth(self):
        return fs.join(self.earth_path, self.simulation_prefix + ".ski")

    # -----------------------------------------------------------------

    @property
    def ski_path_faceon(self):
        return fs.join(self.faceon_path, self.simulation_prefix + ".ski")

    # -----------------------------------------------------------------

    @property
    def ski_path_edgeon(self):
        return fs.join(self.edgeon_path, self.simulation_prefix + ".ski")

    # -----------------------------------------------------------------

    @property
    def has_skifile_earth(self):
        return fs.is_file(self.ski_path_earth)

    # -----------------------------------------------------------------

    @property
    def has_skifile_faceon(self):
        return fs.is_file(self.ski_path_faceon)

    # -----------------------------------------------------------------

    @property
    def has_skifile_edgeon(self):
        return fs.is_file(self.ski_path_edgeon)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_earth(self):
        return FrameInstrument.from_projection(self.projection_earth)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_faceon(self):
        return FrameInstrument.from_projection(self.projection_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_edgeon(self):
        return FrameInstrument.from_projection(self.projection_edgeon)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_earth(self):
        return FullInstrument.from_projection(self.projection_earth)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_faceon(self):
        return FullInstrument.from_projection(self.projection_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_edgeon(self):
        return FullInstrument.from_projection(self.projection_edgeon)

    # -----------------------------------------------------------------

    def create_skifile_earth(self):

        """
        This function ...
        :return:
        """

        # Create a ski template
        ski = get_oligochromatic_template()

        # Add the old stellar bulge component
        ski.create_new_stellar_component(component_id=self.name, geometry=self.model, luminosities=[1])

        # Add the instrument
        ski.add_instrument(earth_name, self.frame_instrument_earth)

        # Set the number of photon packages
        ski.setpackages(self.npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.ski_path_earth, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    def create_skifile_faceon(self):

        """
        This function ...
        :return:
        """

        # Create a ski template
        ski = get_oligochromatic_template()

        # Add the old stellar bulge component
        ski.create_new_stellar_component(component_id=self.name, geometry=self.model, luminosities=[1])

        # Add the instrument
        ski.add_instrument(faceon_name, self.frame_instrument_faceon)

        # Set the number of photon packages
        ski.setpackages(self.npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.ski_path_faceon, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    def create_skifile_edgeon(self):

        """
        This function ...
        :return:
        """

        # Create a ski template
        ski = get_oligochromatic_template()

        # Add the old stellar bulge component
        ski.create_new_stellar_component(component_id=self.name, geometry=self.model, luminosities=[1])

        # Add the instrument
        ski.add_instrument(edgeon_name, self.frame_instrument_edgeon)

        # Set the number of photon packages
        ski.setpackages(self.npackages)

        # Remove the dust system
        ski.remove_dust_system()

        # Save the skifile
        ski.saveto(self.ski_path_edgeon, fix=True)

        # Return the skifile
        return ski

    # -----------------------------------------------------------------

    @lazyproperty # lazy because only check and save one time to file (and runned simulations are not static simulations, so output is just a property)
    def earth_projection_output(self):
        output = self.simulation_earth.output
        if not self.has_earth_output_file: output.saveto(self.earth_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection_output(self):
        output = self.simulation_faceon.output
        if not self.has_faceon_output_file: output.saveto(self.faceon_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection_output(self):
        output = self.simulation_edgeon.output
        if not self.has_edgeon_output_file: output.saveto(self.edgeon_output_filepath)
        return output

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_map_path(self):
        return self.earth_projection_output.single_total_images

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_map_path(self):
        return self.faceon_projection_output.single_total_images

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_map_path(self):
        return self.edgeon_projection_output.single_total_images

    # -----------------------------------------------------------------

    @property
    def has_earth_wcs(self):
        return self.earth_wcs is not None

    # -----------------------------------------------------------------

    @property
    def pixelscale_earth(self):
        return self.projection_earth.pixelscale

    # -----------------------------------------------------------------

    @property
    def pixelscale_faceon(self):
        return self.projection_faceon.pixelscale

    # -----------------------------------------------------------------

    @property
    def pixelscale_edgeon(self):
        return self.projection_edgeon.pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def earth(self):

        """
        This function ...
        :return:
        """

        # Get the coordinate system
        if self.has_wcs: wcs = self.wcs
        elif self.has_earth_wcs: wcs = self.earth_wcs
        else: raise ValueError("No WCS info for the earth map")

        # Create the frame
        return Frame.from_file(self.earth_map_path, wcs=wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon(self):
        return Frame.from_file(self.faceon_map_path, distance=self.distance, pixelscale=self.pixelscale_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon(self):
        return Frame.from_file(self.edgeon_map_path, distance=self.distance, pixelscale=self.pixelscale_edgeon)

# -----------------------------------------------------------------
