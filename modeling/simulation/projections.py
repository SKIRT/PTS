#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.simulation.projections Contains the ComponentProjections class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.simulation.simulation import createsimulations
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...magic.core.frame import Frame
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.simulation.execute import run_simulation
from ...core.prep.smile import get_oligochromatic_template
from ...core.tools import introspection
from ..basics.projection import GalaxyProjection
from ..basics.instruments import FrameInstrument, FullInstrument
from ..build.representations.galaxy import create_faceon_projection, create_edgeon_projection, create_projection_from_deprojection
from ..basics.models import DeprojectionModel3D
from ...core.units.parsing import parse_unit as u

# -----------------------------------------------------------------

# Number of photon packages
default_npackages = 5e7

# Instruments/orientations
earth_name = "earth"
faceon_name = "faceon"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

default_scale_heights = 15. # number of times to take the scale height as the vertical radius of the model

# -----------------------------------------------------------------

class ComponentProjections(object):

    """
    This class ...
    """

    def __init__(self, name, model, projection=None, projection_faceon=None, projection_edgeon=None,
                 path=None, earth=True, faceon=True, edgeon=True, npackages=default_npackages,
                 description=None, input_filepaths=None, distance=None, wcs=None, center=None, radial_factor=1,
                 earth_wcs=None):

        """
        This function ...
        :param name: name of the component
        :param model: the geometric model of the component
        :param projection:
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

        # Set the earth projection
        if projection is None: projection = self.create_projection_earth()
        self.projection_earth = projection

        # Set the face-on projection
        if projection_faceon is None: projection_faceon = self.create_projection_faceon(radial_factor=radial_factor)
        self.projection_faceon = projection_faceon

        # Set the edge-on projection
        if projection_edgeon is None: projection_edgeon = self.create_projection_edgeon(radial_factor=radial_factor)
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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

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

        """
        Thisnf unction ...
        :return:
        """

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

    def create_projection_edgeon(self, radial_factor=1):

        """
        This function ...
        :param radial_factor:
        :return:
        """

        return create_edgeon_projection_from_earth_projection(self.projection_earth, self.scaleheight, radial_factor=radial_factor)

    # -----------------------------------------------------------------

    @property
    def has_earth(self):

        """
        This function ...
        :return:
        """

        return self.simulation_earth is not None

    # -----------------------------------------------------------------

    @property
    def has_faceon(self):

        """
        This function ...
        :return:
        """

        return self.simulation_faceon is not None

    # -----------------------------------------------------------------

    @property
    def has_edgeon(self):

        """
        This function ...
        :return:
        """

        return self.simulation_edgeon is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, earth_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.earth_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, faceon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.faceon_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.path, edgeon_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.edgeon_path, "out")

    # -----------------------------------------------------------------

    @property
    def has_earth_simulation(self):

        """
        Thisn function ...
        :return:
        """

        return fs.has_files_in_path(self.earth_out_path, extension="fits")

    # -----------------------------------------------------------------

    @property
    def has_faceon_simulation(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.faceon_out_path, extension="fits")

    # -----------------------------------------------------------------

    @property
    def has_edgeon_simulation(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.edgeon_out_path, extension="fits")

    # -----------------------------------------------------------------

    def get_simulation_earth(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        if self.has_earth_simulation: return createsimulations(self.earth_out_path, single=True)

        # Run the simulation
        else: return self.run_earth_simulation()

    # -----------------------------------------------------------------

    def get_simulation_faceon(self):

        """
        This function ...
        :return:
        """

        # Simulation already performed?
        if self.has_faceon_simulation: return createsimulations(self.faceon_out_path, single=True)

        # Run the simulation
        else: return self.run_faceon_simulation()

    # -----------------------------------------------------------------

    def get_simulation_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

        # Simulation already performed?
        if self.has_edgeon_simulation: return createsimulations(self.edgeon_out_path, single=True)

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
    def ski_path_earth(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.earth_path, self.name + ".ski")

    # -----------------------------------------------------------------

    @property
    def ski_path_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.faceon_path, self.name + ".ski")

    # -----------------------------------------------------------------

    @property
    def ski_path_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.edgeon_path, self.name + ".ski")

    # -----------------------------------------------------------------

    @property
    def has_skifile_earth(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ski_path_earth)

    # -----------------------------------------------------------------

    @property
    def has_skifile_faceon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ski_path_faceon)

    # -----------------------------------------------------------------

    @property
    def has_skifile_edgeon(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.ski_path_edgeon)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_earth(self):

        """
        This function ...
        :return:
        """

        return FrameInstrument.from_projection(self.projection_earth)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_faceon(self):

        """
        This function ...
        :return:
        """

        return FrameInstrument.from_projection(self.projection_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_instrument_edgeon(self):

        """
        This function ...
        :return:
        """

        return FrameInstrument.from_projection(self.projection_edgeon)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_earth(self):

        """
        This function ...
        :return:
        """

        return FullInstrument.from_projection(self.projection_earth)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_faceon(self):

        """
        This function ...
        :return:
        """

        return FullInstrument.from_projection(self.projection_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_instrument_edgeon(self):

        """
        This function ...
        :return:
        """

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
        #add_new_stellar_component(ski, bulge_component_name, self.old_bulge_component)
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
        #add_new_stellar_component(ski, bulge_component_name, self.old_bulge_component)
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
        #add_new_stellar_component(ski, bulge_component_name, self.old_bulge_component)
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

    @lazyproperty
    def old_bulge_earth_projection_output(self):

        """
        This function ...
        :return:
        """

        return self.simulation_earth.output

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_faceon_projection_output(self):

        """
        This function ...
        :return:
        """

        return self.simulation_faceon.output

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_edgeon_projection_output(self):

        """
        This function ...
        :return:
        """

        return self.simulation_edgeon.output

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_earth_map_path(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_earth_projection_output.single_total_images

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_faceon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_faceon_projection_output.single_total_images

    # -----------------------------------------------------------------

    @lazyproperty
    def old_bulge_edgeon_map_path(self):

        """
        This function ...
        :return:
        """

        return self.old_bulge_edgeon_projection_output.single_total_images

    # -----------------------------------------------------------------

    @property
    def has_earth_wcs(self):

        """
        This function ...
        :return:
        """

        return self.earth_wcs is not None

    # -----------------------------------------------------------------

    @property
    def pixelscale_earth(self):

        """
        This function ...
        :return:
        """

        return self.projection_earth.pixelscale

    # -----------------------------------------------------------------

    @property
    def pixelscale_faceon(self):

        """
        This function ...
        :return:
        """

        return self.projection_faceon.pixelscale

    # -----------------------------------------------------------------

    @property
    def pixelscale_edgeon(self):

        """
        Thisf unction ...
        :return:
        """

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
        return Frame.from_file(self.old_bulge_earth_map_path, wcs=wcs, distance=self.distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_bulge_faceon_map_path, distance=self.distance, pixelscale=self.pixelscale_faceon)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon(self):

        """
        This function ...
        :return:
        """

        return Frame.from_file(self.old_bulge_edgeon_map_path, distance=self.distance, pixelscale=self.pixelscale_edgeon)

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

def create_edgeon_projection_from_earth_projection(earth_projection, scaleheight, radial_factor=1):

    """
    This function ...
    :param earth_projection:
    :param scaleheight:
    :param radial_factor:
    :return:
    """

    # Get pixelscale
    #physical_pixelscale = earth_projection.physical_pixelscale
    physical_pixelscale = earth_projection.physical_pixelscale.average

    # Determine extent in the radial and in the vertical direction
    radial_extent = max(earth_projection.field_x, earth_projection.field_y)

    # Determine the z extent
    z_extent = 2. * scaleheight * default_scale_heights

    # Determine number of pixels
    nx = int(round(radial_extent / physical_pixelscale)) * radial_factor
    nz = int(round(z_extent / physical_pixelscale))

    # Create and return
    return create_edgeon_projection(nx, nz, physical_pixelscale, earth_projection.distance)

# -----------------------------------------------------------------
