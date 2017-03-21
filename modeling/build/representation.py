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

# Import astronomical modules
from astropy.units import dimensionless_angles
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools import tables
from ...core.tools import filesystem as fs
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ...core.tools.logging import log
from ..build.component import get_stellar_component_names, get_dust_component_names, load_stellar_component, load_dust_component
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.basics.configuration import prompt_string
from ...core.basics.quantity import represent_quantity
from ...core.simulation.grids import load_grid

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

class RepresentationBuilder(BuildComponent):
    
    """
    This class...
    """

    def __init__(self, config=None, interactive=False):

        """
        The constructor ...
        :param config:
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(RepresentationBuilder, self).__init__(config, interactive)

        # The model definition
        self.definition = None

        # The representation
        self.representation = None

        # The deprojections
        self.deprojections = dict()

        # The reference deprojection
        self.reference_deprojection = None

        # The projections
        self.projections = dict()

        # The instruments
        self.instruments = dict()

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

        # 3. Prompt for the resolution of this representation
        self.prompt_resolution()

        # 4. Create the projections
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Create the dust grids
        self.create_dust_grid()

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
        super(RepresentationBuilder, self).setup(**kwargs)

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
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.representation.model_name

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

            # Load the component
            component = load_stellar_component(self.config.path, self.model_name, name)

            # Set deprojection
            if "deprojection" in component:

                # Get title
                title = component.parameters.title

                # Add to the dictionary of deprojections
                self.deprojections[(name, title)] = component.deprojection

            # Check if this is a new component, add geometry, SED and normalization all at once
            if "geometry" in component.parameters:

                # Get title
                title = component.parameters.title

                # Check whether this is a read FITS geometry
                geometry_type = component.parameters.geometry
                if geometry_type != "ReadFitsGeometry": continue

                # Get properties for each of the three classes
                geometry_properties = component.properties["geometry"]

                # Get the path of the input map
                filepath = geometry_properties["filename"]

                # Get the scale height
                scale_height = geometry_properties["axialScale"]

                # Create the deprojection
                wcs = CoordinateSystem.from_file(filepath)
                deprojection = self.create_deprojection_for_wcs(wcs, filepath, scale_height)

                # Add to the dictionary
                self.deprojections[(name, title)] = deprojection

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

            # Load the component
            component = load_dust_component(self.config.path, self.model_name, name)

            # Set deprojection
            if "deprojection" in component:

                # Get title
                title = component.parameters.title

                # Add to the dictionary of deprojections
                self.deprojections[(name, title)] = component.deprojection

            # Check if this is a new dust component, add geometry, mix and normalization all at once
            if "geometry" in component.parameters:

                # Get title
                title = component.parameters.title

                # Check whether this is a read FITS geometry
                geometry_type = component.parameters.geometry
                if geometry_type != "ReadFitsGeometry": continue

                # Get properties for each of the three classes
                geometry_properties = component.properties["geometry"]

                # Get the path of the input map
                filepath = geometry_properties["filename"]

                # Get the scale height
                scale_height = geometry_properties["axialScale"]

                # Create the deprojection
                wcs = CoordinateSystem.from_file(filepath)
                deprojection = self.create_deprojection_for_wcs(wcs, filepath, scale_height)

                # Add to the dictionary
                self.deprojections[(name, title)] = deprojection

    # -----------------------------------------------------------------

    def prompt_resolution(self):

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
        answer = prompt_string("reference map", "input map to use as the reference for the resolution of the model representation", choices=options, default=lowest_pixelscale_name)

        # Set the reference deprojection
        answer_title = None
        for name, title in self.deprojections:
            if name == answer:
                answer_title = title
                break
        self.reference_deprojection = self.deprojections[(answer, answer_title)]

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projection systems ...")

        # Create the 'earth' projection system
        azimuth = 0.0
        self.projections["earth"] = GalaxyProjection.from_deprojection(self.reference_deprojection, self.galaxy_distance, azimuth)

        # Create the face-on projection system
        self.projections["faceon"] = FaceOnProjection.from_deprojection(self.reference_deprojection, self.galaxy_distance)

        # Create the edge-on projection system
        self.projections["edgeon"] = EdgeOnProjection.from_deprojection(self.reference_deprojection, self.galaxy_distance)

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

        # TODO: implement

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the projections
        self.write_projections()

        # 2. Write the instruments
        self.write_instruments()

        # 3. Write the dust grids
        self.write_dust_grid()

        # 4. Write the representations table
        self.write_table()

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

    def write_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grid ...")

        # Loop over the grids
        index = 0
        for grid in self.dg_generator.grids:

            # Determine the path to the grid
            path = fs.join(self.fitting_run.dust_grids_path, str(index) + ".dg")

            # Save the dust grid
            grid.saveto(path)

            # Increment the index
            index += 1

        # Write the dust grids table
        tables.write(self.dg_generator.table, self.fitting_run.dust_grids_table_path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing the representations table ...")

        table = self.representations_table
        table.add_entry(self.representation_name, self.model_name)
        table.save()

# -----------------------------------------------------------------
