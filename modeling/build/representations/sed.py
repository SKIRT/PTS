#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.representations.sed Contains the SEDRepresentationBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.tools import filesystem as fs
from ...basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ....core.basics.log import log
from ..representation import Representation
from ...component.sed import get_ski_template
from ...basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from .base import RepresentationBuilderBase

# -----------------------------------------------------------------

class SEDRepresentationBuilder(RepresentationBuilderBase):
    
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
        super(SEDRepresentationBuilder, self).__init__(*args, **kwargs)

        # The ski file template
        self.ski = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the ski file
        self.load_ski()

        # 3. Load the instrument
        self.load_projection()

        # 4. Create projections
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Create the dust grids
        self.create_dust_grid()

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
        super(SEDRepresentationBuilder, self).setup(**kwargs)

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

    @property
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return self.projections["earth"]

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the ski template ...")

        # Load ski template
        self.ski = get_ski_template(self.config.path)

    # -----------------------------------------------------------------

    def load_projection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the galaxy projection from the ski file ...")

        # Get the instrument names
        names = self.ski.get_instrument_names()

        if len(names) == 0: raise ValueError("Ski file template does not contain any instruments")
        if len(names) > 1: raise ValueError("Ski file contains more than one instrument")

        # Load the instrument object
        instrument = self.ski.get_instrument_object(names[0])

        # Get the size of the dust grid
        x_range = self.ski.get_dust_grid_x_range()
        y_range = self.ski.get_dust_grid_y_range()
        z_range = self.ski.get_dust_grid_z_range()
        max_field = max(x_range.span, y_range.span, z_range.span)

        # Create projection
        projection = GalaxyProjection.from_instrument(instrument, default_pixels_x=200, default_pixels_y=200, default_field_x=max_field, default_field_y=max_field)

        # Set the projection
        self.projections["earth"] = projection

    # -----------------------------------------------------------------

    def create_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the projections ...")

        # Create face on projection
        self.projections["faceon"] = FaceOnProjection.from_projection(self.earth_projection)

        # Create edge on projection
        self.projections["edgeon"] = EdgeOnProjection.from_projection(self.earth_projection)

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

        # Get the dust grid
        self.dust_grid = self.ski.get_dust_grid_object()

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

        # 1. Write the instruments
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

        # Save the dust grid
        self.dust_grid.saveto(self.representation.dust_grid_path)

    # -----------------------------------------------------------------

    def write_table(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Writing the representations table ...")

        # Add entry and save
        table = self.representations_table
        table.add_entry(self.representation_name, self.model_name)
        table.save()

# -----------------------------------------------------------------
