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

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ....core.tools import tables
from ....core.tools import filesystem as fs
from ...basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument
from ....core.tools.logging import log
from ....core.prep.wavelengthgrids import WavelengthGridGenerator
from ....core.prep.dustgrids import DustGridGenerator
from ....core.basics.range import RealRange, QuantityRange
from ...component.galaxy import GalaxyModelingComponent
from ....core.basics.unit import parse_unit as u
from ...build.component import get_stellar_component_names, get_dust_component_names, load_stellar_component, load_dust_component
from ..basics.models import DeprojectionModel3D
from ..basics.projection import EdgeOnProjection, FaceOnProjection, GalaxyProjection
from ...magic.basics.coordinatesystem import CoordinateSystem

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

        # The deprojections
        self.deprojections = dict()

        # The projections
        self.projections = dict()

        # The instruments
        self.instruments = dict()

        # The wavelength grid and dust grid generators
        self.wg_generator = None
        self.dg_generator = None

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

        # 3. Create the wavelength grids
        self.create_wavelength_grids()

        # 4. Create the projections
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Create the dust grids
        self.create_dust_grids()

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
        super(RepresentationBuilder, self).setup(self, **kwargs)

        # Create the model definition
        self.definition = self.get_model_definition(self.config.model_name)

        # Create a WavelengthGridGenerator
        self.wg_generator = WavelengthGridGenerator()

        # Create the DustGridGenerator
        self.dg_generator = DustGridGenerator()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.model_name

    # -----------------------------------------------------------------

    def load_deprojections(self):

        """
        This function ...
        :return:
        """

        # Load stellar deprojections
        self.load_stellar_deprojections()

        # load dust deprojections
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

    def create_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grids ...")

        # Fixed wavelengths (always in the grid)
        fixed = [self.i1_filter.pivotwavelength(), self.fuv_filter.pivotwavelength()]

        # Set options
        self.wg_generator.config.show = False
        self.wg_generator.config.write = False

        # Generate the wavelength grids
        self.wg_generator.run(npoints_range=self.config.wg.npoints_range, ngrids=self.config.wg.ngrids,
                              fixed=fixed, add_emission_lines=self.config.wg.add_emission_lines,
                              min_wavelength=self.config.wg.min_wavelength, max_wavelength=self.config.wg.max_wavelength,
                              filters=self.fitting_run.fitting_filters)

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
        self.projections["earth"] = GalaxyProjection.from_wcs(self.reference_wcs, self.galaxy_properties.center,
                                                              self.galaxy_properties.distance,
                                                              self.galaxy_properties.inclination, azimuth, self.disk_pa)

        # Create the face-on projection system
        self.projections["faceon"] = FaceOnProjection.from_wcs(self.reference_wcs, self.galaxy_properties.center,
                                                               self.galaxy_properties.distance)

        # Create the edge-on projection system
        self.projections["edgeon"] = EdgeOnProjection.from_wcs(self.reference_wcs, self.galaxy_properties.center,
                                                               self.galaxy_properties.distance)

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

    def create_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the grids ...")

        # Calculate the major radius of the truncation ellipse in physical coordinates (pc)
        semimajor_angular = self.truncation_ellipse.semimajor  # semimajor axis length of the sky ellipse
        radius_physical = (semimajor_angular * self.galaxy_properties.distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        distance = self.galaxy_properties.distance
        pixelscale_angular = self.reference_wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # BINTREE: (smallest_cell_pixels, min_level, max_mass_fraction)
        # Low-resolution: 10., 6, 1e-5
        # High-resolution: 0.5, 9, 0.5e-6

        # OCTTREE:
        # Low-resolution: 10., 2, 1e-5
        # High-resolution: 0.5, 3, 0.5e-6

        # Because we (currently) can't position the grid exactly as the 2D pixels (rotation etc.),
        # take half of the pixel size to avoid too much interpolation
        min_scale = self.config.dg.scale_range.min * pixelscale
        max_scale = self.config.dg.scale_range.max * pixelscale
        scale_range = QuantityRange(min_scale, max_scale, invert=True)

        # The range of the max mass fraction
        mass_fraction_range = RealRange(self.config.dg.mass_fraction_range.min, self.config.dg.mass_fraction_range.max, invert=True) # must be inverted

        # Set fixed grid properties
        self.dg_generator.grid_type = self.config.dg.grid_type # set grid type
        self.dg_generator.x_radius = radius_physical
        self.dg_generator.y_radius = radius_physical
        self.dg_generator.z_radius = 3. * u("kpc")

        # Set options
        self.dg_generator.show = False
        self.dg_generator.write = False

        # Generate the dust grids
        self.dg_generator.run(scale_range=scale_range, level_range=self.config.dg.level_range, mass_fraction_range=mass_fraction_range, ngrids=10)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the instruments
        self.write_instruments()

        # 6. Write the wavelength grids
        self.write_wavelength_grids()

        # 7. Write the dust grids
        self.write_dust_grids()

        # Write the representations table
        self.write_table()

    # -----------------------------------------------------------------

    def write_instruments(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the SED, frame and simple instruments ...")

        # Write the SED instrument
        self.instruments["SED"].saveto(self.fitting_run.sed_instrument_path)

        # Write the frame instrument
        self.instruments["frame"].saveto(self.fitting_run.frame_instrument_path)

        # Write the simple instrument
        self.instruments["simple"].saveto(self.fitting_run.simple_instrument_path)

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Loop over the grids
        index = 0
        for grid in self.wg_generator.grids:

            # Determine the path to the grid
            path = fs.join(self.fitting_run.wavelength_grids_path, str(index) + ".txt")

            # Save the wavelength grid
            grid.to_skirt_input(path)

            # Increment the index
            index += 1

        # Write the wavelength grids table
        tables.write(self.wg_generator.table, self.fitting_run.wavelength_grids_table_path)

    # -----------------------------------------------------------------

    def write_dust_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust grids ...")

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

        table = self.representations_table
        table.add_entry()

# -----------------------------------------------------------------
