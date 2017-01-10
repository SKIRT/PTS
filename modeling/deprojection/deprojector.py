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

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from .component import DeprojectionComponent
from ..basics.models import DeprojectionModel3D
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.execute import SkirtExec
from ...core.tools import introspection
from ...core.simulation.skifile import LabeledSkiFile
from ..basics.instruments import FrameInstrument
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

template_ski_path = fs.join(introspection.pts_dat_dir("modeling"), "ski", "labeled_template.ski")

# -----------------------------------------------------------------

class Deprojector(DeprojectionComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(Deprojector, self).__init__(config)

        # The SKIRT execution environment
        self.skirt = SkirtExec()

        # The deprojection models for the different components
        self.deprojections = dict()

        # The ski files
        self.ski_files = dict()

        # The ski file paths
        self.ski_paths = dict()

        self.output_paths = dict()

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # Create the deprojection models
        self.create_deprojection_models()

        # Create the ski files
        self.create_ski_files()

        # Write
        self.write()

        # Simulate
        self.simulate()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(Deprojector, self).setup()

        self.ski_paths["old stars"] = fs.join(self.deprojection_path, "old_stars.ski")
        self.ski_paths["young stars"] = fs.join(self.deprojection_path, "young_stars.ski")
        self.ski_paths["ionizing stars"] = fs.join(self.deprojection_path, "ionizing_stars.ski")
        self.ski_paths["dust"] = fs.join(self.deprojection_path, "dust.ski")

        self.output_paths["old stars"] = fs.create_directory_in(self.deprojection_path, "old stars")
        self.output_paths["young stars"] = fs.create_directory_in(self.deprojection_path, "young stars")
        self.output_paths["ionizing stars"] = fs.create_directory_in(self.deprojection_path, "ionizing stars")
        self.output_paths["dust"] = fs.create_directory_in(self.deprojection_path, "dust")

    # -----------------------------------------------------------------

    def create_deprojection_models(self):

        """
        This function ...
        :return:
        """

        # Generate default deprojection model
        deprojection = self.create_default_deprojection_model()

        # Create the ...
        self.create_old_deprojection_model(deprojection)
        self.create_young_deprojection_model(deprojection)
        self.create_ionizing_deprojection_model(deprojection)
        self.create_dust_deprojection_model(deprojection)

    # -----------------------------------------------------------------

    def create_default_deprojection_model(self):

        """
        This function ...
        :return:
        """

        # Not specified
        filename = None
        hz = None

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        pa = self.earth_projection.position_angle

        # Get the center pixel
        pixel_center = self.galaxy_properties.center.to_pixel(self.reference_wcs)
        xc = pixel_center.x
        yc = pixel_center.y

        # Get the pixelscale in physical units
        pixelscale_angular = self.reference_wcs.average_pixelscale * Unit("pix")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the number of x and y pixels
        x_size = self.reference_wcs.xsize
        y_size = self.reference_wcs.ysize

        # Create the deprojection model
        deprojection = DeprojectionModel3D(filename, pixelscale, pa, inclination, x_size, y_size, xc, yc, hz)

        # Return the default deprojection model
        return deprojection

    # -----------------------------------------------------------------

    def create_old_deprojection_model(self, default):

        """
        This function ...
        :return:
        """

        # Does not really matter
        #scale_height = self.disk2d_model.scalelength / 8.26
        scale_height = 100. * Unit("pc")

        # Set the parameters of the evolved stellar component
        deprojection = default.copy()
        deprojection.filename = self.old_stellar_map_filename
        deprojection.scale_height = scale_height
        self.deprojections["old stars"] = deprojection

    # -----------------------------------------------------------------

    def create_young_deprojection_model(self, default):

        """
        This function ...
        :return:
        """

        # Does not really matter
        scale_height = 100. * Unit("pc")

        # Set the parameters of the young stellar component
        deprojection = default.copy()
        deprojection.filename = self.young_stellar_map_filename
        deprojection.scale_height = scale_height
        self.deprojections["young stars"] = deprojection

    # -----------------------------------------------------------------

    def create_ionizing_deprojection_model(self, default):

        """
        This function ...
        :return:
        """

        # Does not really matter
        scale_height = 100. * Unit("pc")

        # Set the parameters of the ionizing stellar component
        deprojection = default.copy()
        deprojection.filename = self.ionizing_stellar_map_filename
        deprojection.scale_height = scale_height
        self.deprojections["ionizing stars"] = deprojection

    # -----------------------------------------------------------------

    def create_dust_deprojection_model(self, default):

        """
        This function ...
        :param default:
        :return:
        """

        # Does not really matter
        scale_height = 100. * Unit("pc")

        # Set the parameters of the ionizing stellar component
        deprojection = default.copy()
        deprojection.filename = self.dust_map_filename
        deprojection.scale_height = scale_height
        self.deprojections["dust"] = deprojection

    # -----------------------------------------------------------------

    def create_ski_files(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski files ...")

        # Create instrument
        #self.instrument = FrameInstrument.from_projection(self.earth_projection)
        self.instrument = FrameInstrument.from_projection(self.faceon_projection)

        # Load ski template
        self.ski_template = LabeledSkiFile(template_ski_path)

        # Convert to oligochromatic simulation
        self.ski_template.to_oligochromatic(1. * Unit("micron"))

        # Remove the dust system
        self.ski_template.remove_dust_system()

        # Set number of packages per wavelength
        self.ski_template.setpackages(1e6)

        # Add one instrument
        self.ski_template.remove_all_instruments()
        self.ski_template.add_instrument("faceon", self.instrument)



        # Old

        old_ski = self.ski_template.copy()
        old_ski.remove_stellar_components_except("Evolved stellar disk")
        #old_ski.remove_dust_system()
        old_ski.set_stellar_component_geometry("Evolved stellar disk", self.deprojections["old stars"])
        self.ski_files["old stars"] = old_ski

        young_ski = self.ski_template.copy()
        young_ski.remove_stellar_components_except("Young stars")
        young_ski.set_stellar_component_geometry("Young stars", self.deprojections["young stars"])
        self.ski_files["young stars"] = young_ski


        ionizing_ski = self.ski_template.copy()
        ionizing_ski.remove_stellar_components_except("Ionizing stars")
        ionizing_ski.set_stellar_component_geometry("Ionizing stars", self.deprojections["ionizing stars"])
        self.ski_files["ionizing stars"] = ionizing_ski

        dust_ski = self.ski_template.copy()
        dust_ski.remove_stellar_components_except("Ionizing stars")
        dust_ski.set_stellar_component_geometry("Ionizing stars", self.deprojections["dust"])
        self.ski_files["dust"] = dust_ski

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        self.write_deprojections()

        self.write_ski_files()

    # -----------------------------------------------------------------

    def write_deprojections(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_ski_files(self):

        """
        This function ...
        :return:
        """

        # Write the ...
        for name in self.ski_files:

            self.ski_files[name].saveto(self.ski_paths[name])

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting by launching SKIRT simulations ...")

        for name in self.ski_files:

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

# -----------------------------------------------------------------
