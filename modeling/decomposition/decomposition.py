#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.decomposition Contains the GalaxyDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .component import DecompositionComponent
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from ...core.simulation.skifile import SkiFile
from ...core.simulation.arguments import SkirtArguments
from ...core.simulation.execute import SkirtExec
from ...magic.basics.vector import Extent, Position
from ...magic.basics.skygeometry import SkyEllipse
from ...magic.basics.skyregion import SkyRegion
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem
from ..basics.models import SersicModel3D, ExponentialDiskModel3D
from ..basics.instruments import SimpleInstrument
from ...magic.misc.kernels import AnianoKernels
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from .s4g import S4GDecomposer
from .fitting import FittingDecomposer

# -----------------------------------------------------------------

# The path to the template ski files directory
template_path = fs.join(introspection.pts_dat_dir("modeling"), "ski")

# -----------------------------------------------------------------

class GalaxyDecomposer(DecompositionComponent):
    
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
        super(GalaxyDecomposer, self).__init__(config)

        # The SKIRT execution context
        self.skirt = SkirtExec()

        # The 2D components
        self.components = None

        # The position angle of the disk of the galaxy (used as the position angle of the galaxy)
        self.disk_pa = None

        # The bulge and disk model
        self.bulge = None
        self.disk = None

        # The bulge and disk image
        self.bulge_image = None
        self.disk_image = None
        self.model_image = None

        # The projection systems
        self.projections = dict()

        # The instruments
        self.instruments = dict()

        # The PSF (of the reference image) for convolution with the simulated images
        self.psf = None

        # Paths to ...
        self.images_bulge2d_path = None
        self.images_bulge_path = None
        self.images_disk_path = None
        self.images_model_path = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Get the decomposition parameters
        self.decompose()

        # 3. Create the 3D models (deproject the 2D models)
        self.create_models()

        # 4. Create the projection systems
        self.create_projections()

        # 5. Create the instruments
        self.create_instruments()

        # 6. Simulate the bulge and disk images
        self.create_images()

        # 7. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyDecomposer, self).setup()

        # TEMP: provide a cfg file for this class
        self.config.bulge_packages = 1e7
        self.config.disk_packages = 1e8

        # Load the PSF
        aniano = AnianoKernels()
        self.psf = aniano.get_psf(self.reference_filter)

        # Create the directory to simulate the bulge (2D method)
        self.images_bulge2d_path = fs.create_directory_in(self.components_images_path, "bulge2D")

        # Create the directory to simulate the bulge (3D)
        self.images_bulge_path = fs.create_directory_in(self.components_images_path, "bulge")

        # Create the directory to simulate the disk (3D)
        self.images_disk_path = fs.create_directory_in(self.components_images_path, "disk")

        # Create the directory to simulate the model (3D bulge + disk)
        self.images_model_path = fs.create_directory_in(self.components_images_path, "model")

    # -----------------------------------------------------------------

    def decompose(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the decomposition parameters ...")

        self.decompose_s4g()

        #self.decompose_fit()

    # -----------------------------------------------------------------

    def decompose_s4g(self):

        """
        This function ...
        :return:
        """

        # Create ...
        decomposer = S4GDecomposer()
        #parameters = FittingDecompositionParameters()

        # Run ...
        decomposer.run()

        # Add the models
        self.components = decomposer.components

        # Set the disk position angle
        self.disk_pa = self.components["disk"].position_angle

    # -----------------------------------------------------------------

    def decompose_fit(self):

        """
        This function ...
        :return:
        """

        decomposer = FittingDecomposer()

        decomposer.run()

        # Add the components
        #self.components = decomposer.components

    # -----------------------------------------------------------------

    def create_models(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the 3D bulge and disk models ...")

        # Create the bulge model
        self.create_bulge_model()

        # Create the disk model
        self.create_disk_model()

    # -----------------------------------------------------------------

    def create_bulge_model(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the bulge model ...")

        # Create a Sersic model for the bulge
        self.bulge = SersicModel3D.from_2d(self.components["bulge"], self.galaxy_properties.inclination, self.disk_pa)

    # -----------------------------------------------------------------

    def create_disk_model(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the disk model ...")

        # Create an exponential disk model for the disk
        self.disk = ExponentialDiskModel3D.from_2d(self.components["disk"], self.galaxy_properties.inclination, self.disk_pa)

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

        # Loop over the projection systems
        for name in self.projections:

            # Create the instrument from the projection system
            self.instruments[name] = SimpleInstrument.from_projection(self.projections[name])

    # -----------------------------------------------------------------

    def create_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the images of the bulge, disk and bulge+disk model ...")

        # Simulate the stellar bulge without deprojection
        self.simulate_bulge2d()

        # Simulate the stellar bulge
        self.simulate_bulge()

        # Simulate the stellar disk
        self.simulate_disk()

        # Simulate the bulge + disk model
        self.simulate_model()

    # -----------------------------------------------------------------

    def simulate_bulge2d(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = fs.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.bulge_packages)

        # Change the ski file parameters
        # component_id, index, radius, y_flattening=1, z_flattening=1
        ski.set_stellar_component_sersic_geometry(0, self.components["bulge"].index, self.components["bulge"].effective_radius, y_flattening=self.components["bulge"].axial_ratio)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Create the instrument
        distance = self.galaxy_properties.distance
        inclination = 0.0
        azimuth = Angle(90., "deg")
        #position_angle = self.parameters.bulge.PA + Angle(90., "deg") # + 90° because we can only do y_flattening and not x_flattening
        position_angle = self.components["bulge"].position_angle
        pixels_x = self.reference_wcs.xsize
        pixels_y = self.reference_wcs.ysize
        pixel_center = self.galaxy_properties.center.to_pixel(self.reference_wcs)
        center = Position(0.5*pixels_x - pixel_center.x - 0.5, 0.5*pixels_y - pixel_center.y - 0.5)
        center_x = center.x * Unit("pix")
        center_y = center.y * Unit("pix")
        center_x = (center_x * self.reference_wcs.pixelscale.x.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * self.reference_wcs.pixelscale.y.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        field_x_angular = self.reference_wcs.pixelscale.x.to("deg/pix") * pixels_x * Unit("pix")
        field_y_angular = self.reference_wcs.pixelscale.y.to("deg/pix") * pixels_y * Unit("pix")
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())
        fake = SimpleInstrument(distance, inclination, azimuth, position_angle, field_x_physical, field_y_physical, pixels_x, pixels_y, center_x, center_y)

        # Add the instrument
        ski.add_instrument("earth", fake)

        # Determine the path to the ski file
        ski_path = fs.join(self.images_bulge2d_path, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.create_directory_in(self.images_bulge2d_path, "out")

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        bulge_image_path = fs.join(out_path, "bulge_earth_total.fits")

        # Check if the output contains the "bulge_earth_total.fits" file
        if not fs.is_file(bulge_image_path): raise RuntimeError("Something went wrong with the simple bulge simulation: output FITS file missing")

    # -----------------------------------------------------------------

    def simulate_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = fs.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.bulge_packages)

        # Set the bulge geometry
        ski.set_stellar_component_geometry(0, self.bulge)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.images_bulge_path, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory and create it
        out_path = fs.create_directory_in(self.images_bulge_path, "out")

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        bulge_image_path = fs.join(out_path, "bulge_earth_total.fits")

        # Check if the output contains the "bulge_earth_total.fits" file
        if not fs.is_file(bulge_image_path):
            raise RuntimeError("Something went wrong with the bulge simulation: output FITS file missing")

        # Open the bulge image
        self.bulge_image = Frame.from_file(bulge_image_path)

        # Set the coordinate system of the bulge image
        self.bulge_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the bulge image to the bulge flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.components["bulge"].fluxdensity
        self.bulge_image *= fluxdensity.to("Jy").value / np.sum(self.bulge_image)
        self.bulge_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the bulge image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.bulge_image = self.bulge_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def simulate_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the disk image ...")

        # Load the disk ski file template
        disk_template_path = fs.join(template_path, "disk.ski")
        ski = SkiFile(disk_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.disk_packages)

        # Change the ski file parameters
        ski.set_stellar_component_geometry(0, self.disk)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.images_disk_path, "disk.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory and create it
        out_path = fs.create_directory_in(self.images_disk_path, "out")

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the disk simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        disk_image_path = fs.join(out_path, "disk_earth_total.fits")

        # Check if the output contains the "disk_earth_total.fits" file
        if not fs.is_file(disk_image_path):
            raise RuntimeError("Something went wrong with the disk simulation: output FITS file missing")

        # Open the disk image
        self.disk_image = Frame.from_file(disk_image_path)

        # Set the coordinate system of the disk image
        self.disk_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the disk image to the disk flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.components["disk"].fluxdensity
        self.disk_image *= fluxdensity.to("Jy").value / np.sum(self.disk_image)
        self.disk_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the disk image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.disk_image = self.disk_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def simulate_model(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge+disk model image ...")

        # Load the disk ski file template
        disk_template_path = fs.join(template_path, "model.ski")
        ski = SkiFile(disk_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.disk_packages)

        # Change the ski file parameters
        ski.set_stellar_component_geometry(0, self.disk)
        ski.set_stellar_component_geometry(1, self.bulge)

        # Set the luminosities of the two components
        ski.set_stellar_component_luminosities(0, [self.components["disk"].rel_contribution])
        ski.set_stellar_component_luminosities(1, [self.components["bulge"].rel_contribution])

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add the instruments
        for name in self.instruments: ski.add_instrument(name, self.instruments[name])

        # Determine the path to the ski file
        ski_path = fs.join(self.images_model_path, "model.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory and create it
        out_path = fs.create_directory_in(self.images_model_path, "out")

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the disk+bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug() else True)

        # Determine the path to the output FITS file
        model_image_path = fs.join(out_path, "model_earth_total.fits")

        # Check if the output contains the "model_earth_total.fits" file
        if not fs.is_file(model_image_path):
            raise RuntimeError("Something went wrong with the disk+bulge simulation: output FITS file missing")

        # Open the model image
        self.model_image = Frame.from_file(model_image_path)

        # Set the coordinate system of the model image
        self.model_image.wcs = self.reference_wcs

        # Debugging
        log.debug("Rescaling the model image to the bulge+disk flux density at 3.6 micron ...")

        # Rescale to the 3.6um flux density
        fluxdensity = self.components["bulge"].fluxdensity + self.components["disk"].fluxdensity # sum of bulge and disk component flux density
        self.model_image *= fluxdensity.to("Jy").value / np.sum(self.model_image)
        self.model_image.unit = "Jy"

        # Debugging
        log.debug("Convolving the model image to the PACS 160 resolution ...")

        # Convolve the frame to the PACS 160 resolution
        self.model_image = self.model_image.convolved(self.psf)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the model descriptions
        self.write_models()

        # Write out the final bulge and disk images
        self.write_images()

        # Write the projection systems
        self.write_projections()

        # Write out the disk ellipse
        self.write_disk_ellipse()

    # -----------------------------------------------------------------

    def write_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the models ...")

        # Write the disk model
        self.disk.save(self.disk_model_path)

        # Write the bulge model
        self.bulge.save(self.bulge_model_path)

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Determine the path to the bulge image and save it
        final_bulge_path = fs.join(self.components_images_path, "bulge.fits")
        self.bulge_image.save(final_bulge_path)

        # Determine the path to the disk image and save it
        final_disk_path = fs.join(self.components_images_path, "disk.fits")
        self.disk_image.save(final_disk_path)

        # Determine the path to the model image and save it
        final_model_path = fs.join(self.components_images_path, "model.fits")
        self.model_image.save(final_model_path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections["earth"].save(self.earth_projection_path)

        # Write the edgeon projection system
        self.projections["edgeon"].save(self.edgeon_projection_path)

        # Write the faceon projection system
        self.projections["faceon"].save(self.faceon_projection_path)

    # -----------------------------------------------------------------

    def write_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing regions file with disk ellipse ...")

        minor = (1.0 - self.parameters.ellipticity) * self.parameters.major_arcsec

        # Ellipse radius
        radius = Extent(self.parameters.major_arcsec, minor)

        # Create sky ellipse
        sky_ellipse = SkyEllipse(self.parameters.center, radius, self.parameters.position_angle)

        # Create region
        region = SkyRegion()
        region.append(sky_ellipse)
        region_path = fs.join(self.components_path, "disk.reg")
        region.save(region_path)

# -----------------------------------------------------------------
