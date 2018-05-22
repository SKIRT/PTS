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
from ...core.basics.log import log
from ...core.simulation.skifile import SkiFile
from ...magic.basics.vector import Position
from ...magic.basics.stretch import SkyStretch
from ...magic.region.ellipse import SkyEllipseRegion
from ...magic.region.list import SkyRegionList
from ..basics.models import SersicModel3D, ExponentialDiskModel3D
from ..basics.instruments import SimpleInstrument
from ...magic.convolution.aniano import AnianoKernels
from ..basics.projection import GalaxyProjection, FaceOnProjection, EdgeOnProjection
from .s4g import S4GDecomposer
#from .fitting import FittingDecomposer
#from .imfit import ImfitDecomposer
from ...core.launch.launcher import SingleImageSKIRTLauncher
from ...magic.core.frame import Frame
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

template_path = introspection.pts_modeling_ski_templates_path

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

# Slightly varying values
kregel_ratio = 8.21
degeyter_ratio = 8.26
mosenkov_ratio = 9.06

# Values and descriptions
scalelength_scaleheight_ratios = dict()
scalelength_scaleheight_ratios[degeyter_ratio] = "De Geyter et al. (2014)"
scalelength_scaleheight_ratios[kregel_ratio] = "Kregel et al. (2002)"
scalelength_scaleheight_ratios[mosenkov_ratio] = "Mosenkov et al. (2015)"

# -----------------------------------------------------------------

q_kregel = 1./kregel_ratio
q_degeyter = 1./degeyter_ratio
q_mosenkov = 1./mosenkov_ratio

# -----------------------------------------------------------------

def get_logq0_mosenkov(hubble_stage):

    """
    Thisn function ...
    :param hubble_stage:
    :return:
    """

    # Below or equal to 7
    if hubble_stage <= 7: return -0.026 * float(hubble_stage) - 0.774

    # Above 7
    else: return 0.107 * float(hubble_stage) - 1.705

# -----------------------------------------------------------------

def axial_ratio_to_inclination_mosenkov(ratio, hubble_stage):

    """
    This function ...
    :param ratio:
    :param hubble_stage:
    :return:
    """

    # Check that axial ratio is smaller than one!
    if ratio >= 1.: raise ValueError("Axial ratio must be ratio of minor axis to major axis length")

    # From Mosenkov et al., 2017 (Appendix)

    # Calculate logq and logq0
    logq = np.log10(ratio)
    logq0 = get_logq0_mosenkov(hubble_stage)

    # Calculate numerator and denominator of the formula
    numerator = 1. - 10**(2. * logq)
    denominator = 1. - 10**(2. * logq0)

    # Calculate the inclination angle
    sin2i = numerator / denominator
    inclination_radians = np.arcsin(np.sqrt(sin2i))
    inclination = Angle(inclination_radians, unit="rad")
    return inclination.to("deg")

# -----------------------------------------------------------------

def ellipticity_to_axial_ratio(ellipticity):

    """
    This function ...
    :param ellipticity:
    :return:
    """

    return 1. - ellipticity

# -----------------------------------------------------------------

def axial_ratio_to_ellipticity(axial_ratio):

    """
    This function ...
    :param axial_ratio:
    :return:
    """

    return 1. - axial_ratio

# -----------------------------------------------------------------

def ellipticity_to_inclination_mosenkov(ellipticity, hubble_stage):

    """
    This function ...
    :param ellipticity:
    :param hubble_stage:
    :return:
    """

    axial_ratio = ellipticity_to_axial_ratio(ellipticity)
    return axial_ratio_to_inclination_mosenkov(axial_ratio, hubble_stage)

# -----------------------------------------------------------------

def axial_ratio_to_inclination_with_intrinsic(ratio, intrinsic_ratio):

    """
    This function ...
    :param ratio:
    :param intrinsic_ratio:
    :return:
    """

    # Check ratios
    if ratio >= 1.: raise ValueError("Axial ratio must be ratio of minor axis to major axis length")
    if intrinsic_ratio >= 1: raise ValueError("Intrinsic axial ratio must be ratio of scale height (smaller) to scale length (larger)")

    # Calculate logq and logq0
    logq = np.log10(ratio)
    logq0 = np.log10(intrinsic_ratio)

    # Calculate numerator and denominator of the formula
    numerator = 1. - 10 ** (2. * logq)
    denominator = 1. - 10 ** (2. * logq0)

    # Calculate the inclination angle
    sin2i = numerator / denominator
    inclination_radians = np.arcsin(np.sqrt(sin2i))
    inclination = Angle(inclination_radians, unit="rad")
    return inclination.to("deg")

# -----------------------------------------------------------------

class GalaxyDecomposer(DecompositionComponent):
    
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
        super(GalaxyDecomposer, self).__init__(*args, **kwargs)

        # The WCS
        self.wcs = None

        # The SKIRT launching environment
        self.launcher = SingleImageSKIRTLauncher()

        # The 2D components
        self.components = None

        # The bulge and disk model
        self.bulge = None
        self.disk = None

        # The bulge and disk image
        self.bulge2d_image = None
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

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Get the decomposition parameters
        if not self.has_components: self.decompose()

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

    @property
    def has_components(self):

        """
        This function ...
        :return:
        """

        return self.disk is not None and self.bulge is not None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyDecomposer, self).setup(**kwargs)

        # Get provided models
        self.disk = kwargs.pop("disk", None)
        self.bulge = kwargs.pop("bulge", None)

        # Check the method and filter
        if self.config.method == "s4g" and self.config.filter != "IRAC I1": raise ValueError("When using the S4G method, the filter can either be 'IRAC I1' or 'IRAC I2'")
        #in [parse_filter("IRAC I1"), parse_filter("IRAC I2")]

        # Set the WCS for the filter
        self.wcs = self.wcs_for_filter(self.config.filter)

        # TEMP: provide a cfg file for this class
        self.config.bulge_packages = 1e7
        self.config.disk_packages = 1e8

        # Load the PSF kernel and prepare
        aniano = AnianoKernels()
        self.psf = aniano.get_psf(self.config.filter)
        self.psf.prepare_for(self.wcs)

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

        # Use the S4G database
        if self.config.method == "s4g": self.decompose_s4g()

        # Fit using python
        elif self.config.method == "fit": self.decompose_fit()

        # Fit using Imfit
        elif self.config.method == "imfit": self.decompose_imfit()

        # Invalid
        else: raise ValueError("Invalid option for 'method': " + self.config.method)

    # -----------------------------------------------------------------

    def decompose_s4g(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the decomposition parameters from the S4G database ...")

        # Create ...
        decomposer = S4GDecomposer()
        #parameters = FittingDecompositionParameters()

        # Run ...
        decomposer.run()

        # Add the models
        self.components = decomposer.components

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_pa(self):

        """
        This function ...
        :return:
        """

        return self.components["disk"].position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_q(self):

        """
        This function ...
        :return:
        """

        return self.components["disk"].axial_ratio

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_scalelength(self):

        """
        This function ...
        :return:
        """

        return self.components["disk"].scalelength

    # -----------------------------------------------------------------

    @property
    def intrinsic_ratio(self):

        """
        Thisf unction ...
        :return:
        """

        return 1. / self.config.scalelength_to_scaleheight

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_inclination(self):

        """
        This function ...
        :return:
        """

        #return axial_ratio_to_inclination_mosenkov(self.disk_q, self.hubble_stage)
        return axial_ratio_to_inclination_with_intrinsic(self.disk_q, intrinsic_ratio=self.intrinsic_ratio)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_ellipticity(self):

        """
        This function ...
        :return:
        """

        return axial_ratio_to_ellipticity(self.disk_q)

    # -----------------------------------------------------------------

    def decompose_fit(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

        # Create the decomposer
        #decomposer = FittingDecomposer()

        # Run the decomposition
        #decomposer.run()

        # Add the components
        #self.components = decomposer.components

    # -----------------------------------------------------------------

    def decompose_imfit(self):

        """
        This function ...
        :return:
        """

        raise NotImplementedError("Not implemented yet")

        # Create the decomposer
        #decomposer = ImfitDecomposer()

        # Run the decomposer
        #decomposer.run()

    # -----------------------------------------------------------------

    def create_models(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the 3D bulge and disk models ...")

        # Create the bulge model
        if self.bulge is None: self.create_bulge_model()

        # Create the disk model
        if self.disk is None: self.create_disk_model()

    # -----------------------------------------------------------------

    def create_bulge_model(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the bulge model ...")

        # Create a Sersic model for the bulge
        self.bulge = SersicModel3D.from_2d(self.components["bulge"], self.disk_inclination, self.disk_pa, azimuth_or_tilt=self.config.bulge_deprojection_method)

    # -----------------------------------------------------------------

    def create_disk_model(self):

        """
        :return:
        """

        # Inform the user
        log.info("Creating the disk model ...")

        # Create an exponential disk model for the disk
        self.disk = ExponentialDiskModel3D.from_2d(self.components["disk"], self.disk_inclination, self.disk_pa)

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
        self.projections["earth"] = GalaxyProjection.from_wcs(self.wcs, self.galaxy_properties.center, self.galaxy_properties.distance, self.disk_inclination, azimuth, self.disk_pa)

        # Create the face-on projection system
        self.projections["faceon"] = FaceOnProjection.from_wcs(self.wcs, self.galaxy_properties.center, self.galaxy_properties.distance)

        # Create the edge-on projection system
        self.projections["edgeon"] = EdgeOnProjection.from_wcs(self.wcs, self.galaxy_properties.center, self.galaxy_properties.distance)

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

    @property
    def simulated_bulge2d_image_path(self):

        """
        This function ...
        :return:
        """

        prefix = "bulge"

        # Determine the name of the SKIRT output FITS file
        fits_name = prefix + "_" + instrument_name + "_total.fits"

        # Determine the path to the output FITS file
        out_path = fs.join(self.images_bulge2d_path, "out")
        fits_path = fs.join(out_path, fits_name)

        # Return the path
        return fits_path

    # -----------------------------------------------------------------

    @property
    def simulated_bulge_image_path(self):

        """
        This function ...
        :return:
        """

        prefix = "bulge"

        # Determine the name of the SKIRT output FITS file
        fits_name = prefix + "_" + instrument_name + "_total.fits"

        # Determine the path to the output FITS file
        out_path = fs.join(self.images_bulge_path, "out")
        fits_path = fs.join(out_path, fits_name)

        # Return the path
        return fits_path

    # -----------------------------------------------------------------

    @property
    def simulated_disk_image_path(self):

        """
        This function ...
        :return:
        """

        prefix = "disk"

        # Determine the name of the SKIRT output FITS file
        fits_name = prefix + "_" + instrument_name + "_total.fits"

        # Determine the path to the output FITS file
        out_path = fs.join(self.images_disk_path, "out")
        fits_path = fs.join(out_path, fits_name)

        # Return the path
        return fits_path

    # -----------------------------------------------------------------

    @property
    def simulated_model_image_path(self):

        """
        Thisf unction ...
        :return:
        """

        prefix = "model"

        # Determine the name of the SKIRT output FITS file
        fits_name = prefix + "_" + instrument_name + "_total.fits"

        # Determine the path to the output FITS file
        out_path = fs.join(self.images_model_path, "out")
        fits_path = fs.join(out_path, fits_name)

        # Return the path
        return fits_path

    # -----------------------------------------------------------------

    @property
    def has_simulated_bulge2d_image(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_bulge2d_image_path)

    # -----------------------------------------------------------------

    @property
    def has_simulated_bulge_image(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_bulge_image_path)

    # -----------------------------------------------------------------

    @property
    def has_simulated_disk_image(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_disk_image_path)

    # -----------------------------------------------------------------

    @property
    def has_simulated_model_image(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.simulated_model_image_path)

    # -----------------------------------------------------------------

    def create_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the images of the bulge, disk and bulge+disk model ...")

        # Simulate the stellar bulge without deprojection
        if self.has_simulated_bulge2d_image: self.load_bulge2d()
        else: self.simulate_bulge2d()

        # Simulate the stellar bulge
        if self.has_simulated_bulge_image: self.load_bulge()
        else: self.simulate_bulge()

        # Simulate the stellar disk
        if self.has_simulated_disk_image: self.load_disk()
        else: self.simulate_disk()

        # Simulate the bulge + disk model
        if self.has_simulated_model_image: self.load_model()
        else: self.simulate_model()

    # -----------------------------------------------------------------

    def load_bulge2d(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated bulge-2D image ...")

        # Open the simulated frame
        simulated_frame = Frame.from_file(self.simulated_bulge2d_image_path)
        fluxdensity = self.components["bulge"].fluxdensity

        # Set the coordinate system of the disk image
        simulated_frame.wcs = self.wcs

        # Debugging
        log.debug("Rescaling the bulge image to a flux density of " + str(fluxdensity) + " ...")

        # Rescale to the 3.6um flux density
        simulated_frame.normalize(to=fluxdensity)

        # Debugging
        log.debug("Convolving the bulge image to the resolution of the " + str(self.config.filter) + " filter ...")

        # Convolve the frame
        simulated_frame.convolve(self.psf)

        # Set
        self.bulge2d_image = simulated_frame

    # -----------------------------------------------------------------

    def load_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated bulge image ...")

        # Open the simulated frame
        simulated_frame = Frame.from_file(self.simulated_bulge_image_path)
        fluxdensity = self.components["bulge"].fluxdensity

        # Set the coordinate system of the disk image
        simulated_frame.wcs = self.wcs

        # Debugging
        log.debug("Rescaling the bulge image to a flux density of " + str(fluxdensity) + " ...")

        # Rescale to the 3.6um flux density
        simulated_frame.normalize(to=fluxdensity)

        # Debugging
        log.debug("Convolving the bulge image to the resolution of the " + str(self.config.filter) +  " filter ...")

        # Convolve the frame
        simulated_frame.convolve(self.psf)

        # Set
        self.bulge_image = simulated_frame

    # -----------------------------------------------------------------

    def load_disk(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated disk image ...")

        # Open the simulated frame
        simulated_frame = Frame.from_file(self.simulated_disk_image_path)
        fluxdensity = self.components["disk"].fluxdensity

        # Set the coordinate system of the disk image
        simulated_frame.wcs = self.wcs

        # Debugging
        log.debug("Rescaling the disk image to a flux density of " + str(fluxdensity) + " ...")

        # Rescale to the 3.6um flux density
        simulated_frame.normalize(to=fluxdensity)

        # Debugging
        log.debug("Convolving the disk image to the resolution of the " + str(self.config.filter) + " filter ...")

        # Convolve the frame
        simulated_frame.convolve(self.psf)

        # Set
        self.disk_image = simulated_frame

    # -----------------------------------------------------------------

    def load_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the simulated model image ...")

        # Open the simulated frame
        simulated_frame = Frame.from_file(self.simulated_model_image_path)
        fluxdensity = self.components["bulge"].fluxdensity + self.components["disk"].fluxdensity

        # Set the coordinate system of the disk image
        simulated_frame.wcs = self.wcs

        # Debugging
        log.debug("Rescaling the model image to a flux density of " + str(fluxdensity) + " ...")

        # Rescale to the 3.6um flux density
        simulated_frame.normalize(to=fluxdensity)

        # Debugging
        log.debug("Convolving the model image to the resolution of the " + str(self.config.filter) + " filter ...")

        # Convolve the frame
        simulated_frame.convolve(self.psf)

        # Set
        self.model_image = simulated_frame

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
        inclination = Angle(0.0, "deg")
        #inclination = 0.0 # doesn't matter, also works! (thanks to intelligent parse_quantity)
        azimuth = Angle(90., "deg")
        #position_angle = self.parameters.bulge.PA + Angle(90., "deg") # + 90° because we can only do y_flattening and not x_flattening
        position_angle = self.components["bulge"].position_angle
        pixels_x = self.wcs.xsize
        pixels_y = self.wcs.ysize
        pixel_center = self.galaxy_properties.center.to_pixel(self.wcs)
        center = Position(0.5*pixels_x - pixel_center.x - 0.5, 0.5*pixels_y - pixel_center.y - 0.5)
        center_x = center.x
        center_y = center.y
        center_x = (center_x * self.wcs.pixelscale.x.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * self.wcs.pixelscale.y.to("deg") * distance).to("pc", equivalencies=dimensionless_angles())
        field_x_angular = self.wcs.pixelscale.x.to("deg") * pixels_x
        field_y_angular = self.wcs.pixelscale.y.to("deg") * pixels_y
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())
        fake = SimpleInstrument(distance=distance, inclination=inclination, azimuth=azimuth, position_angle=position_angle,
                                field_x=field_x_physical, field_y=field_y_physical, pixels_x=pixels_x, pixels_y=pixels_y,
                                center_x=center_x, center_y=center_y)

        # Add the instrument
        ski.add_instrument(instrument_name, fake)

        # Determine the path to the ski file
        ski_path = fs.join(self.images_bulge2d_path, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = fs.create_directory_in(self.images_bulge2d_path, "out")

        # Inform the user
        log.info("Running the bulge 2D simulation ...")

        # Simulate the bulge image
        fluxdensity = self.components["bulge"].fluxdensity
        self.bulge2d_image = self.launcher.run(ski_path, out_path, self.wcs, fluxdensity, self.psf)

        # Check WCS
        if self.bulge2d_image.wcs != self.wcs: raise RuntimeError("Something went wrong setting the coordinate system")

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

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Simulate the bulge image
        fluxdensity = self.components["bulge"].fluxdensity
        self.bulge_image = self.launcher.run(ski_path, out_path, self.wcs, fluxdensity, self.psf, instrument_name=instrument_name)

        # Check WCS
        if self.bulge_image.wcs != self.wcs: raise RuntimeError("Something went wrong setting the coordinate system")

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

        # Inform the user
        log.info("Running the disk simulation ...")

        # Simulate the disk image
        fluxdensity = self.components["disk"].fluxdensity
        self.disk_image = self.launcher.run(ski_path, out_path, self.wcs, fluxdensity, self.psf, instrument_name=instrument_name)

        # Check WCS
        if self.disk_image.wcs != self.wcs: raise RuntimeError("Something went wrong setting the coordinate system")

    # -----------------------------------------------------------------

    def simulate_model(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge+disk model image ...")

        # Load the disk ski file template
        model_template_path = fs.join(template_path, "model.ski")
        ski = SkiFile(model_template_path)

        # Set the number of photon packages
        ski.setpackages(self.config.disk_packages)

        # Change the ski file parameters
        ski.set_stellar_component_geometry(0, self.disk)
        ski.set_stellar_component_geometry(1, self.bulge)

        #print("disk", [self.components["disk"].rel_contribution])
        #print("bulge", [self.components["bulge"].rel_contribution])

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

        # Inform the user
        log.info("Running the model simulation ...")

        # Simulate the model image
        fluxdensity = self.components["bulge"].fluxdensity + self.components["disk"].fluxdensity  # sum of bulge and disk component flux density
        self.model_image = self.launcher.run(ski_path, out_path, self.wcs, fluxdensity, self.psf, instrument_name=instrument_name)

        # Check WCS
        if self.model_image.wcs != self.wcs: raise RuntimeError("Something went wrong setting the coordinate system")

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
        self.disk.saveto(self.disk_model_path)

        # Write the bulge model
        self.bulge.saveto(self.bulge_model_path)

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Write bulge 2D image
        self.write_bulge2d_image()

        # Write bulge image
        self.write_bulge_image()

        # Write disk image
        self.write_disk_image()

        # Write model image
        self.write_model_image()

    # -----------------------------------------------------------------

    def write_bulge2d_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge 2D image ...")

        # Determine the path to the bulge 2D image and save it
        bulge_2d_path = fs.join(self.bulge2d_image_path)
        self.bulge2d_image.saveto(bulge_2d_path)

    # -----------------------------------------------------------------

    def write_bulge_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the bulge image ...")

        # Determine the path to the bulge image and save it
        self.bulge_image.saveto(self.bulge_image_path)

    # -----------------------------------------------------------------

    def write_disk_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the disk image ...")

        # Determine the path to the disk image and save it
        self.disk_image.saveto(self.disk_image_path)

    # -----------------------------------------------------------------

    def write_model_image(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model image ...")

        # Determine the path to the model image and save it
        self.model_image.saveto(self.model_image_path)

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the projection systems ...")

        # Write the earth projection system
        self.projections["earth"].saveto(self.earth_projection_path)

        # Write the edgeon projection system
        self.projections["edgeon"].saveto(self.edgeon_projection_path)

        # Write the faceon projection system
        self.projections["faceon"].saveto(self.faceon_projection_path)

    # -----------------------------------------------------------------

    def write_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing regions file with disk ellipse ...")

        # Determine minor axis length: NO GALAXY PROPERTIES NOT ENTIRELY CONSISTENT WITH DECOMPOSITION!!
        #minor = (1.0 - self.galaxy_properties.ellipticity) * self.galaxy_properties.major_arcsec
        minor = self.disk_q * self.galaxy_properties.major_arcsec

        # Ellipse radius
        radius = SkyStretch(self.galaxy_properties.major_arcsec, minor)

        # Create sky ellipse
        sky_ellipse = SkyEllipseRegion(self.galaxy_properties.center, radius, self.disk_pa)

        # Create region
        region = SkyRegionList()
        region.append(sky_ellipse)
        region.saveto(self.disk_region_path)

# -----------------------------------------------------------------
