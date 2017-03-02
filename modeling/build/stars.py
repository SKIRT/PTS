#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.stars Contains the StarsBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math

# Import the relevant PTS classes and modules
from .component import BuildComponent
from ...core.tools.logging import log
from ...core.basics.configuration import ConfigurationDefinition
from ...core.basics.configuration import InteractiveConfigurationSetter, prompt_proceed
from ...core.basics.unit import parse_unit as u
from ..core.mappings import Mappings
from ...core.prep.smile import SKIRTSmileSchema

# -----------------------------------------------------------------

class StarsBuilder(BuildComponent):
    
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
        super(StarsBuilder, self).__init__(config, interactive)

        # The parameters
        self.parameters = dict()

        # The stellar components
        self.components = dict()

        # The deprojections
        self.deprojections = dict()

        # The SKIRT smile schema
        self.smile = None

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Build bulge
        if self.config.bulge: self.build_bulge()

        # 3. Set the evolved stellar disk component
        if self.config.old: self.build_old()

        # 4. Set the young stellar component
        if self.config.young: self.build_young()

        # 5. Set the ionizing stellar component
        if self.config.ionizing: self.build_ionizing()

        # 6. Build additional stellar components
        if self.config.additional: self.build_additional()

        # 7. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(StarsBuilder, self).setup()

        # Create the SKIRT smile schema
        self.smile = SKIRTSmileSchema()

    # -----------------------------------------------------------------

    def build_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the old stellar bulge component ...")

        # Get the parameters
        self.get_bulge_parameters()

        # Load the model
        self.load_bulge_model()

    # -----------------------------------------------------------------

    def get_bulge_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the bulge component ...")

        # Like M31
        bulge_template = "BruzualCharlot"
        bulge_age = 12
        bulge_metallicity = 0.03

        # Get the flux density of the bulge
        fluxdensity = self.bulge2d_model.fluxdensity

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=bulge_template, choices=[bulge_template])
        definition.add_optional("age", "positive_real", "age in Gyr", default=bulge_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=bulge_metallicity)
        definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=fluxdensity)

        # Prompt for the values
        setter = InteractiveConfigurationSetter("bulge")
        config = setter.run(definition)

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)

        # Set the luminosity
        config.luminosity = luminosity

        # Set the bulge parameters
        self.parameters["bulge"] = config

    # -----------------------------------------------------------------

    def load_bulge_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the bulge model ...")

        # Set the parameters of the bulge
        #self.ski_template.set_stellar_component_geometry("Evolved stellar bulge", self.bulge_model)
        #self.ski_template.set_stellar_component_sed("Evolved stellar bulge", bulge_template, bulge_age, bulge_metallicity)  # SED
        ## self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1_filter.centerwavelength() * u("micron"))

    # -----------------------------------------------------------------

    def build_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the old stellar disk component ...")

        # Get the parameters
        self.get_old_parameters()

        # Load the map
        self.load_old_map()

        # Create the deprojection
        self.create_deprojection_old()

    # -----------------------------------------------------------------

    def get_old_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the old stellar disk component ...")

        # Like M31
        disk_template = "BruzualCharlot"
        disk_age = 8
        # disk_metallicity = 0.02
        disk_metallicity = 0.03

        # Get the scale height
        scale_height = self.disk2d_model.scalelength / 8.26  # De Geyter et al. 2014
        bulge_fluxdensity = self.bulge2d_model.fluxdensity

        # Get the 3.6 micron flux density with the bulge subtracted
        fluxdensity = self.observed_flux(self.i1_filter, unit="Jy") - bulge_fluxdensity

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=disk_template, choices=[disk_template])
        definition.add_optional("age", "positive_real", "age in Gyr", default=disk_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=disk_metallicity)
        definition.add_optional("scale_height", "quantity", "scale height", default=scale_height)
        definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=fluxdensity)

        # Prompt for the values
        setter = InteractiveConfigurationSetter("old stellar disk")
        config = setter.run(definition)

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)

        # Set the luminosity
        config.luminosity = luminosity

        # Set the parameters
        self.parameters["old"] = config

        ## NEW: SET FIXED PARAMETERS
        #self.fixed["metallicity"] = disk_metallicity
        #self.fixed["old_scaleheight"] = scale_height
        #self.fixed["i1_old"] = luminosity
        ##

        # Get the spectral luminosity in solar units
        # luminosity = luminosity.to(self.sun_i1).value

        # Set the parameters of the evolved stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.old_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["old stars"] = deprojection

        # Adjust the ski file
        #self.ski_template.set_stellar_component_geometry("Evolved stellar disk", deprojection)
        #self.ski_template.set_stellar_component_sed("Evolved stellar disk", disk_template, disk_age, disk_metallicity)  # SED
        # self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1) # normalization by band
        #self.ski_template.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1_filter.centerwavelength() * u("micron"))

    # -----------------------------------------------------------------

    def load_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of old stars ...")



    # -----------------------------------------------------------------

    def create_deprojection_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the old stellar disk ...")

        filename = None
        hz = None

        # Get the galaxy distance, the inclination and position angle
        distance = self.galaxy_properties.distance
        inclination = self.galaxy_properties.inclination
        pa = self.earth_projection.position_angle

        ## NEW: SET FIXED PARAMETERS
        #self.fixed["distance"] = distance
        #self.fixed["inclination"] = inclination
        #self.fixed["position_angle"] = pa
        ##

        # Get the center pixel
        pixel_center = self.galaxy_properties.center.to_pixel(self.reference_wcs)
        xc = pixel_center.x
        yc = pixel_center.y

        # Get the pixelscale in physical units
        pixelscale_angular = self.reference_wcs.average_pixelscale.to("deg")  # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the number of x and y pixels
        x_size = self.reference_wcs.xsize
        y_size = self.reference_wcs.ysize

        # Create the deprojection model
        deprojection = DeprojectionModel3D(filename, pixelscale, pa, inclination, x_size, y_size, xc, yc, hz)

    # -----------------------------------------------------------------

    def build_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the young stellar component ...")

        # Get the parameters
        self.get_young_parameters()

        # Load the map
        self.load_young_map()

        # Create deprojection
        self.create_deprojection_young()

    # -----------------------------------------------------------------

    def get_young_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the young stellar component ...")

        # Like M31
        young_template = "BruzualCharlot"
        young_age = 0.1
        # young_metallicity = 0.02
        young_metallicity = 0.03

        # Get the scale height
        # scale_height = 150 * Unit("pc") # first models
        scale_height = 100. * u("pc")  # M51

        # Get the FUV flux density
        fluxdensity = 2. * self.observed_flux(self.fuv_filter, unit="Jy")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=young_template, choices=[young_template])
        definition.add_optional("age", "positive_real", "age in Gyr", default=young_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=young_metallicity)
        definition.add_optional("scale_height", "quantity", "scale height", default=scale_height)
        definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=fluxdensity)

        # Prompt for the values
        setter = InteractiveConfigurationSetter("young stellar disk")
        config = setter.run(definition)

        # Convert the flux density into a spectral luminosity
        luminosity = fluxdensity_to_luminosity(fluxdensity, self.fuv_filter.pivot, self.galaxy_properties.distance)

        # Set the luminosity
        config.luminosity = luminosity

        # Set the parameters
        self.parameters["young"] = config

        # Get the spectral luminosity in solar units
        # luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        ## NEW: SET FIXED PARAMETERS
        #self.fixed["young_scaleheight"] = scale_height
        ##

        # Set the parameters of the young stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.young_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["young stars"] = deprojection

        # Adjust the ski file
        #self.ski_template.set_stellar_component_geometry("Young stars", deprojection)
        #self.ski_template.set_stellar_component_sed("Young stars", young_template, young_age, young_metallicity)  # SED
        # self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv) # normalization by band
        # self.ski_template.set_stellar_component_luminosity("Young stars", luminosity, self.fuv_filter.centerwavelength() * u("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        #self.ski_template.set_stellar_component_normalization_wavelength("Young stars", self.fuv_filter.centerwavelength() * u("micron"))
        #self.ski_template.set_labeled_value("fuv_young", luminosity)

    # -----------------------------------------------------------------

    def load_young_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_deprojection_young(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def build_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the ionizing stellar component ...")

        # Set the parameters
        self.set_ionizing_parameters()

        # Load map
        self.load_ionizing_map()

        # Create deprojection model
        self.create_deprojection_ionizing()

    # -----------------------------------------------------------------

    def set_ionizing_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the ionizing stellar component ...")

        # Like M51 and M31
        # ionizing_metallicity = 0.02
        ionizing_metallicity = 0.03  # XU KONG et al. 2000
        ionizing_compactness = 6.
        ionizing_pressure = 1e12 * u("K/m3")
        ionizing_covering_factor = 0.2

        # Get the scale height
        # scale_height = 150 * Unit("pc") # first models
        scale_height = 100. * u("pc")  # M51

        ## NEW: SET FIXED PARAMETERS
        #self.fixed["ionizing_scaleheight"] = scale_height
        #self.fixed["sfr_compactness"] = ionizing_compactness
        #self.fixed["sfr_covering"] = ionizing_covering_factor
        #self.fixed["sfr_pressure"] = ionizing_pressure
        ##

        # Convert the SFR into a FUV luminosity
        sfr = 0.8  # The star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)

        # Create definition
        definition = ConfigurationDefinition()
        #definition.add_optional("template", "string", "template SED family", default=young_template, choices=[young_template])
        #definition.add_optional("age", "positive_real", "age in Gyr", default=young_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=ionizing_metallicity)
        definition.add_optional("compactness", "positive_real", "compactness", default=ionizing_compactness)
        definition.add_optional("pressure", "quantity", "pressure", default=ionizing_pressure)
        definition.add_optional("covering_factor", "positive_real", "covering factor", default=ionizing_covering_factor)
        definition.add_optional("scale_height", "quantity", "scale height", default=scale_height)
        #definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=fluxdensity)
        definition.add_optional("sfr", "positive_real", "SFR", default=sfr)

        # Prompt for the values
        setter = InteractiveConfigurationSetter("ionizing stellar disk")
        config = setter.run(definition)

        # Generate Mappings template for the specified parameters
        mappings = Mappings(ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor, sfr)
        # luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        # Get the spectral luminosity at the FUV wavelength
        luminosity = mappings.luminosity_at(self.fuv_filter.pivot)

        # Set the luminosity
        config.luminosity = luminosity

        # Set the parameters
        self.parameters["ionizing"] = config

        # Set the parameters of the ionizing stellar component
        #deprojection = self.deprojection.copy()
        #deprojection.filename = self.ionizing_stellar_map_filename
        #deprojection.scale_height = scale_height
        #self.deprojections["ionizing stars"] = deprojection

        # Adjust the ski file
        #self.ski_template.set_stellar_component_geometry("Ionizing stars", deprojection)
        #self.ski_template.set_stellar_component_mappingssed("Ionizing stars", ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor)  # SED
        # self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv) # normalization by band
        # self.ski_template.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv_filter.centerwavelength() * Unit("micron"))

        # SET NORMALIZATION (IS FREE PARAMETER)
        #self.ski_template.set_stellar_component_normalization_wavelength("Ionizing stars", self.fuv_filter.centerwavelength() * u("micron"))
        #self.ski_template.set_labeled_value("fuv_ionizing", luminosity)  # keep label

    # -----------------------------------------------------------------

    def load_ionizing_map(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def create_deprojection_ionizing(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def build_additional(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building additional stellar components ...")

        # Proceed?
        while prompt_proceed():

            # Set parameters
            name = self.set_additional_parameters()

            # Set properties
            normalization_parameters = self.set_normalization_properties(name)

            # Set geometry properties
            geometry_parameters = self.set_geometry_properties(name)

    # -----------------------------------------------------------------

    def set_additional_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring an additional stellar component ...")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_required("name", "string", "name for this stellar component")
        definition.add_optional("description", "string", "description for the component")
        definition.add_optional("geometry", "string", "SKIRT base geometry for the component", self.smile.concrete_geometries)
        definition.add_optional("normalization", "string", "normalization for the component", self.smile.concrete_stellar_normalizations)

        # Prompt for settings
        setter = InteractiveConfigurationSetter("additional stellar component", add_cwd=False, add_logging=False)
        config = setter.run(definition)

        # Set the parameters
        self.parameters[config.name] = config

        # Return the name of the new component
        return config.name

    # -----------------------------------------------------------------

    def set_normalization_properties(self, name):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the normalization of stellar component '" + name + "' ...")

        # Get the selected type of normalization
        normalization_type = self.parameters[name].normalization

        # Get parameters for this simulation item
        parameters, children = self.smile.prompt_parameters_for_type(normalization_type)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def set_geometry_properties(self, name):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the geometry of stellar component '" + name + "' ...")

        # Get the selected type of geometry
        geometry_type = self.parameters[name].geometry

        # Get parameters for this simulation item
        parameters, children = self.smile.prompt_parameters_for_type(geometry_type)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

# -----------------------------------------------------------------

# The speed of light
speed_of_light = constants.c

# -----------------------------------------------------------------

def spectral_factor_hz_to_micron(wavelength):

    """
    This function ...
    :param wavelength:
    :return:
    """

    wavelength_unit = "micron"
    frequency_unit = "Hz"

    # Convert string units to Unit objects
    if isinstance(wavelength_unit, basestring): wavelength_unit = u(wavelength_unit)
    if isinstance(frequency_unit, basestring): frequency_unit = u(frequency_unit)

    conversion_factor_unit = wavelength_unit / frequency_unit

    # Calculate the conversion factor
    factor = (wavelength ** 2 / speed_of_light).to(conversion_factor_unit).value
    return 1. / factor

# -----------------------------------------------------------------

def fluxdensity_to_luminosity(fluxdensity, wavelength, distance):

    """
    This function ...
    :param fluxdensity:
    :param wavelength:
    :param distance:
    :return:
    """

    luminosity = (fluxdensity * 4. * math.pi * distance ** 2.).to("W/Hz")

    # 3 ways:
    #luminosity_ = luminosity.to("W/micron", equivalencies=spectral_density(wavelength)) # does not work
    luminosity_ = (speed_of_light * luminosity / wavelength**2).to("W/micron")
    luminosity = luminosity.to("W/Hz").value * spectral_factor_hz_to_micron(wavelength) * u("W/micron")
    #print(luminosity_, luminosity) # is OK!

    return luminosity

# -----------------------------------------------------------------
