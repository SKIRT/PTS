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
import numpy as np

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.basics.configuration import ConfigurationDefinition
from ...core.basics.configuration import InteractiveConfigurationSetter, prompt_proceed, prompt_string
from ...core.basics.unit import parse_unit as u
from ..core.mappings import Mappings
from .general import GeneralBuilder
from .component import model_map_filename
from ...core.tools import filesystem as fs
from ...magic.core.frame import Frame
from ..maps.component import get_old_stars_maps_path, get_young_stars_maps_path, get_ionizing_stars_maps_path

# -----------------------------------------------------------------

# Define titles for the different fixed components
titles = dict()
titles["bulge"] = "Evolved stellar bulge"
titles["old"] = "Evolved stellar disk"
titles["young"] = "Young stars"
titles["ionizing"] = "Ionizing stars"

# -----------------------------------------------------------------

class StarsBuilder(GeneralBuilder):
    
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
        bulge_age = 10
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
        luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)
        luminosity = config.fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_properties.distance)
        assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Set the luminosity
        config.filter = str(self.i1_filter)
        config.luminosity = luminosity

        # Set the title
        config.title = titles["bulge"]

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

        self.models["bulge"] = None

    # -----------------------------------------------------------------

    def build_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the old stellar disk component ...")

        # 1. Get the parameters
        self.get_old_parameters()

        # 2. Load the map
        self.load_old_map()

        # 3. Create the deprojection
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
        luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)
        luminosity = config.fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_properties.distance)
        assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Set the luminosity
        config.filter = str(self.i1_filter)
        config.luminosity = luminosity

        # Set title
        config.title = titles["old"]

        # Set the parameters
        self.parameters["old"] = config

    # -----------------------------------------------------------------

    def load_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

        # Get the path to the old stellar maps directory
        directory_path = get_old_stars_maps_path(self.config.path)

        # Get the present filenames
        names = fs.files_in_path(directory_path, extension="fits", returns="name")

        # Ask for the old stellar map
        name = prompt_string("old_map", "old stellar disk map to use for this model", choices=names)
        path = fs.join(directory_path, name + ".fits")

        # Set the map
        self.maps["old"] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def create_deprojection_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the old stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.maps["old"], model_map_filename, self.parameters["old"].scale_height)

        # Set the deprojection model
        self.deprojections["old"] = deprojection

    # -----------------------------------------------------------------

    def build_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the young stellar component ...")

        # 1. Get the parameters
        self.get_young_parameters()

        # 2. Load the map
        self.load_young_map()

        # 3. Create deprojection
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
        luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.fuv_filter.pivot, self.galaxy_properties.distance)
        luminosity = config.fluxdensity.to("W/micron", fltr=self.fuv_filter, distance=self.galaxy_properties.distance)
        assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Set the luminosi
        config.filter = str(self.fuv_filter)
        config.luminosity = luminosity

        # Set the title
        config.title = titles["young"]

        # Set the parameters
        self.parameters["young"] = config

    # -----------------------------------------------------------------

    def load_young_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of young stars ...")

        # Get the path to the young stellar maps directory
        directory_path = get_young_stars_maps_path(self.config.path)

        # Get the present filenames
        names = fs.files_in_path(directory_path, extension="fits", returns="name")

        # Ask for the young stellar map
        name = prompt_string("young_map", "young stellar disk map to use for this model", choices=names)
        path = fs.join(directory_path, name + ".fits")

        # Set the map
        self.maps["young"] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def create_deprojection_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the young stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.maps["young"], model_map_filename, self.parameters["young"].scale_height)

        # Set the deprojection model
        self.deprojections["young"] = deprojection

    # -----------------------------------------------------------------

    def build_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the ionizing stellar component ...")

        # 1. Set the parameters
        self.set_ionizing_parameters()

        # 2. Load map
        self.load_ionizing_map()

        # 3. Create deprojection model
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
        config.filter = str(self.fuv_filter)
        config.luminosity = luminosity

        # Set title
        config.title = titles["ionizing"]

        # Set the parameters
        self.parameters["ionizing"] = config

    # -----------------------------------------------------------------

    def load_ionizing_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of ionizing stars ...")

        # Get the path to the ionizing stellar maps directory
        directory_path = get_ionizing_stars_maps_path(self.config.path)

        # Get the present filenames
        names = fs.files_in_path(directory_path, extension="fits", returns="name")

        # Ask for the ionizing stellar map
        name = prompt_string("ionizing_map", "ionizing stellar disk map to use for this model", choices=names)
        path = fs.join(directory_path, name + ".fits")

        # Set the map
        self.maps["ionizing"] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def create_deprojection_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the ionizing stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.maps["ionizing"], model_map_filename, self.parameters["ionizing"].scale_height)

        # Set the deprojection model
        self.deprojections["ionizing"] = deprojection

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

            # Set SED properties
            sed_parameters = self.set_sed_properties(name)

            # Set properties
            normalization_parameters = self.set_normalization_properties(name)

            # Set geometry properties
            geometry_parameters = self.set_geometry_properties(name)

            # Set the properties
            properties = {"sed": sed_parameters, "normalization": normalization_parameters, "geometry": geometry_parameters}

            # Load map if necessary
            if self.parameters[name].geometry == "ReadFitsGeometry": self.load_additional_map(name, geometry_parameters)

            # Add the properties
            self.properties[name] = properties

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
        definition.add_optional("title", "string", "short description for this component")
        definition.add_optional("geometry", "string", "SKIRT base geometry for the component", self.smile.concrete_geometries)
        definition.add_optional("sed", "string", "SED template for the component", self.smile.concrete_stellar_seds)
        definition.add_optional("normalization", "string", "normalization for the component", self.smile.concrete_stellar_normalizations)

        # Prompt for settings
        setter = InteractiveConfigurationSetter("additional stellar component", add_cwd=False, add_logging=False)
        config = setter.run(definition)

        # Check the name
        if config.name in self.parameters: raise ValueError("You cannot use this name for the stellar component: already in use")

        # Set the parameters
        self.parameters[config.name] = config

        # Return the name of the new component
        return config.name

    # -----------------------------------------------------------------

    def set_sed_properties(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Configuring the SED template of stellar component '" + name + "' ...")

        # Get the selected type of SED
        sed_type = self.parameters[name].sed

        # Get parameters for this simulation item
        parameters = self.smile.prompt_parameters_for_type(sed_type, merge=True)

        # Return the parameters
        return parameters

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
        parameters = self.smile.prompt_parameters_for_type(normalization_type, merge=True)

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
        parameters = self.smile.prompt_parameters_for_type(geometry_type, merge=True)

        # Return the parameters
        return parameters

    # -----------------------------------------------------------------

    def load_additional_map(self, name, parameters):

        """
        This function ...
        :param name:
        :param parameters:
        :return:
        """

        # Inform the user
        log.info("Loading the input map for stellar component '" + name + "' ...")

        # Get the absolute file path to the map
        path = fs.absolute_path(parameters["filename"])

        # Load the map
        self.maps[name] = Frame.from_file(path)

        # Change the filename in the geometry parameters
        parameters["filename"] = model_map_filename

# -----------------------------------------------------------------

# Import astronomical modules
from astropy import constants

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
