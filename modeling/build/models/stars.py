#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.models.stars Contains the StarsBuilder class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import numpy as np

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configuration import ConfigurationDefinition
from ....core.basics.configuration import InteractiveConfigurationSetter, prompt_proceed, prompt_string, prompt_yn, prompt_filepath
from ....core.units.parsing import parse_unit as u
from ...core.mappings import Mappings
from .general import GeneralBuilder
from ..suite import model_map_filename
from ....core.tools import filesystem as fs
from ....magic.core.frame import Frame
from ...component.galaxy import GalaxyModelingComponent
from ....core.tools import types
from ....magic.tools import extinction
from ....core.tools.utils import lazyproperty

# -----------------------------------------------------------------

basic_old_map_name = "old_disk"
basic_young_map_name = "young"
basic_ionizing_map_name = "ionizing"
basic_stellar_map_names = [basic_old_map_name, basic_young_map_name, basic_ionizing_map_name]

# -----------------------------------------------------------------

bulge_component_name = "bulge"
old_component_name = "old"
young_component_name = "young"
ionizing_component_name = "ionizing"
basic_stellar_component_names = [bulge_component_name, old_component_name, young_component_name, ionizing_component_name]

# -----------------------------------------------------------------

# Define titles for the different fixed components
titles = dict()
titles[bulge_component_name] = "Evolved stellar bulge"
titles[old_component_name] = "Evolved stellar disk"
titles[young_component_name] = "Young stars"
titles[ionizing_component_name] = "Ionizing stars"

# -----------------------------------------------------------------

component_name_for_map_name = dict()
component_name_for_map_name["old_disk"] = old_component_name
component_name_for_map_name["young"] = young_component_name
component_name_for_map_name["ionizing"] = ionizing_component_name

# -----------------------------------------------------------------

class StarsBuilder(GeneralBuilder, GalaxyModelingComponent):
    
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
        #super(StarsBuilder, self).__init__(*args, **kwargs)
        GeneralBuilder.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # Name of the maps from the maps selection
        self.old_map_name = None
        self.young_map_name = None
        self.ionizing_map_name = None

        # The Mappings template for the ionizing stellar component
        self.ionizing_mappings = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Build bulge
        if self.config.bulge: self.build_bulge()

        # 3. Set the evolved stellar disk component
        if self.config.old: self.build_old()

        # 4. Set the ionizing stellar component
        if self.config.ionizing: self.build_ionizing()

        # 5. Set the young stellar component
        if self.config.young: self.build_young()

        # 6. Build additional stellar components
        if self.config.additional: self.build_additional()

        # 7. Write
        self.write()

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        #super(StarsBuilder, self).setup()
        GeneralBuilder.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_flux(self):

        """
        This function ...
        :return: 
        """

        # Determine unattenuated flux
        factor = extinction.observed_to_intrinsic_factor(self.config.fuv_attenuation)
        return self.observed_flux(self.fuv_filter, unit="Jy") * factor

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_fuv_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_fuv_flux.to("W/micron", wavelength=self.fuv_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_ionizing_fuv_luminosity(self):

        """
        This function ...
        :return:
        """

        # Get the spectral luminosity at the FUV wavelength
        return self.ionizing_mappings.luminosity_at(self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_ionizing_fuv_flux(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_ionizing_fuv_luminosity.to("Jy", wavelength=self.fuv_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_ionizing_fuv_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_ionizing_fuv_luminosity.to("Lsun", density=True, density_strict=True, wavelength=self.fuv_wavelength)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_young_fuv_luminosity(self):

        """
        This function ...
        :return:
        """

        # Checks
        if self.intrinsic_ionizing_fuv_luminosity >= self.intrinsic_fuv_luminosity: raise ValueError("Cannot determine the initial normalization of young and ionizing component: intrinsic FUV luminosity of ionizing stars based on SFR is larger than the total unattenuated FUV luminosity")
        if self.intrinsic_ionizing_fuv_luminosity / self.intrinsic_fuv_luminosity > 0.5: log.warning("The contribution of ionizing stars to the intrinsic FUV luminosity is more than 50%")
        if self.intrinsic_ionizing_fuv_luminosity / self.intrinsic_fuv_luminosity < 0.1: log.warning("The contribution of ionizing stars to the intrinsic FUV luminosity is less than 10%")

        # Return the difference fo the total unattenuated FUV luminosity and the intrinsic FUV luminosity of the ionizing stars
        return self.intrinsic_fuv_luminosity - self.intrinsic_ionizing_fuv_luminosity

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_young_fuv_flux(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_young_fuv_luminosity.to("Jy", wavelength=self.fuv_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def intrinsic_young_fuv_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.intrinsic_young_fuv_luminosity.to("Lsun", density=True, density_strict=True, wavelength=self.fuv_wavelength)

    # -----------------------------------------------------------------

    def build_bulge(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the old stellar bulge component ...")

        # Get the parameters
        self.set_bulge_parameters()

        # Load the model
        self.load_bulge_model()

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_fluxdensity(self):

        """
        This function ...
        :return: 
        """

        # Get the flux density of the bulge
        fluxdensity = self.bulge2d_model.fluxdensity
        return fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.bulge_fluxdensity.to("W/micron", wavelength=self.i1_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.bulge_luminosity.to("Lsun", density=True, density_strict=True, wavelength=self.i1_wavelength)

    # -----------------------------------------------------------------

    def set_bulge_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the bulge component ...")

        # Check defaults
        if self.config.default_old_bulge_template is None: raise ValueError("Default old stellar bulge template cannot be undefined")
        if self.config.default_old_bulge_age is None: raise ValueError("Default bulge age cannot be undefined")
        if self.config.default_old_bulge_metallicity is None: raise ValueError("Default bulge metallicity cannot be undefined")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=self.config.default_old_bulge_template)
        definition.add_optional("age", "positive_real", "age in Gyr", default=self.config.default_old_bulge_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=self.config.default_old_bulge_metallicity)
        #definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=self.bulge_fluxdensity)
        definition.add_optional("neutral_luminosity", "photometric_density_quantity", "intrinsic neutral luminosity density at IRAC 3.6 micron wavelength", default=self.bulge_neutral_luminosity)

        # Get the parameters
        label = bulge_component_name
        config = self.get_parameters(label, definition)

        # Convert the flux density into a spectral luminosity
        #luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)
        #luminosity = config.fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_properties.distance)
        #assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Convert neutral luminosity to spectral luminosity
        luminosity = config.neutral_luminosity.to("W/micron", wavelength=self.i1_wavelength)

        # Set the luminosity
        config.filter = str(self.i1_filter)
        config.wavelength = self.i1_wavelength
        config.luminosity = luminosity

        # Calculate flux density
        config.fluxdensity = luminosity.to("Jy", wavelength=self.i1_wavelength, distance=self.galaxy_distance)

        # Set the title
        config.title = titles[bulge_component_name]

        # Set the bulge parameters
        self.parameters[bulge_component_name] = config

    # -----------------------------------------------------------------

    def load_bulge_model(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the bulge model ...")

        # Load the bulge model
        self.models[bulge_component_name] = self.bulge_model

    # -----------------------------------------------------------------

    def build_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the old stellar disk component ...")

        # 1. Get the parameters
        self.set_old_parameters()

        # 2. Load the map
        self.load_old_map()

        # 3. Create the deprojection
        self.create_deprojection_old()

    # -----------------------------------------------------------------

    @lazyproperty
    def old_scaleheight(self):

        """
        This function ...
        :return: 
        """

        scale_height = self.disk2d_model.scalelength / self.config.scalelength_to_scaleheight
        return scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def old_fluxdensity(self):

        """
        This function ...
        :return: 
        """

        # Get the flux
        bulge_fluxdensity = self.bulge2d_model.fluxdensity

        # Get the 3.6 micron flux density with the bulge subtracted
        fluxdensity = self.observed_flux(self.i1_filter, unit="Jy") - bulge_fluxdensity

        # Return the flux density
        return fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def old_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.old_fluxdensity.to("W/micron", wavelength=self.i1_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_neutral_luminosity(self):

        """
        This function ...
        :return:
        """

        return self.old_fluxdensity.to("Lsun", density=True, density_strict=True, wavelength=self.i1_wavelength, distance=self.galaxy_distance)

    # -----------------------------------------------------------------

    def set_old_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the old stellar disk component ...")

        # Check defaults
        if self.config.default_old_disk_template is None: raise ValueError("Default old stellar disk template cannot be undefined")
        if self.config.default_old_disk_age is None: raise ValueError("Default old stellar disk age cannot be undefined")
        if self.config.default_old_disk_metallicity is None: raise ValueError("Default old stellar disk metallicity cannot be undefined")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=self.config.default_old_disk_template)
        definition.add_optional("age", "positive_real", "age in Gyr", default=self.config.default_old_disk_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=self.config.default_old_disk_metallicity)
        definition.add_optional("scale_height", "quantity", "scale height", default=self.old_scaleheight)
        #definition.add_optional("fluxdensity", "photometric_quantity", "flux density", default=self.old_fluxdensity)
        definition.add_optional("neutral_luminosity", "photometric_density_quantity", "intrinsic neutral luminosity density at IRAC 3.6 micron wavelength", default=self.old_neutral_luminosity)

        # Set the label
        label = "old stellar disk"

        # Get the parameters
        config = self.get_parameters(label, definition)

        # Convert the flux density into a spectral luminosity
        #luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.i1_filter.pivot, self.galaxy_properties.distance)
        #luminosity = config.fluxdensity.to("W/micron", fltr=self.i1_filter, distance=self.galaxy_properties.distance)
        #assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Convert neutral luminosity to spectral luminosity
        luminosity = config.neutral_luminosity.to("W/micron", wavelength=self.i1_wavelength)

        # Set the luminosity
        config.filter = str(self.i1_filter)
        config.wavelength = self.i1_wavelength
        config.luminosity = luminosity

        # Set flux density
        config.fluxdensity = luminosity.to("Jy", wavelength=self.i1_wavelength, distance=self.galaxy_distance)

        # Set title
        config.title = titles[old_component_name]

        # Set the parameters
        self.parameters[old_component_name] = config

    # -----------------------------------------------------------------

    def load_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

        # Ask whether a custom map should be used
        custom = prompt_yn("custom", "use a custom map for the old stellar disk (instead of one of those created in the modeling pipeline)", default=False)

        # Load a custom old stellar disk map
        if custom: self.load_custom_old_map()

        # Load an old stellar disk map from modeling pipeline
        else: self.load_modeling_old_map()

    # -----------------------------------------------------------------

    def load_custom_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a custom map of the old stars ...")

        # Get the path
        path = prompt_filepath("filepath", "custom old stellar disk map path")

        # Load the map
        self.maps[old_component_name] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_modeling_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of old stars from the modeling pipeline ...")

        # Get the map names
        names = self.static_maps_selection.old_map_names

        # Ask for the old stellar map
        name = prompt_string("old_map", "old stellar disk map to use for this model", choices=names)

        # Set the name of the old stellar map
        self.old_map_name = name

        # Get the filepath
        filepath = self.static_maps_selection.old_map_paths[name]

        # Set the map
        self.maps[old_component_name] = Frame.from_file(filepath)

    # -----------------------------------------------------------------

    def create_deprojection_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the old stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.galaxy_properties, self.disk_position_angle,
                                                        self.maps[old_component_name], model_map_filename,
                                                        self.parameters[old_component_name].scale_height,
                                                        inclination=self.disk_inclination)

        # Set the deprojection model
        self.deprojections[old_component_name] = deprojection

    # -----------------------------------------------------------------

    def build_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building the young stellar component ...")

        # 1. Get the parameters
        self.set_young_parameters()

        # 2. Load the map
        self.load_young_map()

        # 3. Create deprojection
        self.create_deprojection_young()

    # -----------------------------------------------------------------

    @lazyproperty
    def young_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.config.young_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    def set_young_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the young stellar component ...")

        # Check defaults
        if self.config.default_young_template is None: raise ValueError("The default young stellar template cannot be undefined")
        if self.config.default_young_age is None: raise ValueError("The default young stellar age cannot be undefined")
        if self.config.default_young_metallicity is None: raise ValueError("The default young stellar metallicity cannot be undefined")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("template", "string", "template SED family", default=self.config.default_young_template)
        definition.add_optional("age", "positive_real", "age in Gyr", default=self.config.default_young_age)
        definition.add_optional("metallicity", "positive_real", "metallicity", default=self.config.default_young_metallicity)
        definition.add_optional("scale_height", "quantity", "scale height", default=self.young_scaleheight)
        definition.add_optional("neutral_luminosity", "photometric_density_quantity", "intrinsic neutral luminosity density at FUV wavelength", default=self.intrinsic_young_fuv_luminosity)

        # Set label
        label = "young stellar disk"

        # Get the parameters
        config = self.get_parameters(label, definition)

        # Convert neutral luminosity to spectral luminosity
        luminosity = config.neutral_luminosity.to("W/micron", wavelength=self.fuv_wavelength)

        # Set the luminosity
        config.filter = str(self.fuv_filter)
        config.wavelength = self.fuv_wavelength
        config.luminosity = luminosity

        # Set flux density
        config.fluxdensity = luminosity.to("Jy", wavelength=self.fuv_wavelength, distance=self.galaxy_distance)

        # Set the title
        config.title = titles[young_component_name]

        # Set the parameters
        self.parameters[young_component_name] = config

    # -----------------------------------------------------------------

    def load_young_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of young stars ...")

        # Ask whether a custom map should be used
        custom = prompt_yn("custom", "use a custom map for the young stars (instead of one of those created in the modeling pipeline)", default=False)

        # Load a custom young stellar disk map
        if custom: self.load_custom_young_map()

        # Load a young stellar disk map from modeling pipeline
        else: self.load_modeling_young_map()

    # -----------------------------------------------------------------

    def load_custom_young_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a custom map for the young stars ...")

        # Get the path
        path = prompt_filepath("filepath", "custom young stellar disk map path")

        # Load the map
        self.maps[young_component_name] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_modeling_young_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a young stellar disk map from the modeling pipeline ...")

        # Get the map names
        names = self.static_maps_selection.young_map_names

        # Ask for the young stellar map
        name = prompt_string("young_map", "young stellar disk map to use for this model", choices=names)

        # Set the young map name
        self.young_map_name = name

        # Get the path
        filepath = self.static_maps_selection.young_map_paths[name]

        # Set the map
        self.maps[young_component_name] = Frame.from_file(filepath)

    # -----------------------------------------------------------------

    def create_deprojection_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the young stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.galaxy_properties, self.disk_position_angle,
                                                        self.maps[young_component_name], model_map_filename,
                                                        self.parameters[young_component_name].scale_height,
                                                        inclination=self.disk_inclination)

        # Set the deprojection model
        self.deprojections[young_component_name] = deprojection

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

    @lazyproperty
    def ionizing_scaleheight(self):

        """
        This function ...
        :return: 
        """

        return self.config.ionizing_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    def set_ionizing_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Configuring the ionizing stellar component ...")

        # Check defaults
        if self.config.default_ionizing_metallicity is None: raise ValueError("The default ionizing stellar metallicity cannot be undefined")
        if self.config.default_ionizing_compactness is None: raise ValueError("The default ionizing stellar compactness cannot be undefined")
        if self.config.default_ionizing_pressure is None: raise ValueError("The default ionizing stellar pressure cannot be undefined")
        if self.config.default_covering_factor is None: raise ValueError("The default ionizing covering factor cannot be undefined")
        if self.config.default_sfr is None: raise ValueError("The default sSFR cannot be undefined")

        # Create definition
        definition = ConfigurationDefinition()
        definition.add_optional("metallicity", "positive_real", "metallicity", default=self.config.default_ionizing_metallicity)
        definition.add_optional("compactness", "positive_real", "compactness", default=self.config.default_ionizing_compactness)
        definition.add_optional("pressure", "quantity", "pressure", default=self.config.default_ionizing_pressure)
        definition.add_optional("covering_factor", "positive_real", "covering factor", default=self.config.default_covering_factor)
        definition.add_optional("scale_height", "quantity", "scale height", default=self.ionizing_scaleheight)
        definition.add_optional("sfr", "positive_real", "SFR", default=self.config.default_sfr)

        # Set label
        label = "ionizing stellar disk"

        # Get the parameters for this component
        config = self.get_parameters(label, definition)

        # Get the parameters
        metallicity = config.metallicity
        compactness = config.compactness
        pressure = config.pressure
        covering_factor = config.covering_factor
        sfr = config.sfr

        # Generate Mappings template for the specified parameters
        self.ionizing_mappings = Mappings(metallicity, compactness, pressure, covering_factor, sfr)
        # luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        # Get the spectral luminosity at the FUV wavelength
        luminosity = self.ionizing_mappings.luminosity_at(self.fuv_wavelength)

        # Set the luminosity
        config.filter = str(self.fuv_filter)
        config.wavelength = self.fuv_wavelength
        config.luminosity = luminosity

        # Set neutral luminosity and flux density
        config.neutral_luminosity = luminosity.to("Lsun", density=True, density_strict=True, wavelength=self.fuv_wavelength)
        config.fluxdensity = luminosity.to("Jy", wavelength=self.fuv_wavelength, distance=self.galaxy_distance)

        # Set title
        config.title = titles[ionizing_component_name]

        # Set the parameters
        self.parameters[ionizing_component_name] = config

    # -----------------------------------------------------------------

    def load_ionizing_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the map of ionizing stars ...")

        # Ask whether a custom map should be used
        custom = prompt_yn("custom", "use a custom map for the ionizing stellar disk (instead of one of those created in the modeling pipeline)", default=False)

        # Load a custom ionizing stellar disk map
        if custom: self.load_custom_ionizing_map()

        # Load a ionizing stellar disk map from modeling pipeline
        else: self.load_modeling_ionizing_map()

    # -----------------------------------------------------------------

    def load_custom_ionizing_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a custom map for the ionizing stellar disk ...")

        # Get the path
        path = prompt_filepath("filepath", "custom ionizing stellar disk map path")

        # Load the map
        self.maps[ionizing_component_name] = Frame.from_file(path)

    # -----------------------------------------------------------------

    def load_modeling_ionizing_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading a ionizing stellar disk map from the modeling pipeline ...")

        # Get map names
        names = self.static_maps_selection.ionizing_map_names

        # Ask for the ionizing stellar map
        name = prompt_string("ionizing_map", "ionizing stellar disk map to use for this model", choices=names)

        # Set the name of the ionizing stellar map
        self.ionizing_map_name = name

        # Set the path
        filepath = self.static_maps_selection.ionizing_map_paths[name]

        # Set the map
        self.maps[ionizing_component_name] = Frame.from_file(filepath)

    # -----------------------------------------------------------------

    def create_deprojection_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the deprojection model for the ionizing stellar disk ...")

        # Create the deprojection model
        deprojection = self.create_deprojection_for_map(self.galaxy_properties, self.disk_position_angle,
                                                        self.maps[ionizing_component_name], model_map_filename,
                                                        self.parameters[ionizing_component_name].scale_height,
                                                        inclination=self.disk_inclination)

        # Set the deprojection model
        self.deprojections[ionizing_component_name] = deprojection

    # -----------------------------------------------------------------

    def build_additional(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Building additional stellar components ...")

        # Proceed?
        while prompt_proceed("build additional stellar components?"):

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
        config = setter.run(definition, prompt_optional=True)

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

    @property
    def bulge_path(self):

        """
        This function ...
        :return:
        """

        component_name = bulge_component_name
        return self.paths[component_name]

    # -----------------------------------------------------------------

    @property
    def old_stars_path(self):

        """
        This function ...
        :return:
        """

        component_name = old_component_name
        return self.paths[component_name]

    # -----------------------------------------------------------------

    @property
    def young_stars_path(self):

        """
        This function ...
        :return:
        """

        component_name = young_component_name
        return self.paths[component_name]

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_path(self):

        """
        This function ...
        :return:
        """

        component_name = ionizing_component_name
        return self.paths[component_name]

    # -----------------------------------------------------------------

    @property
    def additional_names(self):

        """
        This function ...
        :return:
        """

        for name in self.component_names:
            if name in basic_stellar_component_names: continue
            else: yield name

    # -----------------------------------------------------------------

    @property
    def old_stars_map_path(self):

        """
        This function ...
        :return: 
        """

        component_name = component_name_for_map_name[basic_old_map_name]
        return self.map_paths[component_name]

    # -----------------------------------------------------------------

    @property
    def young_stars_map_path(self):

        """
        Thisf unction ...
        :return: 
        """

        component_name = component_name_for_map_name[basic_young_map_name]
        return self.map_paths[component_name]

    # -----------------------------------------------------------------

    @property
    def ionizing_stars_map_path(self):

        """
        This function ...
        :return: 
        """

        component_name = component_name_for_map_name[basic_ionizing_map_name]
        return self.map_paths[component_name]

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
    if types.is_string_type(wavelength_unit): wavelength_unit = u(wavelength_unit)
    if types.is_string_type(frequency_unit): frequency_unit = u(frequency_unit)

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
