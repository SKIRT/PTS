#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.build.adapter Contains the GalaxyModelAdapter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ..component.galaxy import GalaxyModelingComponent
from ...core.basics.log import log
from ...dustpedia.core.database import get_cigale_parameters
from .component import BuildComponent
from ...core.tools.utils import lazyproperty
from ...core.tools import formatting as fmt
from .models.galaxy import show_component
from .models.galaxy import metallicities, solar_metallicity
from ..core.mappings import Mappings
from .models.stars import bulge_component_name, old_component_name, young_component_name, ionizing_component_name, basic_stellar_component_names
from .models.dust import disk_component_name, basic_dust_component_names
from ...core.tools import sequences
from ...magic.tools import extinction

# -----------------------------------------------------------------

stellar = "stellar"
dust = "dust"

# -----------------------------------------------------------------

class GalaxyModelAdapter(BuildComponent, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructors of the base classes
        BuildComponent.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The model definition
        self.definition = None

        # Global SED fitting parameters
        self.parameters = None

        # The stellar components
        self.stellar_components = OrderedDict()

        # The dust components
        self.dust_components = OrderedDict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get parameters
        self.get_parameters()

        # Load stellar components
        if self.stellar: self.load_stellar()

        # Load dust components
        if self.dust: self.load_dust()

        # 5. Adapt stellar
        if self.stellar: self.adapt_stellar()

        # 6. Adapt dust
        if self.dust: self.adapt_dust()

        # 7. Write
        self.write()

        # 8. Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup functions of the base classes
        BuildComponent.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

        # Load the model
        self.definition = self.get_model_definition(self.model_name)

        # Check configuration
        if self.config.component_name is not None:

            # Stellar component
            if self.only_stellar:

                stellar_component_name = self.config.component_name
                if not self.definition.is_stellar_component(stellar_component_name): raise ValueError("Stellar component '" + stellar_component_name + "' does not exist")

            # Dust component
            elif self.only_dust:

                dust_component_name = self.config.component_name
                if not self.definition.is_dust_component(dust_component_name): raise ValueError("Dust component '" + dust_component_name + "' does not exist")

            # Invalid
            else: raise ValueError("Cannot specifiy component name when both 'stellar' and 'dust' are used")

    # -----------------------------------------------------------------

    @property
    def stellar(self):

        """
        This function ...
        :return:
        """

        return stellar in self.config.dust_or_stellar

    # -----------------------------------------------------------------

    @property
    def only_stellar(self):

        """
        This function ...
        :return:
        """

        return self.config.dust_or_stellar == [stellar]

    # -----------------------------------------------------------------

    @lazyproperty
    def single_stellar_component(self):

        """
        This function ...
        :return:
        """

        if not self.only_stellar: raise ValueError("Not only stellar")
        if self.config.component_name is None: raise ValueError("Not a single stellar component")
        return self.definition.get_stellar_component(self.config.component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_component_names(self):

        """
        This function ...
        :return:
        """

        if self.only_dust: raise ValueError("Only dust")
        if self.config.component_name is None: return self.definition.stellar_component_names
        elif not self.only_stellar: raise ValueError("Not only stellar")
        else: return [self.config.component_name]

    # -----------------------------------------------------------------

    @property
    def dust(self):

        """
        This function ...
        :return:
        """

        return dust in self.config.dust_or_stellar

    # -----------------------------------------------------------------

    @property
    def only_dust(self):

        """
        This function ...
        :return:
        """

        return self.config.dust_or_stellar == [dust]

    # -----------------------------------------------------------------

    @lazyproperty
    def single_dust_component(self):

        """
        This function ...
        :return:
        """

        if not self.only_dust: raise ValueError("Not only dust")
        if self.config.component_name is None: raise ValueError("Not a single dust component")
        return self.definition.get_dust_component(self.config.component_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component_names(self):

        """
        This function ...
        :return:
        """

        if self.only_stellar: raise ValueError("Only stellar")
        if self.config.component_name is None: return self.definition.dust_component_names
        elif not self.only_dust: raise ValueError("Not only dust")
        else: return [self.config.component_name]

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the global SED fitting results ...")

        # Get the parameters from the DustPedia database
        self.parameters = get_cigale_parameters(self.ngc_name_nospaces)

    # -----------------------------------------------------------------

    @property
    def dust_mass(self):

        """
        This function ...
        :return:
        """

        return self.parameters["dust_mass"]

    # -----------------------------------------------------------------

    @property
    def sfr(self):

        """
        This function ...
        :return:
        """

        return self.parameters["sfr"]

    # -----------------------------------------------------------------

    @property
    def sfr_msun_per_year(self):

        """
        Thisf unction ...
        :return:
        """

        return self.sfr.to("Msun/yr").value

    # -----------------------------------------------------------------

    @property
    def fuv_attenuation(self):

        """
        This function ...
        :return:
        """

        return self.parameters["fuv_attenuation"]

    # -----------------------------------------------------------------

    @property
    def metallicity(self):

        """
        Thisf unction ...
        :return:
        """

        if self.config.metallicity is not None: return self.config.metallicity
        elif self.ngc_name_nospaces in metallicities: return metallicities[self.ngc_name_nospaces]
        else: return solar_metallicity

    # -----------------------------------------------------------------

    def load_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading stellar components ...")

        # Loop over the names
        for name in self.stellar_component_names:

            # Load the component
            component = self.definition.get_stellar_component(name)

            # Add the component
            self.stellar_components[name] = component

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading dust components ...")

        # Loop over the names
        for name in self.dust_component_names:

            # Load the component
            component = self.definition.get_dust_component(name)

            # Add the component
            self.dust_components[name] = component

    # -----------------------------------------------------------------

    @property
    def has_bulge(self):

        """
        This function ...
        :return:
        """

        return bulge_component_name in self.stellar_component_names

    # -----------------------------------------------------------------

    @property
    def has_old(self):

        """
        This function ...
        :return:
        """

        return old_component_name in self.stellar_component_names

    # -----------------------------------------------------------------

    @property
    def has_young(self):

        """
        This function ...
        :return:
        """

        return young_component_name in self.stellar_component_names

    # -----------------------------------------------------------------

    @property
    def has_ionizing(self):

        """
        This ufnction ...
        :return:
        """

        return ionizing_component_name in self.stellar_component_names

    # -----------------------------------------------------------------

    @property
    def has_extra_stellar(self):

        """
        This function ...
        :return:
        """

        return sequences.has_other(self.stellar_component_names, basic_stellar_component_names)

    # -----------------------------------------------------------------

    def adapt_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting stellar components ...")

        # Bulge
        if self.has_bulge: self.adapt_bulge()

        # Old
        if self.has_old: self.adapt_old()

        # Young
        if self.has_young: self.adapt_young()

        # Ionizing
        if self.has_ionizing: self.adapt_ionizing()

        # Extra
        if self.has_extra_stellar: self.adapt_extra_stellar()

    # -----------------------------------------------------------------

    def adapt_bulge(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def adapt_old(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def adapt_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the young stellar component ...")

        fuv_attenuation = self.fuv_attenuation

        fuv_ionizing_contribution = 0.3
        factor = extinction.observed_to_intrinsic_factor(fuv_attenuation)
        unattenuated_fuv_flux = self.observed_flux(self.fuv_filter, unit="Jy") * factor

        unattenuated_fuv_flux_young_stars = fuv_ionizing_contribution * unattenuated_fuv_flux

        # Get the FUV flux density
        fluxdensity = unattenuated_fuv_flux_young_stars

    # -----------------------------------------------------------------

    def adapt_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the ionizing stellar component ...")

        metallicity = self.metallicity
        compactness = self.config.default_ionizing_compactness
        pressure = self.config.default_ionizing_pressure
        covering_factor = self.config.default_covering_factor
        sfr = self.sfr_msun_per_year

        print("SFR", sfr)

        # Generate Mappings template for the specified parameters
        mappings = Mappings(metallicity, compactness, pressure, covering_factor, sfr)
        # luminosity = luminosity.to(self.sun_fuv).value # for normalization by band

        # Get the spectral luminosity at the FUV wavelength
        luminosity = mappings.luminosity_at(self.fuv_filter.pivot)

        print(luminosity)

        print(luminosity.to("W", density=True, wavelength=self.fuv_filter.pivot))

        # Set the luminosity
        #config.filter = str(self.fuv_filter)
        #config.luminosity = luminosity

        # Set title
        #config.title = titles[ionizing_component_name]

        # Set the parameters
        #self.parameters[ionizing_component_name] = config

    # -----------------------------------------------------------------

    def adapt_extra_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting extra stellar components ...")

    # -----------------------------------------------------------------

    @property
    def has_disk(self):

        """
        This function ...
        :return:
        """

        return disk_component_name in self.dust_component_names

    # -----------------------------------------------------------------

    @property
    def has_extra_dust(self):

        """
        This function ...
        :return:
        """

        return sequences.has_other(self.dust_component_names, basic_dust_component_names)

    # -----------------------------------------------------------------

    def adapt_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting dust components ...")

        # Dust disk
        if self.has_disk: self.adapt_disk()

        # Extra
        if self.has_extra_dust: self.adapt_extra_dust()

    # -----------------------------------------------------------------

    def adapt_disk(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting dust disk component ...")

    # -----------------------------------------------------------------

    def adapt_extra_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting extra dust components ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write stellar components
        if self.stellar: self.write_stellar()

        # Write dust components
        if self.dust: self.write_dust()

    # -----------------------------------------------------------------

    def write_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing stellar components ...")

    # -----------------------------------------------------------------

    def write_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing dust components ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model composition ...")

        # Show stellar
        self.show_stellar()

        # Show dust
        self.show_dust()

    # -----------------------------------------------------------------

    def show_stellar(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model stellar components ...")

        print("")
        print(fmt.green + fmt.underlined + "STELLAR COMPONENTS" + fmt.reset)
        print("")

        # Loop over the stellar components
        for name in self.definition.stellar_component_names:

            # Show name
            print("  " + fmt.magenta + name.upper() + fmt.reset + ": ")
            print("")

            # Get path
            path = self.definition.get_stellar_component_path(name)
            map_path = None

            # Show
            show_component(path, map_path, line_prefix="    ")
            print("")

    # -----------------------------------------------------------------

    def show_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the model dust components ...")

        print("")
        print(fmt.green + fmt.underlined + "DUST COMPONENTS" + fmt.reset)
        print("")

        # Loop over the dust components
        for name in self.definition.dust_component_names:

            # Show name
            print("  " + fmt.magenta + name.upper() + fmt.reset + ": ")
            print("")

            # Get path
            path = self.definition.get_dust_component_path(name)
            map_path = None

            # Show
            show_component(path, map_path, line_prefix="    ")
            print("")

# -----------------------------------------------------------------
