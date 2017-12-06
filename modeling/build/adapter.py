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
from ...core.tools.stringify import tostr

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

        # The representations
        self.model_representations = OrderedDict()

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

        # 3. Load stellar components
        if self.stellar: self.load_stellar()

        # 4. Load dust components
        if self.dust: self.load_dust()

        # 5. Adapt stellar
        if self.stellar: self.adapt_stellar()

        # 6. Adapt dust
        if self.dust: self.adapt_dust()

        # 7. Load representations
        if self.representations: self.load_representations()

        # 8. Adapt representations
        if self.representations: self.adapt_representations()

        # 9. Write
        self.write()

        # 10. Show
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

        # Check
        if self.config.no_components and self.config.representations is None: raise ValueError("If not adapting components, select representations to adapt")
        if self.config.no_components and self.config.component_name is not None: raise ValueError("Not adapting components but component name is specified")

        # Set
        if self.config.no_components: self.config.dust_or_stellar = []

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

    @property
    def representations(self):

        """
        This function ...
        :return:
        """

        return self.config.representations is not None

    # -----------------------------------------------------------------

    @property
    def representation_names(self):

        """
        This function ...
        :return:
        """

        return self.config.representations

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
    def bulge(self):

        """
        This funtion ...
        :return:
        """

        return self.stellar_components[bulge_component_name]

    # -----------------------------------------------------------------

    @property
    def bulge_parameters(self):

        """
        This function ...
        :return:
        """

        return self.bulge.parameters

    # -----------------------------------------------------------------

    @property
    def bulge_model(self):

        """
        This function ...
        :return:
        """

        return self.bulge.model

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
    def old(self):

        """
        This function ...
        :return:
        """

        return self.stellar_components[old_component_name]

    # -----------------------------------------------------------------

    @property
    def old_parameters(self):

        """
        This function ...
        :return:
        """

        return self.old.parameters

    # -----------------------------------------------------------------

    @property
    def old_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.old.deprojection

    # -----------------------------------------------------------------

    @property
    def old_map_path(self):

        """
        Thisfunction ...
        :return:
        """

        return self.old.map_path

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
    def young(self):

        """
        Thisn function ...
        :return:
        """

        return self.stellar_components[young_component_name]

    # -----------------------------------------------------------------

    @property
    def young_parameters(self):

        """
        This function ...
        :return:
        """

        return self.young.parameters

    # -----------------------------------------------------------------

    @property
    def young_deprojection(self):

        """
        Thisf unction ...
        :return:
        """

        return self.young.deprojection

    # -----------------------------------------------------------------

    @property
    def young_map_path(self):

        """
        This function ...
        :return:
        """

        return self.young.map_path

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
    def ionizing(self):

        """
        This function ...
        :return:
        """

        return self.stellar_components[ionizing_component_name]

    # -----------------------------------------------------------------

    @property
    def ionizing_parameters(self):

        """
        This function ...
        :return:
        """

        return self.ionizing.parameters

    # -----------------------------------------------------------------

    @property
    def ionizing_deprojection(self):

        """
        Thisnf unction ...
        :return:
        """

        return self.ionizing.deprojection

    # -----------------------------------------------------------------

    @property
    def ionizing_map_path(self):

        """
        This function ...
        :return:
        """

        return self.ionizing.map_path

    # -----------------------------------------------------------------

    @property
    def has_extra_stellar(self):

        """
        This function ...
        :return:
        """

        return sequences.has_other(self.stellar_component_names, basic_stellar_component_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_stellar_component_names(self):

        """
        This function ...
        :return:
        """

        return sequences.get_other(self.stellar_component_names, basic_stellar_component_names)

    # -----------------------------------------------------------------

    def get_extra_stellar(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.stellar_components[name]

    # -----------------------------------------------------------------

    def get_extra_stellar_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.stellar_components[name].parameters

    # -----------------------------------------------------------------

    def get_extra_stellar_properties(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.stellar_components[name].properties

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

        # Inform the user
        log.info("Adapting the old stellar bulge component ...")

        for label in self.bulge_parameters: print(label, tostr(self.bulge_parameters[label]))
        print(self.bulge_model)

    # -----------------------------------------------------------------

    def adapt_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the old stellar disk component ...")

        for label in self.old_parameters: print(label, tostr(self.old_parameters[label]))
        print(self.old_deprojection)
        print(self.old_map_path)

    # -----------------------------------------------------------------

    def adapt_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the young stellar component ...")

        #for label in self.young_parameters: print(label, tostr(self.young_parameters[label]))
        #print(self.young_deprojection)
        #print(self.young_map_path)

        fuv_attenuation = self.fuv_attenuation

        #fuv_ionizing_contribution = 0.3
        factor = extinction.observed_to_intrinsic_factor(fuv_attenuation)
        unattenuated_fuv_flux = self.observed_flux(self.fuv_filter, unit="Jy") * factor

        #unattenuated_fuv_flux_young_stars = fuv_ionizing_contribution * unattenuated_fuv_flux

        print(unattenuated_fuv_flux)
        #print(self.young_parameters["fluxdensity"])

        unattenuated_fuv_luminosity = unattenuated_fuv_flux.to("W/micron", fltr=self.fuv_filter, distance=self.galaxy_distance)

        print(unattenuated_fuv_luminosity)

        # Get the FUV flux density
        #fluxdensity = unattenuated_fuv_flux_young_stars

        # Convert the flux density into a spectral luminosity
        #luminosity_manual = fluxdensity_to_luminosity(config.fluxdensity, self.fuv_filter.pivot, self.galaxy_properties.distance)
        #luminosity = config.fluxdensity.to("W/micron", fltr=self.fuv_filter, distance=self.galaxy_properties.distance)
        ##assert np.isclose(luminosity_manual.to("W/micron").value, luminosity.to("W/micron").value)

        # Set the luminosi
        #config.filter = str(self.fuv_filter)
        #config.luminosity = luminosity

        # Set the title
        #config.title = titles[young_component_name]

    # -----------------------------------------------------------------

    def adapt_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting the ionizing stellar component ...")

        for label in self.ionizing_parameters: print(label, tostr(self.ionizing_parameters[label]))
        print(self.ionizing_deprojection)
        print(self.ionizing_map_path)

        # Get relevant parameters
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

        for name in self.extra_stellar_component_names:

            parameters = self.get_extra_stellar_parameters(name)
            properties = self.get_extra_stellar_properties(name)

            for label in parameters: print(label, tostr(parameters[label]))
            print(properties)

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
    def disk(self):

        """
        This function ...
        :return:
        """

        return self.dust_components[disk_component_name]

    # -----------------------------------------------------------------

    @property
    def disk_parameters(self):

        """
        This function ...
        :return:
        """

        return self.disk.parameters

    # -----------------------------------------------------------------

    @property
    def disk_deprojection(self):

        """
        This function ...
        :return:
        """

        return self.disk.deprojection

    # -----------------------------------------------------------------

    @property
    def disk_map_path(self):

        """
        This function ...
        :return:
        """

        return self.disk.map_path

    # -----------------------------------------------------------------

    @property
    def has_extra_dust(self):

        """
        This function ...
        :return:
        """

        return sequences.has_other(self.dust_component_names, basic_dust_component_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def extra_dust_component_names(self):

        """
        This function ...
        :return:
        """

        return sequences.get_other(self.dust_component_names, basic_dust_component_names)

    # -----------------------------------------------------------------

    def get_extra_dust(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.dust_components[name]

    # -----------------------------------------------------------------

    def get_extra_dust_parameters(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_components[name].parameters

    # -----------------------------------------------------------------

    def get_extra_dust_properties(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.dust_components[name].properties

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

        for label in self.disk_parameters: print(label, tostr(self.disk_parameters[label]))
        print(self.disk_deprojection)
        print(self.disk_map_path)

    # -----------------------------------------------------------------

    def adapt_extra_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting extra dust components ...")

        for name in self.extra_dust_component_names:

            parameters = self.get_extra_dust_parameters(name)
            properties = self.get_extra_dust_properties(name)

            for label in parameters: print(label, tostr(parameters[label]))
            print(properties)

    # -----------------------------------------------------------------

    def load_representations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading model representations ...")

        # Loop over the representation names
        for name in self.representation_names:

            # Debugging
            log.debug("Loading the '" + name + "' model representation ...")

            # Load the representation
            representation = self.suite.get_representation(name)

            # Add
            self.model_representations[name] = representation

    # -----------------------------------------------------------------

    def get_projections(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        return self.model_representations[representation_name].get_projections()

    # -----------------------------------------------------------------

    def get_instruments(self, representation_name):

        """
        This function ...
        :param representation_name:
        :return:
        """

        return self.model_representations[representation_name].get_instruments()

    # -----------------------------------------------------------------

    def get_dust_grid(self, representation_name):

        """
        This function ....
        :param representation_name:
        :return:
        """

        return self.model_representations[representation_name].dust_grid

    # -----------------------------------------------------------------

    def adapt_representations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting model representations ...")

        # Projections
        self.adapt_projections()

        # Instruments
        self.adapt_instruments()

        # Dust grids
        self.adapt_dust_grids()

    # -----------------------------------------------------------------

    def adapt_projections(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adapting projections ...")

        # Loop over the representations
        for name in self.representation_names:

            # Debugging
            log.debug("Adapting the '" + name + "' projections ...")

            # Get projections
            projections = self.get_projections(name)

            # Loop over the projections
            for label in projections:

                # Get the projection
                projection = projections[label]

                # Propt
                projection.prompt_properties(recursive=True, contains=self.config.contains, not_contains=self.config.not_contains,
                                             exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                             startswith=self.config.startswith, endswith=self.config.endswith)

                # Save

    # -----------------------------------------------------------------

    def adapt_instruments(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Adapting instruments ...")

        # Loop over the representations
        for name in self.representation_names:

            # Debugging
            log.debug("Adapting the '" + name + "' instruments ...")

            # Get instruments
            instruments = self.get_instruments(name)

            # Loop over the instruments
            for label in instruments:

                # Get the instrument
                instrument = instruments[label]

                # Prompt
                instrument.prompt_properties(recursive=True, contains=self.config.contains, not_contains=self.config.not_contains,
                                             exact_name=self.config.exact_name, exact_not_name=self.config.exact_not_name,
                                             startswith=self.config.startswith, endswith=self.config.endswith)

                # Save

    # -----------------------------------------------------------------

    def adapt_dust_grids(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Adapting dust grids ...")

        # Loop over the representations
        for name in self.representation_names:

            # Debugging
            log.debug("Adapting the '" + name + "' dust grid ...")

            # Get dust grid
            grid = self.get_dust_grid(name)

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
