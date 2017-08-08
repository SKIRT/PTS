#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.summarize_ski Summarize the contents of a ski file.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import Quantity

# Import the relevant PTS classes and modules
from ..simulation.skifile import SkiFile, LabeledSkiFile
from ..tools import formatting as fmt
from ..tools import stringify
from ..filter.filter import Filter
from ..units.stringify import represent_quantity
from ..basics.configurable import Configurable
from ..basics.log import log

# -----------------------------------------------------------------

class SkiSummarizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This fucntion ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(SkiSummarizer, self).__init__(*args, **kwargs)

        # The ski file
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

        # Find ski file(s)?
        if self.ski is None: self.find()

        # Show
        if self.config.show: self.show()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SkiSummarizer, self).setup(**kwargs)

        # Get the ski file
        if "ski" in kwargs: self.ski = kwargs.pop("ski")
        elif self.config.ski is not None: self.ski = LabeledSkiFile(self.config.ski)

    # -----------------------------------------------------------------

    def find(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Finding ski files ...")

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

        print("")
        print(fmt.yellow + "GENERAL:" + fmt.reset)
        print("")

        if self.ski.oligochromatic(): print(" - oligochromatic simulation")
        elif self.ski.panchromatic(): print(" - panchromatic simulation")

        print("")
        print(fmt.yellow + "WAVELENGTH GRID:" + fmt.reset)
        print("")

        print(" - grid type: " + self.ski.get_wavelength_grid().tag)
        if self.ski.oligochromatic(): print(" - wavelengths: " + ",".join(self.ski.wavelengths()))
        else:
            print(" - number of wavelengths: " + str(self.ski.nwavelengths()))
            print(" - minimum wavelength: " + represent_quantity(self.ski.min_wavelength))
            print(" - maximum wavelength: " + represent_quantity(self.ski.max_wavelength))

        print("")
        print(fmt.yellow + "INSTRUMENTS:" + fmt.reset)
        print("")

        # Loop over the instruments
        for name in self.ski.get_instrument_names():

            # Get the instrument
            instrument = self.ski.get_instrument_object(name)

            print(fmt.red + fmt.underlined + name + fmt.reset)
            print("")

            properties = instrument.get_properties()
            for label in properties:
                print(" - " + label + ": " + stringify.stringify(properties[label])[1])
            print("")

        print(fmt.yellow + "STELLAR COMPONENTS:" + fmt.reset)
        print("")

        # Loop over the stellar components
        for name in self.ski.get_stellar_component_ids():

            print(fmt.red + fmt.underlined + str(name) + fmt.reset)
            print("")

            properties = self.ski.get_stellar_component_properties(name)
            #print(properties)

            # Get normalization
            lum, wavelength_or_filter = self.ski.get_stellar_component_luminosity(name)

            # Print
            if isinstance(wavelength_or_filter, Filter):
                print(" - filter: " + str(wavelength_or_filter))
                print(" - (neutral) spectral luminosity: " + represent_quantity(lum, scientific=True))
            elif isinstance(wavelength_or_filter, Quantity):
                print(" - wavelength: " + represent_quantity(wavelength_or_filter, scientific=True))
                print(" - spectral luminosity: " + represent_quantity(lum, scientific=True))
            elif wavelength_or_filter is None:
                print(" - bolometric luminosity: " + represent_quantity(lum, scientific=True))
            else: raise ValueError("Unrecognized filter or wavelength: " + str(wavelength_or_filter))

            #print(" - geometry: " + ski.get_stellar_component_geometry(name).tag)

            print(" - geometry: " + " > ".join(self.ski.get_stellar_component_geometry_hierarchy_names(name)))
            print(" - sed: " + self.ski.get_stellar_component_sed(name).tag)

            print("")

        print(fmt.yellow + "DUST COMPONENTS:" + fmt.reset)
        print("")

        # Loop over the dust components
        for name in self.ski.get_dust_component_ids():

            print(fmt.red + fmt.underlined + str(name) + fmt.reset)
            print("")

            mass = self.ski.get_dust_component_mass(name)
            print(" - mass: " + represent_quantity(mass, scientific=True))
            print(" - geometry: " + " > ".join(self.ski.get_dust_component_geometry_hierarchy_names(name)))
            print(" - mix: " + self.ski.get_dust_component_mix(name).tag)
            if self.ski.transient_dust_emissivity: print(" - emissivity: transient")
            if self.ski.grey_body_dust_emissivity: print(" - emissivity: grey body")

            print("")

        print(fmt.yellow + "DUST GRID:" + fmt.reset)
        print("")

        print(" - grid type: " + self.ski.get_dust_grid().tag)

        grid = self.ski.get_dust_grid_object()
        properties = grid.get_properties()
        for label in properties:
            print(" - " + label + ": " + stringify.stringify(properties[label], scientific=True)[1])
        print("")

# -----------------------------------------------------------------
