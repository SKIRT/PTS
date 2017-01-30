#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************


# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

# Import standard modules
from collections import defaultdict
from textwrap import wrap

# Import astronomical modules
from astropy.units import spectral

# Import the relevant PTS classes and modules
from .broad import BroadBandFilter
from .broad import identifiers as broad_identifiers
from .broad import generate_aliases as broad_generate_aliases
from .narrow import NarrowBandFilter, defines_wavelength
from .narrow import identifiers as narrow_identifiers
from .narrow import generate_aliases as narrow_generate_aliases
from ..tools import formatting as fmt
from pts.core.tools import stringify
from pts.core.plot.transmission import TransmissionPlotter
from pts.core.data.transmission import TransmissionCurve
from pts.core.basics.quantity import represent_quantity
from ..basics.configurable import Configurable
from ..tools.logging import log
from ..tools import parsing

# -----------------------------------------------------------------

class FilterShower(Configurable):

    """
    This class ...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        """

        # Call the constructor of the base class
        super(FilterShower, self).__init__(config)

        # Categorized broad and narrow band filters
        self.broad = defaultdict(list)
        self.narrow = defaultdict(list)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. call the setup fucntion
        self.setup()

        # 2. Categorize
        self.categorize()

        # 2. Show
        if self.config.show: self.show()

        # 3. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base clas
        super(FilterShower, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def categorize(self):

        """
        This function ...
        :return:
        """

        # Categorize the broad bands
        self.categorize_broad()

        # Categorize the narrow bands
        self.categorize_narrow()

    # -----------------------------------------------------------------

    def categorize_broad(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing broad band filters ...")

        # Categorize
        for spec in broad_identifiers:

            identifier = broad_identifiers[spec]
            if "instruments" in identifier:
                if "observatories" in identifier:
                    self.broad[identifier.observatories[0] + " " + identifier.instruments[0]].append(spec)
                else: self.broad[identifier.instruments[0]].append(spec)
            elif "observatories" in identifier:
                self.broad[identifier.observatories[0]].append(spec)
            elif "system" in identifier:
                self.broad[identifier.system].append(spec)
            else: self.broad[spec].append(spec)

    # -----------------------------------------------------------------

    def categorize_narrow(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing narrow band filters ...")

        # Categorize
        for spec in narrow_identifiers:

            #identifiers = narrow_identifiers[spec]
            self.narrow[spec].append(spec)

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the filters ...")

        print("")
        print("NARROW BAND FILTERS")
        print("")

        self.show_narrow()

        print("")
        print("BROAD BAND FILTERS")
        print("")

        self.show_broad()

    # -----------------------------------------------------------------

    def show_narrow(self):

        """
        This function ...
        :return:
        """

        for label in sorted(self.narrow.keys(), key=lambda x: narrow_identifiers.keys().index(self.narrow[x][0])):

            print(fmt.yellow + fmt.bold + label + fmt.reset)
            if not self.config.short: print("")

            # Loop over the filters
            for spec in self.narrow[label]:

                identifier = narrow_identifiers[spec]

                if defines_wavelength(identifier):

                    # Load the filter
                    fltr = NarrowBandFilter(spec)
                    wavelength = fltr.wavelength
                    wavelength_range = None

                else:

                    wavelength = None
                    if "wavelength_range" in identifier: wavelength_range = parsing.quantity_range(identifier.wavelength_range)
                    elif "frequency_range" in identifier: wavelength_range = parsing.quantity_range(identifier.frequency_range).to("micron", equivalencies=spectral())
                    else: raise ValueError("Wavelength range or frequency range is not defined for filter spec " + spec)

                print("   " + fmt.green + fmt.bold + spec + fmt.reset)

                if not self.config.short:

                    print("")

                    if wavelength is not None: print("     - Wavelength: " + represent_quantity(wavelength))
                    if wavelength_range is not None: print("    - Wavelength range: " + str(wavelength_range))

                    print("")

                if self.config.aliases:

                    print("   " + "\n     ".join(wrap(stringify.stringify(list(narrow_generate_aliases(narrow_identifiers[spec])))[1].replace(",", " :: "), 100)))
                    print("")

    # -----------------------------------------------------------------

    def show_broad(self):

        """
        This function ...
        :return:
        """

        # Loop over the labels
        for label in sorted(self.broad.keys(), key=lambda x: broad_identifiers.keys().index(self.broad[x][0])):

            print(fmt.yellow + fmt.bold + label + fmt.reset)
            if not self.config.short: print("")

            # Loop over the filters
            for spec in self.broad[label]:

                # Load the filter
                fltr = BroadBandFilter(spec)

                print("   " + fmt.green + fmt.bold + spec + fmt.reset)
                if not self.config.short:

                    print("")

                    print("    - Minimum wavelength: " + represent_quantity(fltr.min))
                    print("    - Maximum wavelength: " + represent_quantity(fltr.max))
                    print("    - Mean wavelength: " + represent_quantity(fltr.mean))
                    if fltr.effective is not None: print("    - Effective wavelength: " + represent_quantity(fltr.effective))
                    print("    - Pivot wavelength: " + represent_quantity(fltr.pivot))
                    if fltr.effective is not None: print("    - Effective bandwidth: " + represent_quantity(fltr.bandwidth))

                    print("")

                if self.config.aliases:

                    print("     " + "\n     ".join(wrap(stringify.stringify(list(broad_generate_aliases(broad_identifiers[spec])))[1].replace(",", " :: "), 100)))
                    print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        plotter = TransmissionPlotter()

        # Loop over the filters
        for label in self.broad:

            # Loop over the filters
            for spec in self.broad[label]:

                fltr = BroadBandFilter(spec)

                curve = TransmissionCurve.from_filter(fltr)
                plotter.add_transmission_curve(curve, spec)

        # Loop over the narrow band filters
        for label in self.narrow:

            # Loop over the filters
            for spec in self.narrow[label]:

                fltr = NarrowBandFilter(spec)

                # Loop over the filters
                plotter.add_wavelength(fltr.wavelength)

        # Run the plotter
        plotter.run()

# -----------------------------------------------------------------
