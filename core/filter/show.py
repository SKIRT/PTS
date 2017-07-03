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

# Import the relevant PTS classes and modules
from .broad import BroadBandFilter
from .broad import identifiers as broad_identifiers
from .broad import generate_aliases as broad_generate_aliases
from .broad import categorize_filters as categorize_broad
from .narrow import NarrowBandFilter, defines_wavelength, wavelengths_for_spec, wavelength_range_for_spec
from .narrow import identifiers as narrow_identifiers
from .narrow import generate_aliases as narrow_generate_aliases
from .narrow import categorize_filters as categorize_narrow
from ..tools import formatting as fmt
from pts.core.tools import stringify
from pts.core.plot.transmission import TransmissionPlotter
from pts.core.data.transmission import TransmissionCurve
from pts.core.units.stringify import represent_quantity
from ..basics.configurable import Configurable
from ..tools.logging import log
from ..tools import parsing
from .broad import categorized_filters_sorted_labels as sorted_broad_labels
from .narrow import categorized_filters_sorted_labels as sorted_narrow_labels

# -----------------------------------------------------------------

class FilterShower(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FilterShower, self).__init__(*args, **kwargs)

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
        self.setup(**kwargs)

        # 2. Categorize
        self.categorize()

        # 3. Show
        if self.config.show: self.show()

        # 4. Plot
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

        # Inform the user
        log.info("Categorizing ...")

        # Categorize the narrow bands
        if "narrow" in self.config.types: self.categorize_narrow()

        # Categorize the broad bands
        if "broad" in self.config.types: self.categorize_broad()

    # -----------------------------------------------------------------

    def categorize_broad(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing broad band filters ...")

        # If filters should match with a certain string
        if self.config.matching is not None:

            filters = parsing.lazy_broad_band_filter_list(self.config.matching)
            #self.broad[self.config.matching] = filters
            for fltr in filters: self.broad[fltr.spec].append(fltr.spec)

        # ALl broad band filters
        else: self.broad = categorize_broad()

    # -----------------------------------------------------------------

    def categorize_narrow(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Categorizing narrow band filters ...")

        # If filters should match with a certain string
        if self.config.matching is not None:

            filters = parsing.lazy_narrow_band_filter_list(self.config.matching)
            #self.narrow[self.config.matching] = filters
            for fltr in filters: self.narrow[fltr.spec].append(fltr.spec)

        # All narrow band filters
        else: self.narrow = categorize_narrow()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing the filters ...")

        # Show narrow band filters
        if "narrow" in self.config.types: self.show_narrow()

        # Show broad band filters
        if "broad" in self.config.types: self.show_broad()

    # -----------------------------------------------------------------

    def show_narrow(self):

        """
        This function ...
        :return:
        """

        print("")
        print(fmt.bold + fmt.blue + "NARROW BAND FILTERS" + fmt.reset)
        print("")

        #for label in sorted(self.narrow.keys(), key=lambda x: narrow_identifiers.keys().index(self.narrow[x][0])):
        for label in sorted_narrow_labels(self.narrow):

            print(fmt.yellow + fmt.bold + label + fmt.reset)
            if not self.config.short: print("")

            # Loop over the filters
            for spec in self.narrow[label]:

                identifier = narrow_identifiers[spec]

                if defines_wavelength(spec):

                    # Load the filter
                    fltr = NarrowBandFilter(spec)
                    wavelength = fltr.wavelength
                    wavelength_range = None

                else:

                    wavelength = None
                    wavelength_range = wavelength_range_for_spec(spec)

                print("   " + fmt.green + fmt.bold + spec + fmt.reset)

                if not self.config.short:

                    print("")

                    if wavelength is not None: print("     - Wavelength: " + represent_quantity(wavelength))
                    if wavelength_range is not None: print("    - Wavelength range: " + str(wavelength_range))

                    print("")

                if self.config.aliases:

                    print("    " + stringify.stringify_list_fancy(list(narrow_generate_aliases(narrow_identifiers[spec])), 100, " :: ", "    ")[1])
                    print("")

    # -----------------------------------------------------------------

    def show_broad(self):

        """
        This function ...
        :return:
        """

        print("")
        print(fmt.bold + fmt.blue + "BROAD BAND FILTERS" + fmt.reset)
        print("")

        # Loop over the labels
        for label in sorted_broad_labels(self.broad):

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

                    print("    " + stringify.stringify_list_fancy(list(broad_generate_aliases(broad_identifiers[spec])), 100, " :: ", "    ")[1])
                    print("")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the transmission curves and narrow bands ...")

        # Create the transmission plotter
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

                # Add the wavelengths
                for wavelength_label, wavelength in wavelengths_for_spec(spec).items(): plotter.add_wavelength(wavelength, label=wavelength_label)

        # Run the plotter
        plotter.run()

# -----------------------------------------------------------------
