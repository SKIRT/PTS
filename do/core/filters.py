#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.core.filters Show available filters, their aliases and properties.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict
from textwrap import wrap

# Import the relevant PTS classes and modules
from pts.core.basics.filter import Filter, identifiers, generate_aliases
from pts.core.tools import formatting as fmt
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.tools import stringify
from pts.core.plot.transmission import TransmissionPlotter
from pts.core.data.transmission import TransmissionCurve
from pts.core.basics.quantity import represent_quantity

# -----------------------------------------------------------------

# Create config
definition = ConfigurationDefinition()
definition.add_flag("short", "short: only show the filter names", letter="s")
definition.add_flag("aliases", "show aliases", letter="a")
definition.add_flag("categorize", "categorize per instrument/observatory/system", True, letter="c")
definition.add_flag("plot", "plot the filter transmission curves", letter="p")
setter = ArgumentConfigurationSetter("filters")
config = setter.run(definition)

# -----------------------------------------------------------------

categorized = defaultdict(list)

# Categorize
for spec in identifiers:

    identifier = identifiers[spec]
    if "instruments" in identifier:
        if "observatories" in identifier: categorized[identifier.observatories[0] + " " + identifier.instruments[0]].append(spec)
        else: categorized[identifier.instruments[0]].append(spec)
    elif "observatories" in identifier: categorized[identifier.observatories[0]].append(spec)
    elif "system" in identifier: categorized[identifier.system].append(spec)
    else: categorized[spec].append(spec)

print("")

if config.plot: plotter = TransmissionPlotter()
else: plotter = None

# Loop over the labels
for label in sorted(categorized.keys(), key=lambda x: identifiers.keys().index(categorized[x][0])):

    print(fmt.yellow + fmt.bold + label + fmt.reset)
    if not config.short: print("")

    # Loop over the filters
    for spec in categorized[label]:

        # Load the filter
        fltr = Filter(spec)

        print("   " + fmt.green + fmt.bold + spec + fmt.reset)
        if not config.short:

            print("")

            print("    - Minimum wavelength: " + represent_quantity(fltr.min))
            print("    - Maximum wavelength: " + represent_quantity(fltr.max))
            print("    - Mean wavelength: " + represent_quantity(fltr.mean))
            if fltr.effective is not None: print("    - Effective wavelength: " + represent_quantity(fltr.effective))
            print("    - Pivot wavelength: " + represent_quantity(fltr.pivot))
            if fltr.effective is not None: print("    - Effective bandwidth: " + represent_quantity(fltr.bandwidth))

            print("")

        if config.aliases:

            print("     " + "\n     ".join(wrap(stringify.stringify(list(generate_aliases(identifiers[spec])))[1].replace(",", " :: "), 100)))
            print("")

        if config.plot:
            curve = TransmissionCurve.from_filter(fltr)
            plotter.add_transmission_curve(curve, spec)

# Plot
if config.plot: plotter.run()

# -----------------------------------------------------------------
