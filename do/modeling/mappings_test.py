#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.mappings_test This ...

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.filter import Filter
from pts.core.simulation.execute import SkirtExec
from pts.modeling.core.mappings import Mappings
from pts.modeling.core.sed import IntrinsicSED
from pts.core.simulation.skifile import LabeledSkiFile

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# The galaxy name
#definition.add_required("galaxy", "string", "galaxy name")

# Get configuration
setter = ArgumentConfigurationSetter("mappings_test")
config = setter.run(definition)

# -----------------------------------------------------------------

# Determine the log file path
logfile_path = fs.join(fs.cwd(), "log", time.unique_name("log") + ".txt") if config.report else None

# Determine the log level
level = "DEBUG" if config.debug else "INFO"

# Initialize the logger
log = logging.setup_log(level=level, path=logfile_path)
log.start("Starting mappings_test ...")

# -----------------------------------------------------------------

# Modeling path
modeling_path = fs.cwd()

# Ski template path
ski_template_path = fs.join(modeling_path, "fit", "template.ski")

# Open the ski template
ski = LabeledSkiFile(ski_template_path)

values = ski.get_labeled_values()

# -----------------------------------------------------------------

print(values)

metallicity = values["metallicity"]
compactness = values["sfr_compactness"]
pressure = values["sfr_pressure"]
covering_factor = values["sfr_covering"]
sfr = None

print(metallicity)
print(compactness)
print(pressure)
print(covering_factor)
print(sfr)

exit()

# -----------------------------------------------------------------

# Create the MAPPINGS template
mappings = Mappings(metallicity, compactness, pressure, covering_factor, sfr)

# -----------------------------------------------------------------

# Perform the SKIRT simulation
simulation = SkirtExec().execute(ski_path, brief=True, inpath=outpath, outpath=outpath)[0]

# Load the fluxes, convert them to luminosities in erg/s
sedpath = simulation.seddatpaths()[0]
lambdav, fv = np.loadtxt(sedpath, usecols=(0, 1), unpack=True)
lambdav = simulation.convert(lambdav, to_unit='micron', quantity='wavelength')
lambdaLlambdav = simulation.luminosityforflux(fv, simulation.instrumentdistance(unit='m'), distance_unit='m', luminositydensity_unit='W/micron', wavelength=lambdav) * lambdav * 1e7

# Create the SED
sed = IntrinsicSED.from_luminosities(lambdav, lambdaLlambdav, luminosity_unit="erg/s")

# -----------------------------------------------------------------
