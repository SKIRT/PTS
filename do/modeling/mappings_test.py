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
import math
import numpy as np

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from pts.core.tools import logging, time
from pts.core.tools import filesystem as fs
from pts.core.basics.configuration import ConfigurationDefinition, ArgumentConfigurationSetter
from pts.core.basics.filter import Filter
from pts.core.simulation.execute import SkirtExec
from pts.modeling.core.mappings import Mappings
from pts.core.data.sed import IntrinsicSED, SED
from pts.core.simulation.skifile import LabeledSkiFile
from pts.core.basics.map import Map
from pts.core.basics.configuration import load_mapping
from pts.modeling.basics.properties import GalaxyProperties
from pts.modeling.fitting.initialization import spectral_factor_hz_to_micron
from pts.core.basics.range import QuantityRange
from pts.core.plot.sed import SEDPlotter
from pts.core.simulation.wavelengthgrid import WavelengthGrid
from pts.core.plot.seds import plotseds
from pts.modeling.basics.projection import GalaxyProjection
from pts.modeling.basics.instruments import FullInstrument, FrameInstrument, SimpleInstrument

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

fuv_filter = Filter.from_string("FUV")
fuv_wavelength = fuv_filter.pivot

# -----------------------------------------------------------------

# Modeling path
modeling_path = fs.cwd()

# Path to the fixed parameters file
fixed_path = fs.join(modeling_path, "fit", "fixed.dat")

# Load the fixed parameters map
fixed = Map()
with open(fixed_path, 'r') as f: load_mapping(f, fixed)

# Get parameters
metallicity = fixed["metallicity"]
compactness = fixed["sfr_compactness"]
pressure = fixed["sfr_pressure"]
covering_factor = fixed["sfr_covering"]

# -----------------------------------------------------------------

test_path = fs.create_directory_in(modeling_path, "test")

# -----------------------------------------------------------------

# Ski template path
ski_template_path = fs.join(modeling_path, "fit", "template.ski")

# Open the ski template
ski = LabeledSkiFile(ski_template_path)

# Get initial values
initial_values = ski.get_labeled_values()

# -----------------------------------------------------------------

# Initial guess for the star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)
sfr = 0.8 # in Msun yr-1

# Create the MAPPINGS template
mappings_initial = Mappings(metallicity, compactness, pressure, covering_factor, sfr)

# Luminosity for FUV (wavelength)
lum = mappings_initial.luminosity_at(fuv_wavelength)
lum2 = mappings_initial.luminosity_for_filter(fuv_filter)

log.info("For SFR of 0.8, the luminosity at FUV wavelength is " + str(lum))
log.info("When convolved over the FUV filter, the luminosity is " + str(lum2))

# -----------------------------------------------------------------

mappings_normalized = Mappings(metallicity, compactness, pressure, covering_factor, 1.0)

# Luminosity for FUV (wavelength)
lum = mappings_normalized.luminosity_at(fuv_wavelength)
lum2 = mappings_normalized.luminosity_for_filter(fuv_filter)

log.info("For SFR of 1.0, the luminosity at FUV wavelength is " + str(lum))
log.info("When convolved over the FUV filter, the luminosity is " + str(lum2))

# -----------------------------------------------------------------

# Load the galaxy properties
prop_path = fs.join(modeling_path, "data", "properties.dat")
prop = GalaxyProperties.from_file(prop_path)

# -----------------------------------------------------------------

denominator = 4. * math.pi * prop.distance**2.

# -----------------------------------------------------------------

sed = mappings_initial.sed
wavelengths = sed.wavelengths(unit="micron")
luminosities = sed.luminosities(unit="W/micron")

fluxes = []

# Convert the spectral luminosities from W / micron to W / Hz
for i in range(len(wavelengths)):

    wavelength = wavelengths[i]
    luminosity = luminosities[i]

    new_luminosity = luminosity.to("W/micron").value / spectral_factor_hz_to_micron(wavelength) * Unit("W/Hz")

    #new_luminosities.append(new_luminosity)

    # Calculate flux density in Jansky
    flux = (new_luminosity / denominator).to("Jy")

    # Add flux density
    fluxes.append(flux)

# Create an SED
flux_sed = SED()

for i in range(len(wavelengths)):

    flux_sed.add_entry(wavelengths[i], fluxes[i])

#plotter = SEDPlotter()

#plotter.add_modeled_sed(flux_sed, "mappings", residuals=False)

min_wavelength = 0.1 #* Unit("micron")
max_wavelength = 1000. #* Unit("micron")
min_flux = 1. #* Unit("Jy")
max_flux = 100. #* Unit("Jy")

min_wavelength = 0.08322738 #np.min(flux_sed.table["Wavelength"])
max_wavelength = 0.776812 # np.max(flux_sed.table["Wavelength"])
min_flux = 1e-8 # np.min(flux_sed.table["Flux"])
max_flux = 1e-4 # np.max(flux_sed.table["Flux"])

#plotter.run(min_wavelength=min_wavelength, max_wavelength=max_wavelength, min_flux=min_flux, max_flux=max_flux)

#exit()

# -----------------------------------------------------------------

sim_wavelengths_range = QuantityRange(0.1, 1000., "micron")
sim_wavelengths = sim_wavelengths_range.log(150, as_list=True)

# Add FUV wavelength
sim_wavelengths.append(fuv_wavelength)

sim_wavelengths = sorted(sim_wavelengths)

wavelengths_array = np.array([wav.to("micron").value for wav in sim_wavelengths])

wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths_array, "micron")

# -----------------------------------------------------------------

# Remove all stellar components except the ionizing stars
ski.remove_stellar_components_except("Ionizing stars")

# Convert to oligochromatic simulation
#ski.to_oligochromatic(sim_wavelengths)

# Remove the dust system
ski.remove_dust_system()

# Set number of packages per wavelength
ski.setpackages(1e5)

# -----------------------------------------------------------------

earth_projection_path = fs.join(modeling_path, "components", "projections", "earth.proj")
earth_projection = GalaxyProjection.from_file(earth_projection_path)
#full_instrument = FullInstrument.from_projection(earth_projection)
simple_instrument = SimpleInstrument.from_projection(earth_projection)

faceon_projection_path = fs.join(modeling_path, "components", "projections", "faceon.proj")
faceon_projection = GalaxyProjection.from_file(faceon_projection_path)

faceon_instrument = FrameInstrument.from_projection(faceon_projection)

ski.remove_all_instruments()
ski.add_instrument("earth", simple_instrument)
ski.add_instrument("faceon", faceon_instrument)

# -----------------------------------------------------------------

# Save the ski template in the test directory
ski_path = fs.join(test_path, "mappings.ski")
ski.saveto(ski_path)

# Output path
out_path = fs.create_directory_in(test_path, "out")

# Input path
in_path = fs.create_directory_in(test_path, "in")

# Copy file i
ionizing_stars_path = fs.join(modeling_path, "maps", "ionizing_stars.fits")
fs.copy_file(ionizing_stars_path, in_path)

# Save the wavelength grid
wavelength_grid_path = fs.join(in_path, "wavelengths.txt")
wavelength_grid.to_skirt_input(wavelength_grid_path)

# -----------------------------------------------------------------

# Perform the SKIRT simulation
simulation = SkirtExec().execute(ski_path, inpath=in_path, outpath=out_path)[0]

# Load the fluxes, convert them to luminosities in erg/s
sedpath = simulation.seddatpaths()[0]
lambdas, lambda_flambda = np.loadtxt(sedpath, usecols=(0, 1), unpack=True)
#lambdav = simulation.convert(lambdav, to_unit='micron', quantity='wavelength')
#lambdaLlambdav = simulation.luminosityforflux(fv, simulation.instrumentdistance(unit='m'), distance_unit='m', luminositydensity_unit='W/micron', wavelength=lambdav) * lambdav * 1e7

# Create the SED
#sed = IntrinsicSED.from_luminosities(lambdav, lambdaLlambdav, luminosity_unit="erg/s")

# Divide by wavelength in micron
flambda = lambda_flambda / lambdas

# Unit is now W / [ m2 * micron]
# convert to W / [m2 * Hz]

fnus = []

for i in range(len(flambda)):

    fnu = flambda[i] / spectral_factor_hz_to_micron(lambdas[i] * Unit("micron")) * Unit("W / (m2 * Hz)")
    fnu = fnu.to("Jy").value
    fnus.append(fnu)

# Create the SED
sed = SED()

for i in range(len(lambdas)):
    sed.add_entry(lambdas[i] * Unit("micron"), fnus[i] * Unit("Jy"))

plotter = SEDPlotter()
plotter.add_modeled_sed(sed, "...", residuals=False)

#min_wavelength = None
#max_wavelength = None
#min_flux = None
#max_flux = None

#plotter.run(min_wavelength=min_wavelength, max_wavelength=max_wavelength, min_flux=min_flux, max_flux=max_flux)

plotseds(simulation)

# -----------------------------------------------------------------
