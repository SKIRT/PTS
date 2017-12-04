#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.sfr_to_lum Convert a SFR to a luminosity in a certain band.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from pts.modeling.preparation.unitconversion import ab_to_jansky
from pts.core.units.parsing import parse_unit as u
from pts.core.filter.filter import parse_filter
from pts.modeling.core.mappings import Mappings
from pts.core.basics.log import log
from pts.core.tools.stringify import tostr
from pts.core.basics.range import RealRange
from pts.core.prep.smile import get_oligochromatic_template
from pts.modeling.misc.playground import MappingsPlayground
from pts.core.plot.sed import SEDPlotter
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.data.sun import Sun

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()
definition.add_required("filter", "filter", "filter in which to calculate the luminosity")
definition.add_positional_optional("sfr", "quantity", "star formation rate", "1. Msun/yr", convert_default=True)
definition.add_positional_optional("metallicity", "positive_real", "metallicity", 0.02)
definition.add_positional_optional("compactness", "positive_real", "compactness", 6)
definition.add_positional_optional("pressure", "quantity", "pressure", 1e12 * u("K/m3"))
definition.add_positional_optional("covering_factor", "positive_real", "covering factor", 0.2) # fPDR

# Get the arguments
config = parse_arguments("sfr_to_lum", definition, "Convert a SFR to a luminosity in a certain band")

# -----------------------------------------------------------------

fltr = config.filter
fltr_wavelength = fltr.wavelength

sun = Sun()

solar_neutral_density = sun.luminosity_for_wavelength(fltr_wavelength, unit="W", density=True)
solar_wavelength_density = sun.luminosity_for_wavelength(fltr_wavelength, unit="W/micron")
log.info("solar in neutral density: " + tostr(solar_neutral_density))
log.info("solar in wavelength density: " + tostr(solar_wavelength_density))

print("")

# -----------------------------------------------------------------

# Séba:
# 3.53e14 Lsun^FUV MAPPINGS 1 SFR
# 7.7e34W if Lsun(FUV) is 1.425e21 W

# FOR FUV:

seba_lum_solar = 3.53e14 * u("Lsun") # Lsun as neutral luminosity density (lambda * Lambda) put into Lsun or Lsun_FUV?
#seba_lum = 7.7e34 * u("W/micron") # WRONG
#seba_solar = 1.425e21 * u("W") # Where does this come from?
seba_solar = 1.425e21 * u("W/micron")
# I GET 1.57e21 W/micron

seba_lum = seba_lum_solar.value * seba_solar

log.info("seba_lum_solar (FUV): " + tostr(seba_lum_solar))
log.info("seba_solar (FUV): " + tostr(seba_solar))
log.info("seba_lum (FUV): " + tostr(seba_lum))

# -----------------------------------------------------------------

# AB magnitude sun: 16.42

#fuv = parse_filter("FUV")
#fuv_wavelength = fuv.wavelength

#jansky = ab_to_jansky(16.42) * u("Jy")
#print(jansky)
#solar_distance = 149597870700 * u("m")
#wmicron = jansky.to("W/micron", wavelength=fuv.pivot, distance=solar_distance)
#print(wmicron)
#watts = wmicron.to("W", density=True, wavelength=fuv.pivot)
#print(watts)

#SFR: 0.351

#sfr = 0.351 * u("Msun/yr")

# 1.63849785244e+42 W / micron
# 2.51522448305e+41 W

#metallicity = 0.03
#compactness = 6 # logC
#pressure = 1e12 * u("K/m3")
#covering_factor = 0.2 # fPDR

#logp = np.log10(pressure.to("K/cm3").value)
#print(logp)

#minlogp_value = 4.
#maxlogp_value = 8.

# Q_CLASSINFO("MinValue", "1e10 K/m3")
# Q_CLASSINFO("MaxValue", "1e14 K/m3")
# Q_CLASSINFO("Default", "1e11 K/m3")

#min_pressure = 1e10 * u("K/m3")
#max_pressure = 1e14 * u("K/m3")

#minlogp = np.log10(min_pressure.to("K/cm3").value)
#maxlogp = np.log10(max_pressure.to("K/cm3").value)

#print(minlogp_value, minlogp)
#print(maxlogp_value, maxlogp)

#logp: 4 -> 8

# LIMITS:
# Zrel = max(Zrel,0.05);
#Zrel = min(Zrel,2.0-1e-8);
#logC = max(logC,4.0);
#logC = min(logC,6.5-1e-8);
#logp = max(logp,4.0);
#logp = min(logp,8.0-1e-8);

#sfr_scalar = sfr.value

# -----------------------------------------------------------------

# Mappings SED
#mappings = Mappings(metallicity, compactness, pressure, covering_factor, sfr_scalar)

# -----------------------------------------------------------------

# Luminosity for FUV (wavelength)
#lum = mappings.luminosity_at(fuv_wavelength)
#lum2 = mappings.luminosity_at(fuv_wavelength, interpolate=False)
#lum3 = mappings.luminosity_for_filter(fuv)

#log.info("For SFR of " + tostr(sfr) + ", the luminosity at FUV wavelength is " + str(lum))
#log.info("Without interpolation: " + str(lum2))
#log.info("When convolved over the FUV filter, the luminosity is " + str(lum3))

# -----------------------------------------------------------------

#mappings_normalized = Mappings(metallicity, compactness, pressure, covering_factor, 1.0)

# Luminosity for FUV (wavelength)
#lum_one = mappings_normalized.luminosity_at(fuv_wavelength)
#lum_one2 = mappings_normalized.luminosity_at(fuv_wavelength, interpolate=False)
#lum_one3 = mappings_normalized.luminosity_for_filter(fuv)

#log.info("For SFR of 1.0, the luminosity at FUV wavelength is " + str(lum_one))
#log.info("Without interpolation: " + str(lum_one2))
#log.info("When convolved over the FUV filter, the luminosity is " + str(lum_one3))

# -----------------------------------------------------------------

#sim_wavelengths_range = QuantityRange(0.1, 1000., "micron")
#sim_wavelengths = sim_wavelengths_range.log(150, as_list=True)
# Add FUV wavelength
#sim_wavelengths.append(fuv_wavelength)
#sim_wavelengths = sorted(sim_wavelengths)
#wavelengths_array = np.array([wav.to("micron").value for wav in sim_wavelengths])

#sim_wavelengths = RealRange(0.05, 10, "micron").log(150, as_list=True)
#sim_wavelengths.append(fuv_wavelength.to("micron").value)
#wavelength_grid = WavelengthGrid.from_wavelengths(wavelengths_array, "micron", sort=True)

# -----------------------------------------------------------------

print("")

sfr_scalar = config.sfr.to("Msun/yr").value

#ski = get_oligochromatic_template()

# Remove all stellar components except the ionizing stars
#ski.remove_stellar_components_except("Ionizing stars")

# Remove the dust system
#ski.remove_dust_system()

# Set number of packages per wavelength
#ski.setpackages(1e5)

# Perform the SKIRT simulation
#simulation = SkirtExec().execute(ski_path, inpath=in_path, outpath=out_path)[0]

# -----------------------------------------------------------------

print("")

playground = MappingsPlayground()

logp = np.log10(config.pressure.to("K/cm3").value)

sed_skirt = playground.simulate_sed(logp, sfr_scalar, config.metallicity, config.compactness, config.covering_factor)

#plotter = SEDPlotter()
#plotter.config.unit = u("W", density=True)
#plotter.add_sed(sed, "simulation")
##plotter.add_sed(sed2, "simulation2")
#plotter.run()

lum_skirt = sed_skirt.photometry_at(fltr_wavelength, unit="W/micron")
lum_skirt2 = sed_skirt.photometry_at(fltr_wavelength, unit="W/micron", interpolate=False)
#lum3 = sed.luminosity_for_filter(fuv)

print("")
print("2 different methods:")

print("")
print("USING SKIRT:")
print("")

log.info("Luminosity: " + tostr(lum_skirt))
log.info("No interpolation: " + tostr(lum_skirt2))

# Convert to neutral
lum_skirt_neutral = lum_skirt.to("W", density=True, wavelength=fltr_wavelength)
lum_skirt_solar = lum_skirt.to("Lsun", density=True, wavelength=fltr_wavelength)

print("")
log.info("Luminosity in neutral units: " + tostr(lum_skirt_neutral))
log.info("Luminosity in solar units: " + tostr(lum_skirt_solar))

#test = lum_skirt_neutral.value * fltr_wavelength.value
#print(test)

# -----------------------------------------------------------------

# # 1 SFR:
# sed_one = playground.simulate_sed(logp, 1., metallicity, compactness, covering_factor)
#
# lum_skirt = sed_one.photometry_at(fuv_wavelength, unit="W/micron")
# lum_skirt2 = sed_one.photometry_at(fuv_wavelength, unit="W/micron", interpolate=False)
#
# log.info("Luminosity [SFR = 1]: " + tostr(lum_skirt))
# log.info("No interpolation [SFR = 1]: " + tostr(lum_skirt2))
#
# # Convert to neutral
# lum_skirt_neutral = lum_skirt.to("W", density=True, wavelength=fuv_wavelength)
# lum_skirt_solar = lum_skirt.to("Lsun", density=True, wavelength=fuv_wavelength)
#
# log.info("Luminosity in neutral units [SFR = 1]: " + tostr(lum_skirt_neutral))
# log.info("Luminosity in solar units [SFR = 1]: " + tostr(lum_skirt_solar))

# -----------------------------------------------------------------

print("")
print("USING PTS:")
print("")

# Mappings SED
mappings = Mappings(config.metallicity, config.compactness, config.pressure, config.covering_factor, sfr_scalar)
sed = mappings.sed

lum = sed.photometry_at(fltr_wavelength, unit="W/micron")
lum2 = sed.photometry_at(fltr_wavelength, unit="W/micron", interpolate=False)

log.info("Luminosity: " + tostr(lum))
log.info("No interpolation: " + tostr(lum2))

# Convert to neutral
lum_neutral = lum.to("W", density=True, wavelength=fltr_wavelength)
lum_solar = lum.to("Lsun", density=True, wavelength=fltr_wavelength)

print("")
log.info("Luminosity in neutral units: " + tostr(lum_neutral))
log.info("Luminosity in solar units: " + tostr(lum_solar))

# Luminosity for FUV (wavelength)
#lum = mappings.luminosity_at(fuv_wavelength)
#lum2 = mappings.luminosity_at(fuv_wavelength, interpolate=False)
#lum3 = mappings.luminosity_for_filter(fuv)

# -----------------------------------------------------------------
