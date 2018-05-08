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
import inspect
import numpy as np

# Import the relevant PTS classes and modules
from pts.core.units.parsing import parse_unit as u
from pts.modeling.core.mappings import Mappings
from pts.core.basics.log import log
from pts.core.tools.stringify import tostr
from pts.modeling.misc.playground import MappingsPlayground
from pts.core.plot.sed import SEDPlotter
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.data.sun import Sun
from pts.core.tools import filesystem as fs
from pts.core.data.sed import SED
from pts.core.simulation.wavelengthgrid import WavelengthGrid

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# Filter
definition.add_required("filter", "filter", "filter in which to calculate the luminosity")

# SFR
definition.add_positional_optional("sfr", "quantity", "star formation rate", "1. Msun/yr", convert_default=True)

# MAPPINGS parameters
definition.add_optional("metallicity", "positive_real", "metallicity", 0.02)
definition.add_optional("compactness", "positive_real", "compactness", 6)
definition.add_optional("pressure", "quantity", "pressure", 1e12 * u("K/m3"))
definition.add_optional("covering_factor", "positive_real", "covering factor", 0.2) # fPDR

# Flags
definition.add_flag("plot", "plot the SEDs", False)
definition.add_flag("skirt", "use SKIRT", True)
definition.add_flag("pts", "use PTS", True)
definition.add_flag("sampled", "use SKIRT luminosities already sampled on a wavelength grid", True)
definition.add_flag("only_skirt", "only SKIRT", False)
definition.add_flag("only_pts", "only PTS", False)
definition.add_flag("only_sampled", "only sampled", False)

# Output path
definition.add_optional("output", "directory_path", "output path")

# Get the configuration
config = parse_arguments("sfr_to_lum", definition, "Convert a SFR to a luminosity in a certain band")

# -----------------------------------------------------------------

# Check
if config.only_skirt:
    config.pts = config.sampled = False
    config.skirt = True
elif config.only_pts:
    config.skirt = config.sampled = False
    config.pts = True
elif config.only_sampled:
    config.skirt = config.pts = False
    config.sampled = True

# -----------------------------------------------------------------

fltr = config.filter
fltr_wavelength = fltr.wavelength

sun = Sun()

print("")
solar_neutral_density = sun.luminosity_for_wavelength(fltr_wavelength, unit="W", density=True)
solar_wavelength_density = sun.luminosity_for_wavelength(fltr_wavelength, unit="W/micron")
log.info("solar in neutral density: " + tostr(solar_neutral_density))
log.info("solar in wavelength density: " + tostr(solar_wavelength_density))
log.info("bolometric solar luminosity: " + tostr(sun.total_luminosity()))
print("")

# -----------------------------------------------------------------

sfr_scalar = config.sfr.to("Msun/yr").value

# -----------------------------------------------------------------

def show_luminosities(sed):

    """
    This function ...
    :param sed:
    :return:
    """

    # Get spectral luminosity density
    lum = sed.photometry_at(fltr_wavelength, unit="W/micron")
    lum2 = sed.photometry_at(fltr_wavelength, unit="W/micron", interpolate=False)

    #
    log.info("Luminosity: " + tostr(lum))
    log.info("No interpolation: " + tostr(lum2))

    # Convert to solar SPECTRAL luminosity DENSITY at wavelength
    lum_spectral_solar = lum.to("W/micron").value / solar_wavelength_density.to("W/micron").value

    # Convert to neutral
    lum_neutral = lum.to("W", density=True, wavelength=fltr_wavelength)
    lum_solar = lum.to("Lsun", density=True, wavelength=fltr_wavelength)

    # Neutral and solar
    log.info("Luminosity in spectral solar units: " + tostr(lum_spectral_solar) + " Lsun_" + fltr.band)
    log.info("Luminosity in neutral units: " + tostr(lum_neutral))
    log.info("Luminosity in solar units: " + tostr(lum_solar))

# -----------------------------------------------------------------

# Calculate from file with sampled luminosituies
if config.sampled:

    # Determine the path to this directory
    this_filepath = fs.absolute_or_in_cwd(inspect.getfile(inspect.currentframe()))
    directory_path = fs.directory_of(this_filepath)

    # Determine the filepath
    filepath = fs.join(directory_path, "MappingsTemplate.dat")

    # Get wavelength grid, calculate luminosities in W
    wg = WavelengthGrid.from_text_file(filepath, "m")
    deltas = wg.deltas(asarray=True, unit="m")
    lums = np.loadtxt(filepath, skiprows=0, unpack=True, usecols=(1))
    lums *= sfr_scalar # CORRECT FOR SFR

    # Construct the SED
    spectrallums = lums / deltas
    sampled_sed = SED.from_arrays(wg.wavelengths(asarray=True, unit="m"), spectrallums, "m", "W/m")

    print("FROM SAMPLED LUMINOSITIES:")
    print("")

    # Show
    show_luminosities(sampled_sed)

    print("")

# Don't calculate
else: sampled_sed = None

# -----------------------------------------------------------------

# Calculate with PTS
if config.pts:

    # Mappings SED
    mappings = Mappings(config.metallicity, config.compactness, config.pressure, config.covering_factor, sfr_scalar)
    sed = mappings.sed

    print("USING PTS:")
    print("")

    # Show
    show_luminosities(sed)

    print("")

# Don't calculate with PTS
else: sed = None

# -----------------------------------------------------------------

# Calculate with SKIRT
if config.skirt:

    logp = np.log10(config.pressure.to("K/cm3").value)

    # Get the SED
    playground = MappingsPlayground()
    sed_skirt = playground.simulate_sed(logp, sfr_scalar, config.metallicity, config.compactness, config.covering_factor, output_path=config.output)

    lum_skirt = sed_skirt.photometry_at(fltr_wavelength, unit="W/micron")
    lum_skirt2 = sed_skirt.photometry_at(fltr_wavelength, unit="W/micron", interpolate=False)

    print("USING SKIRT:")
    print("")

    # Show
    show_luminosities(sed_skirt)
    print("")

    # EXTRA CHECK:
    if config.output is not None:

        # Determine path to SKIRT output file
        luminosities_path = fs.join(config.output, "oneparticle_hii_luminosities.dat")

        # Load the SED
        sed_luminosities_skirt = SED.from_skirt(luminosities_path)

        # Get the deltas
        deltas = sed_luminosities_skirt.wavelength_deltas(unit="micron", asarray=True)

        # Get the spectral luminosities
        luminosities = sed_luminosities_skirt.photometry(unit="W", asarray=True)
        spectral_luminosities = luminosities / deltas

        # Create SED with wavelength luminosity densities
        wavelengths = sed_luminosities_skirt.wavelengths(unit="micron", asarray=True)
        extra_sed = SED.from_arrays(wavelengths, spectral_luminosities, wavelength_unit="micron", photometry_unit="W/micron")

        print("EXTRA:")
        print("")

        # Show
        show_luminosities(extra_sed)

        print("")

# Don't calculate with SKIRT
else: sed_skirt = None

# -----------------------------------------------------------------

# Plot?
if config.plot:

    plotter = SEDPlotter()
    plotter.config.unit = u("W/micron", density=True)
    if sed is not None: plotter.add_sed(sed, "PTS")
    if sed_skirt is not None: plotter.add_sed(sed_skirt, "SKIRT")
    if sampled_sed is not None: plotter.add_sed(sampled_sed, "Sampled")
    plotter.run()

# -----------------------------------------------------------------

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
