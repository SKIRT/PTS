#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.plot_galaxies Plot the positions of the galaxies in the DustPedia database.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase, get_account
from pts.dustpedia.core.photometry import DustPediaPhotometry
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.services.s4g import get_galaxy_names, has_galaxy
from pts.core.basics.map import Map
from pts.core.basics.log import log
from pts.core.basics.table import SmartTable
from pts.magic.tools.catalogs import get_galaxy_info, get_galaxy_s4g_one_component_info
from pts.core.tools import strings
from pts.magic.tools.colours import calculate_colour
from pts.core.filter.filter import parse_filter

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Maximum number of galaxies
definition.add_optional("ngalaxies", "positive_integer", "max number of galaxies")

# Get configuration
config = parse_arguments("plot_galaxies", definition)

# -----------------------------------------------------------------

# Create the database instance
database = DustPediaDatabase()
#username, password = get_account()
#database.login(username, password)

# -----------------------------------------------------------------

sample = DustPediaSample()
galaxy_names = sample.get_names()
ngalaxies = len(galaxy_names)

# Photometry
photometry = DustPediaPhotometry()

# -----------------------------------------------------------------

s4g_names = get_galaxy_names()

# -----------------------------------------------------------------

galaxies = []

# -----------------------------------------------------------------

galaxy_catalog_names = ["IC", "UGC", "MCG", "CGCG", "VCC", "NPM1", "ESO", "NGC", "PGC", "UGCA", "IRAS", "KUG", "KIG", "MRC", "MRK", "DUKST", "GIN"]

# -----------------------------------------------------------------

# Initialize filter objects
fuv_filter = parse_filter("GALEX FUV")
nuv_filter = parse_filter("GALEX NUV")
u_filter = parse_filter("SDSS u")
g_filter = parse_filter("SDSS g")
r_filter = parse_filter("SDSS r")
i_filter = parse_filter("SDSS i")
z_filter = parse_filter("SDSS z")
j_filter = parse_filter("2MASS J")
h_filter = parse_filter("2MASS H")
k_filter = parse_filter("2MASS Ks")
w1_filter = parse_filter("WISE W1")
w2_filter = parse_filter("WISE W2")
w3_filter = parse_filter("WISE W3")
w4_filter = parse_filter("WISE W4")
iras12_filter = parse_filter("IRAS 12mu")
iras25_filter = parse_filter("IRAS 25mu")
iras60_filter = parse_filter("IRAS 60mu")
iras100_filter = parse_filter("IRAS 100mu")
i1_filter = parse_filter("IRAC I1")
i2_filter = parse_filter("IRAC I2")
i3_filter = parse_filter("IRAC I3")
i4_filter = parse_filter("IRAC I4")
mips24_filter = parse_filter("MIPS 24mu")
mips70_filter = parse_filter("MIPS 70mu")
mips160_filter = parse_filter("MIPS 160mu")
pblue_filter = parse_filter("Pacs blue")
pgreen_filter = parse_filter("Pacs green")
pred_filter = parse_filter("Pacs red")
psw_filter = parse_filter("SPIRE PSW")
pmw_filter = parse_filter("SPIRE PMW")
plw_filter = parse_filter("SPIRE PLW")
hfi_350_filter = parse_filter("HFI 857")
hfi_550_filter = parse_filter("HFI 545")
hfi_850_filter = parse_filter("HFI 353")
hfi_1380_filter = parse_filter("HFI 217")
hfi_2100_filter = parse_filter("HFI 143")
hfi_3000_filter = parse_filter("HFI 100")

# -----------------------------------------------------------------

# Loop over the names
for galaxy_index, galaxy_name in enumerate(galaxy_names):

    # Inform the user
    log.info("Processing galaxy '" + galaxy_name + "' (" + str(galaxy_index+1) + " out of " + str(ngalaxies) + ") ...")

    # Get info
    info = database.get_galaxy_info(galaxy_name)

    # Get photometry
    observed_sed = photometry.get_sed(galaxy_name, add_iras=True, add_planck=True, exact_name=True)

    # Initialize fluxes
    fuv_flux = None
    nuv_flux = None
    u_flux = None
    g_flux = None
    r_flux = None
    i_flux = None
    z_flux = None
    j_flux = None
    h_flux = None
    k_flux = None
    w1_flux = None
    w2_flux = None
    w3_flux = None
    w4_flux = None
    iras12_flux = None
    iras25_flux = None
    iras60_flux = None
    iras100_flux = None
    i1_flux = None
    i2_flux = None
    i3_flux = None
    i4_flux = None
    mips24_flux = None
    mips70_flux = None
    mips160_flux = None
    pblue_flux = None
    pgreen_flux = None
    pred_flux = None
    psw_flux = None
    pmw_flux = None
    plw_flux = None
    hfi_350_flux = None
    hfi_550_flux = None
    hfi_850_flux = None
    hfi_1380_flux = None
    hfi_2100_flux = None
    hfi_3000_flux = None

    # Loop over the entries in the observed SED
    for index in range(len(observed_sed)):

        # Get the filter
        fltr = observed_sed.get_filter(index)
        #filter_name = str(fltr)

        # Get the flux
        flux = observed_sed.get_photometry(index, add_unit=True)

        # Set correct flux
        if fltr == fuv_filter: fuv_flux = flux
        elif fltr == nuv_filter: nuv_flux = flux
        elif fltr == u_filter: u_flux = flux
        elif fltr == g_filter: g_flux = flux
        elif fltr == r_filter: r_flux = flux
        elif fltr == i_filter: i_flux = flux
        elif fltr == z_filter: z_flux = flux
        elif fltr == j_filter: j_flux = flux
        elif fltr == h_filter: h_flux = flux
        elif fltr == k_filter: k_flux = flux
        elif fltr == w1_filter: w1_flux = flux
        elif fltr == w2_filter: w2_flux = flux
        elif fltr == w3_filter: w3_flux = flux
        elif fltr == w4_filter: w4_flux = flux
        elif fltr == iras12_filter: iras12_flux = flux
        elif fltr == iras25_filter: iras25_flux = flux
        elif fltr == iras60_filter: iras60_flux = flux
        elif fltr == iras100_filter: iras100_flux = flux
        elif fltr == i1_filter: i1_flux = flux
        elif fltr == i2_filter: i2_flux = flux
        elif fltr == i3_filter: i3_flux = flux
        elif fltr == i4_filter: i4_flux = flux
        elif fltr == mips24_filter: mips24_flux = flux
        elif fltr == mips70_filter: mips70_flux = flux
        elif fltr == mips160_filter: mips160_flux = flux
        elif fltr == pblue_filter: pblue_flux = flux
        elif fltr == pgreen_filter: pgreen_flux = flux
        elif fltr == pred_filter: pred_flux = flux
        elif fltr == psw_filter: psw_flux = flux
        elif fltr == pmw_filter: pmw_flux = flux
        elif fltr == plw_filter: plw_flux = flux
        elif fltr == hfi_350_filter: hfi_350_flux = flux
        elif fltr == hfi_550_filter: hfi_550_flux = flux
        elif fltr == hfi_850_filter: hfi_850_flux = flux
        elif fltr == hfi_1380_filter: hfi_1380_flux = flux
        elif fltr == hfi_2100_filter: hfi_2100_flux = flux
        elif fltr == hfi_3000_filter: hfi_3000_flux = flux
        else: log.warning("Unknown filter: " + str(fltr))

    #print(observed_sed)
    #continue

    # CALCULATE COLOURS

    # FUV-NUV
    if fuv_flux is not None and nuv_flux is not None: fuv_nuv = calculate_colour(fuv_flux, nuv_flux)
    else: fuv_nuv = None

    # FUV-H
    if fuv_flux is not None and h_flux is not None: fuv_h = calculate_colour(fuv_flux, h_flux)
    else: fuv_h = None

    # FUV-J
    if fuv_flux is not None and j_flux is not None: fuv_j = calculate_colour(fuv_flux, j_flux)
    else: fuv_j = None

    # FUV-K
    if fuv_flux is not None and k_flux is not None: fuv_k = calculate_colour(fuv_flux, k_flux)
    else: fuv_k = None

    # FUV-u
    if fuv_flux is not None and u_flux is not None: fuv_u = calculate_colour(fuv_flux, u_flux)
    else: fuv_u = None

    # FUV-g
    if fuv_flux is not None and g_flux is not None: fuv_g = calculate_colour(fuv_flux, g_flux)
    else: fuv_g = None

    # FUV-r
    if fuv_flux is not None and r_flux is not None: fuv_r = calculate_colour(fuv_flux, r_flux)
    else: fuv_r = None

    # FUV-i
    if fuv_flux is not None and i_flux is not None: fuv_i = calculate_colour(fuv_flux, i_flux)
    else: fuv_i = None

    # FUV-z
    if fuv_flux is not None and z_flux is not None: fuv_z = calculate_colour(fuv_flux, z_flux)
    else: fuv_z = None

    # FUV-3.5
    if fuv_flux is not None and i1_flux is not None: fuv_mu3 = calculate_colour(fuv_flux, i1_flux)
    elif fuv_flux is not None and w1_flux is not None: fuv_mu3 = calculate_colour(fuv_flux, w1_flux)
    else: fuv_mu3 = None

    # FUV-4.5
    if fuv_flux is not None and i2_flux is not None: fuv_mu4 = calculate_colour(fuv_flux, i2_flux)
    elif fuv_flux is not None and w2_flux is not None: fuv_mu4 = calculate_colour(fuv_flux, w2_flux)
    else: fuv_mu4 = None

    # NUV-H
    if nuv_flux is not None and h_flux is not None: nuv_h = calculate_colour(nuv_flux, h_flux)
    else: nuv_h = None

    # NUV-J
    if nuv_flux is not None and j_flux is not None: nuv_j = calculate_colour(nuv_flux, j_flux)
    else: nuv_j = None

    # NUV-K
    if nuv_flux is not None and k_flux is not None: nuv_k = calculate_colour(nuv_flux, k_flux)
    else: nuv_k = None

    # NUV-u
    if nuv_flux is not None and u_flux is not None: nuv_u = calculate_colour(nuv_flux, u_flux)
    else: nuv_u = None

    # NUV-g
    if nuv_flux is not None and g_flux is not None: nuv_g = calculate_colour(nuv_flux, g_flux)
    else: nuv_g = None

    # NUV-r
    if nuv_flux is not None and r_flux is not None: nuv_r = calculate_colour(nuv_flux, r_flux)
    else: nuv_r = None

    # NUV-i
    if nuv_flux is not None and i_flux is not None: nuv_i = calculate_colour(nuv_flux, i_flux)
    else: nuv_i = None

    # NUV-z
    if nuv_flux is not None and z_flux is not None: nuv_z = calculate_colour(nuv_flux, z_flux)
    else: nuv_z = None

    # NUV-3.5
    if nuv_flux is not None and i1_flux is not None: nuv_mu3 = calculate_colour(nuv_flux, i1_flux)
    elif nuv_flux is not None and w1_flux is not None: nuv_mu3 = calculate_colour(nuv_flux, w1_flux)
    else: nuv_mu3 = None

    # NUV-4.5
    if nuv_flux is not None and i2_flux is not None: nuv_mu4 = calculate_colour(nuv_flux, i2_flux)
    elif nuv_flux is not None and w2_flux is not None: nuv_mu4 = calculate_colour(nuv_flux, w2_flux)
    else: nuv_mu4 = None

    # 25-70
    if mips24_flux is not None and pblue_flux is not None: mu25_mu70 = calculate_colour(mips24_flux, pblue_flux)
    elif mips24_flux is not None and mips70_flux is not None: mu25_mu70 = calculate_colour(mips24_flux, mips70_flux)
    elif iras25_flux is not None and pblue_flux is not None: mu25_mu70 = calculate_colour(iras25_flux, pblue_flux)
    elif iras25_flux is not None and mips70_flux is not None: mu25_mu70 = calculate_colour(iras25_flux, mips70_flux)
    else: mu25_mu70 = None

    # 25-60
    if iras25_flux is not None and iras60_flux is not None: mu25_mu70 = calculate_colour(iras25_flux, iras60_flux)
    elif mips24_flux is not None and iras60_flux is not None: mu25_mu70 = calculate_colour(mips24_flux, iras60_flux)
    else: mu25_mu70 = None

    # 60-100
    if iras60_flux is not None and pgreen_flux is not None: mu60_mu100 = calculate_colour(iras60_flux, pgreen_flux)
    elif iras60_flux is not None and iras100_flux is not None: mu60_mu100 = calculate_colour(iras60_flux, iras100_flux)
    else: mu60_mu100 = None

    # 70-100
    if pblue_flux is not None and pgreen_flux is not None: mu70_mu100 = calculate_colour(pblue_flux, pgreen_flux)
    elif pblue_flux is not None and iras100_flux is not None: mu70_mu100 = calculate_colour(pblue_flux, iras100_flux)
    elif mips70_flux is not None and pgreen_flux is not None: mu70_mu100 = calculate_colour(mips70_flux, pgreen_flux)
    elif mips70_flux is not None and iras100_flux is not None: mu70_mu100 = calculate_colour(mips70_flux, iras100_flux)
    else: mu70_mu100 = None

    # 100-160
    if pgreen_flux is not None and pred_flux is not None: mu100_mu160 = calculate_colour(pgreen_flux, pred_flux)
    elif pgreen_flux is not None and mips160_flux is not None: mu100_mu160 = calculate_colour(pgreen_flux, mips160_flux)
    elif iras100_flux is not None and pred_flux is not None: mu100_mu160 = calculate_colour(iras100_flux, pred_flux)
    elif iras100_flux is not None and mips160_flux is not None: mu100_mu160 = calculate_colour(iras100_flux, mips160_flux)
    else: mu100_mu160 = None

    # 160-250
    if pred_flux is not None and psw_flux is not None: mu160_mu250 = calculate_colour(pred_flux, psw_flux)
    elif mips160_flux is not None and psw_flux is not None: mu160_mu250 = calculate_colour(mips160_flux, psw_flux)
    else: mu160_mu250 = None

    # 250-350
    if psw_flux is not None and pmw_flux is not None: mu250_mu350 = calculate_colour(psw_flux, pmw_flux)
    elif psw_flux is not None and hfi_350_flux is not None: mu250_mu350 = calculate_colour(psw_flux, hfi_350_flux)
    else: mu250_mu350 = None

    # 350-500
    if pmw_flux is not None and plw_flux is not None: mu350_mu500 = calculate_colour(pmw_flux, plw_flux)
    elif hfi_350_flux is not None and plw_flux is not None: mu350_mu500 = calculate_colour(hfi_350_flux, plw_flux)
    else: mu350_mu500 = None

    # 350-550
    if pmw_flux is not None and hfi_550_flux is not None: mu350_mu550 = calculate_colour(pmw_flux, hfi_550_flux)
    elif hfi_350_flux is not None and hfi_550_flux is not None: mu350_mu550 = calculate_colour(hfi_350_flux, hfi_550_flux)
    else: mu350_mu550 = None

    # 500-850
    if plw_flux is not None and hfi_850_flux is not None: mu500_mu850 = calculate_colour(plw_flux, hfi_850_flux)
    else: mu500_mu850 = None

    # 550-850
    if hfi_550_flux is not None and hfi_850_flux is not None: mu550_mu850 = calculate_colour(hfi_550_flux, hfi_850_flux)
    else: mu550_mu850 = None

    # 850-1380
    if hfi_850_flux is not None and hfi_1380_flux is not None: mu850_mu1380 = calculate_colour(hfi_850_flux, hfi_1380_flux)
    else: mu850_mu1380 = None

    # 1380-2100
    if hfi_1380_flux is not None and hfi_2100_flux is not None: mu1380_mu2100 = calculate_colour(hfi_1380_flux, hfi_2100_flux)
    else: mu1380_mu2100 = None

    # 2100-3000
    if hfi_2100_flux is not None and hfi_3000_flux is not None: mu2100_mu3000 = calculate_colour(hfi_2100_flux, hfi_3000_flux)
    else: mu2100_mu3000 = None

    # Get properties
    gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa = get_galaxy_info(galaxy_name, info.position)
    #print(gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa)
    #continue

    common_name = None
    if gal_name is not None and not strings.startswith_any(gal_name, galaxy_catalog_names): common_name = gal_name
    #if "ESO" not in gal_name and "NGC" not in gal_name and "PGC" not in gal_name and not gal_name.startswith("IC") and "UGCA" not in gal_name:
    #    common_name = gal_name

    #print(common_name)
    #print(gal_names)
    #continue
    #print(gal_distance)

    # Convert fluxes to luminosities
    if gal_distance is not None:
        fuv_lum = fuv_flux.to("W/micron", distance=gal_distance, filter=fuv_filter) if fuv_flux is not None else None
        nuv_lum = nuv_flux.to("W/micron", distance=gal_distance, filter=nuv_filter) if nuv_flux is not None else None
        u_lum = u_flux.to("W/micron", distance=gal_distance, filter=u_filter) if u_flux is not None else None
        g_lum = g_flux.to("W/micron", distance=gal_distance, filter=g_filter) if g_flux is not None else None
        r_lum = r_flux.to("W/micron", distance=gal_distance, filter=r_filter) if r_flux is not None else None
        i_lum = i_flux.to("W/micron", distance=gal_distance, filter=i_filter) if i_flux is not None else None
        z_lum = z_flux.to("W/micron", distance=gal_distance, filter=z_filter) if z_flux is not None else None
        j_lum = j_flux.to("W/micron", distance=gal_distance, filter=j_filter) if j_flux is not None else None
        h_lum = h_flux.to("W/micron", distance=gal_distance, filter=h_filter) if h_flux is not None else None
        k_lum = k_flux.to("W/micron", distance=gal_distance, filter=k_filter) if k_flux is not None else None
        w1_lum = w1_flux.to("W/micron", distance=gal_distance, filter=w1_filter) if w1_flux is not None else None
        w2_lum = w2_flux.to("W/micron", distance=gal_distance, filter=w2_filter) if w2_flux is not None else None
        w3_lum = w3_flux.to("W/micron", distance=gal_distance, filter=w3_filter) if w3_flux is not None else None
        w4_lum = w4_flux.to("W/micron", distance=gal_distance, filter=w4_filter) if w4_flux is not None else None
        iras12_lum = iras12_flux.to("W/micron", distance=gal_distance, filter=iras12_filter) if iras12_flux is not None else None
        iras25_lum = iras25_flux.to("W/micron", distance=gal_distance, filter=iras25_filter) if iras25_flux is not None else None
        iras60_lum = iras60_flux.to("W/micron", distance=gal_distance, filter=iras60_filter) if iras60_flux is not None else None
        iras100_lum = iras100_flux.to("W/micron", distance=gal_distance, filter=iras100_filter) if iras100_flux is not None else None
        i1_lum = i1_flux.to("W/micron", distance=gal_distance, filter=i1_filter) if i1_flux is not None else None
        i2_lum = i2_flux.to("W/micron", distance=gal_distance, filter=i2_filter) if i2_flux is not None else None
        i3_lum = i3_flux.to("W/micron", distance=gal_distance, filter=i3_filter) if i3_flux is not None else None
        i4_lum = i4_flux.to("W/micron", distance=gal_distance, filter=i4_filter) if i4_flux is not None else None
        mips24_lum = mips24_flux.to("W/micron", distance=gal_distance, filter=mips24_filter) if mips24_flux is not None else None
        mips70_lum = mips70_flux.to("W/micron", distance=gal_distance, filter=mips70_filter) if mips70_flux is not None else None
        mips160_lum = mips160_flux.to("W/micron", distance=gal_distance, filter=mips160_filter) if mips160_flux is not None else None
        pblue_lum = pblue_flux.to("W/micron", distance=gal_distance, filter=pblue_filter) if pblue_flux is not None else None
        pgreen_lum = pgreen_flux.to("W/micron", distance=gal_distance, filter=pgreen_filter) if pgreen_flux is not None else None
        pred_lum = pred_flux.to("W/micron", distance=gal_distance, filter=pred_filter) if pred_flux is not None else None
        psw_lum = psw_flux.to("W/micron", distance=gal_distance, filter=psw_filter) if psw_flux is not None else None
        pmw_lum = pmw_flux.to("W/micron", distance=gal_distance, filter=pmw_filter) if pmw_flux is not None else None
        plw_lum = plw_flux.to("W/micron", distance=gal_distance, filter=plw_filter) if plw_flux is not None else None
        hfi_350_lum = hfi_350_flux.to("W/micron", distance=gal_distance, filter=hfi_350_filter) if hfi_350_flux is not None else None
        hfi_550_lum = hfi_550_flux.to("W/micron", distance=gal_distance, filter=hfi_550_filter) if hfi_550_flux is not None else None
        hfi_850_lum = hfi_850_flux.to("W/micron", distance=gal_distance, filter=hfi_850_filter) if hfi_850_flux is not None else None
        hfi_1380_lum = hfi_1380_flux.to("W/micron", distance=gal_distance, filter=hfi_1380_filter) if hfi_1380_flux is not None else None
        hfi_2100_lum = hfi_2100_flux.to("W/micron", distance=gal_distance, filter=hfi_2100_filter) if hfi_2100_flux is not None else None
        hfi_3000_lum = hfi_3000_flux.to("W/micron", distance=gal_distance, filter=hfi_3000_filter) if hfi_3000_flux is not None else None
    else: fuv_lum = nuv_lum = u_lum = g_lum = r_lum = i_lum = z_lum = j_lum = h_lum = k_lum = w1_lum = w2_lum =w3_lum = w4_lum = iras12_lum = iras25_lum = iras60_lum = iras100_lum = i1_lum = i2_lum = i3_lum = i4_lum = mips24_lum = mips70_lum = mips160_lum = pblue_lum = pgreen_lum = pred_lum = psw_lum = pmw_lum = plw_lum =hfi_350_lum = hfi_550_lum = hfi_850_lum = hfi_1380_lum = hfi_2100_lum = hfi_3000_lum = None

    # Initialize the parameters
    sfr = None
    sfr_error = None
    dust_mass = None
    dust_mass_error = None
    dust_luminosity = None
    dust_luminosity_error = None
    dust_temperature = None
    dust_temperature_error = None
    stellar_mass = None
    stellar_mass_error = None
    stellar_luminosity = None
    stellar_luminosity_error = None

    position_angle = None
    ellipticity = None
    sersic_index = None
    effective_radius = None

    # Has S4G decomposition
    has_s4g = galaxy_name in s4g_names

    if has_s4g:
        s4g_name, position_angle, ellipticity, sersic_index, effective_radius, magnitude = get_galaxy_s4g_one_component_info(galaxy_name)
        #print(position_angle, ellipticity, sersic_index, effective_radius)

    #continue

    # There are results from CIGALE fits
    if database.has_cigale_parameters(galaxy_name):

        # Get the parameters
        parameters = database.get_cigale_parameters(galaxy_name)

        # Dust parameters
        dust_mass = parameters["dust_mass"]
        dust_mass_error = parameters["dust_mass_error"]
        dust_luminosity = parameters["dust_luminosity"]
        dust_luminosity_error = parameters["dust_luminosity_error"]
        fuv_attenuation = parameters["fuv_attenuation"]
        fuv_attenuation_error = parameters["fuv_attenuation_error"]

        # Stellar parameters
        sfr = parameters["sfr"]
        sfr_error = parameters["sfr_error"]
        stellar_mass = parameters["stellar_mass"]
        stellar_mass_error = parameters["stellar_mass_error"]
        stellar_luminosity = parameters["stellar_luminosity"]
        stellar_luminosity_error = parameters["stellar_luminosity_error"]

    # There are results from black body fits
    elif database.has_dust_black_body_table_parameters(galaxy_name):

        # Get parameters
        parameters = database.get_dust_black_body_table_parameters(galaxy_name)

        # Dust parameters
        dust_mass = parameters["dust_mass"]
        dust_mass_error = parameters["dust_mass_error"]
        dust_luminosity = parameters["dust_luminosity"]
        dust_luminosity_error = parameters["dust_luminosity_error"]
        dust_temperature = parameters["dust_temperature"]
        dust_temperature_error = parameters["dust_temperature_error"]

    else: log.warning("No model parameters available for galaxy '" + galaxy_name + "'")

    # Check presence of data
    has_galex = database.has_galex(galaxy_name)
    has_sdss = database.has_sdss(galaxy_name)
    has_2mass = database.has_2mass(galaxy_name)
    has_spitzer = database.has_spitzer(galaxy_name)
    has_wise = database.has_wise(galaxy_name)
    has_pacs = database.has_pacs(galaxy_name)
    has_spire = database.has_spire(galaxy_name)
    has_planck = database.has_planck(galaxy_name)

    # Set galaxy properties
    properties = Map()

    # Set basic info
    properties.name = info.name
    properties.common_name = common_name
    properties.names = "  ".join(gal_names) if gal_names is not None and len(gal_names) > 0 else None
    properties.ra = info.position.ra
    properties.dec = info.position.dec
    properties.stage = info.stage
    properties.type = info.type
    properties.velocity = info.velocity
    properties.d25 = info.d25
    properties.inclination = info.inclination

    # Set other properties
    properties.distance = gal_distance
    properties.redshift = gal_redshift
    properties.has_s4g = has_s4g
    properties.position_angle = position_angle
    properties.ellipticity = ellipticity
    properties.sersic_index = sersic_index
    properties.effective_radius = effective_radius

    # Set fluxes
    properties.galex_fuv = fuv_lum
    properties.galex_nuv = nuv_lum
    properties.sdss_u = u_lum
    properties.sdss_g = g_lum
    properties.sdss_r = r_lum
    properties.sdss_i = i_lum
    properties.sdss_z = z_lum
    properties["2mass_j"] = j_lum
    properties["2mass_h"] = h_lum
    properties["2mass_k"] = k_lum
    properties.wise_w1 = w1_lum
    properties.wise_w2 = w2_lum
    properties.wise_w3 = w3_lum
    properties.wise_w4 = w4_lum
    properties.iras_12 = iras12_lum
    properties.iras25 = iras25_lum
    properties.iras60 = iras60_lum
    properties.iras100 = iras100_lum
    properties.i1 = i1_lum
    properties.i2 = i2_lum
    properties.i3 = i3_lum
    properties.i4 = i4_lum
    properties.mips24 = mips24_lum
    properties.mips70 = mips70_lum
    properties.mips160 = mips160_lum
    properties.pblue_flux = pblue_lum
    properties.pgreen_flux = pgreen_lum
    properties.pred_flux = pred_lum
    properties.psw_flux = psw_lum
    properties.pmw_flux = pmw_lum
    properties.plw_flux = plw_lum
    properties.hfi_350_flux = hfi_350_lum
    properties.hfi_550_flux = hfi_550_lum
    properties.hfi_850_flux = hfi_850_lum
    properties.hfi_1380_flux = hfi_1380_lum
    properties.hfi_2100_flux = hfi_2100_lum
    properties.hfi_3000_flux = hfi_3000_lum

    properties.fuv_nuv = fuv_nuv
    properties.fuv_h = fuv_h
    properties.fuv_j = fuv_j
    properties.fuv_k = fuv_k
    properties.fuv_u = fuv_u
    properties.fuv_g = fuv_g
    properties.fuv_r = fuv_r
    properties.fuv_i = fuv_i
    properties.fuv_z = fuv_z
    properties.fuv_mu3 = fuv_mu3
    properties.fuv_mu4 = fuv_mu4
    properties.nuv_h = nuv_h
    properties.nuv_j = nuv_j
    properties.nuv_k = nuv_k
    properties.nuv_u = nuv_u
    properties.nuv_g = nuv_g
    properties.nuv_r = nuv_r
    properties.nuv_i = nuv_i
    properties.nuv_z = nuv_z
    properties.nuv_mu3 = nuv_mu3
    properties.nuv_mu4 = nuv_mu4
    properties.mu25_mu70 = mu25_mu70
    properties.mu25_mu70 = mu25_mu70
    properties.mu60_mu100 = mu60_mu100
    properties.mu70_mu100 = mu70_mu100
    properties.mu100_mu160 = mu100_mu160
    properties.mu160_mu250 = mu160_mu250
    properties.mu250_mu350 = mu250_mu350
    properties.mu350_mu500 = mu350_mu500
    properties.mu350_mu550 = mu350_mu550
    properties.mu500_mu850 = mu500_mu850
    properties.mu550_mu850 = mu550_mu850
    properties.mu850_mu1380 = mu850_mu1380
    properties.mu1380_mu2100 = mu1380_mu2100
    properties.mu2100_mu3000 = mu2100_mu3000

    # Set model properties
    properties.sfr = sfr
    properties.sfr_error = sfr_error
    properties.dust_mass = dust_mass
    properties.dust_mass_error = dust_mass_error
    properties.dust_luminosity = dust_luminosity
    properties.dust_luminosity_error = dust_luminosity_error
    properties.dust_temperature = dust_temperature
    properties.dust_temperature_error = dust_temperature_error
    properties.stellar_mass = stellar_mass
    properties.stellar_mass_error = stellar_mass_error
    properties.stellar_luminosity = stellar_luminosity
    properties.stellar_luminosity_error = stellar_luminosity_error

    # Set flags
    properties.has_galex = has_galex
    properties.has_sdss = has_sdss
    properties.has_2mass = has_2mass
    properties.has_spitzer = has_spitzer
    properties.has_wise = has_wise
    properties.has_pacs = has_pacs
    properties.has_spire = has_spire
    properties.has_planck = has_planck

    # Add properties
    galaxies.append(properties)

    # Limit of number of galaxies
    if config.ngalaxies is not None and len(galaxies) == config.ngalaxies: break

# -----------------------------------------------------------------

# Create table
table = SmartTable.from_dictionaries(*galaxies, first="name", ignore_none=True)

# -----------------------------------------------------------------

# Save the table
filename = "galaxies"
table.saveto_pts(filename)

# -----------------------------------------------------------------
