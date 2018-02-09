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

# Import standard modules
import numpy as np
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.dustpedia.core.database import DustPediaDatabase
from pts.dustpedia.core.photometry import DustPediaPhotometry
from pts.dustpedia.core.sample import DustPediaSample
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.magic.services.s4g import get_galaxy_names
from pts.core.basics.log import log
from pts.core.basics.table import SmartTable
from pts.magic.tools.catalogs import get_galaxy_info, get_galaxy_s4g_one_component_info
from pts.core.tools import strings
from pts.magic.tools.colours import calculate_colour
from pts.core.filter.filter import parse_filter
from pts.core.tools import sequences
from pts.core.tools import tables

# -----------------------------------------------------------------

# Basic info
name_name = "Name"
ra_name = "RA"
dec_name = "DEC"
stage_name = "Stage"
type_name = "Type"
velocity_name = "Velocity"
d25_name = "D25"
inclination_name = "Inclination"

old_name_name = "name"
old_ra_name = "ra"
old_dec_name = "dec"
old_stage_name = "stage"
old_type_name = "type"
old_velocity_name = "velocity"
old_d25_name = "d25"
old_inclination_name = "inclination"

# Properties
common_name_name = "Common name"
names_name = "Names"
distance_name = "Distance"
redshift_name = "Redshift"

old_common_name_name = "common_name"
old_names_name = "names"
old_distance_name = "distance"
old_redshift_name = "redshift"

# Morphology
has_s4g_name = "Has S4G"
position_angle_name = "Position angle"
ellipticity_name = "Ellipticity"
sersic_index_name = "Sersic index"
effective_radius_name = "Effective radius"

old_has_s4g_name = "has_s4g"
old_position_angle_name = "position_angle"
old_ellipticity_name = "ellipticity"
old_sersic_index_name = "sersic_index"
old_effective_radius_name = "effective_radius"

# Set luminosities
galex_fuv_name = "GALEX FUV"
galex_nuv_name = "GALEX NUV"
sdss_u_name = "SDSS u"
sdss_g_name = "SDSS g"
sdss_r_name = "SDSS r"
sdss_i_name = "SDSS i"
sdss_z_name = "SDSS z"
twomassj_name = "2MASS J"
twomassh_name = "2MASS H"
twomassk_name = "2MASS Ks"
wise_w1_name = "WISE W1"
wise_w2_name = "WISE W2"
wise_w3_name = "WISE W3"
wise_w4_name = "WISE W4"
iras12_name = "IRAS 12 micron"
iras25_name = "IRAS 25 micron"
iras60_name = "IRAS 60 micron"
iras100_name = "IRAS 100 micron"
i1_name = "IRAC I1"
i2_name = "IRAC I2"
i3_name = "IRAC I3"
i4_name = "IRAC I4"
mips24_name = "MIPS 24 micron"
mips70_name = "MIPS 70 micron"
mips160_name = "MIPS 160 micron"
pblue_name = "Pacs blue"
pgreen_name = "Pacs green"
pred_name = "Pacs red"
psw_name = "SPIRE PSW"
pmw_name = "SPIRE PMW"
plw_name = "SPIRE PLW"
hfi_350_name = "Planck 350 micron"
hfi_550_name = "Planck 550 micron"
hfi_850_name = "Planck 850 micron"
hfi_1380_name = "Planck 1380 micron"
hfi_2100_name = "Planck 2100 micron"
hfi_3000_name = "Planck 3000 micron"

old_galex_fuv_name = "galex_fuv"
old_galex_nuv_name = "galex_nuv"
old_sdss_u_name = "sdss_u"
old_sdss_g_name = "sdss_g"
old_sdss_r_name = "sdss_r"
old_sdss_i_name = "sdss_i"
old_sdss_z_name = "sdss_z"
old_twomassj_name = "2mass_j"
old_twomassh_name = "2mass_h"
old_twomassk_name = "2mass_k"
old_wise_w1_name = "wise_w1"
old_wise_w2_name = "wise_w2"
old_wise_w3_name = "wise_w3"
old_wise_w4_name = "wise_w4"
old_iras12_name = "iras12"
old_iras25_name = "iras25"
old_iras60_name = "iras60"
old_iras100_name = "iras100"
old_i1_name = "i1"
old_i2_name = "i2"
old_i3_name = "i3"
old_i4_name = "i4"
old_mips24_name = "mips24"
old_mips70_name = "mips70"
old_mips160_name = "mips160"
old_pblue_name = "pblue"
old_pgreen_name = "pgreen"
old_pred_name = "pred"
old_psw_name = "psw"
old_pmw_name = "pmw"
old_plw_name = "plw"
old_hfi_350_name = "hfi_350"
old_hfi_550_name = "hfi_550"
old_hfi_850_name = "hfi_850"
old_hfi_1380_name = "hfi_1380"
old_hfi_2100_name = "hfi_2100"
old_hfi_3000_name = "hfi_3000"

# Colours
fuv_nuv_name = "FUV-NUV"
fuv_h_name = "FUV-H"
fuv_j_name = "FUV-J"
fuv_k_name = "FUV-Ks"
fuv_u_name = "FUV-u"
fuv_g_name = "FUV-g"
fuv_r_name = "FUV_r"
fuv_i_name = "FUV_i"
fuv_z_name = "FUV_z"
fuv_mu3_name = "FUV - 3.5 micron"
fuv_mu4_name = "FUV - 4.5 micron"
nuv_h_name = "NUV-H"
nuv_j_name = "NUV-J"
nuv_k_name = "NUV-Ks"
nuv_u_name = "NUV-u"
nuv_g_name = "NUV-g"
nuv_r_name = "NUV-r"
nuv_i_name = "NUV-i"
nuv_z_name = "NUV-z"
nuv_mu3_name = "NUV - 3.5 micron"
nuv_mu4_name = "NUV - 4.5 micron"
mu25_mu70_name = "25 micron - 70 micron"
mu25_mu60_name = "25 micron - 60 micron"
mu60_mu100_name = "60 micron - 100 micron"
mu70_mu100_name = "70 micron - 100 micron"
mu100_mu160_name = "100 micron - 160 micron"
mu160_mu250_name = "160 micron - 250 micron"
mu250_mu350_name = "250 micron - 350 micron"
mu350_mu500_name = "350 micron - 500 micron"
mu350_mu550_name = "350 micron - 550 micron"
mu500_mu850_name = "500 micron - 850 micron"
mu550_mu850_name = "550 micron - 850 micron"
mu850_mu1380_name = "850 micron - 1380 micron"
mu1380_mu2100_name = "1380 micron - 2100 micron"
mu2100_mu3000_name = "2100 micron - 3000 micron"

old_fuv_nuv_name = "fuv_nuv"
old_fuv_h_name = "fuv_h"
old_fuv_j_name = "fuv_j"
old_fuv_k_name = "fuv_k"
old_fuv_u_name = "fuv_u"
old_fuv_g_name = "fuv_g"
old_fuv_r_name = "fuv_r"
old_fuv_i_name = "fuv_i"
old_fuv_z_name = "fuv_z"
old_fuv_mu3_name = "fuv_mu3"
old_fuv_mu4_name = "fuv_mu4"
old_nuv_h_name = "nuv_h"
old_nuv_j_name = "nuv_j"
old_nuv_k_name = "nuv_k"
old_nuv_u_name = "nuv_u"
old_nuv_g_name = "nuv_g"
old_nuv_r_name = "nuv_r"
old_nuv_i_name = "nuv_i"
old_nuv_z_name = "nuv_z"
old_nuv_mu3_name = "nuv_mu3"
old_nuv_mu4_name = "nuv_mu4"
old_mu25_mu70_name = "mu25_mu70"
old_mu25_mu60_name = "mu25_mu60"
old_mu60_mu100_name = "mu60_mu100"
old_mu70_mu100_name = "mu70_mu100"
old_mu100_mu160_name = "mu100_mu160"
old_mu160_mu250_name = "mu160_mu250"
old_mu250_mu350_name = "mu250_mu350"
old_mu350_mu500_name = "mu350_mu500"
old_mu350_mu550_name = "mu350_mu550"
old_mu500_mu850_name = "mu500_mu850"
old_mu550_mu850_name = "mu550_mu850"
old_mu850_mu1380_name = "mu850_mu1380"
old_mu1380_mu2100_name = "mu1380_mu2100"
old_mu2100_mu3000_name = "mu2100_mu3000"

# Models
sfr_name = "SFR"
sfr_error_name = "SFR error"
dust_mass_name = "Dust mass"
dust_mass_error_name = "Dust mass error"
dust_luminosity_name = "Dust luminosity"
dust_luminosity_error_name = "Dust luminosity error"
dust_temperature_name = "Dust temperature"
dust_temperature_error_name = "Dust temperature error"
stellar_mass_name = "Stellar mass"
stellar_mass_error_name = "Stellar mass error"
stellar_luminosity_name = "Stellar luminosity"
stellar_luminosity_error_name = "Stellar luminosity error"

old_sfr_name = "sfr"
old_sfr_error_name = "sfr_error"
old_dust_mass_name = "dust_mass"
old_dust_mass_error_name = "dust_mass_error"
old_dust_luminosity_name = "dust_luminosity"
old_dust_luminosity_error_name = "dust_luminosity_error"
old_dust_temperature_name = "dust_temperature"
old_dust_temperature_error_name = "dust_temperature_error"
old_stellar_mass_name = "stellar_mass"
old_stellar_mass_error_name = "stellar_mass_error"
old_stellar_luminosity_name = "stellar_luminosity"
old_stellar_luminosity_error_name = "stellar_luminosity_error"

# Presence
has_galex_name = "Has GALEX"
has_sdss_name = "Has SDSS"
has_2mass_name = "Has 2MASS"
has_irac_name = "Has IRAC"
has_mips_name = "Has MIPS"
has_wise_name = "Has WISE"
has_pacs_name = "Has Pacs"
has_spire_name = "Has SPIRE"
has_planck_name = "Has Planck"

old_has_galex_name = "has_galex"
old_has_sdss_name = "has_sdss"
old_has_2mass_name = "has_2mass"
old_has_irac_name = "has_irac"
old_has_mips_name = "has_mips"
old_has_wise_name = "has_wise"
old_has_pacs_name = "has_pacs"
old_has_spire_name = "has_spire"
old_has_planck_name = "has_planck"

# -----------------------------------------------------------------

info_names = [name_name, ra_name, dec_name, stage_name, type_name, velocity_name, d25_name,
                inclination_name]

properties_names = [common_name_name, names_name, distance_name, redshift_name]

morphology_names = [has_s4g_name, position_angle_name, ellipticity_name, sersic_index_name, effective_radius_name]

luminosities_names = [galex_fuv_name, galex_nuv_name, sdss_u_name, sdss_g_name, sdss_r_name, sdss_i_name, sdss_z_name,
                twomassj_name, twomassh_name, twomassk_name, wise_w1_name, wise_w2_name, wise_w3_name, wise_w4_name, iras12_name,
                iras25_name, iras60_name, iras100_name, i1_name, i2_name, i3_name, i4_name, mips24_name, mips70_name, mips160_name,
                pblue_name, pgreen_name, pred_name, psw_name, pmw_name, plw_name, hfi_350_name, hfi_550_name, hfi_850_name,
                hfi_1380_name, hfi_2100_name, hfi_3000_name]

colours_names = [fuv_nuv_name, fuv_h_name, fuv_j_name, fuv_k_name, fuv_u_name,
                fuv_g_name, fuv_r_name, fuv_i_name, fuv_z_name, fuv_mu3_name, fuv_mu4_name, nuv_h_name, nuv_j_name, nuv_k_name,
                nuv_u_name, nuv_g_name, nuv_r_name, nuv_i_name, nuv_z_name, nuv_mu3_name, nuv_mu4_name, mu25_mu70_name,
                mu25_mu60_name, mu60_mu100_name, mu70_mu100_name, mu100_mu160_name, mu160_mu250_name, mu250_mu350_name,
                mu350_mu500_name, mu350_mu550_name, mu500_mu850_name, mu550_mu850_name, mu850_mu1380_name, mu1380_mu2100_name,
                mu2100_mu3000_name]

models_names = [sfr_name, sfr_error_name, dust_mass_name, dust_mass_error_name, dust_luminosity_name,
                dust_luminosity_error_name, dust_temperature_name, dust_temperature_error_name, stellar_mass_name,
                stellar_mass_error_name, stellar_luminosity_name, stellar_luminosity_error_name]

presence_names = [has_galex_name, has_sdss_name, has_2mass_name, has_irac_name, has_mips_name, has_wise_name, has_pacs_name, has_spire_name, has_planck_name]

# -----------------------------------------------------------------

column_names = info_names + properties_names + morphology_names + luminosities_names + colours_names + models_names + presence_names

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Maximum number of galaxies
definition.add_optional("ngalaxies", "positive_integer", "max number of galaxies")
definition.add_optional("filename", "string", "galaxy table filename", "galaxies")
definition.add_flag("random", "make random subset of ngalaxies galaxies")

# Galaxy names
definition.add_optional("names", "string_list", "galaxy names")
definition.add_optional("extend", "file_path", "extend existing table file")

# Adjust existing table
definition.add_optional("adjust", "file_path", "adjust existing table file")
definition.add_optional("columns", "string_list", "column names to adjust", choices=column_names)

# Get configuration
config = parse_arguments("make_galaxy_table", definition)

# -----------------------------------------------------------------

# Create the database instance
database = DustPediaDatabase()
#username, password = get_account()
#database.login(username, password)

# -----------------------------------------------------------------

current_table = None

if config.names is not None:

    if config.extend is not None: raise ValueError("Cannot specify names and extend")
    galaxy_names = config.names

else:

    sample = DustPediaSample()
    galaxy_names = sample.get_names()

    # Check from current table
    if config.extend is not None:

        current_table = SmartTable.from_file(config.extend)
        current_galaxy_names = list(current_table["name"])
        galaxy_names = sequences.elements_not_in_other(galaxy_names, current_galaxy_names)

    # Make subset
    if config.ngalaxies is not None:

        if config.random: galaxy_names = sequences.random_subset(galaxy_names, config.ngalaxies, avoid_duplication=True)
        else: galaxy_names = sequences.get_first_values(galaxy_names, config.ngalaxies)

# -----------------------------------------------------------------

ngalaxies = len(galaxy_names)

# -----------------------------------------------------------------

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

# Check
if config.adjust is not None and config.columns is None: raise ValueError("Column names to adjust must be specified")

# Load table to extend
if config.extend is not None: table = current_table
elif config.adjust is not None: table = SmartTable.from_file(config.adjust)
else: table = None

# -----------------------------------------------------------------

def get_fluxes(galaxy_name):

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

    # Return the fluxes
    return fuv_flux, nuv_flux, u_flux, g_flux, r_flux, i_flux, z_flux, j_flux, h_flux, k_flux, w1_flux, w2_flux, \
    w3_flux, w4_flux, iras12_flux, iras25_flux, iras60_flux, iras100_flux, i1_flux, i2_flux, i3_flux, i4_flux, \
    mips24_flux, mips70_flux, mips160_flux, pblue_flux, pgreen_flux, pred_flux, psw_flux, pmw_flux, plw_flux, \
    hfi_350_flux, hfi_550_flux, hfi_850_flux, hfi_1380_flux, hfi_2100_flux, hfi_3000_flux

# -----------------------------------------------------------------

def get_colours(fuv_flux, nuv_flux, u_flux, g_flux, r_flux, i_flux, z_flux, j_flux, h_flux, k_flux, w1_flux, w2_flux,
                iras25_flux, iras60_flux, iras100_flux, i1_flux, i2_flux, mips24_flux, mips70_flux, mips160_flux,
                pblue_flux, pgreen_flux, pred_flux, psw_flux, pmw_flux, plw_flux, hfi_350_flux, hfi_550_flux,
                hfi_850_flux, hfi_1380_flux, hfi_2100_flux, hfi_3000_flux):

    """
    This function ...
    :return:
    """

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
    if iras25_flux is not None and iras60_flux is not None: mu25_mu60 = calculate_colour(iras25_flux, iras60_flux)
    elif mips24_flux is not None and iras60_flux is not None: mu25_mu60 = calculate_colour(mips24_flux, iras60_flux)
    else: mu25_mu60 = None

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

    # Check for NaN
    if fuv_nuv is not None and np.isnan(fuv_nuv): fuv_nuv = None
    if fuv_h is not None and np.isnan(fuv_h): fuv_h = None
    if fuv_j is not None and np.isnan(fuv_j): fuv_j = None
    if fuv_k is not None and np.isnan(fuv_k): fuv_k = None
    if fuv_u is not None and np.isnan(fuv_u): fuv_u = None
    if fuv_g is not None and np.isnan(fuv_g): fuv_g = None
    if fuv_r is not None and np.isnan(fuv_r): fuv_r = None
    if fuv_i is not None and np.isnan(fuv_i): fuv_i = None
    if fuv_z is not None and np.isnan(fuv_z): fuv_z = None
    if fuv_mu3 is not None and np.isnan(fuv_mu3): fuv_mu3 = None
    if fuv_mu4 is not None and np.isnan(fuv_mu4): fuv_mu4 = None
    if nuv_h is not None and np.isnan(nuv_h): nuv_h = None
    if nuv_j is not None and np.isnan(nuv_j): nuv_j = None
    if nuv_k is not None and np.isnan(nuv_k): nuv_k = None
    if nuv_u is not None and np.isnan(nuv_u): nuv_u = None
    if nuv_g is not None and np.isnan(nuv_g): nuv_g = None
    if nuv_r is not None and np.isnan(nuv_r): nuv_r = None
    if nuv_i is not None and np.isnan(nuv_i): nuv_i = None
    if nuv_z is not None and np.isnan(nuv_z): nuv_z = None
    if nuv_mu3 is not None and np.isnan(nuv_mu3): nuv_mu3 = None
    if nuv_mu4 is not None and np.isnan(nuv_mu4): nuv_mu4 = None
    if mu25_mu70 is not None and np.isnan(mu25_mu70): mu25_mu70 = None
    if mu25_mu60 is not None and np.isnan(mu25_mu60): mu25_mu60 = None
    if mu60_mu100 is not None and np.isnan(mu60_mu100): mu60_mu100 = None
    if mu70_mu100 is not None and np.isnan(mu70_mu100): mu70_mu100 = None
    if mu100_mu160 is not None and np.isnan(mu100_mu160): mu100_mu160 = None
    if mu160_mu250 is not None and np.isnan(mu160_mu250): mu160_mu250 = None
    if mu250_mu350 is not None and np.isnan(mu250_mu350): mu250_mu350 = None
    if mu350_mu500 is not None and np.isnan(mu350_mu500): mu350_mu500 = None
    if mu350_mu550 is not None and np.isnan(mu350_mu550): mu350_mu550 = None
    if mu500_mu850 is not None and np.isnan(mu500_mu850): mu500_mu850 = None
    if mu550_mu850 is not None and np.isnan(mu550_mu850): mu550_mu850 = None
    if mu850_mu1380 is not None and np.isnan(mu850_mu1380): mu850_mu1380 = None
    if mu1380_mu2100 is not None and np.isnan(mu1380_mu2100): mu1380_mu2100 = None
    if mu2100_mu3000 is not None and np.isnan(mu2100_mu3000): mu2100_mu3000 = None

    # Return the colours
    return fuv_nuv, fuv_h, fuv_j, fuv_k, fuv_u, fuv_g, fuv_r, fuv_i, fuv_z, fuv_mu3, fuv_mu4, nuv_h, nuv_j, nuv_k, nuv_u, \
    nuv_g, nuv_r, nuv_i, nuv_z, nuv_mu3, nuv_mu4, mu25_mu70, mu25_mu60, mu60_mu100, mu70_mu100, mu100_mu160, \
    mu160_mu250, mu250_mu350, mu350_mu500, mu350_mu550, mu500_mu850, mu550_mu850, mu850_mu1380, mu1380_mu2100, mu2100_mu3000

# -----------------------------------------------------------------

def calculate_luminosities(fuv_flux, nuv_flux, u_flux, g_flux, r_flux, i_flux, z_flux, j_flux, h_flux, k_flux, w1_flux,
                           w2_flux, w3_flux, w4_flux, iras12_flux, iras25_flux, iras60_flux, iras100_flux, i1_flux,
                           i2_flux, i3_flux, i4_flux, mips24_flux, mips70_flux, mips160_flux, pblue_flux, pgreen_flux,
                           pred_flux, psw_flux, pmw_flux, plw_flux, hfi_350_flux, hfi_550_flux, hfi_850_flux,
                           hfi_1380_flux, hfi_2100_flux, hfi_3000_flux, gal_distance):

    """
    This function ...
    :return:
    """

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

    # No distance: luminosities cannot be calculated
    else: fuv_lum = nuv_lum = u_lum = g_lum = r_lum = i_lum = z_lum = j_lum = h_lum = k_lum = w1_lum = w2_lum =w3_lum = w4_lum = iras12_lum = iras25_lum = iras60_lum = iras100_lum = i1_lum = i2_lum = i3_lum = i4_lum = mips24_lum = mips70_lum = mips160_lum = pblue_lum = pgreen_lum = pred_lum = psw_lum = pmw_lum = plw_lum =hfi_350_lum = hfi_550_lum = hfi_850_lum = hfi_1380_lum = hfi_2100_lum = hfi_3000_lum = None

    # Return
    return fuv_lum, nuv_lum, u_lum, g_lum, r_lum, i_lum, z_lum, j_lum, h_lum, k_lum, w1_lum, w2_lum, w3_lum, w4_lum, \
    iras12_lum, iras25_lum, iras60_lum, iras100_lum, i1_lum, i2_lum, i3_lum, i4_lum, mips24_lum, mips70_lum, \
    mips160_lum, pblue_lum, pgreen_lum, pred_lum, psw_lum, pmw_lum, plw_lum, hfi_350_lum, hfi_550_lum, hfi_850_lum, \
    hfi_1380_lum, hfi_2100_lum, hfi_3000_lum

# -----------------------------------------------------------------

def get_model_parameters(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

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

    # There are results from CIGALE fits
    if database.has_cigale_parameters(galaxy_name):

        # Get the parameters
        parameters = database.get_cigale_parameters(galaxy_name)

        # Dust parameters
        dust_mass = parameters["dust_mass"]
        dust_mass_error = parameters["dust_mass_error"]
        dust_luminosity = parameters["dust_luminosity"]
        dust_luminosity_error = parameters["dust_luminosity_error"]
        #fuv_attenuation = parameters["fuv_attenuation"]
        #fuv_attenuation_error = parameters["fuv_attenuation_error"]

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

    # Return the parameters
    return sfr, sfr_error, dust_mass, dust_mass_error, dust_luminosity, dust_luminosity_error, dust_temperature, \
    dust_temperature_error, stellar_mass, stellar_mass_error, stellar_luminosity, stellar_luminosity_error

# -----------------------------------------------------------------

def get_presence(galaxy_name):

    """
    This function ...
    :param galaxy_name:
    :return:
    """

    #observatories = database.get_observatories(galaxy_name)
    #print("OBS:", observatories)

    #filters = database.get_image_filters_per_observatory(galaxy_name)
    #print(filters)

    # Check presence of data
    has_galex = database.has_galex(galaxy_name)
    has_sdss = database.has_sdss(galaxy_name)
    has_2mass = database.has_2mass(galaxy_name)
    has_irac = database.has_irac(galaxy_name)
    has_mips = database.has_mips(galaxy_name)
    has_wise = database.has_wise(galaxy_name)
    has_pacs = database.has_pacs(galaxy_name)
    has_spire = database.has_spire(galaxy_name)
    has_planck = database.has_planck(galaxy_name)

    print("GALEX:", has_galex)
    print("SDSS:", has_sdss)
    print("2MASS:", has_2mass)
    print("IRAC:", has_irac)
    print("MIPS:", has_mips)
    print("WISE:", has_wise)
    print("Pacs:", has_pacs)
    print("SPIRE:", has_spire)
    print("Planck:", has_planck)

    # Return
    return has_galex, has_sdss, has_2mass, has_irac, has_mips, has_wise, has_pacs, has_spire, has_planck

# -----------------------------------------------------------------

# Set adjust flags
if config.adjust is not None:

    adjust_any = True
    adjust_info = sequences.contains_any(info_names, config.columns)
    adjust_properties = sequences.contains_any(properties_names, config.columns)
    adjust_morphology = sequences.contains_any(morphology_names, config.columns)
    adjust_luminosities = sequences.contains_any(luminosities_names, config.columns)
    adjust_colours = sequences.contains_any(colours_names, config.columns)
    adjust_models = sequences.contains_any(models_names, config.columns)
    adjust_presence = sequences.contains_any(presence_names, config.columns)

# Nothing to adjust
else: adjust_any = adjust_info = adjust_properties = adjust_morphology = adjust_luminosities = adjust_colours = adjust_models = adjust_presence = False

# -----------------------------------------------------------------

# Set flags
needs_info = not adjust_any or adjust_info
needs_properties = not adjust_any or adjust_properties
needs_morphology = not adjust_any or adjust_morphology
needs_luminosities = not adjust_any or adjust_luminosities
needs_colours = not adjust_any or adjust_colours
needs_models = not adjust_any or adjust_models
needs_presence = not adjust_any or adjust_presence

# -----------------------------------------------------------------

# Loop over the names
for galaxy_index, galaxy_name in enumerate(galaxy_names):

    # Inform the user
    log.info("Processing galaxy '" + galaxy_name + "' (" + str(galaxy_index+1) + " out of " + str(ngalaxies) + ") ...")

    # Has galaxy on the database?
    if not database.has_galaxy(galaxy_name): continue

    # Set galaxy properties
    properties = OrderedDict()

    # Get row index
    if config.adjust is not None: row_index = tables.find_index(table, galaxy_name, name_name)
    else: row_index = None

    # Get info
    if needs_info:

        # Get info
        info = database.get_galaxy_info(galaxy_name)

        # Check
        if info is None:
            log.warning("Could not find data for galaxy '" + galaxy_name + "' on the database")
            continue

        # Check name
        if info.name != galaxy_name: raise RuntimeError("Something went wrong")

        # Adjust current table?
        if adjust_info:
            for column_name in config.columns:
                if column_name == name_name: table[column_name][row_index] = info.name
                elif column_name == ra_name: table[column_name][row_index] = info.position.ra
                elif column_name == dec_name: table[column_name][row_index] = info.position.dec
                elif column_name == stage_name: table[column_name][row_index] = info.stage
                elif column_name == type_name: table[column_name][row_index] = info.type
                elif column_name == velocity_name: table[column_name][row_index] = info.velocity
                elif column_name == d25_name: table[column_name][row_index] = info.d25
                elif column_name == inclination_name: table[column_name][row_index] = info.inclination

        # Set basic info
        properties[name_name] = info.name
        properties[ra_name] = info.position.ra
        properties[dec_name] = info.position.dec
        properties[stage_name] = info.stage
        properties[type_name] = info.type
        properties[velocity_name] = info.velocity
        properties[d25_name] = info.d25
        properties[inclination_name] = info.inclination

    # No info
    else: info = None

    # Get properties
    if needs_properties:

        # Get properties
        position = info.position if info is not None else None
        gal_name, position, gal_redshift, gal_type, gal_names, gal_distance, gal_inclination, gal_d25, gal_major, gal_minor, gal_pa = get_galaxy_info(galaxy_name, position)

        # Set commmon name
        if gal_name is not None and not strings.startswith_any(gal_name, galaxy_catalog_names): common_name = gal_name
        else: common_name = None

        # Adjust current table?
        if adjust_properties:
            for column_name in config.columns:
                if column_name == common_name_name: table[column_name][row_index] = common_name
                elif column_name == names_name: table[column_name][row_index] = "  ".join(gal_names) if gal_names is not None and len(gal_names) > 0 else None
                elif column_name == distance_name: table[column_name][row_index] = gal_distance
                elif column_name == redshift_name: table[column_name][row_index] = gal_redshift

        # Set properties
        properties[common_name_name] = common_name
        properties[names_name] = "  ".join(gal_names) if gal_names is not None and len(gal_names) > 0 else None
        properties[distance_name] = gal_distance
        properties[redshift_name] = gal_redshift

    # Get morphology
    if needs_morphology:

        # Has S4G decomposition
        has_s4g = galaxy_name in s4g_names
        if has_s4g: s4g_name, position_angle, ellipticity, sersic_index, effective_radius, magnitude = get_galaxy_s4g_one_component_info(galaxy_name)
        else: position_angle = ellipticity = sersic_index = effective_radius = None

        # Adjust current table?
        if adjust_morphology:
            for column_name in config.columns:
                if column_name == has_s4g_name: table[column_name][row_index] = has_s4g
                elif column_name == position_angle_name: table[column_name][row_index] = position_angle
                elif column_name == ellipticity_name: table[column_name][row_index] = ellipticity
                elif column_name == sersic_index_name: table[column_name][row_index] = sersic_index
                elif column_name == effective_radius_name: table[column_name][row_index] = effective_radius

        # Set morphology properties
        properties[has_s4g_name] = has_s4g
        properties[position_angle_name] = position_angle
        properties[ellipticity_name] = ellipticity
        properties[sersic_index_name] = sersic_index
        properties[effective_radius_name] = effective_radius

    # Needs luminosities
    if needs_luminosities or needs_colours:

        # Get the observed fluxes
        fuv_flux, nuv_flux, u_flux, g_flux, r_flux, i_flux, z_flux, j_flux, h_flux, k_flux, w1_flux, w2_flux, \
        w3_flux, w4_flux, iras12_flux, iras25_flux, iras60_flux, iras100_flux, i1_flux, i2_flux, i3_flux, i4_flux, \
        mips24_flux, mips70_flux, mips160_flux, pblue_flux, pgreen_flux, pred_flux, psw_flux, pmw_flux, plw_flux, \
        hfi_350_flux, hfi_550_flux, hfi_850_flux, hfi_1380_flux, hfi_2100_flux, hfi_3000_flux = get_fluxes(galaxy_name)

        # Calculate luminosities
        fuv_lum, nuv_lum, u_lum, g_lum, r_lum, i_lum, z_lum, j_lum, h_lum, k_lum, w1_lum, w2_lum, w3_lum, w4_lum, \
        iras12_lum, iras25_lum, iras60_lum, iras100_lum, i1_lum, i2_lum, i3_lum, i4_lum, mips24_lum, mips70_lum, \
        mips160_lum, pblue_lum, pgreen_lum, pred_lum, psw_lum, pmw_lum, plw_lum, hfi_350_lum, hfi_550_lum, hfi_850_lum, \
        hfi_1380_lum, hfi_2100_lum, hfi_3000_lum = calculate_luminosities(fuv_flux, nuv_flux, u_flux, g_flux, r_flux,
                                                                          i_flux,
                                                                          z_flux, j_flux, h_flux, k_flux, w1_flux,
                                                                          w2_flux, w3_flux, w4_flux, iras12_flux,
                                                                          iras25_flux, iras60_flux, iras100_flux,
                                                                          i1_flux,
                                                                          i2_flux, i3_flux, i4_flux, mips24_flux,
                                                                          mips70_flux, mips160_flux, pblue_flux,
                                                                          pgreen_flux,
                                                                          pred_flux, psw_flux, pmw_flux, plw_flux,
                                                                          hfi_350_flux, hfi_550_flux, hfi_850_flux,
                                                                          hfi_1380_flux, hfi_2100_flux, hfi_3000_flux,
                                                                          gal_distance)

        # Adjust current table?
        if adjust_luminosities:
            for column_name in config.columns:
                if column_name == galex_fuv_name: table[column_name][row_index] = fuv_lum
                elif column_name == galex_nuv_name: table[column_name][row_index] = nuv_lum
                elif column_name == sdss_u_name: table[column_name][row_index] = u_lum
                elif column_name == sdss_g_name: table[column_name][row_index] = g_lum
                elif column_name == sdss_r_name: table[column_name][row_index] = r_lum
                elif column_name == sdss_i_name: table[column_name][row_index] = i_lum
                elif column_name == sdss_z_name: table[column_name][row_index] = z_lum
                elif column_name == twomassj_name: table[column_name][row_index] = j_lum
                elif column_name == twomassh_name: table[column_name][row_index] = h_lum
                elif column_name == twomassk_name: table[column_name][row_index] = k_lum
                elif column_name == wise_w1_name: table[column_name][row_index] = w1_lum
                elif column_name == wise_w2_name: table[column_name][row_index] = w2_lum
                elif column_name == wise_w3_name: table[column_name][row_index] = w3_lum
                elif column_name == wise_w4_name: table[column_name][row_index] = w4_lum
                elif column_name == iras12_name: table[column_name][row_index] = iras12_lum
                elif column_name == iras25_name: table[column_name][row_index] = iras25_lum
                elif column_name == iras60_name: table[column_name][row_index] = iras60_lum
                elif column_name == iras100_name: table[column_name][row_index] = iras100_lum
                elif column_name == i1_name: table[column_name][row_index] = i1_lum
                elif column_name == i2_name: table[column_name][row_index] = i2_lum
                elif column_name == i3_name: table[column_name][row_index] = i3_lum
                elif column_name == i4_name: table[column_name][row_index] = i4_lum
                elif column_name == mips24_name: table[column_name][row_index] = mips24_lum
                elif column_name == mips70_name: table[column_name][row_index] = mips70_lum
                elif column_name == mips160_name: table[column_name][row_index] = mips160_lum
                elif column_name == pblue_name: table[column_name][row_index] = pblue_lum
                elif column_name == pgreen_name: table[column_name][row_index] = pgreen_lum
                elif column_name == pred_name: table[column_name][row_index] = pred_lum
                elif column_name == psw_name: table[column_name][row_index] = psw_lum
                elif column_name == pmw_name: table[column_name][row_index] = pmw_lum
                elif column_name == plw_name: table[column_name][row_index] = plw_lum
                elif column_name == hfi_350_name: table[column_name][row_index] = hfi_350_lum
                elif column_name == hfi_550_name: table[column_name][row_index] = hfi_550_lum
                elif column_name == hfi_850_name: table[column_name][row_index] = hfi_850_lum
                elif column_name == hfi_1380_name: table[column_name][row_index] = hfi_1380_lum
                elif column_name == hfi_2100_name: table[column_name][row_index] = hfi_2100_lum
                elif column_name == hfi_3000_name: table[column_name][row_index] = hfi_3000_lum

        # Set luminosities
        properties[galex_fuv_name] = fuv_lum
        properties[galex_nuv_name] = nuv_lum
        properties[sdss_u_name] = u_lum
        properties[sdss_g_name] = g_lum
        properties[sdss_r_name] = r_lum
        properties[sdss_i_name] = i_lum
        properties[sdss_z_name] = z_lum
        properties[twomassj_name] = j_lum
        properties[twomassh_name] = h_lum
        properties[twomassk_name] = k_lum
        properties[wise_w1_name] = w1_lum
        properties[wise_w2_name] = w2_lum
        properties[wise_w3_name] = w3_lum
        properties[wise_w4_name] = w4_lum
        properties[iras12_name] = iras12_lum
        properties[iras25_name] = iras25_lum
        properties[iras60_name] = iras60_lum
        properties[iras100_name] = iras100_lum
        properties[i1_name] = i1_lum
        properties[i2_name] = i2_lum
        properties[i3_name] = i3_lum
        properties[i4_name] = i4_lum
        properties[mips24_name] = mips24_lum
        properties[mips70_name] = mips70_lum
        properties[mips160_name] = mips160_lum
        properties[pblue_name] = pblue_lum
        properties[pgreen_name] = pgreen_lum
        properties[pred_name] = pred_lum
        properties[psw_name] = psw_lum
        properties[pmw_name] = pmw_lum
        properties[plw_name] = plw_lum
        properties[hfi_350_name] = hfi_350_lum
        properties[hfi_550_name] = hfi_550_lum
        properties[hfi_850_name] = hfi_850_lum
        properties[hfi_1380_name] = hfi_1380_lum
        properties[hfi_2100_name] = hfi_2100_lum
        properties[hfi_3000_name] = hfi_3000_lum

    # Needs colours
    if needs_colours:

        # CALCULATE COLOURS
        fuv_nuv, fuv_h, fuv_j, fuv_k, fuv_u, fuv_g, fuv_r, fuv_i, fuv_z, fuv_mu3, fuv_mu4, nuv_h, nuv_j, nuv_k, nuv_u, \
        nuv_g, nuv_r, nuv_i, nuv_z, nuv_mu3, nuv_mu4, mu25_mu70, mu25_mu60, mu60_mu100, mu70_mu100, mu100_mu160, \
        mu160_mu250, mu250_mu350, mu350_mu500, mu350_mu550, mu500_mu850, mu550_mu850, mu850_mu1380, mu1380_mu2100, \
        mu2100_mu3000 = get_colours(fuv_flux, nuv_flux, u_flux, g_flux, r_flux, i_flux, z_flux, j_flux, h_flux, k_flux, w1_flux, w2_flux,
                                    iras25_flux, iras60_flux, iras100_flux, i1_flux, i2_flux, mips24_flux, mips70_flux, mips160_flux,
                                    pblue_flux, pgreen_flux, pred_flux, psw_flux, pmw_flux, plw_flux, hfi_350_flux, hfi_550_flux,
                                    hfi_850_flux, hfi_1380_flux, hfi_2100_flux, hfi_3000_flux)

        # Adjust current table?
        if adjust_colours:
            for column_name in config.columns:
                if column_name == fuv_nuv_name: table[column_name][row_index] = fuv_nuv
                elif column_name == fuv_h_name: table[column_name][row_index] = fuv_h
                elif column_name == fuv_j_name: table[column_name][row_index] = fuv_j
                elif column_name == fuv_k_name: table[column_name][row_index] = fuv_k
                elif column_name == fuv_u_name: table[column_name][row_index] = fuv_u
                elif column_name == fuv_g_name: table[column_name][row_index] = fuv_g
                elif column_name == fuv_r_name: table[column_name][row_index] = fuv_r
                elif column_name == fuv_i_name: table[column_name][row_index] = fuv_i
                elif column_name == fuv_z_name: table[column_name][row_index] = fuv_z
                elif column_name == fuv_mu3_name: table[column_name][row_index] = fuv_mu3
                elif column_name == fuv_mu4_name: table[column_name][row_index] = fuv_mu4
                elif column_name == nuv_h_name: table[column_name][row_index] = nuv_h
                elif column_name == nuv_j_name: table[column_name][row_index] = nuv_j
                elif column_name == nuv_k_name: table[column_name][row_index] = nuv_k
                elif column_name == nuv_u_name: table[column_name][row_index] = nuv_u
                elif column_name == nuv_g_name: table[column_name][row_index] = nuv_g
                elif column_name == nuv_r_name: table[column_name][row_index] = nuv_r
                elif column_name == nuv_i_name: table[column_name][row_index] = nuv_i
                elif column_name == nuv_z_name: table[column_name][row_index] = nuv_z
                elif column_name == nuv_mu3_name: table[column_name][row_index] = nuv_mu3
                elif column_name == nuv_mu4_name: table[column_name][row_index] = nuv_mu4
                elif column_name == mu25_mu70_name: table[column_name][row_index] = mu25_mu70
                elif column_name == mu25_mu60_name: table[column_name][row_index] = mu25_mu60
                elif column_name == mu60_mu100_name: table[column_name][row_index] = mu60_mu100
                elif column_name == mu70_mu100_name: table[column_name][row_index] = mu70_mu100
                elif column_name == mu100_mu160_name: table[column_name][row_index] = mu100_mu160
                elif column_name == mu160_mu250_name: table[column_name][row_index] = mu160_mu250
                elif column_name == mu250_mu350_name: table[column_name][row_index] = mu250_mu350
                elif column_name == mu350_mu500_name: table[column_name][row_index] = mu350_mu500
                elif column_name == mu350_mu550_name: table[column_name][row_index] = mu350_mu550
                elif column_name == mu500_mu850_name: table[column_name][row_index] = mu500_mu850
                elif column_name == mu550_mu850_name: table[column_name][row_index] = mu550_mu850
                elif column_name == mu850_mu1380_name: table[column_name][row_index] = mu850_mu1380
                elif column_name == mu1380_mu2100_name: table[column_name][row_index] = mu1380_mu2100
                elif column_name == mu2100_mu3000_name: table[column_name][row_index] = mu2100_mu3000

        # Set colours
        properties[fuv_nuv_name] = fuv_nuv
        properties[fuv_h_name] = fuv_h
        properties[fuv_j_name] = fuv_j
        properties[fuv_k_name] = fuv_k
        properties[fuv_u_name] = fuv_u
        properties[fuv_g_name] = fuv_g
        properties[fuv_r_name] = fuv_r
        properties[fuv_i_name] = fuv_i
        properties[fuv_z_name] = fuv_z
        properties[fuv_mu3_name] = fuv_mu3
        properties[fuv_mu4_name] = fuv_mu4
        properties[nuv_h_name] = nuv_h
        properties[nuv_j_name] = nuv_j
        properties[nuv_k_name] = nuv_k
        properties[nuv_u_name] = nuv_u
        properties[nuv_g_name] = nuv_g
        properties[nuv_r_name] = nuv_r
        properties[nuv_i_name] = nuv_i
        properties[nuv_z_name] = nuv_z
        properties[nuv_mu3_name] = nuv_mu3
        properties[nuv_mu4_name] = nuv_mu4
        properties[mu25_mu70_name] = mu25_mu70
        properties[mu25_mu60_name] = mu25_mu60
        properties[mu60_mu100_name] = mu60_mu100
        properties[mu70_mu100_name] = mu70_mu100
        properties[mu100_mu160_name] = mu100_mu160
        properties[mu160_mu250_name] = mu160_mu250
        properties[mu250_mu350_name] = mu250_mu350
        properties[mu350_mu500_name] = mu350_mu500
        properties[mu350_mu550_name] = mu350_mu550
        properties[mu500_mu850_name] = mu500_mu850
        properties[mu550_mu850_name] = mu550_mu850
        properties[mu850_mu1380_name] = mu850_mu1380
        properties[mu1380_mu2100_name] = mu1380_mu2100
        properties[mu2100_mu3000_name] = mu2100_mu3000

    # Needs model parameters
    if needs_models:

        # Get model parameters
        sfr, sfr_error, dust_mass, dust_mass_error, dust_luminosity, dust_luminosity_error, dust_temperature, \
        dust_temperature_error, stellar_mass, stellar_mass_error, stellar_luminosity, stellar_luminosity_error = get_model_parameters(galaxy_name)

        # Adjust current table?
        if adjust_models:
            for column_name in config.columns:
                if column_name == sfr_name: table[column_name][row_index] = sfr
                elif column_name == sfr_error_name: table[column_name][row_index] = sfr_error
                elif column_name == dust_mass_name: table[column_name][row_index] = dust_mass
                elif column_name == dust_mass_error_name: table[column_name][row_index] = dust_mass_error
                elif column_name == dust_luminosity_name: table[column_name][row_index] = dust_luminosity
                elif column_name == dust_luminosity_error_name: table[column_name][row_index] = dust_luminosity_error
                elif column_name == dust_temperature_name: table[column_name][row_index] = dust_temperature
                elif column_name == dust_temperature_error_name: table[column_name][row_index] = dust_temperature_error
                elif column_name == stellar_mass_name: table[column_name][row_index] = stellar_mass
                elif column_name == stellar_mass_error_name: table[column_name][row_index] = stellar_mass_error
                elif column_name == stellar_luminosity_name: table[column_name][row_index] = stellar_luminosity
                elif column_name == stellar_luminosity_error_name: table[column_name][row_index] = stellar_luminosity_error

        # Set model properties
        properties[sfr_name] = sfr
        properties[sfr_error_name] = sfr_error
        properties[dust_mass_name] = dust_mass
        properties[dust_mass_error_name] = dust_mass_error
        properties[dust_luminosity_name] = dust_luminosity
        properties[dust_luminosity_error_name] = dust_luminosity_error
        properties[dust_temperature_name] = dust_temperature
        properties[dust_temperature_error_name] = dust_temperature_error
        properties[stellar_mass_name] = stellar_mass
        properties[stellar_mass_error_name] = stellar_mass_error
        properties[stellar_luminosity_name] = stellar_luminosity
        properties[stellar_luminosity_error_name] = stellar_luminosity_error

    # Needs presence flags
    if needs_presence:

        # Get the presence of data for different observatories
        has_galex, has_sdss, has_2mass, has_irac, has_mips, has_wise, has_pacs, has_spire, has_planck = get_presence(galaxy_name)

        # Adjust current table?
        if adjust_models:
            for column_name in config.columns:
                if column_name == has_galex_name: table[column_name][row_index] = has_galex
                elif column_name == has_sdss_name: table[column_name][row_index] = has_sdss
                elif column_name == has_2mass_name: table[column_name][row_index] = has_2mass
                #elif column_name == has_spitzer_name: table[column_name][row_index] = has_spitzer
                elif column_name == has_irac_name: table[column_name][row_index] = has_irac
                elif column_name == has_mips_name: table[column_name][row_index] = has_mips
                elif column_name == has_wise_name: table[column_name][row_index] = has_wise
                elif column_name == has_pacs_name: table[column_name][row_index] = has_pacs
                elif column_name == has_spire_name: table[column_name][row_index] = has_spire
                elif column_name == has_planck_name: table[column_name][row_index] = has_planck

        # Set flags
        properties[has_galex_name] = has_galex
        properties[has_sdss_name] = has_sdss
        properties[has_2mass_name] = has_2mass
        #properties[has_spitzer_name] = has_spitzer
        properties[has_irac_name] = has_irac
        properties[has_mips_name] = has_mips
        properties[has_wise_name] = has_wise
        properties[has_pacs_name] = has_pacs
        properties[has_spire_name] = has_spire
        properties[has_planck_name] = has_planck

    # Add properties
    galaxies.append(properties)

    # Extend?
    if config.extend is not None: table.add_row_from_dict(properties)

# -----------------------------------------------------------------

# Create table
if config.extend is None and config.adjust is None: table = SmartTable.from_dictionaries(*galaxies, first="name", ignore_none=True)

# -----------------------------------------------------------------

# Save the table
table.saveto_pts(config.filename)

# -----------------------------------------------------------------
