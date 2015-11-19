#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.plotresults Plot histograms and scaling relations for a collection of EAGLE SKIRT-runs.
#
# The facilities in this module serve to plot histograms and scaling relations for the results in a set of EAGLE
# SKIRT-runs that have been previously collected in a single data file.

# -----------------------------------------------------------------

import pickle
import os.path
import numpy as np
from scipy.optimize import curve_fit

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

from pts.skirtunits import SkirtUnits
from pts.filter import Filter
import eagle.config as config

# -----------------------------------------------------------------

## This dictionary contains plot axis type definitions for use in combination with EAGLE SKIRT-run result collections.
# The dictionary keys function as axis type identifiers. Each value specifies a plot axis as a tuple containing
# a human-readable label and a callable that returns the axis value for a given SKIRT-run result. The callable is
# executed in a context that has a global variable corresponding to each item in a SKIRT-run info file.
# If an instrument is specified on a higher level in the plotting function, the instrument name is removed from the
# global variable names.
axistypes = {

    # dust grid properties
    'Ngasparts': ( r"$N_{\mathrm{particles},\mathrm{gas}}/10^3$", lambda: setup_particles_cold_gas/1e3 ),
    'Ncells': ( r"$N_\mathrm{cells}/10^6$", lambda: setup_cells_dust_grid/1e6 ),
    'taumax': ( r"$\tau_\mathrm{V,max}$", lambda: setup_optical_depth_maximum ),
    'tau90': ( r"$\tau_\mathrm{V,90}$", lambda: setup_optical_depth_percentile90 ),
    'dusterror': ( r"$\mathrm{1-(M_\mathrm{grid}/M_\mathrm{dust})\,[\%]}$",
        lambda: 100*divide_if_positive(setup_mass_dust-setup_mass_dust_grid,setup_mass_dust) ),

    # intrinsic properties
    'logMstar': ( r"$\log_{10}(M_*)\,[M_\odot]$", lambda: log_if_positive(original_mass_stars) ),
    'logMdust': ( r"$\log_{10}(M_\mathrm{dust})\,[M_\odot]$", lambda: log_if_positive(setup_mass_dust) ),
    'logMdust/Mstar': ( r"$\log_{10}(M_\mathrm{dust}/M_*)$", lambda: log_divide_if_positive(setup_mass_dust,original_mass_stars) ),
    'logMhii': ( r"$\log_{10}(M_\mathrm{HII})\,[M_\odot]$", lambda: log_if_positive(setup_mass_hii_regions) ),
    'fracMhii.fromgas': ( r"$M_{\mathrm{HII},\mathrm{from gas}}/M_{\mathrm{HII},\mathrm{total}}$",
        lambda: divide_if_positive(exported_mass_hii_regions_from_gas,exported_mass_hii_regions) ),
    'logMdust+hii': ( r"$\log_{10}(M_\mathrm{dust}+\frac{1}{100}M_\mathrm{HII})\,[M_\odot]$",
        lambda: log_if_positive(setup_mass_dust+0.01*setup_mass_hii_regions) ),

    'logLtot': ( r"$\log_{10}(L_\mathrm{tot})\,[L_\odot]$", lambda: log_if_positive(setup_luminosity_stars+setup_luminosity_hii_regions) ),
    'logLhii': ( r"$\log_{10}(L_\mathrm{HII})\,[L_\odot]$", lambda: log_if_positive(setup_luminosity_hii_regions) ),
    'Zgas': ( r"$Z_\mathrm{gas}$", lambda: divide_if_positive(setup_mass_metallic_gas,setup_mass_cold_gas) ),
    'fdust': ( r"$f_\mathrm{dust}$", lambda: divide_if_positive(setup_mass_dust,setup_mass_metallic_gas) ),
    'Mgas/Mdust': ( r"$M_\mathrm{gas}/M_\mathrm{dust}$", lambda: divide_if_positive(setup_mass_cold_gas,setup_mass_dust) ),
    'fracMgas': ( r"$M_\mathrm{gas}/(M_*+M_\mathrm{gas})$", lambda: divide_if_positive(setup_mass_cold_gas,original_mass_stars+setup_mass_cold_gas) ),
    'logM/L': ( r"$\log_{10}(M_*/L_\mathrm{tot})\,[M_\odot/L_\odot]$",
        lambda: log_divide_if_positive(original_mass_stars,setup_luminosity_stars+setup_luminosity_hii_regions) ),
    'Mgas/Mhii': ( r"$M_\mathrm{gas}/M_\mathrm{HII}$", lambda: divide_if_positive(setup_mass_cold_gas,setup_mass_hii_regions) ),

    # magnitudes and colors
    'g': ( r"$M_\mathrm{r}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_g ),
    'r': ( r"$M_\mathrm{r}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_r ),
    'i': ( r"$M_\mathrm{i}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_i ),
    'g-r': ( r"$\mathrm{g}-\mathrm{r}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_g - instr_magnitude_sdss_r ),
    'g-i': ( r"$\mathrm{g}-\mathrm{i}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_g - instr_magnitude_sdss_i ),
    'i-H': ( r"$\mathrm{i}-\mathrm{H}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_i - instr_magnitude_2mass_h ),
    'i-H.zib': ( r"$\mathrm{i}-\mathrm{H}\,[\mathrm{mag}]$", lambda: instr_magnitude_sdss_i - instr_magnitude_2mass_h + 1.39 ),
    'NUV-r': ( r"$\mathrm{NUV}-\mathrm{r}\,[\mathrm{mag}]$", lambda: instr_magnitude_galex_nuv - instr_magnitude_sdss_r ),

    # flux densities (Jy)
    'fmax': ( r"$f_{\nu,\mathrm{max}}\,[\mathrm{kJy}]$",
        lambda: np.maximum(instr_xy_fluxdensity_maximum,instr_xz_fluxdensity_maximum,instr_yz_fluxdensity_maximum)/1e3 ),

    # ratios of flux densities using unlimited continuum fluxes (Jy/Jy)
    'logf250/f500': ( r"$\log_{10}(f_{250}/f_{500})$",
        lambda: log_divide_if_positive(instr_fluxdensity_spire_psw_continuum,instr_fluxdensity_spire_plw_continuum) ),
    'logf250/fNUV': ( r"$f_{250}/f_\mathrm{NUV}$",
        lambda: log_divide_if_positive(instr_fluxdensity_spire_psw_continuum,instr_fluxdensity_galex_nuv) ),
    'f250/f350': ( r"$f_{250}/f_{350}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_psw_continuum,instr_fluxdensity_spire_pmw_continuum) ),
    'f250/f500': ( r"$f_{250}/f_{500}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_psw_continuum,instr_fluxdensity_spire_plw_continuum) ),
    'f350/f500': ( r"$f_{350}/f_{500}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_pmw_continuum,instr_fluxdensity_spire_plw_continuum) ),

    # ratios of flux densities using limited continuum fluxes (Jy/Jy)
    'f250/f350.lim': ( r"$f_{250,\mathrm{lim}}/f_{350,\mathrm{lim}}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_psw_limited,instr_fluxdensity_spire_pmw_limited) ),
    'f250/f500.lim': ( r"$f_{250,\mathrm{lim}}/f_{500,\mathrm{lim}}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_psw_limited,instr_fluxdensity_spire_plw_limited) ),
    'f350/f500.lim': ( r"$f_{350,\mathrm{lim}}/f_{500,\mathrm{lim}}$",
        lambda: divide_if_positive(instr_fluxdensity_spire_pmw_limited,instr_fluxdensity_spire_plw_limited) ),

    # luminosities in specific bands
    'logLk': ( r"$\log_{10}(L_\mathrm{K})\,[L_{\odot,\mathrm{K}}]$",
        lambda: log_if_positive(units.luminosityforflux(instr_fluxdensity_2mass_k,setup_distance_instrument,'W/Hz')/LsunK) ),
    'logL250': ( r"$\log_{10}(L_{250})\,[\mathrm{W}/\mathrm{Hz}]$",
        lambda: log_if_positive(units.luminosityforflux(instr_fluxdensity_spire_psw_continuum,setup_distance_instrument,'W/Hz')) ),
    'logLdust': ( r"$\log_{10}(L_{dust})\,[L_\odot]$",
        lambda: log_if_positive(units.luminosityforflux(instr_fluxdensity_uniform_8_1000,setup_distance_instrument,'W/micron',
                                                        wavelength=np.sqrt(8*1000))*(1000-8)/Lsun) ),
    'logM/Lh': ( r"$\log_{10}(M_*/L_\mathrm{H})\,[M_\odot/L_{\odot,\mathrm{H}}]$",
        lambda: log_divide_if_positive(original_mass_stars,units.luminosityforflux(instr_fluxdensity_2mass_h,setup_distance_instrument,'W/Hz')*LsunH) ),

    # other ratios
    'logMdust/f350/D2' : ( r"$\log_{10}(M_\mathrm{dust}/(f_{350}D^2))\,[\mathrm{kg}\,\mathrm{W}^{-1}\,\mathrm{Hz}]$",
        lambda: log_divide_if_positive(setup_mass_dust*Msun,instr_fluxdensity_spire_pmw_continuum*1e-26*(setup_distance_instrument*pc)**2) ),

    # observationally derived mass properties
    'logMstar.zib': ( r"$\log_{10}(M_{*,\mathrm{zib}})\,[M_\odot]$", lambda: log_stellar_mass_as_zibetti() ),
    'logMdust.fit': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit}})\,[M_\odot]$", lambda: log_dust_mass_from_grey_body_fit("continuum") ),
    'logMdust.fit.lim': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit},\mathrm{lim}})\,[M_\odot]$", lambda: log_dust_mass_from_grey_body_fit("limited") ),
    'logMdust.cort': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{cort}})\,[M_\odot]$", lambda: log_dust_mass_as_cortese() ),
    'logMdust.grid': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{grid}})\,[M_\odot]$", lambda: log_dust_mass_from_grid_temperature() ),
    'logMdust.fit/Mstar.zib': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit}}/M_{*,\mathrm{zib}})$",
        lambda: log_dust_mass_from_grey_body_fit() - log_stellar_mass_as_zibetti() ),
    'logMdust.cort/Mstar.zib': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{cort}}/M_{*,\mathrm{zib}})$",
        lambda: log_dust_mass_as_cortese() - log_stellar_mass_as_zibetti() ),

    # dust temperature
    'Tdust.fit': ( r"$T_{\mathrm{dust},\mathrm{fit}}\,[\mathrm{K}]$", lambda: dust_temperature_from_grey_body_fit("continuum") ),
    'Tdust.fit.lim': ( r"$T_{\mathrm{dust},\mathrm{fit},\mathrm{lim}}\,[\mathrm{K}]$", lambda: dust_temperature_from_grey_body_fit("limited") ),
    'Tdust.grid': ( r"$T_{\mathrm{dust},\mathrm{grid}}\,[\mathrm{K}]$", lambda: probe_average_temperature_dust ),
}

# -----------------------------------------------------------------

# some globals used in the axis type definitions
units = SkirtUnits('stellar','frequency')       # distance in pc, flux density in Jy
pc = units.convert(1., from_unit='pc', to_unit='m')
Msun = units.convert(1., from_unit='Msun', to_unit='kg')
Lsun = units.convert(1., from_unit='Lsun', to_unit='W')
LsunK = 10**((34.1-5.19)/2.5)  # solar luminosity in K band expressed in W/Hz  (AB magnitude is 5.19)
LsunH = 10**((34.1-4.71)/2.5)  # solar luminosity in H band expressed in W/Hz  (AB magnitude is 4.71)
Zsun = 0.0127

c = 2.99792458e8;
h = 6.62606957e-34;
k = 1.3806488e-23;

# -----------------------------------------------------------------

# some generic functions used in the axis type definitions

# return log10(x), or smallest other result for x=0
def log_if_positive(x):
    result = np.zeros_like(x)
    positive = x>0
    result[positive] = np.log10(x[positive])
    result[~positive] = result[positive].min()
    return result

# return x/y, or zero for y<=0
def divide_if_positive(x,y):
    result = np.zeros_like(x)
    result[y>0] = x[y>0] / y[y>0]
    return result

# return log10(x/y), or smallest other result for x<=0 or y<=0
def log_divide_if_positive(x,y):
    result = np.zeros_like(x)
    result[y>0] = x[y>0] / y[y>0]
    positive = result>0
    result[positive] = np.log10(result[positive])
    result[~positive] = result[positive].min()
    return result

# -----------------------------------------------------------------

# functions used in the axis type definitions to derive mass properties from observations

# stellar mass according to Zibetti et al 2009, table B1, using color g-i and i-band luminosity
# returns log10 of stellar mass in solar units
def log_stellar_mass_as_zibetti():
    color = instr_magnitude_sdss_g - instr_magnitude_sdss_i   # AB color g - i
    logUpsi = -0.963 + 1.032*color    # Upsilon in solar units (coefficients a_i and b_i for color g-i in table B1)
    logLi = (4.54 - instr_magnitude_sdss_i) / 2.5   # solar AB magnitude in i band is 4.54
    return logUpsi + logLi

# returns black body emission B_nu(nu,T)
def Bnu(nu,T):
    return (2*h*nu**3/c**2) / (np.exp(h*nu/k/T) - 1)

# returns grey body flux in Jy for wavelength(s) specified in micron,
# for given temperature T (K) and dust mass M (Msun), using global distance and beta=2
def greybody(wave, T, M):
    beta = 2
    kappa350 = 0.192                                        # m2/kg (Cortese)
    nu350 = c/350e-6                                        # Hz
    nu = c / (wave * 1e-6)                                  # Hz
    kappa = kappa350 * (nu/nu350)**beta                     # m2/kg
    D = setup_distance_instrument[0] * pc                   # m      !!! assumes that distance is same for all galaxies in the collection
    flux = M * Msun * kappa * Bnu(nu,T) / D**2              # W/m2/Hz
    return flux * 1e26                                      # Jy

# returns the Herschel 160, 250, 350, 500 data points as a tuple (N, waves, fluxes, sigmas)
# with shapes (), (4,) (4, N) (4,) where N is the number of galaxies in the collection
def herschel_data(fluxtype):
    waves = np.array( [ Filter(fs).pivotwavelength() for fs in ("Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW")] )
    fluxstring = '''( instr_fluxdensity_pacs_red_{0}, instr_fluxdensity_spire_psw_{0},
                      instr_fluxdensity_spire_pmw_{0}, instr_fluxdensity_spire_plw_{0} )'''.format(fluxtype)
    fluxes = np.array(eval(fluxstring))
    sigmas = np.array(( 3,1,1,3 ))      # pacs is less sensitive; longer wavelength fluxes are harder to measure
    return fluxes.shape[1], waves, fluxes, sigmas

# returns temperature T (K) for best fit with Herschel 160, 250, 350, 500 data points
# using beta and kappa as set by the greybody() function above
def dust_temperature_from_grey_body_fit(fluxtype="continuum"):
    N, waves, fluxes, sigmas = herschel_data(fluxtype)
    T = np.zeros(N)
    for i in range(N):
        if fluxes[-1,i]>0:
            (T[i],M),dummy = curve_fit(greybody, waves, fluxes[:,i], p0=(17,1e7), sigma=sigmas, absolute_sigma=False, maxfev=5000)
    return T

# dust mass M (Msun) for best fit with Herschel 160, 250, 350, 500 data points
# using beta and kappa as set by the greybody() function above
# returns log10 of dust mass in solar units
def log_dust_mass_from_grey_body_fit(fluxtype="continuum"):
    N, waves, fluxes, sigmas = herschel_data(fluxtype)
    M = np.zeros(N)
    for i in range(N):
        if fluxes[-1,i]>0:
            (T,M[i]),dummy = curve_fit(greybody, waves, fluxes[:,i], p0=(17,1e7), sigma=sigmas, absolute_sigma=False, maxfev=5000)
    return log_if_positive(M)

# dust mass according to Cortese et al 2012, appendix B, using beta=2 for extended sources
# returns log10 of dust mass in solar units
def log_dust_mass_as_cortese():
    x = log_divide_if_positive(instr_fluxdensity_spire_psw_continuum,instr_fluxdensity_spire_plw_continuum)
    logMFD = 16.880 - 1.559*x + 0.160*x**2 - 0.079*x**3 - 0.363*x**4
    logD = np.log10(setup_distance_instrument/1e6)
    logF = log_if_positive(instr_fluxdensity_spire_pmw_continuum)
    logDust = logMFD + 2*logD + logF - 11.32
    #logDust += np.log10( 0.192 / 0.330 )     # kappa at 350 micron assumed in Cortese vs actual for Zubko
    return logDust

# dust mass based on the dust temperature probed in the dust grid and the 350 micron flux
# returns log10 of dust mass in solar units
def log_dust_mass_from_grid_temperature():
    nu350 = c / 350e-6                                      # Hz
    kappa350 = 0.192                                        # m2/kg (Cortese)
    #kappa350 = 0.330                                        # m2/kg  (Zubko)
    f350 = instr_fluxdensity_spire_pmw_continuum * 1e-26    # W/m2
    D = setup_distance_instrument * pc                      # m
    T = probe_average_temperature_dust                      # K
    T[T<1] = 1
    return log_divide_if_positive(f350*D*D, kappa350 * Bnu(nu350,T) * Msun)

# -----------------------------------------------------------------

## An instance of the Collection class represents the contents of a particular collection of EAGLE SKIRT-run results,
# so that the data is ready for plotting. A Collection object has two public properties:
# - \em name: a name identifying the collection
# - \em info: a dictionary of info dictionaries, keyed on instrument name
class Collection:

    ## The constructor loads the contents of the specified collection so that it is ready for plotting.
    # The collection name should \em not include the directory (which is taken from eagle.conf) nor the
    # postfix "_info_collection.dat".
    def __init__(self, collectionname, collectionlabel):
        self.name = collectionlabel
        self.info = {}

        # load the collection
        infilepath = os.path.join(config.collections_path, collectionname+"_info_collection.dat")
        infile = open(infilepath, "r")
        self.info['any'] = pickle.load(infile)
        infile.close()

        # construct filtered dicts for each instrument name
        names = set([ key.split("_")[1] for key in filter(lambda key: key.startswith("instr_"), self.info['any'].keys()) ])
        for name in names:
            self.info[name] = {}
            for key,value in self.info['any'].iteritems():
                if key.startswith("instr_"):
                    segments = key.split("_")
                    if segments[1]==name:
                        segments.pop(1)
                        cleankey = "_".join(segments)
                        self.info[name][cleankey] = value
                else:
                    self.info[name][key] = value

# -----------------------------------------------------------------

## This function produces a one-page pdf file with one or more plots for the specified SKIRT-run collections.
# It expects the following arguments:
# - collections: a sequence of Collection instances
# - plotname: name of the output plot \em not including the directory, nor the filename extension
# - plotdefs: sequence of plot definitions; each item is a dictionary specifying a single plot as described below
# - pagesize: a 2-tuple specifying the size of the complete page in inch; default is A4 format
# - layout: a 2-tuple specifying the number of columns and rows in the layout of the plots; default is 2 by 3
#           (the layout must accomodate all items in the plotdefs sequence)
# - title: title of the plot; default is the value of plotname; specify empty string to omit title
#
# The following table describes the key-value pairs in a plot definition dictionary.
#
#| Key | Presence | Description of Value
#|-----|----------|---------------------
#| x   | required | one of the axis type identifiers in the \em axistypes dictionary
#| y   | required | one of the axis type identifiers in the \em axistypes dictionary, or 'hist' for a histogram
#| instr | optional | the name of the instrument for which to plot data for both x and y axes; defaults to 'any'
#| xinstr | optional | the name of the instrument for the x axis; defaults to the value of \em instr
#| yinstr | optional, used only if y!='hist' | the name of the instrument for the y axis; defaults to the value of \em instr
#| bins | optional, used only if y=='hist' | the number of bins in a histogram; defaults to 10
#| log | optional, used only if y=='hist' | True for histogram on log scale, False for linear scale (the default)
#| xmin | optional | the minimum x value shown; default is smallest x value
#| xmax | optional | the maximum x value shown; default is largest x value
#| ymin | optional | the minimum y value shown; default is smallest y value
#| ymax | optional | the maximum y value shown; default is largest y value
#| diag | optional | if present and True, a dashed diagonal is drawn from (xmin,ymin) to (xmax,ymax)
#
def plotresults(collections, plotname, plotdefs, layout=(2,3), pagesize=(8.268,11.693), title=None):

    # setup the figure
    figure = plt.figure(figsize=pagesize)
    figure.subplots_adjust(wspace=0.15, hspace=0.23,
                           left=0.08, right=0.97, top=0.93, bottom=0.1)
    colors = ('r', 'g', 'b', 'm', 'c', 'y')

    # add figure title
    if title==None: title=plotname
    if len(title)>0: plt.suptitle(title)

    # loop over the plots
    plotindex = 0
    for plotdef in plotdefs:
        # start the appropriate subplot
        plotindex += 1
        ax = plt.subplot(layout[1], layout[0], plotindex)

        # extract the main specifications from the plot definition
        xaxis = plotdef['x']
        yaxis = plotdef['y']
        instr = plotdef.get('instr','any')
        xinstr = plotdef.get('xinstr',instr)
        yinstr = plotdef.get('yinstr',instr)

        # for a regular relation plot...
        if yaxis!='hist':
            # get the specifications from the axis type definitions
            xlabel,xvalue = axistypes[xaxis]
            ylabel,yvalue = axistypes[yaxis]

            # loop over the collections
            for collection,color in zip(collections,colors):
                # setup the x and y values for each axes,
                # after loading the statistics for the appropriate instrument as global variables
                # that can be used in the callables that setup the x and y values for each axis
                globals().update(collection.info[xinstr])
                x = xvalue()
                globals().update(collection.info[yinstr])
                y = yvalue()

                # plot the relation
                plt.scatter(x, y, marker='o', s=10, alpha=0.5, edgecolors='k', linewidths=(1,), facecolors=color)

                # fit a line through the data and plot it
                xmin = plotdef.get('xmin', x.min())
                xmax = plotdef.get('xmax', x.max())
                if xmin>xmax: xmin,xmax = xmax,xmin
                ymin = plotdef.get('ymin', y.min())
                ymax = plotdef.get('ymax', y.max())
                if ymin>ymax: ymin,ymax = ymax,ymin
                mask = (x>=xmin) & (x<=xmax) & (y>=ymin) & (y<=ymax)
                if np.any(mask):
                    rico, y0 = np.polyfit(x[mask], y[mask], 1)
                    x1 = xmin
                    x2 = xmax
                    y1 = y0 + rico*x1
                    y2 = y0 + rico*x2
                    plt.plot([x1,x2], [y1,y2], color=color, label=collection.name)

            # if requested, plot a dashed diagonal
            if plotdef.get('diag', False):
                plt.plot([xmin,xmax], [ymin,ymax], color='k', ls='dashed', alpha=0.7)

        # for a histogram...
        else:
            # get the histogram options
            bins = plotdef.get('bins', 10)
            log = plotdef.get('log', False)

            # get the specifications from the x-axis type definition
            xlabel,xvalue = axistypes[xaxis]

            # setup the y-axis label (force the x-axis instrument to 'any' to avoid changes to the label)
            ylabel = r"$\log_{10}(N_\mathrm{galaxies})$" if log else r"$N_\mathrm{galaxies}$"
            yinstr = 'any'

            # loop over the collections
            for collection,color in zip(collections,colors):
                # setup the x values in the same way as for a regular plot
                globals().update(collection.info[xinstr])
                x = xvalue()

                # the plt.hist() function does not support square axes with mixed linear/log scale;
                # so, compute the histogram
                xmin = plotdef.get('xmin', x.min())
                xmax = plotdef.get('xmax', x.max())
                counts,binedges = np.histogram(x, bins=bins, range=(xmin,xmax))
                if log:
                    counts[counts<1] = 1
                    counts = np.log10(counts)

                # and, plot the histogram
                xpoints = np.zeros(2*len(binedges))
                ypoints = np.zeros(2*len(binedges))
                xpoints[0::2] = binedges
                xpoints[1::2] = binedges
                ypoints[1:-1:2] = counts
                ypoints[2::2] = counts
                plt.plot(xpoints, ypoints, ls='solid', color=color, label=collection.name)

        # set the data limits, if requested
        plt.xlim( xmin=plotdef.get('xmin'), xmax=plotdef.get('xmax') )
        plt.ylim( ymin=plotdef.get('ymin'), ymax=plotdef.get('ymax') )

        # make the plot axes square
        ax.set_aspect(1./ax.get_data_ratio())

        # include instrument names in axis labels if relevant
        if xinstr != 'any': xlabel += r"$\;\triangleright\mathrm{"+xinstr+"}$"
        if yinstr != 'any': ylabel += r"$\;\triangleright\mathrm{"+yinstr+"}$"

        # add axis labels
        plt.xlabel(xlabel, fontsize='medium')
        plt.ylabel(ylabel, fontsize='medium')

        # fine-tune the tick label size
        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize('x-small')

        # add a legend
        plt.legend(loc='best', prop={'size':'small'})
        plt.grid(True)

    # save and close the figure
    plotfilepath = os.path.join(config.plots_path, plotname+".pdf")
    #plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.2)
    plt.savefig(plotfilepath)
    plt.close()
    print "Created results plot file", plotfilepath

# -----------------------------------------------------------------
