#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.eagle.plotresults Plot histograms and scaling relations for a collection of EAGLE SKIRT-runs.
#
# The facilities in this module serve to plot histograms and scaling relations for the results in a set of EAGLE
# SKIRT-runs that have been previously collected in a single data file.

# -----------------------------------------------------------------

import os.path
import numpy as np

import matplotlib.pyplot as plt

from .collections import CollectionData, log_if_positive, divide_if_positive, log_divide_if_positive
from ..core.simulation.units import SkirtUnits
from ..core.filter.broad import BroadBandFilter
from . import config

# -----------------------------------------------------------------

# some globals used in the axis type definitions
units = SkirtUnits('stellar','frequency')       # distance in pc, flux density in Jy
pc = units.convert(1., from_unit='pc', to_unit='m')
Msun = units.convert(1., from_unit='Msun', to_unit='kg')
Lsun = units.convert(1., from_unit='Lsun', to_unit='W')
LsunK = 10**((34.1-5.19)/2.5)  # solar luminosity in K band expressed in W/Hz  (AB magnitude is 5.19)
LsunH = 10**((34.1-4.71)/2.5)  # solar luminosity in H band expressed in W/Hz  (AB magnitude is 4.71)
c = 2.99792458e8               # speed of light in m/s

# global to hold instance of CollectionData, assigned in plotresults() and used in axis type definitions
cd = None

# -----------------------------------------------------------------

## This dictionary contains plot axis type definitions for use in combination with EAGLE SKIRT-run result collections.
# The dictionary keys function as axis type identifiers. Each value specifies a plot axis as a tuple containing
# a human-readable label and a callable that returns the axis value for a given SKIRT-run result. The callable is
# executed in a context that has a global variable "cd" holding the relevant CollectionData instance.
# If an instrument is specified on a higher level in the plotting function, the instrument name is stripped from the
# property names in the CollectionData instance.
axistypes = {

    # dust grid properties
    'Ngasparts': ( r"$N_{\mathrm{particles},\mathrm{gas}}/10^3$", lambda: cd.setup_particles_cold_gas/1e3 ),
    'Ncells': ( r"$N_\mathrm{cells}/10^6$", lambda: cd.setup_cells_dust_grid/1e6 ),
    'taumax': ( r"$\tau_\mathrm{V,max}$", lambda: cd.setup_optical_depth_maximum ),
    'tau90': ( r"$\tau_\mathrm{V,90}$", lambda: cd.setup_optical_depth_percentile90 ),
    'dusterror': ( r"$\mathrm{|1-(M_\mathrm{grid}/M_\mathrm{dust})|\,[\%]}$",
        lambda: 100*divide_if_positive(np.abs(cd.setup_mass_dust-cd.setup_mass_dust_grid),cd.setup_mass_dust) ),

    # intrinsic properties
    'logMstar': ( r"$\log_{10}(M_*)\,[M_\odot]$", lambda: log_if_positive(cd.original_mass_stars) ),
    'logMdust': ( r"$\log_{10}(M_\mathrm{dust})\,[M_\odot]$", lambda: log_if_positive(cd.exported_mass_dust) ),
    'logMdust/Mstar': ( r"$\log_{10}(M_\mathrm{dust}/M_*)$", lambda: log_divide_if_positive(cd.exported_mass_dust,cd.original_mass_stars) ),
    'logMhii': ( r"$\log_{10}(M_\mathrm{HII})\,[M_\odot]$", lambda: log_if_positive(cd.exported_mass_hii_regions) ),
    'fracMhii.fromgas': ( r"$M_{\mathrm{HII},\mathrm{from gas}}/M_{\mathrm{HII},\mathrm{total}}$",
        lambda: divide_if_positive(cd.exported_mass_hii_regions_from_gas,cd.exported_mass_hii_regions) ),
    'logLtot': ( r"$\log_{10}(L_\mathrm{tot})\,[L_\odot]$", lambda: log_if_positive(cd.setup_luminosity_stars+cd.setup_luminosity_hii_regions) ),
    'logLhii': ( r"$\log_{10}(L_\mathrm{HII})\,[L_\odot]$", lambda: log_if_positive(cd.setup_luminosity_hii_regions) ),
    'Zgas': ( r"$Z_\mathrm{gas}$", lambda: divide_if_positive(cd.exported_mass_metallic_gas-cd.exported_mass_negative_metallic_gas,
                                                              cd.exported_mass_cold_gas-cd.exported_mass_negative_cold_gas) ),
    'fdust': ( r"$f_\mathrm{dust}$", lambda: divide_if_positive(cd.exported_mass_dust,cd.exported_mass_metallic_gas) ),
    'Mgas/Mdust': ( r"$M_\mathrm{gas}/M_\mathrm{dust}$", lambda: divide_if_positive(cd.exported_mass_cold_gas,cd.exported_mass_dust) ),
    'fracMgas': ( r"$M_\mathrm{gas}/(M_*+M_\mathrm{gas})$",
        lambda: divide_if_positive(cd.exported_mass_cold_gas,cd.original_mass_stars+cd.exported_mass_cold_gas) ),
    'logM/L': ( r"$\log_{10}(M_*/L_\mathrm{tot})\,[M_\odot/L_\odot]$",
        lambda: log_divide_if_positive(cd.original_mass_stars,cd.setup_luminosity_stars+cd.setup_luminosity_hii_regions) ),
    'Mgas/Mhii': ( r"$M_\mathrm{gas}/M_\mathrm{HII}$", lambda: divide_if_positive(cd.exported_mass_cold_gas,cd.exported_mass_hii_regions) ),

    # magnitudes and colors
    'g': ( r"$M_\mathrm{r}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_g ),
    'r': ( r"$M_\mathrm{r}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_r ),
    'i': ( r"$M_\mathrm{i}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_i ),
    'g-r': ( r"$\mathrm{g}-\mathrm{r}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_g - cd.instr_magnitude_sdss_r ),
    'g-i': ( r"$\mathrm{g}-\mathrm{i}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_g - cd.instr_magnitude_sdss_i ),
    'i-H': ( r"$\mathrm{i}-\mathrm{H}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_i - cd.instr_magnitude_2mass_h ),
    'i-H.zib': ( r"$\mathrm{i}-\mathrm{H}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_sdss_i - cd.instr_magnitude_2mass_h + 1.39 ),
    'NUV-r': ( r"$\mathrm{NUV}-\mathrm{r}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_galex_nuv - cd.instr_magnitude_sdss_r ),
    'K': ( r"$M_\mathrm{K}\,[\mathrm{mag}]$", lambda: cd.instr_magnitude_2mass_k ),

    # flux densities (Jy)
    'fmax': ( r"$f_{\nu,\mathrm{max}}\,[\mathrm{kJy}]$",
        lambda: np.maximum(cd.instr_xy_fluxdensity_maximum,cd.instr_xz_fluxdensity_maximum,cd.instr_yz_fluxdensity_maximum)/1e3 ),

    # ratios of flux densities (Jy/Jy)
    'logf250/f500': ( r"$\log_{10}(f_{250}/f_{500})$",
        lambda: log_divide_if_positive(cd.instr_fluxdensity_spire_psw_limited,cd.instr_fluxdensity_spire_plw_limited) ),
    'logf250/fNUV': ( r"$f_{250}/f_\mathrm{NUV}$",
        lambda: log_divide_if_positive(cd.instr_fluxdensity_spire_psw_limited,cd.instr_fluxdensity_galex_nuv) ),
    'f250/f350': ( r"$f_{250}/f_{350}$",
        lambda: divide_if_positive(cd.instr_fluxdensity_spire_psw_limited,cd.instr_fluxdensity_spire_pmw_limited) ),
    'f250/f500': ( r"$f_{250}/f_{500}$",
        lambda: divide_if_positive(cd.instr_fluxdensity_spire_psw_limited,cd.instr_fluxdensity_spire_plw_limited) ),
    'f350/f500': ( r"$f_{350}/f_{500}$",
        lambda: divide_if_positive(cd.instr_fluxdensity_spire_pmw_limited,cd.instr_fluxdensity_spire_plw_limited) ),

    # luminosities in specific bands
    'logLk': ( r"$\log_{10}(L_\mathrm{K})\,[L_{\odot,\mathrm{K}}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_2mass_k,cd.setup_distance_instrument,'W/Hz')/LsunK) ),
    'logL24': ( r"$\log_{10}(L_{24})\,[\mathrm{W}/\mathrm{Hz}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_mips_24,cd.setup_distance_instrument,'W/Hz')) ),
    'logL250': ( r"$\log_{10}(L_{250})\,[\mathrm{W}/\mathrm{Hz}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_spire_psw_limited,cd.setup_distance_instrument,'W/Hz')) ),
    'logLdust': ( r"$\log_{10}(L_{dust})\,[L_\odot]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_uniform_8_1000,cd.setup_distance_instrument,'W/micron',
                                                        wavelength=np.sqrt(8*1000))*(1000-8)/Lsun) ),
    'logM/Lh': ( r"$\log_{10}(M_*/L_\mathrm{H})\,[M_\odot/L_{\odot,\mathrm{H}}]$",
        lambda: log_divide_if_positive(cd.original_mass_stars,units.luminosityforflux(cd.instr_fluxdensity_2mass_h,cd.setup_distance_instrument,'W/Hz')*LsunH) ),

    # other ratios
    'logMdust/f350/D2' : ( r"$\log_{10}(M_\mathrm{dust}/(f_{350}D^2))\,[\mathrm{kg}\,\mathrm{W}^{-1}\,\mathrm{Hz}]$",
        lambda: log_divide_if_positive(cd.exported_mass_dust*Msun,cd.instr_fluxdensity_spire_pmw_limited*1e-26*(cd.setup_distance_instrument*pc)**2) ),

    # observationally derived mass properties
    'logMstar.zib': ( r"$\log_{10}(M_{*,\mathrm{zib}})\,[M_\odot]$", lambda: cd.log_stellar_mass_as_zibetti() ),
    'logMdust.fit.unlim': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit},\mathrm{unlim}})\,[M_\odot]$",
        lambda: cd.log_dust_mass_from_grey_body_fit("continuum") ),
    'logMdust.fit': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit}})\,[M_\odot]$",
        lambda: cd.log_dust_mass_from_grey_body_fit("limited") ),
    'logMdust.hii.fit': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{HII},\mathrm{fit}})\,[M_\odot]$",
        lambda: log_if_positive(cd.dust_temperature_and_mass_from_grey_body_fit("limited")[1] * cd.dust_fraction_in_hii_regions()) ),
    'logMdust.other.fit': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{other},\mathrm{fit}})\,[M_\odot]$",
        lambda: log_if_positive(cd.dust_temperature_and_mass_from_grey_body_fit("limited")[1] * (1 - cd.dust_fraction_in_hii_regions())) ),
    'Mdust.hii.fit/Mdust.fit': ( r"$M_{\mathrm{dust},\mathrm{hii},\mathrm{fit}}/M_{\mathrm{dust},\mathrm{fit}}$",
        lambda: cd.dust_fraction_in_hii_regions() ),
    'Mdust.hii.fit/Mhii': ( r"$M_{\mathrm{dust},\mathrm{hii},\mathrm{fit}}/M_{\mathrm{HII}}$",
            lambda: divide_if_positive(cd.dust_temperature_and_mass_from_grey_body_fit("limited")[1] * cd.dust_fraction_in_hii_regions(), cd.exported_mass_hii_regions) ),
    'logMdust.cort': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{cort}})\,[M_\odot]$", lambda: cd.log_dust_mass_as_cortese() ),
    'logMdust.grid': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{grid}})\,[M_\odot]$", lambda: cd.log_dust_mass_from_grid_temperature() ),
    'logMdust.fit/Mstar.zib': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{fit}}/M_{*,\mathrm{zib}})$",
        lambda: cd.log_dust_mass_from_grey_body_fit("limited") - cd.log_stellar_mass_as_zibetti() ),
    'logMdust.cort/Mstar.zib': ( r"$\log_{10}(M_{\mathrm{dust},\mathrm{cort}}/M_{*,\mathrm{zib}})$",
        lambda: cd.log_dust_mass_as_cortese() - cd.log_stellar_mass_as_zibetti() ),

    # dust temperature
    'Tdust.fit.unlim': ( r"$T_{\mathrm{dust},\mathrm{fit},\mathrm{unlim}}\,[\mathrm{K}]$",
        lambda: cd.dust_temperature_from_grey_body_fit("continuum") ),
    'Tdust.fit': ( r"$T_{\mathrm{dust},\mathrm{fit}}\,[\mathrm{K}]$",
        lambda: cd.dust_temperature_from_grey_body_fit("limited") ),
    'Tdust.grid': ( r"$T_{\mathrm{dust},\mathrm{grid}}\,[\mathrm{K}]$", lambda: cd.probe_average_temperature_dust ),

    # public EAGLE database fields (to be provided from column text files)
    'K.db': ( r"$M_\mathrm{K,db}\,[\mathrm{mag}]$", lambda: cd.public_database_magnitude_ukidss_k ),
    'logLk.db': ( r"$\log_{10}(L_\mathrm{K,db})\,[L_{\odot,\mathrm{K}}]$", lambda: cd.public_database_logluminosity_ukidss_k ),
    'logMstar.db': ( r"$\log_{10}(M_{*,\mathrm{db}})\,[M_\odot]$", lambda: log_if_positive(cd.public_database_mass_stars) ),
    'logSFR.db': ( r"$\log_{10}(\mathrm{SFR}_\mathrm{db})\,[M_\odot\,\mathrm{year}^{-1}]$",
        lambda: log_if_positive(cd.public_database_mass_stars*cd.public_database_specific_star_formation_rate) ),
    'logSSFR.db': ( r"$\log_{10}(\mathrm{SFR}_\mathrm{db}/M_{*,\mathrm{db}})\,[\mathrm{year}^{-1}]$",
        lambda: log_if_positive(cd.public_database_specific_star_formation_rate) ),
    'Zgas.db': ( r"$Z_\mathrm{sfgas,db}$", lambda: cd.public_database_star_forming_gas_metallicity ),

    # Star-formation-rate predictions from observations (see Kennicutt-Evans 2012 table 1)
    'logSFR.NUV': ( r"$\log_{10}(\mathrm{SFR}_\mathrm{NUV})\,[M_\odot\,\mathrm{year}^{-1}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_galex_nuv,cd.setup_distance_instrument,'erg/s/Hz') \
                  * c / (1e-6*BroadBandFilter("GALEX.NUV").pivotwavelength())) - 43.17 ),
    'logSFR.24': ( r"$\log_{10}(\mathrm{SFR}_{24\mu\mathrm{m}})\,[M_\odot\,\mathrm{year}^{-1}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_mips_24,cd.setup_distance_instrument,'erg/s/Hz') \
                * c / (1e-6*BroadBandFilter("MIPS.24").pivotwavelength())) - 42.69 ),
    'logSFR.TIR': ( r"$\log_{10}(\mathrm{SFR}_\mathrm{TIR})\,[M_\odot\,\mathrm{year}^{-1}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_uniform_3_1100,cd.setup_distance_instrument,'W/micron',
                    wavelength=np.sqrt(3*1100))*(1100-3)*1e7) - 43.41 ),
    'logSFR.Halpha': ( r"$\log_{10}(\mathrm{SFR}_{\mathrm{H}\alpha})\,[M_\odot\,\mathrm{year}^{-1}]$",
        lambda: log_if_positive(units.luminosityforflux(cd.instr_fluxdensity_halpha,cd.setup_distance_instrument,'erg/s/Hz') \
                * c / 0.6565e-6) - 41.27 ),

}

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
#| instr | optional | the name of the instrument for which to plot data for both x and y axes; defaults to None
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
    global cd       # allow assigning to the global variable that holds the current CollectionData instance
    np.seterr(invalid='ignore')

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
        instr = plotdef.get('instr',None)
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
                cd = CollectionData(collection, instrument=xinstr)
                x = xvalue()
                cd = CollectionData(collection, instrument=yinstr)
                y = yvalue()

                # create a mask that excludes invalid data (i.e. NaN for one of the axes)
                valid = ~(np.isnan(x) | np.isnan(y))

                # plot the relation
                plt.scatter(x[valid], y[valid], marker='o', s=10, alpha=0.5, edgecolors='k', linewidths=(1,), facecolors=color)

                # fit a line through the data and plot it
                xmin = plotdef.get('xmin', x[valid].min())
                xmax = plotdef.get('xmax', x[valid].max())
                if xmin>xmax: xmin,xmax = xmax,xmin
                ymin = plotdef.get('ymin', y[valid].min())
                ymax = plotdef.get('ymax', y[valid].max())
                if ymin>ymax: ymin,ymax = ymax,ymin
                valid = valid & (x>=xmin) & (x<=xmax) & (y>=ymin) & (y<=ymax)
                if np.any(valid):
                    rico, y0 = np.polyfit(x[valid], y[valid], 1)
                    x1 = xmin
                    x2 = xmax
                    y1 = y0 + rico*x1
                    y2 = y0 + rico*x2
                    plt.plot([x1,x2], [y1,y2], color=color, label=collection.label())

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

            # setup the y-axis label (force the x-axis instrument to None to avoid changes to the label)
            ylabel = r"$\log_{10}(N_\mathrm{galaxies})$" if log else r"$N_\mathrm{galaxies}$"
            yinstr = None

            # loop over the collections
            for collection,color in zip(collections,colors):
                # setup the x values in the same way as for a regular plot
                cd = CollectionData(collection, instrument=xinstr)
                x = xvalue()

                # create a mask that excludes invalid data (i.e. NaN)
                valid = ~np.isnan(x)

                # the plt.hist() function does not support square axes with mixed linear/log scale;
                # so, compute the histogram
                xmin = plotdef.get('xmin', x[valid].min())
                xmax = plotdef.get('xmax', x[valid].max())
                counts,binedges = np.histogram(x[valid], bins=bins, range=(xmin,xmax))
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
                plt.plot(xpoints, ypoints, ls='solid', color=color, label=collection.label())

        # set the data limits, if requested
        plt.xlim( xmin=plotdef.get('xmin'), xmax=plotdef.get('xmax') )
        plt.ylim( ymin=plotdef.get('ymin'), ymax=plotdef.get('ymax') )

        # make the plot axes square
        ax.set_aspect(1./ax.get_data_ratio())

        # include instrument names in axis labels if relevant
        if xinstr is not None: xlabel += r"$\;\triangleright\mathrm{"+xinstr+"}$"
        if yinstr is not None: ylabel += r"$\;\triangleright\mathrm{"+yinstr+"}$"

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
