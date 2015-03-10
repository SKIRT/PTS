#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.plotresults Plot histograms and scaling relations for a collection of EAGLE SKIRT-runs.
#
# The class in this module serves to plot histograms and scaling relations for the results in a set of EAGLE
# SKIRT-runs that have been previously collected in a single data file.

# -----------------------------------------------------------------

import pickle
import os.path
import numpy as np

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

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
    'Ncells': ( r"$N_\mathrm{cells}/10^6$", lambda: setup_cells_dust_grid/1e6 ),
    'taumax': ( r"$\tau_\mathrm{max}$", lambda: setup_optical_depth_maximum ),
    'tau90': ( r"$\tau_\mathrm{90}$", lambda: setup_optical_depth_percentile90 ),

    # intrinsic properties
    'logMstar': ( r"$\log_{10}(M_*\,[M_\odot])$", lambda: np.log10(original_mass_stars) ),
    'logMdust': ( r"$\log_{10}(M_\mathrm{dust}\,[M_\odot])$", lambda: np.log10(setup_mass_dust) ),
    'logMdust/Mstar': ( r"$\log_{10}(M_\mathrm{dust}/M_*)$", lambda: np.log10(setup_mass_dust/original_mass_stars) ),

    # magnitudes and colors
    "g" : ( r"$M_\mathrm{r}$", lambda: instr_magnitude_sdss_g ),
    "r" : ( r"$M_\mathrm{r}$", lambda: instr_magnitude_sdss_r ),
    'g-r' : ( r"$\mathrm{g}-\mathrm{r}$", lambda: instr_magnitude_sdss_g - instr_magnitude_sdss_r ),
    'NUV-r' : ( r"$\mathrm{NUV}-\mathrm{r}$", lambda: instr_magnitude_galex_nuv - instr_magnitude_sdss_r ),
}

# -----------------------------------------------------------------

class Collection:
    # ---------- Constructing -------------------------------------

    ## The constructor loads the contents of the specified collection so that it is ready for plotting.
    # The collection name should \em not include the directory (which is taken from eagle.conf) nor the
    # postfix "_info_collection.dat".
    def __init__(self, collectionname):
        self._collectionname = collectionname
        self._info = {}

        # load the collection
        infilepath = os.path.join(config.collections_path, collectionname+"_info_collection.dat")
        infile = open(infilepath, "r")
        self._info['any'] = pickle.load(infile)
        infile.close()

        # construct filtered dicts for each instrument name
        names = set([ key.split("_")[1] for key in filter(lambda key: key.startswith("instr_"), self._info['any'].keys()) ])
        for name in names:
            self._info[name] = {}
            for key,value in self._info['any'].iteritems():
                if key.startswith("instr_"):
                    segments = key.split("_")
                    if segments[1]==name:
                        segments.pop(1)
                        cleankey = "_".join(segments)
                        self._info[name][cleankey] = value
                else:
                    self._info[name][key] = value

    # ---------- Plotting -----------------------------------------

    ## This function produces a one-page pdf file with one or more plots for the collection of SKIRT-runs.
    # It expects the following arguments:
    # - plotname: name of the output plot \em not including the directory, collection name, nor filename extension
    # - plotdefs: sequence of plot definitions; each item is a dictionary specifying a single plot as described below
    # - pagesize: a 2-tuple specifying the size of the complete page in inch; default is A4 format
    # - layout: a 2-tuple specifying the number of columns and rows in the layout of the plots; default is 2 by 3
    #           (the layout must accomodate all items in the plotdefs sequence)
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
    #| xmax| optional | the maximum x value shown; default is largest x value
    #| ymin | optional | the minimum y value shown; default is smallest y value
    #| ymax| optional | the maximum y value shown; default is largest y value
    #
    def plotresults(self, plotname, plotdefs, layout=(2,3), pagesize=(8.268,11.693)):

        # setup the figure
        figure = plt.figure(figsize=pagesize)
        figure.subplots_adjust(wspace=0.15, hspace=0.23,
                               left=0.08, right=0.97, top=0.93, bottom=0.1)

        # add figure title
        plt.suptitle(plotname + " for " + self._collectionname)

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

                # setup the x and y values for each axes,
                # after loading the statistics for the appropriate instrument as global variables
                # that can be used in the callables that setup the x and y values for each axis
                globals().update(self._info[xinstr])
                x = xvalue()
                globals().update(self._info[yinstr])
                y = yvalue()

                # plot the relation
                plt.scatter(x, y, marker='o', s=15, edgecolors='k', facecolors='r')

            # for a histogram...
            else:
                # get the histogram options
                bins = plotdef.get('bins', 10)
                log = plotdef.get('log', False)

                # setup the x-axis as for a relation plot
                xlabel,xvalue = axistypes[xaxis]
                globals().update(self._info[xinstr])
                x = xvalue()

                # setup the y-axis label (force the x-axis instrument to 'any' to avoid changes to the label)
                ylabel = r"$\log_{10}(N_\mathrm{galaxies})$" if log else r"$N_\mathrm{galaxies}$"
                yinstr = 'any'

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
                plt.plot(xpoints, ypoints, ls='solid', color='r')

            # set the data limits, if requested
            plt.xlim( xmin=plotdef.get('xmin'), xmax=plotdef.get('xmax') )
            plt.ylim( ymin=plotdef.get('ymin'), ymax=plotdef.get('ymax') )

            # make the plot axes square
            ax.set_aspect(1./ax.get_data_ratio())

            # include instrument names in axis labels if relevant
            if xinstr != 'any': xlabel += r"$\;\mathrm{'"+xinstr+"'}$"
            if yinstr != 'any': ylabel += r"$\;\mathrm{'"+yinstr+"'}$"

            # add axis labels
            plt.xlabel(xlabel, fontsize='medium')
            plt.ylabel(ylabel, fontsize='medium')

            # fine-tune the tick label size
            for item in (ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize('x-small')

        # save and close the figure
        plotfilepath = os.path.join(config.plots_path, self._collectionname+"_"+plotname+".pdf")
        plt.savefig(plotfilepath)
        plt.close()
        print "Created results plot file", plotfilepath

# -----------------------------------------------------------------
