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

# dictionary of tuples representing a type of plot axis
#  key: axis type identifier
#  value : ( human-readable label, callable returning the axis value for a given simulation )
_axisdefinitions = {

    # intrinsic properties
    'logMstar': ( r"$\log_{10}(M_*\,[M_\odot])$", lambda: np.log10(original_mass_stars) ),
    'logMdust': ( r"$\log_{10}(M_\mathrm{dust}\,[M_\odot])$", lambda: np.log10(setup_mass_dust) ),
    'logMdust/Mstar': ( r"$\log_{10}(M_\mathrm{dust}/M_*)$", lambda: np.log10(setup_mass_dust/original_mass_stars) ),

    # magnitudes and colors
    "r" : ( r"$M_\mathrm{r}$", lambda: result_magnitude_sdss_r ),
    'g-r' : ( r"$\mathrm{g}-\mathrm{r}$", lambda: result_magnitude_sdss_g - result_magnitude_sdss_r ),
    'NUV-r' : ( r"$\mathrm{NUV}-\mathrm{r}$", lambda: result_magnitude_galex_nuv - result_magnitude_sdss_r ),
}

# -----------------------------------------------------------------

class Collection:
    # ---------- Constructing -------------------------------------

    ## The constructor loads the contents of the specified collection so that it is ready for plotting.
    # The collection name should \em not include the directory (which is taken from eagle.conf) nor the
    # postfix "_info_collection.dat".
    def __init__(self, collectionname):
        self._collectionname = collectionname

        # load the collection
        infilepath = os.path.join(config.collections_path, collectionname+"_info_collection.dat")
        infile = open(infilepath, "r")
        info = pickle.load(infile)
        infile.close()

        # construct filtered dicts for each instrument name
        names = set([ key.split("_")[1] for key in filter(lambda key: key.startswith("result_"), info.keys()) ])
        self._info = {}
        for name in names:
            self._info[name] = {}
            for key,value in info.iteritems():
                if key.startswith("result_"):
                    segments = key.split("_")
                    if segments[1]==name:
                        segments.pop(1)
                        cleankey = "_".join(segments)
                        self._info[name][cleankey] = value
                else:
                    self._info[name][key] = value

    # ---------- Plotting -----------------------------------------

    ## This function produces a one-page pdf file with one or more relation plots for the collection of SKIRT-runs.
    # It expects the following arguments:
    # - plotname: name of the output plot \not including the directory, colllection name, nor filename extension
    # - instrument: name of the instrument for which to plot data; must be a valid instrument even if not used
    # - plotdefs: sequence of plot definitions; each item is a 2-tuple of axis type identifiers (key in above dict)
    # - pagesize: a 2-tuple specifying the size of the complete page in inch
    # - layout: a 2-tuple specifying the number of columns and rows in the layout of the plots
    #           (the layout must accomodate all items in the plotdefs sequence)
    def plotrelations(self, plotname, instrument, plotdefs, layout=(2,3), pagesize=(8,12)):

        # setup the figure
        figure = plt.figure(figsize=pagesize)
        figure.subplots_adjust(wspace=0.35, hspace=0.25)

        # load the statistics for the appropriate instrument as global variables
        # that can be used in the callables that setup the x and y values for each axis
        globals().update(self._info[instrument])

        # loop over the plots
        plotindex = 0
        for xaxis,yaxis in plotdefs:
            xlabel,xvalue = _axisdefinitions[xaxis]
            ylabel,yvalue = _axisdefinitions[yaxis]

            # start the appropriate subplot
            plotindex += 1
            plt.subplot(layout[1], layout[0], plotindex)

            # setup the x and y values for each axes
            x = [ xvalue() ]
            y = [ yvalue() ]

            # plot the scaling relation
            plt.scatter(x, y)

            # add axis labels
            plt.xlabel(xlabel, fontsize='large')
            plt.ylabel(ylabel, fontsize='large')

        # save and close the figure
        plotfilepath = os.path.join(config.plots_path, self._collectionname+"_"+plotname+".pdf")
        plt.savefig(plotfilepath, bbox_inches='tight', pad_inches=0.25)
        plt.close()
        print "Created relations plot file", plotfilepath

# -----------------------------------------------------------------
