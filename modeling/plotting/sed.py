#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Import the relevant PTS classes and modules
from ...core.tools.logging import log

import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
line_styles = ['-','--','-.',':']

# -----------------------------------------------------------------

class SEDPlotter(object):
    
    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :return:
        """

        self.models = []
        self.observations = []

    # -----------------------------------------------------------------

    def add_modeled_sed(self, sed):

        """
        This function ...
        :param sed:
        :return:
        """

        # Add the SED to the list
        self.models.append(sed)

    # -----------------------------------------------------------------

    def add_observed_sed(self, sed):

        """
        This function ...
        :param sed:
        :return:
        """

        # Add the SED to the list
        self.observations.append(sed)

    # -----------------------------------------------------------------

    def run(self, output_path, input=None):

        """
        This function ...
        :param output_path:
        :return:
        """

        #if '-ima_name' in sys.argv:
        #  plotfile = str(sys.argv[sys.argv.index('-ima_name')+1])
        #else:
        #  plotfile = 'sed.pdf'


        # survey	band	lambda(nm)	F_nu(Jy)	sigma(Jy)
        #SURVEY,BAND,WAVELENGTH,FLUX,FLUX_ERR = loadtxt(file_with_data, usecols=[0,1,2,3,4],dtype=str, unpack=True,skiprows=1,delimiter='\t')
        #WAVELENGTH = np.array(WAVELENGTH,float)
        #FLUX = np.array(FLUX,float)
        #FLUX_ERR = np.array(FLUX_ERR,float)

        #LABELS = SURVEY

        # plotseds([WAVELENGTH,FLUX,FLUX_ERR],SKIRT_DATA, plotfile, LABELS, figsize=(10,6), xlim=None, ylim=(-5,2))

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Setup the figure
        figure = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(2, 1,height_ratios=[4,1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1],sharex=ax1)



        # Set colours and markers for the flux points
        WAVELENGTH,FLUX,FLUX_ERR = data
        color=iter(cm.rainbow(np.linspace(0,1,len(WAVELENGTH))))
        Labels = list(set(list(labels)))
        filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
        Markers = filled_markers[:len(Labels)]




        # Set axis limits if requested
        if xlim != None: ax1.set_xlim(xlim)
        if ylim != None: ax1.set_ylim(ylim)

        figure.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in figure.axes[:-1]], visible=False)

        # Add axis labels and a legend
        ax2.set_xscale('log')
        ax1.set_xscale('log')
        ax2.set_xlabel(r"Wavelength $\lambda\,[\mu \mathrm{m}]$", fontsize='large')
        ax1.set_ylabel(r"Log $F_\nu$$[Jy]$", fontsize='large')
        ax2.set_ylabel(r"Residuals $[\%]$", fontsize='large')

        # Add the legend
        ax1.legend(numpoints=1,loc=4,frameon=True,ncol=2,fontsize=11)

        # Save the figure
        plt.savefig(plotfile, bbox_inches='tight', pad_inches=0.25)
        plt.show()
        plt.close()

# -----------------------------------------------------------------
