#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## Plot extinction functions for comparison

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

rcParams['font.family'] = 'serif'

from pts.core.data.extinction import CardelliClaytonMathisExtinctionCurve, ODonnellExtinctionCurve
from pts.core.data.extinction import FitzpatrickExtinctionCurve, FitzpatrickMassaExtinctionCurve
from pts.core.data.extinction import CalzettiExtinctionCurve, BattistiExtinctionCurve
from pts.core.basics.range import QuantityRange

# -----------------------------------------------------------------

def extinction_figure(wave, a_lambda, residual_from, residual_lims=(-0.1, 0.4), title_text='$R_V = 3.1$'):

    """
    This function ...
    :param wave:
    :param a_lambda:
    :param residual_from:
    :param residual_lims:
    :param title_text:
    :return:
    """

    names = list(a_lambda.keys())  # consistent ordering between panels
    fig = plt.figure(figsize=(8.5, 6.))

    ax = plt.axes()
    for name in names:
        plt.plot(wave, a_lambda[name], label=name)

    plt.axvline(x=2700., ls=':', c='k')
    plt.axvline(x=3030.3030, ls=':', c='k')
    plt.axvline(x=9090.9091, ls=':', c='k')
    plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
    plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)    

    plt.text(0.65, 0.95, title_text, transform=ax.transAxes, va='top',
             ha='right', size='x-large')

    plt.ylabel('Extinction ($A(\lambda)$ / $A_V$)')

    plt.legend()

    plt.setp(ax.get_xticklabels(), visible=False)

    divider = make_axes_locatable(ax)
    axresid = divider.append_axes("bottom", size=2.0, pad=0.2, sharex=ax)

    for name in names:
        #print(a_lambda[name], a_lambda[residual_from])
        plt.plot(wave, a_lambda[name] - a_lambda[residual_from])

    plt.axvline(x=2700., ls=':', c='k')
    plt.axvline(x=3030.3030, ls=':', c='k')
    plt.axvline(x=9090.9091, ls=':', c='k')
    plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
    plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)

    plt.xlim(wave[0], wave[-1])
    plt.ylim(ymin=residual_lims[0], ymax=residual_lims[1])

    plt.ylabel('residual from ' + residual_from)
    plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')

    ax.set_xscale('log')
    axresid.set_xscale('log')

    plt.tight_layout()

    #return fig

    plt.show()

# -----------------------------------------------------------------

#wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

# Create wavelengths
wavelenth_range = QuantityRange(910., 30000., unit="AA")
wavelengths = wavelenth_range.log(2000)
#print(wavelengths)

wavelengths_scalar = wavelenth_range.log(2000, add_unit=False, as_list=True)

# -----------------------------------------------------------------

a_lambda = {'ccm89': CardelliClaytonMathisExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True),
            'odonnell94': ODonnellExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True),
            'fitzpatrick99': FitzpatrickExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True),
            'fm07': FitzpatrickMassaExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True),
            'calzetti': CalzettiExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True),
            'battisti': BattistiExtinctionCurve(wavelengths=wavelengths).extinctions(asarray=True)}
extinction_figure(wavelengths_scalar, a_lambda, 'fitzpatrick99')

# -----------------------------------------------------------------