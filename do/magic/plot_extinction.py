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

#import extinction

rcParams['font.family'] = 'serif'

from distutils.core import setup
from Cython.Build import cythonize

from setuptools.extension import Extension

# Import the relevant PTS classes and modules
from pts.core.data.attenuation import CalzettiAttenuationCurve, BattistiAttenuationCurve
#from pts.core.data import extinction
from pts.core.tools import introspection
from pts.core.tools import filesystem as fs

extinction_path = fs.join(introspection.pts_subproject_dir("core"), "data", "extinction.pyx")

#setup(
#    ext_modules = cythonize(extinction_path, include_dirs=[np.get_include()])
#)

#print(ext_modules)

import pyximport

import os
#if os.path.exists(extinction_path):
USE_CYTHON = True
#fname = "extinction.pyx"
fname = extinction_path
#else:
#    USE_CYTHON = False
#    fname = "extinction.c"

core_data_directory_path = fs.join(introspection.pts_subproject_dir("core"), "data")

#extinction_module_name = "extinction"
extinction_module_name = "pts.core.data.extinction"

extern_path = fs.join(introspection.pts_subproject_dir("core"), "data", "extern")
bs_c_path = fs.join(extern_path, "bs.c")
bs_h_path = fs.join(extern_path, "bs.h")
bsplines_path = fs.join(extern_path, "bsplines.pxi")

sourcefiles = [fname, bs_c_path]
dependsfiles = [bs_h_path, bsplines_path]
include_dirs = [np.get_include(), extern_path]
extensions = [Extension(extinction_module_name, sourcefiles, include_dirs=include_dirs,
                        depends=dependsfiles, extra_compile_args=['-std=c99'])]

#from os import environ
#cflags = '-I/' + extern_path[1:]
#ldflags = "-L" + extern_path[1:] + " -lcustomlib"
#print(cflags)
#print(ldflags)

#environ['CFLAGS'] = cflags
##environ['LDFLAGS'] = '-Lpath/to/my/custom/lib -lcustomlib'
#environ['LDFLAGS'] = ldflags

#pyximport.install()
# From https://stackoverflow.com/questions/14657375/cython-fatal-error-numpy-arrayobject-h-no-such-file-or-directory
#pyximport.install(setup_args={"script_args":["--compiler=mingw32"],
#                              "include_dirs":np.get_include()}, reload_support=True)
pyximport.install(build_dir=core_data_directory_path, setup_args={"include_dirs":include_dirs}, reload_support=True, pyimport=True)
pyximport.install(setup_args={"include_dirs":include_dirs}, reload_support=True, pyimport=True)
from pts.core.data import extinction

#extra_compile_args = ['-std=c99']
#pyximport.install(setup_args={"include_dirs": include_dirs, "depends": dependsfiles, "extra_compile_args": extra_compile_args})
#from pts.core.data import extinction

#if USE_CYTHON:
#from Cython.Build import cythonize
#modules = cythonize(extensions)

#from pts.core.data import extinction
#extinction = modules[0]
#import pts.core.data.extinction.pyx

#setup(ext_modules=modules)
#from pts.core.data import extinction
#import pts.core.data.extinction

# -----------------------------------------------------------------

def extinction_figure(wave, a_lambda, residual_from, residual_lims=(-0.1, 0.4), title_text='$R_V = 3.1$'):

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

    return fig

# -----------------------------------------------------------------

wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

a_lambda = {'ccm89': extinction.ccm89(wave, 1.0, 3.1),
            'odonnell94': extinction.odonnell94(wave, 1.0, 3.1),
            'fitzpatrick99': extinction.fitzpatrick99(wave, 1.0),
            'fm07': extinction.fm07(wave, 1.0)}
extinction_figure(wave, a_lambda, 'fitzpatrick99')

#a_lambda = {""}

# -----------------------------------------------------------------