#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package eagle.makeinfofile Creating information files with statistics on EAGLE SKIRT-runs.
#
# The main function in this module creates an information file with statistics on the results of an EAGLE SKIRT-run.
# The information file has a simple text format and includes statistics on the data exported from the EAGLE snapshot
# (i.e. the input data for the SKIRT simulation), on the SKIRT setup (such as the dust grid), and on the results
# of the SKIRT simulation (such as fluxes in various bands).

# -----------------------------------------------------------------

import os.path
import numpy as np
import pts.archive as arch
from pts.filter import Filter

# -----------------------------------------------------------------

# private list of standard filters for which integrated fluxes should be calculated
_filterspecs = ( \
                "GALEX.FUV","GALEX.NUV",
                "MCPS.U","MCPS.B","MCPS.V","MCPS.I",
                "SDSS.u","SDSS.g","SDSS.r","SDSS.i","SDSS.z",
                "UKIDSS.Z","UKIDSS.Y","UKIDSS.K","UKIDSS.J","UKIDSS.H",
                "IRAS.12","IRAS.25","IRAS.60","IRAS.100",
                "Pacs.blue","Pacs.green","Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW",
                )

# private list of uniform filters for which integrated fluxes should be calculated;
# for each filter, specify the filter spec and the wavelength range in micron)
_uniformfilters = ( ("Uniform_8_1000",8,1000), )

# private global dictionary to hold filter objects: (key,value) = (filterspec, Filter object)
_filters = { }

## This private function ensures that the filters relevant for this module are loaded into the private global
# dictionary \em _filters, using the filter specification string as a key and the corresponding Filter instance
# as a value. When the function is called for the first time, it creates the Filter objects and stores them
# in the dictionary. Once the dictionary has been loaded, subsequent invocations of the function don't do anything.
def _loadfilters():
    if len(_filters) == 0:
        # standard filters
        for filterspec in _filterspecs:
            _filters[filterspec] = Filter(filterspec)
        # uniform filters
        for filterspec,wavemin,wavemax in _uniformfilters:
            _filters[filterspec] = Filter((wavemin,wavemax))

# -----------------------------------------------------------------

# This function creates an information file with statistics on the results of an EAGLE SKIRT-run.
# The information file has a simple text format and includes statistics on the data exported from the EAGLE snapshot
# (i.e. the input data for the SKIRT simulation), on the SKIRT setup (such as the dust grid), and on the results
# of the SKIRT simulation (such as fluxes in various bands).
#
# The information file is placed in the simulation's output directory, and is named "prefix_info.txt".
#
def makeinfofile(skirtrun):
    simulation = skirtrun.simulation()

    # the info dict that will be saved at the end
    info = { }

    # load statistics on the EAGLE data from the info file in the input folder
    inpath = skirtrun.inpath()
    infile = arch.listdir(inpath, "_info.txt")[0]
    for line in arch.opentext(os.path.join(inpath,infile)):
        if not line.startswith("#"):
            segments = line.split()
            info[segments[0]] = float(segments[2])

    # gather SKIRT setup statistics
    info["setup_mass_dust"] = simulation.dustmass()
    info["setup_mass_dust_grid"] = simulation.dustgridmass()
    info["setup_mass_cold_gas"] = simulation.coldgasmass()
    info["setup_mass_metallic_gas"] = simulation.metallicgasmass()
    info["setup_initial_mass_stars"] = simulation.initialstellarmass()
    info["setup_mass_hii_regions"] = simulation.hiiregionmass()
    info["setup_luminosity_stars"] = simulation.stellarluminosity()
    info["setup_luminosity_hii_regions"] = simulation.hiiregionluminosity()
    info["setup_cells_dust_grid"], info["setup_optical_depth_maximum"], \
        info["setup_optical_depth_percentile90"] = map(float,simulation.dustcellstats())
    distance = simulation.instrumentdistance()
    info["setup_distance_instrument"] = distance

    # load filters and wavelength grid
    _loadfilters()
    wavebincenters, wavebinwidths = simulation.wavelengthbins(unit='micron')

    # gather statistics on fluxes received by each instrument
    for name in simulation.instrumentnames():
        # maximum flux in Jy
        fluxdensities = simulation.fluxdensities(name, unit='Jy')
        info["result_"+name+"_fluxdensity_maximum"] = fluxdensities.max()

        # get flux densities per unit of wavelength because filter.apply() requires this
        fluxdensities = simulation.fluxdensities(name, unit='W/m2/micron')

        # integrated flux density and absolute magnitude for each filter
        for filterspec,filter in _filters.iteritems():
            fluxdensity = filter.apply(wavebincenters, wavebinwidths, fluxdensities)
            fluxdensity = simulation.convert(fluxdensity, from_unit='W/m2/micron', to_unit='Jy',
                                             wavelength=filter.meanwavelength())
            magnitude = simulation.absolutemagnitude(fluxdensity, distance,
                                                     fluxdensity_unit='Jy', distance_unit='pc')
            filtername = filterspec.replace(".","_").lower()
            info["result_"+name+"_fluxdensity_"+filtername] = fluxdensity
            info["result_"+name+"_magnitude_"+filtername] = magnitude

    # save the info file
    infofilepath = simulation.outfilepath("info.txt")
    infofile = open(infofilepath, 'w')
    infofile.write('# Information file for SKIRT-run {}\n'.format(skirtrun.runid()))
    infofile.write('# cells : 1\n')
    infofile.write('# particles : 1\n')
    infofile.write('# luminosity : Lsun\n')
    infofile.write('# mass : Msun\n')
    infofile.write('# distance : pc\n')
    infofile.write('# fluxdensity : Jy\n')
    infofile.write('# magnitude : dex\n')
    maxkeylen = max(map(len,info.keys()))
    for key in sorted(info.keys()):
        valueformat = ".0f" if "_particles_" in key or "_cells_" in key else ".9e"
        infofile.write( ("{0:"+str(maxkeylen)+"} = {1:16"+valueformat+"}\n").format(key, info[key]) )
    infofile.close()

    # report success
    print "Created info file " + infofilepath

# -----------------------------------------------------------------
