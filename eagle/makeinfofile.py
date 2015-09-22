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
import pyfits
import pts.archive as arch
from pts.filter import Filter

# -----------------------------------------------------------------

# private list of standard filters for which integrated fluxes should be calculated
_filterspecs = ( \
                "GALEX.FUV","GALEX.NUV",
                "MCPS.U","MCPS.B","MCPS.V","MCPS.I",
                "SDSS.u","SDSS.g","SDSS.r","SDSS.i","SDSS.z",
                "UKIDSS.Z","UKIDSS.Y","UKIDSS.J","UKIDSS.H","UKIDSS.K",
                "2MASS.J", "2MASS.H", "2MASS.K",
                "IRAS.12","IRAS.25","IRAS.60","IRAS.100",
                "IRAC.I1","IRAC.I2","IRAC.I3","IRAC.I4",
                "MIPS.24", "MIPS.70", "MIPS.160",
                "Pacs.blue","Pacs.green","Pacs.red",
                "SPIRE.PSW","SPIRE.PMW","SPIRE.PLW",
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
    info["skirt_run_id"] = skirtrun.runid()

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
    info["setup_particles_cold_gas"] = simulation.coldgasparticles()
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
    wavelengths = simulation.wavelengths()
    # create a mask that removes the carbon line emission peaks from the dust continuum emission
    cmask = (np.abs(wavelengths-157.5)>3) & (np.abs(wavelengths-360)>20) & (np.abs(wavelengths-600)>20)

    # gather statistics on fluxes received by each instrument
    for name in simulation.instrumentnames():
        # maximum flux in Jy
        fluxdensities = simulation.fluxdensities(name, unit='Jy')
        info["instr_"+name+"_fluxdensity_maximum"] = fluxdensities.max()

        # get flux densities per unit of wavelength because filter.convolve() requires this
        fluxdensities = simulation.fluxdensities(name, unit='W/m2/micron')

        # integrated flux density and absolute magnitude for each filter
        for filterspec,filter in _filters.iteritems():
            fluxdensity = filter.convolve(wavelengths, fluxdensities)
            fluxdensity = simulation.convert(fluxdensity, from_unit='W/m2/micron', to_unit='Jy',
                                             wavelength=filter.pivotwavelength())
            magnitude = simulation.absolutemagnitude(fluxdensity, distance,
                                                     fluxdensity_unit='Jy', distance_unit='pc')
            filtername = filterspec.replace(".","_").lower()
            info["instr_"+name+"_fluxdensity_"+filtername] = fluxdensity
            info["instr_"+name+"_magnitude_"+filtername] = magnitude

        # for the Herschel filters, calculate flux and magnitude excluding the carbon line emission peaks
        for filterspec in ("Pacs.blue","Pacs.green","Pacs.red","SPIRE.PSW","SPIRE.PMW","SPIRE.PLW"):
            filter = _filters[filterspec]
            fluxdensity = filter.convolve(wavelengths[cmask], fluxdensities[cmask])
            fluxdensity = simulation.convert(fluxdensity, from_unit='W/m2/micron', to_unit='Jy',
                                             wavelength=filter.pivotwavelength())
            magnitude = simulation.absolutemagnitude(fluxdensity, distance,
                                                     fluxdensity_unit='Jy', distance_unit='pc')
            filtername = filterspec.replace(".","_").lower()
            info["instr_"+name+"_fluxdensity_"+filtername+"_continuum"] = fluxdensity
            info["instr_"+name+"_magnitude_"+filtername+"_continuum"] = magnitude

    # estimate a representative temperature and corresponding standard deviation
    info["probe_average_temperature_dust"], info["probe_stddev_temperature_dust"] = dusttemperature(simulation)

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
    infofile.write('# temperature : K\n')
    infofile.write('# magnitude : dex\n')
    maxkeylen = max(map(len,info.keys()))
    for key in sorted(info.keys()):
        valueformat = ".0f" if "_particles_" in key or "_cells_" in key or "run_id" in key else ".9e"
        infofile.write( ("{0:"+str(maxkeylen)+"} = {1:16"+valueformat+"}\n").format(key, info[key]) )
    infofile.close()

    # report success
    print "Created info file " + infofilepath

# -----------------------------------------------------------------

# Estimate a representative temperature and corresponding standard deviation.
# The data is extracted from the coordinate plane cuts generated by SKIRT for dust
# temperature and dust density, so the data does not cover the complete volume.

# the mass fractions of the dust populations for the Zubko mix used for EAGLE
# (it would be better to load this from the output of each simulation)
fdustpop = np.array((   1.023490e-31, 1.177482e-31, 1.214114e-31, 1.160066e-31, 1.065389e-31,
                        9.811228e-32, 9.531694e-32, 1.037440e-31, 1.347960e-31, 2.180059e-31,
                        4.209576e-31, 7.831375e-31, 9.716384e-31, 7.040336e-31, 1.910658e-31,
                        1.047810e-31, 1.231537e-31, 1.437921e-31, 1.687585e-31, 1.997078e-31,
                        2.393029e-31, 2.920088e-31, 3.653303e-31, 4.713983e-31, 6.274145e-31,
                        8.495666e-31, 1.132516e-30, 1.425789e-30, 1.515186e-30, 4.563139e-31,
                        3.304054e-32, 3.666504e-32, 4.005689e-32, 4.277583e-32, 4.414807e-32,
                        4.322263e-32, 3.881002e-32, 2.976585e-32, 1.618320e-32, 3.522952e-33,
                        3.304054e-32, 3.666504e-32, 4.005689e-32, 4.277583e-32, 4.414807e-32,
                        4.322263e-32, 3.881002e-32, 2.976585e-32, 1.618320e-32, 3.522952e-33 ))
fdustpop /= fdustpop.sum()

# load and return the temperature cut for the specified coordinate plane, averaged over the dust populations
def loadtempcut(simulation, plane):
    filepath = simulation.outfilepath("ds_temp"+plane+".fits")
    T = pyfits.getdata(arch.openbinary(filepath))
    T = np.average(T, weights=fdustpop, axis=0)
    return T

# return the density cut for the specified coordinate plane
def loaddenscut(simulation, plane):
    # load the density cut
    filepath = simulation.outfilepath("ds_trho"+plane+".fits")
    rho = pyfits.getdata(arch.openbinary(filepath))
    return rho

# return the average temperature and corresponding standard deviation for a coordinate plane
def dusttemperature(simulation):
    T = np.dstack([ loadtempcut(simulation, plane) for plane in ('xy','xz','yz') ])
    rho = np.dstack([ loaddenscut(simulation, plane) for plane in ('xy','xz','yz') ])

    # return the average temperature and corresponding standard deviation
    if rho.sum() > 0:
        Tavg = np.average(T, weights=rho)
        Tstd = np.sqrt( (rho*((T-Tavg)**2)).sum()/rho.sum() )
        return ( Tavg, Tstd )
    else:
        return ( 0, 0 )

# -----------------------------------------------------------------
