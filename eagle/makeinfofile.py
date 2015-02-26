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

    # save the info file
    infofilepath = simulation.outfilepath("info.txt")
    infofile = open(infofilepath, 'w')
    infofile.write('# Information file for SKIRT-run {}\n'.format(skirtrun.runid()))
    infofile.write('# cells : 1\n')
    infofile.write('# particles : 1\n')
    infofile.write('# luminosity : Lsun\n')
    infofile.write('# mass : Msun\n')
    maxkeylen = max(map(len,info.keys()))
    for key in sorted(info.keys()):
        valueformat = ".0f" if "_particles_" in key or "_cells_" in key else ".9e"
        infofile.write( ("{0:"+str(maxkeylen)+"} = {1:15"+valueformat+"}\n").format(key, info[key]) )
    infofile.close()

    # report success
    print "Created info file " + infofilepath

# -----------------------------------------------------------------
