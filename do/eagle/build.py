#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.eagle.build Build visualization files for completed EAGLE SKIRT-runs.
#
# This script builds plots, images, movies, or info files from the results of completed EAGLE SKIRT-runs.
# The resulting files are placed in the corresponding \c vis directories, replacing any previous versions.
#
# The script expects exactly two command-line arguments. The first argument specifies what type of files to build.
# It should be one of the following strings (it is sufficient to specify an unambiguous portion of the string):
#  - bandimages: band-integrated optical and augmented RGB images for each of the frame instruments
#  - densities: a plot with histograms of stellar, gas and dust densities in function of the galaxy radius
#  - fastimages: a fast optical RGB image for each of the frame instruments
#  - greybodyfit: a plot showing a modified blackbody fit to the dust continuum emission and listing the temperature
#  - infofile: a text info file with statistics on the simulation results
#  - particles: a 3D plot of the gas particles indicating star forming, cold and hot gas in different colors
#  - seds: a plot combining the SEDs for all of the instruments
#  - temperature: a plot with a histogram of dust mass versus dust temperature
#  - wavemovie: a wavelength movie for the three frame instruments
#
# The second argument specifies the SKIRT-runs for which the operation should be performed. It is a
# comma-seperated list of run-ids and/or run-id ranges (expressed as two run-ids with a dash in between).
#
# This script does not verify whether the specified SKIRT-runs have been completed, nor does it update
# the SKIRT-runs database in any way.
#

# -----------------------------------------------------------------

# Import standard modules
import sys
import types
import os.path

# Import the relevant PTS classes and modules
from pts.eagle.skirtrun import SkirtRun
from pts.eagle.skirtrun import runids_in_range

# -----------------------------------------------------------------

# get and parse the command-line arguments
if len(sys.argv) != 3:
    raise ValueError("This script expects exactly two command-line arguments: outputtype runidspec")

vistypes = ('gepimages', 'bandimages', 'densities', 'fastimages', 'greybodyfit', 'infofile',
            'particles', 'seds', 'temperature', 'wavemovie')
inputtype = sys.argv[1].lower()
found = False
for knowntype in vistypes:
    if knowntype.startswith(inputtype):
        vistype = knowntype
        found = True
        break
if not found: raise ValueError("Unknown visualization type: " + inputtype)

runidspec = sys.argv[2]
runids = runids_in_range(runidspec)

# construct a list of SKIRT-runs to be processed
skirtruns = [ SkirtRun(runid) for runid in runids ]

# =================================================================

# build GEP-band-integrated gray-scale images
if vistype=='gepimages':
    from pts.core.plot.rgbimages import makeintegratedrgbimages
    from pts.core.filter.broad import BroadBandFilter
    print "Building GEP-band-integrated gray-scale images for {} SKIRT-runs".format(len(skirtruns))
    # List of filter objects
    bandedges = (
    ( 1,  10.00,  11.33),
    ( 2,  11.33,  12.84),
    ( 3,  12.84,  14.56),
    ( 4,  14.56,  16.50),
    ( 5,  16.50,  18.70),
    ( 6,  18.70,  21.19),
    ( 7,  21.19,  24.02),
    ( 8,  24.02,  27.22),
    ( 9,  27.22,  30.85),
    (10,  30.85,  34.96),
    (11,  34.96,  39.62),
    (12,  39.62,  44.90),
    (13,  44.90,  50.89),
    (14,  50.89,  57.68),
    (15,  57.68,  65.37),
    (16,  65.37,  74.08),
    (17,  74.08,  83.96),
    (18,  83.96,  95.16),
    (19,  95.16, 127.00),
    (20, 127.00, 169.00),
    (21, 169.00, 225.00),
    (22, 225.00, 300.00),
    (23, 300.00, 400.00) )
    for skirtrun in skirtruns:
        print "Building GEP-band-integrated gray-scale images for SKIRT-run {}...".format(skirtrun.runid())
        for n,wmin,wmax in bandedges:
            makeintegratedrgbimages(skirtrun.simulation(),
                [ (BroadBandFilter((wmin,wmax)), 1,1,1) ],
                postfix="_gep{:02d}".format(n), output_path=skirtrun.vispath())

# -----------------------------------------------------------------

# build band-integrated images
if vistype=='bandimages':
    from pts.core.plot.rgbimages import makeintegratedrgbimages
    from pts.core.filter.broad import BroadBandFilter
    print "Building band-integrated images for {} SKIRT-runs".format(len(skirtruns))
    filterR = BroadBandFilter('SDSS.i')
    filterG = BroadBandFilter('SDSS.r')
    filterB = BroadBandFilter('SDSS.g')
    filterIR = BroadBandFilter('Pacs.green')
    filterUV = BroadBandFilter('GALEX.FUV')
    for skirtrun in skirtruns:
        print "Building band-integrated images for SKIRT-run {}...".format(skirtrun.runid())
        fmin,fmax = makeintegratedrgbimages(skirtrun.simulation(),
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1) ],
            postfix="_optical", output_path=skirtrun.vispath())
        makeintegratedrgbimages(skirtrun.simulation(),
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1), (filterIR, 0.02,0,0), (filterUV, 0,0,4) ],
            postfix="_augmented", fmin=fmin, fmax=fmax, output_path=skirtrun.vispath())

# -----------------------------------------------------------------

# build density curves
if vistype=='densities':
    from pts.eagle.plotdensitycurves import plotdensitycurves
    print "Building density curves for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building density curves for SKIRT-run {}...".format(skirtrun.runid())
        plotdensitycurves(skirtrun)

# -----------------------------------------------------------------

# build fast optical images
if vistype=='fastimages':
    from pts.core.plot.rgbimages import makergbimages
    print "Building fast optical images for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building fast optical images for SKIRT-run {}...".format(skirtrun.runid())
        makergbimages(skirtrun.simulation(), wavelength_tuples=((0.753,0.617,0.470),), output_path=skirtrun.vispath())

# -----------------------------------------------------------------

# build grey body fits
if vistype=='greybodyfit':
    from pts.eagle.plotgreybodyfit import plotgreybodyfit
    print "Building grey body fits for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building grey body fit for SKIRT-run {}...".format(skirtrun.runid())
        plotgreybodyfit(skirtrun)

# -----------------------------------------------------------------

# build info files
if vistype=='infofile':
    from pts.eagle.database import Database
    from pts.eagle.makeinfofile import makeinfofile
    print "Building info files for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building info file for SKIRT-run {}...".format(skirtrun.runid())
        record = Database().select("runid = ?", (skirtrun.runid(),))[0]
        makeinfofile(skirtrun, record['snaptag'])

# -----------------------------------------------------------------

# build gas particle plots
if vistype=='particles':
    from pts.eagle.plotgasparticles import plotgasparticles
    print "Building gas particle plots for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building gas particle plot for SKIRT-run {}...".format(skirtrun.runid())
        plotgasparticles(skirtrun)

# -----------------------------------------------------------------

# build SED plots for each SKIRT-run
if vistype=='seds':
    from pts.core.plot.seds import plotseds
    print "Building SED plots for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building SED plot for SKIRT-run {}...".format(skirtrun.runid())
        plotseds(skirtrun.simulation(), output_path=skirtrun.vispath())

# -----------------------------------------------------------------

# build temperature histograms
if vistype=='temperature':
    from pts.eagle.plottemperature import plottemperature
    print "Building temperature histograms for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building temperature histogram for SKIRT-run {}...".format(skirtrun.runid())
        plottemperature(skirtrun)

# -----------------------------------------------------------------

# build wavelength movies for each SKIRT-run
if vistype=='wavemovie':
    from pts.core.plot.wavemovie import makewavemovie
    print "Building wavelength movies for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building wavelength movie for SKIRT-run {}...".format(skirtrun.runid())
        makewavemovie(skirtrun.simulation(), output_path=skirtrun.vispath())

print "Done..."

# -----------------------------------------------------------------
