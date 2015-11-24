#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package do.eagle_build Build visualization files for completed EAGLE SKIRT-runs.
#
# This script builds plots, images, movies, or info files from the results of completed EAGLE SKIRT-runs.
# The resulting files are placed in the corresponding \c vis directories, replacing any previous versions.
#
# The script expects exactly two command-line arguments. The first argument specifies what type of files to build.
# It should be one of the following strings (it is sufficient to specify an unambiguous portion of the string):
#  - densities: a plot with histograms of stellar, gas and dust densities in function of the galaxy radius
#  - greybodyfit: a plot showing a modified blackbody fit to the dust continuum emission and listing the temperature
#  - infofile: a text info file with statistics on the simulation results
#  - particles: a 3D plot of the gas particles indicating star forming, cold and hot gas in different colors
#  - rgbimages: an optical RGB image for each of the frame instruments
#  - seds: a plot combining the SEDs for all of the instruments
#  - temperature: a plot with a histogram of dust mass versus dust temperature
#  - wavemovie: a wavelength movie for the three frame instruments
#
# The second argument further specifies the operation that should be performed. It should be one of the following:
#  - the string "update" to perform the build for all completed SKIRT-runs for which the visualization files do not
#    yet exist.
#  - the string "rebuild" to perform the build for any and all all completed SKIRT-runs, even if the visualization
#    files already exist. Previous versions are replaced.
#  - a comma-seperated list of run-ids and/or run-id ranges (expressed as two run-ids with a dash in between)
#    for which the build should be performed. Visualization files will be rebuilt even if they already exist.
#    The script does not verify whether the specified SKIRT-runs have been completed.
#

# -----------------------------------------------------------------

import sys
import types
import os.path
from eagle.database import Database
from eagle.skirtrun import SkirtRun
from eagle.skirtrun import runids_in_range

# -----------------------------------------------------------------

# a list of relevant filename endings for each visualization type
filenames_for_vistype = {
    'densities':    ( "density_curves.pdf", ),
    'greybodyfit':  ( "dust_body_fit.pdf", ),
    'infofile':     ( "info.txt", ),
    'particles':    ( "gas_particles.pdf", ),
    'rgbimages':    ( "xy_total_optical.png", "xz_total_optical.png", "yz_total_optical.png",
                      "xy_total_augmented.png", "xz_total_augmented.png", "yz_total_augmented.png" ),
    'seds':         ( "sed.pdf", ),
    'temperature':  ( "dust_temperature.pdf", ),
    'wavemovie':    ( "wave.mov", ),
}

# -----------------------------------------------------------------

# returns a list of SkirtRun objects corresponding to all completed or archived skirt-runs, in order of run-id,
# optionally omitting any skirt-runs for which all files in the specified sequence exist in the visualization folder
def completed_skirtruns(unless_filenames=None):
    db = Database()
    runids = sorted([ row['runid'] for row in db.select("runstatus='completed' or runstatus='archived'") ])
    runs = [ SkirtRun(runid) for runid in runids ]
    if unless_filenames!=None:
        runs = filter(lambda run: not has_visualization_files(run,unless_filenames), runs)
    db.close()
    return runs

# returns True if all files in the specified sequence exist in the visualization folder of the specified skirt-run
def has_visualization_files(skirtrun, filenames):
    filenames = filenames if isinstance(filenames, (types.TupleType,types.ListType)) else [ filenames ]
    vispath = skirtrun.vispath()
    prefix = skirtrun.prefix()
    for filename in filenames:
        if not os.path.isfile(os.path.join(vispath, prefix+"_"+filename)): return False
    return True

# move any of the files in the specified sequence from the skirt-run output directory to the visualization directory
def move_visualization_files(skirtrun, filenames):
    for visfile in filter(lambda fn: fn.endswith(filenames), os.listdir(skirtrun.outpath())):
        os.rename(os.path.join(skirtrun.outpath(), visfile),
                  os.path.join(skirtrun.vispath(), visfile))

# -----------------------------------------------------------------

# get and parse the command-line arguments
if len(sys.argv) != 3: raise ValueError("This script expects exactly two command-line arguments")

vistype = sys.argv[1].lower()
found = False
for knowntype in filenames_for_vistype.keys():
    if knowntype.startswith(vistype):
        vistype = knowntype
        found = True
        break
if not found: raise ValueError("Unknown visualization type: " + vistype)

buildrange = sys.argv[2].lower()
update =  "update".startswith(buildrange)
rebuild = "rebuild".startswith(buildrange)
runidrange = runids_in_range(buildrange)
if not update and not rebuild and not runidrange: raise ValueError("Unknown build range: " + buildrange)

# -----------------------------------------------------------------

# get a list of relevant filename endings depending on the visualization type
filenames = filenames_for_vistype[vistype]

# construct the list of SKIRT-runs to be processed
if runidrange:
    skirtruns = [ SkirtRun(runid) for runid in runidrange ]
else:
    skirtruns = completed_skirtruns(filenames if update else None)

# =================================================================

# build density curves
if vistype=='densities':
    from eagle.plotdensitycurves import plotdensitycurves
    print "Building density curves for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building density curves for SKIRT-run {}...".format(skirtrun.runid())
        plotdensitycurves(skirtrun)

# -----------------------------------------------------------------

# build grey body fits
if vistype=='greybodyfit':
    from eagle.plotgreybodyfit import plotgreybodyfit
    print "Building grey body fits for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building grey body fit for SKIRT-run {}...".format(skirtrun.runid())
        plotgreybodyfit(skirtrun)

# -----------------------------------------------------------------

# build info files
if vistype=='infofile':
    from eagle.makeinfofile import makeinfofile
    print "Building info files for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building info file for SKIRT-run {}...".format(skirtrun.runid())
        makeinfofile(skirtrun)
        move_visualization_files(skirtrun, filenames)

# -----------------------------------------------------------------

# build gas particle plots
if vistype=='particles':
    from eagle.plotgasparticles import plotgasparticles
    print "Building gas particle plots for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building gas particle plot for SKIRT-run {}...".format(skirtrun.runid())
        plotgasparticles(skirtrun)

# -----------------------------------------------------------------

# build RGB images
if vistype=='rgbimages':
    from misc.makergbimages import makeintegratedrgbimages
    from pts.filter import Filter
    print "Building RGB images for {} SKIRT-runs".format(len(skirtruns))
    filterR = Filter('SDSS.i')
    filterG = Filter('SDSS.r')
    filterB = Filter('SDSS.g')
    filterIR = Filter('Pacs.green')
    filterUV = Filter('GALEX.FUV')
    for skirtrun in skirtruns:
        print "Building RGB images for SKIRT-run {}...".format(skirtrun.runid())
        fmin,fmax = makeintegratedrgbimages(skirtrun.simulation(),
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1) ], postfix="_optical")
        makeintegratedrgbimages(skirtrun.simulation(),
            [ (filterR, 1,0,0), (filterG, 0,1,0), (filterB, 0,0,1), (filterIR, 0.02,0,0), (filterUV, 0,0,4) ],
            postfix="_augmented", fmin=fmin, fmax=fmax)
        move_visualization_files(skirtrun, filenames)

# -----------------------------------------------------------------

# build SED plots for each SKIRT-run
if vistype=='seds':
    from plotting.seds import plotseds
    print "Building SED plots for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building SED plot for SKIRT-run {}...".format(skirtrun.runid())
        plotseds(skirtrun.simulation())
        move_visualization_files(skirtrun, filenames)

# -----------------------------------------------------------------

# build temperature histograms
if vistype=='temperature':
    from eagle.plottemperature import plottemperature
    print "Building temperature histograms for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building temperature histogram for SKIRT-run {}...".format(skirtrun.runid())
        plottemperature(skirtrun)

# -----------------------------------------------------------------

# build wavelength movies for each SKIRT-run
if vistype=='wavemovie':
    from misc.makewavemovie import makewavemovie
    print "Building wavelength movies for {} SKIRT-runs".format(len(skirtruns))
    for skirtrun in skirtruns:
        print "Building wavelength movie for SKIRT-run {}...".format(skirtrun.runid())
        makewavemovie(skirtrun.simulation())
        move_visualization_files(skirtrun, filenames)

print "Done..."

# -----------------------------------------------------------------
