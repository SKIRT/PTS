#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.visualizer Creating visualizations for SKIRT simulation output.
#
# The Visualizer class in this module assists with creating visualizations for SKIRT simulation output.
# A Visualizer instance holds one or more SkirtSimulation objects, and it provides functions for
# automatically creating visualizations (such an SED plot or an RGB image) for all relevant output files
# in all of these simulations.

# -----------------------------------------------------------------

import os
import os.path
import types
import numpy as np
from pts.skirtsimulation import SkirtSimulation
from pts.plotgridpdf import plotgrid
from pts.plotsed import plotsed
from pts.rgbimage import RGBImage
from pts.plotpolarization import plotpolarizationmap
from pts.wavemovie import makewavemovie

# -----------------------------------------------------------------
#  Visualizer class
# -----------------------------------------------------------------

## An instance of the SkirtSimulation class holds one or more SkirtSimulation objects, and it provides functions
# for automatically creating visualizations (such an SED plot or an RGB image) for all relevant output files
# in all of these simulations.
#
class Visualizer:

    ## The constructor initializes the list of SkirtSimulation objects held by the new Visualizer instance.
    # The constructor determines this list depending on the type and value of the \em simulations argument:
    # - missing: all simulations for which there is a log file in the current directory
    # - empty string: all simulations for which there is a log file in the current directory
    # - string containing a slash: all simulations for which there is a log file in the specified directory
    # - nonempty string without a slash: the simulation with the specified prefix in the current directory
    # - simulation object: the simulation represented by the specified object
    # - list of strings and/or simulation objects: all simulations in the listed objects, as defined above
    #
    # If the \em log argument is \c True, a message is printed for each file created by this visualizer.
    # By default, no messages are printed.
    def __init__(self, simulations="", log=False):
        self.log = log
        self.simulations = [ ]
        sourcelist = simulations if isinstance(simulations, (types.TupleType,types.ListType)) else [ simulations ]
        for source in sourcelist:
            if isinstance(source, types.StringTypes):
                if source == "" or "/" in source:
                    dirpath = os.path.realpath(os.path.expanduser(source))
                    logfiles = filter(lambda fn: fn.endswith("_log.txt"), os.listdir(dirpath))
                    for logfile in logfiles:
                        self.simulations.append(SkirtSimulation(prefix=logfile[:-8], outpath=dirpath))
                else:
                    if os.path.exists(source+"_log.txt"):
                        self.simulations.append(SkirtSimulation(prefix=source))
            elif isinstance(source,SkirtSimulation):
                self.simulations.append(source)
            else:
                raise ValueError("Unsupported source type for simulation")

    # This function creates one or more RGB images for each "total.fits" file in the output of the simulations
    # held by this visualizer. If the \em wavelength_tuples argument is missing, a single image is created
    # for each "total.fits" file using the frames loaded by default by the RGBimage constructor. Otherwise,
    # the \em wavelength_tuples argument must contain a sequence of 3-tuples with (R,G,B) wavelengths; each of
    # these tuples causes an image to be created, loading the frames corresponding to the specified wavelengths.
    # The \em from_percentile and \em to_percentile arguments take the percentile values, in range [0,100], used
    # to clip the luminosity values loaded from the fits file.
    # The images are saved in PNG format and are placed next to the original file(s) with the same name
    # but a different extension. If there are multiple images per fits file, a serial number is added.
    def createRGBimages(self, wavelength_tuples=None, from_percentile=30, to_percentile=100):

        # loop over the wavelength tuples
        if wavelength_tuples==None: wavelength_tuples = [ None ]
        for index in range(len(wavelength_tuples)):

            # loop over the simulations
            for simulation in self.simulations:

                # get the frame indices corresponding to the requested wavelengths
                if wavelength_tuples[index] == None: frames = None
                else: frames = simulation.frameindices(wavelength_tuples[index])

                # get a list of output file names, including extension, one for each instrument
                outnames = simulation.totalfitspaths()
                if len(outnames) > 0:

                    # determine the appropriate pixel range for ALL output images for this galaxy
                    ranges = []
                    for outname in outnames:
                        im = RGBImage(outname, frames=frames)
                        ranges += list(im.percentilepixelrange(from_percentile,to_percentile))
                    rmin = min(ranges)
                    rmax = max(ranges)

                    # create an RGB file for each output file
                    for outname in outnames:
                        im = RGBImage(outname, frames=frames)
                        im.setrange(rmin,rmax)
                        im.applylog()
                        im.applycurve()
                        savename = outname[:-5] + (str(index+1) if index > 0 else "") + ".png"
                        im.saveto(savename)
                        if self.log: print "Created RGB image file " + savename

    # This function creates a plot for each "gridxx.dat" file in the output of the simulations held by
    # this visualizer. The plots are saved in PDF format and are placed next to the original file(s) with
    # the same name but a different extension.
    def plotgrids(self):
        for simulation in self.simulations:
            for name in simulation.gridxxdatpaths():
                plotgrid(name, name[:-4] + ".pdf")
                if self.log: print "Created PDF grid plot file " + name[:-4] + ".pdf"

    # This function creates a plot combining the "sed.dat" files in the output of each simulation held by
    # this visualizer. The plots are saved in PDF format and are placed next to the original file(s) with
    # a similar name (omitting the instrument name) but a different extension.
    # The function takes the following arguments:
    # - xlim: the lower and upper limits of the x axis, specified as a 2-tuple; if missing the x axis is auto-scaled
    # - ylim: the lower and upper limits of the y axis, specified as a 2-tuple; if missing the y axis is auto-scaled
    def plotseds(self, xlim=None, ylim=None):
        for simulation in self.simulations:
            sedpaths = simulation.seddatpaths()
            if len(sedpaths) > 0:
                labels = [ path.rsplit("_",2)[1] for path in sedpaths ]
                outpath = sedpaths[0].rsplit("_",2)[0] + "_sed.pdf"
                success = plotsed(sedpaths, outpath, labels, xlim=xlim, ylim=ylim)
                if success and self.log: print "Created PDF grid plot file " + outpath

    # This function creates a polarization map for each set of "prefix_instr_stokes*.fits" files in the output of
    # the simulations held by this visualizer. Each plot is saved in a PDF file named "prefix_instr_stokes.pdf",
    # placed next to the original set of files.
    def plotpolarizationmaps(self):
        for simulation in self.simulations:
            for nametuple in simulation.stokesfitspaths():
                plotname = nametuple[1].replace("Q.fits",".pdf")
                plotpolarizationmap(nametuple, plotname)
                if self.log: print "Created PDF polarization map" + plotname

    # This function creates a movie for the output of each simulation held by this visualizer. The movie combines
    # the SEDs (bottom panel) and the pixel frames (top panel, from left to right) for up to three instruments,
    # running through all wavelengths in the simulation. The movie is placed next to the original file(s) with
    # a similar name (omitting the instrument name) but a different extension.
    # The function takes the following arguments:
    #  - xlim, ylim:  the lower and upper limits of the x/y axis of the SED plot, specified as a 2-tuple;
    #                 if missing the corresponding axis is auto-scaled
    #  - from_percentile, to_percentile: the percentile values, in range [0,100], used to clip the luminosity values
    #                 loaded from the fits files; the default values are 30 and 100 respectively
    def makewavemovie(self, xlim=None, ylim=None, from_percentile=30, to_percentile=100):
        for simulation in self.simulations:
            sedpaths = simulation.seddatpaths()
            if 1 <= len(sedpaths) <= 3:
                fitspaths = simulation.totalfitspaths()
                wavelengths = simulation.wavelengths()
                if len(fitspaths) == len(sedpaths) and len(wavelengths) > 3:
                    labels = [ path.rsplit("_",2)[1] for path in sedpaths ]
                    outpath = sedpaths[0].rsplit("_",2)[0] + "_wave.mov"
                    success = makewavemovie(sedpaths, fitspaths, labels, outpath, wavelengths,
                                            xlim=xlim, ylim=ylim,
                                            from_percentile=from_percentile, to_percentile=to_percentile)
                    if success and self.log: print "Created wave movie file " + outpath

# -----------------------------------------------------------------
