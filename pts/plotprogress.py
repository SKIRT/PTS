#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.plotprogress Plot progress for the various phases of a SKIRT simulation
#
# The function in this module plots the progress in function of time for certain phases of a SKIRT simulation,
# based on the log messages. A seperate PDF plot is created for each of the following phases:
# - shooting photons for stellar emission;
# - calculating dust emission spectra;
# - shooting photons for dust emission.
#

# -----------------------------------------------------------------

# use a non-interactive back-end to generate high-quality vector graphics
import matplotlib
if matplotlib.get_backend().lower() != "pdf": matplotlib.use("pdf")
import matplotlib.pyplot as plt

from datetime import datetime
import pts.archive as arch

# -----------------------------------------------------------------

# This function plots the progress in function of time for certain phases of a SKIRT simulation, based on the log
# messages. The plots are saved in PDF format and are placed next to the original file(s) with a similar name.
# A seperate PDF plot is created for each of the following phases, if present in the simulation:
# - shooting photons for stellar emission ("prefix_progress_stellar_photons.pdf");
# - calculating dust emission spectra ("prefix_progress_dust_spectra.pdf");
# - shooting photons for dust emission ("prefix_progress_dust_photons.pdf");
#
# The dust self-absorption phase, if present, is ignored in the current implementation.
#
# For multi-process (MPI) simulations with verbose logging (i.e. with a separate log file per process),
# the progress for all processes is displayed on the same plot.
#
# The function takes the following arguments:
# - simulation: the SkirtSimulation object representing the simulation to be handled
# - figsize: the horizontal and vertical size of the output figure in inch (!); default is 10 x 6 inch
#
def plotprogress(simulation, figsize=(10,6)):
    _plot_photon_progress(simulation, figsize, 'stellar')
    _plot_spectra_progress(simulation, figsize)
    _plot_photon_progress(simulation, figsize, 'dust')

# -----------------------------------------------------------------

# this private function plots progress of the indicated photon shooting phase ('stellar' or 'dust')
def _plot_photon_progress(simulation, figsize, phase):

    # setup the figure
    figure = plt.figure(figsize=figsize)
    numplots = 0

    # gather and add the data points for each logfile
    for logfile in simulation.logfilepaths():
        percentages = [ ]
        times = [ ]
        triggered = False
        for line in open(logfile, 'r'):
            if not triggered and "Starting the {} emission phase".format(phase) in line:
                triggered = True
                process = line.split()[2][1:5] if "[P" in line else "P000"
            if triggered and " photon packages for " in line:
                startline = line
            if triggered and " photon packages: " in line and "Launched" in line:
                percentages.append( float(line.split()[-1][:-1]) )
                times.append( _timelapse(startline,line) )
            if triggered and " Finished " in line:
                break
        if len(percentages)>1:
            plt.plot(percentages, times, label=process)
            numplots += 1

    # if we actually plotted something, set appropriate limits and labels, and save the figure
    if numplots>0:
        plt.xlim(0,100)
        plt.ylim(0)
        plt.grid('on')
        plt.xlabel("progress (%)", fontsize='large')
        plt.ylabel("time (s)", fontsize='large')
        plt.title("Progress of emitting {} photons".format(phase))
        plt.legend(loc='lower right', ncol=4, prop={'size':8})

        outpath = simulation.outfilepath("progress_{}_photons.pdf".format(phase))
        plt.savefig(outpath, bbox_inches='tight', pad_inches=0.25)
        print "Created PDF progress plot file " + outpath

    # finalize things
    plt.close()

# -----------------------------------------------------------------

# this private function plots progress of the dust spectra calculation
def _plot_spectra_progress(simulation, figsize):

    # setup the figure
    figure = plt.figure(figsize=figsize)
    numplots = 0

    # gather and add the data points for each logfile
    staggered = simulation.staggered()
    logfiles = simulation.logfilepaths()
    numprocesses = len(logfiles)
    for rank in range(numprocesses):
        percentages = [ ]
        times = [ ]
        triggered = False
        for line in open(logfiles[rank], 'r'):
            if not triggered and "Starting the dust emission phase" in line:
                triggered = True
                process = line.split()[2][1:5] if "[P" in line else "P000"
            if triggered and " Calculating dust emission spectra" in line:
                startline = line
            if triggered and " Library entries in use: " in line:
                totalentries = float(line.split()[-1])
                entriesperprocess = totalentries/numprocesses
            if triggered and " Calculating emission for " in line:
                entry = float(line.split()[-1][:-3])
                if staggered:
                    fraction = entry/totalentries
                else:
                    fraction = (entry-rank*entriesperprocess)/entriesperprocess
                percentages.append(100*fraction)
                times.append(_timelapse(startline,line))
            if triggered and " Dust emission spectra calculated" in line:
                break
        if len(percentages)>1:
            plt.plot(percentages, times, label=process)
            numplots += 1

    # if we actually plotted something, set appropriate limits and labels, and save the figure
    if numplots>0:
        plt.xlim(0,100)
        plt.ylim(0)
        plt.grid('on')
        plt.xlabel("progress (%)", fontsize='large')
        plt.ylabel("time (s)", fontsize='large')
        plt.title("Progress of calculating dust spectra")
        plt.legend(loc='lower right', ncol=4, prop={'size':8})

        outpath = simulation.outfilepath("progress_dust_spectra.pdf")
        plt.savefig(outpath, bbox_inches='tight', pad_inches=0.25)
        print "Created PDF progress plot file " + outpath

    # finalize things
    plt.close()

# -----------------------------------------------------------------

# this private helper function returns the datetime object corresponding to the time stamp in a line
def _timestamp(line):
    date, time, dummy = line.split(None, 2)
    day, month, year = date.split('/')
    hour, minute, second = time.split(':')
    second, microsecond = second.split('.')
    return datetime(year=int(year), month=int(month), day=int(day),
                    hour=int(hour), minute=int(minute), second=int(second), microsecond=int(microsecond))

# this private helper function returns the difference in seconds between the time stamps in the two lines
def _timelapse(line1, line2):
    return (_timestamp(line2)-_timestamp(line1)).total_seconds()

# -----------------------------------------------------------------
