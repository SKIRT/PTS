#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.skirtsimulation Representing the files related to a SKIRT simulation
#
# An instance of the SkirtSimulation class in this module represents all input and output files
# related to a single performed SKIRT simulation.

# -----------------------------------------------------------------

import os
import os.path
import types
import numpy as np
from pts.skifile import SkiFile
import pts.archive as arch

# -----------------------------------------------------------------

## This function returns a list of SkirtSimulation objects depending on the type and value of the \em source argument:
# - missing: all simulations for which there is a log file in the current directory
# - empty string: all simulations for which there is a log file in the current directory
# - string containing a slash: all simulations for which there is a log file in the specified directory
# - nonempty string without a slash: the simulation with the specified prefix in the current directory
# - simulation object: the simulation represented by the specified object
# - list of strings and/or simulation objects: all simulations in the listed objects, as defined above
#
def createsimulations(source=""):
    simulations = [ ]
    sourcelist = source if isinstance(source, (types.TupleType,types.ListType)) else [ source ]
    for source in sourcelist:
        if isinstance(source, types.StringTypes):
            if source == "" or "/" in source:
                dirpath = os.path.realpath(os.path.expanduser(source))
                logfiles = arch.listdir(dirpath, "_log.txt")
                for logfile in logfiles:
                    simulations.append(SkirtSimulation(prefix=logfile[:-8], outpath=dirpath))
            else:
                if os.path.exists(source+"_log.txt"):
                    simulations.append(SkirtSimulation(prefix=source))
        elif isinstance(source,SkirtSimulation):
            simulations.append(source)
        else:
            raise ValueError("Unsupported source type for simulation")
    return simulations

# -----------------------------------------------------------------
#  SkirtSimulation class
# -----------------------------------------------------------------

## An instance of the SkirtSimulation class represents all input and output files related to a single performed
# SKIRT simulation. To create an instance of the class, one specifies the name of the ski file (used as prefix
# for all output filenames) plus an input path and an output path. The current implementation only supports
# paths on the local file system; support for remote file systems may be added in the future.
#
# The methods of the class allow retrieving all kinds of high-level information about the simulation results.
#
# The code in this class uses ad-hoc knowledge about SKIRT's data formats and naming schemes; for example:
#  - SKIRT's builtin output filenames are used to locate various output files.
#  - the structure of the SKIRT parameter file (ski file) is used to determine things like the wavelengths used
#    in an oligochromatic simulation, the names of the instruments, or input filenames.
#  - the format of the log file is used to extract success or error messages.
#
# Combining this ad-hoc knowledge in a single class (as much as possible) will ease the pain of updating things
# if and when the SKIRT output schemes change.
#
class SkirtSimulation:

    # -----------------------------------------------------------------

    ## The constructor accepts the following arguments:
    # - prefix: the name of the ski file for which the simulation was performed (without the directory path)
    #   and which has been used as a prefix for all output filenames; if the prefix is empty (or missing) it is
    #   automatically derived from the files in the output directory, assuming the directory contains output from
    #   a single simulation (which means the simulation must already have been run, or at least started).
    # - inpath: the input path of the simulation; the path may be absolute, relative to a user's home folder,
    #   or relative to the current working directory. If there were no input files for the simulation, or if
    #   access to the input files is not needed, the inpath may be missing or empty.
    # - outpath: the output path of the simulation; the path may be absolute, relative to a user's home folder,
    #   or relative to the current working directory. A missing or empty outpath means the current working directory.
    #
    def __init__(self, prefix="", inpath="", outpath=""):
        self._inpath = os.path.realpath(os.path.expanduser(inpath))
        self._outpath = os.path.realpath(os.path.expanduser(outpath))
        self._prefix = prefix
        if self._prefix.endswith(".ski"):
            self._prefix = self._prefix[0:len(self._prefix)-len(".ski")]
        if self._prefix=="":
            logfiles = arch.listdir(self._outpath, "_log.txt")
            if len(logfiles) == 0: raise ValueError("No log file in path: " + self._outpath)
            if len(logfiles) > 1: raise ValueError("Multiple log files in path: " + self._outpath)
            self._prefix = logfiles[0][0:len(logfiles[0])-len("_log.txt")]

        # provide placeholders for caching frequently-used objects
        self._parameters = None
        self._units = None

    ## This function returns the simulation name, used as a prefix for output filenames
    def prefix(self):
        return self._prefix

    ## This function returns the absolute input path of the simulation
    def inpath(self):
        return self._inpath

    ## This function returns the absolute output path of the simulation
    def outpath(self):
        return self._outpath

    ## This function returns the absolute path for a simulation input file, given the file's name
    def infilepath(self, name):
        return os.path.join(self._inpath, name)

    ## This function returns the absolute path for a simulation output file, given the file's partial name
    # (the partial name does not include the prefix and the subsequent underscore).
    def outfilepath(self, partialname):
        return os.path.join(self._outpath, self._prefix + "_" + partialname)

    ## This function returns the absolute path for the simulation log file
    def logfilepath(self):
        return self.outfilepath("log.txt")

    # -----------------------------------------------------------------

    ## This function returns the status of the simulation, based on the contents of its log file, as one of the
    # following strings:
    #  - 'NotStarted': the log file does not exist, so the simulation has not been started
    #  - 'Running': the log file ends without a Finished or Error message, so it must still be running
    #               (unless the process was terminated without leaving an error message!)
    #  - 'Crashed': the log file ends with an Error message
    #  - 'Finished': the log file ends with a proper Finished message
    def status(self):
        logpath = self.logfilepath()

        # handle file no found
        if not arch.isfile(logpath): return "NotStarted"

        # get the last few lines of the file (assume the relevant portion is not longer than 500 characters)
        logfile = arch.openbinary(logpath)
        logfile.seek(0, os.SEEK_END)
        logfile.seek(-min(logfile.tell(),500), os.SEEK_END)
        chunk = logfile.read()
        logfile.close()
        lines = chunk.splitlines()
        last = lines[len(lines)-1] if len(lines)>0 else ""
        lastbutone = lines[len(lines)-2] if len(lines)>1 else ""

        # handle contents of the last lines
        if " Available memory: " in last: last = lastbutone
        if " Finished simulation " + self._prefix in last: return "Finished"
        if " *** Error: " in last: return "Crashed"
        return "Running"

    # -----------------------------------------------------------------

    ## This function returns a SkiFile object representing the parameter file for this simulation.
    def parameters(self):
        if self._parameters==None:
            self._parameters = SkiFile(self.outfilepath("parameters.xml"))
        return self._parameters

    ## This function returns a SkirtUnits object representing the default SKIRT units for this simulation.
    def units(self):
        if self._units==None:
            self._units = self.parameters().units()
        return self._units

    ## This function allows invoking any SkiFile or SkirtUnits function directly on a simulation object. For example,
    # self.instrumentshape() is automatically translated to self.parameters().instrumentshape(); and
    # self.convert() is automatically translated to self.units().convert().
    def __getattr__(self, attrname):
        # if this is not a system attribute
        if not attrname.startswith("__"):
            # attempt delegating to our SkiFile object
            try: return getattr(self.parameters(), attrname)
            except AttributeError: pass
            # attempt delegating to our SkirtUnits object
            try: return getattr(self.units(), attrname)
            except AttributeError: pass
        raise AttributeError("Can't delegate this attribute")

    # -----------------------------------------------------------------

    ## This function returns a list of absolute filepaths for all "total.fits" files produced by the simulation,
    # in the same order as the corresponding instruments occur in the ski file.
    def totalfitspaths(self):
        return [ self.outfilepath(name + "_total.fits") for name in self.parameters().instrumentnames() \
             if arch.isfile(self.outfilepath(name + "_total.fits")) ]

    ## This function returns a list of absolute filepath tuples for all sets of "stokes*.fits" files
    # produced by the simulation, in the same order as the corresponding instruments occur in the ski file.
    # Each returned tuple includes the four files paths corresponding to the components of the Stokes vector,
    # i.e. ("total.fits", "stokesQ.fits", "stokesU.fits", "stokesV.fits"), in that order.
    def stokesfitspaths(self):
        return [ ( self.outfilepath(name + "_total.fits"),
                   self.outfilepath(name + "_stokesQ.fits"),
                   self.outfilepath(name + "_stokesU.fits"),
                   self.outfilepath(name + "_stokesV.fits") ) for name in self.parameters().instrumentnames() \
             if arch.isfile(self.outfilepath(name + "_stokesQ.fits")) ]

    ## This function returns a list of absolute filepaths for all "sed.dat" files produced by the simulation,
    # in the same order as the corresponding instruments occur in the ski file.
    def seddatpaths(self):
        return [ self.outfilepath(name + "_sed.dat") for name in self.parameters().instrumentnames() \
             if arch.isfile(self.outfilepath(name + "_sed.dat")) ]

    ## This function returns a list of absolute filepaths for all "gridxx.dat" files produced by the simulation,
    # in the order xy, xz, yz, xyz.
    def gridxxdatpaths(self):
        return [ self.outfilepath(candidate) for candidate in \
                                ("ds_gridxy.dat", "ds_gridxz.dat", "ds_gridyz.dat", "ds_gridxyz.dat") \
                 if arch.isfile(self.outfilepath(candidate)) ]

    ## This function returns a list of the wavelengths used by the simulation, in micron, if available.
    # For an oligochromatic simulation, the wavelengths are obtained from the ski file.
    # For a panchromatic simulation, the wavelengths are read from the "wavelength.dat" file optionally
    # written by the WavelengthGrid class, or from one of the "sed.dat" files written by instruments.
    # If none of these files is present, the function raises an error.
    def wavelengths(self):
        # first try the ski file (for oligochromatic simulations)
        result = self.parameters().wavelengths()
        if len(result) > 0: return result

        # if that fails, try an SED data file or the wavelengths data file
        sedpaths = self.seddatpaths()
        if len(sedpaths) > 0:
            filepath = sedpaths[0]
        else:
            filepath = self.outfilepath("wavelengths.dat")
        if arch.isfile(filepath):
            result = np.loadtxt(arch.opentext(filepath), usecols=(0,))
            if len(result) > 0:
                return list( self.units().convert(result, to_unit='micron', quantity='wavelength') )

        # if everything fails, raise an error
        raise ValueError("Can't determine wavelengths for simulation")

    ## This function returns a list of the frame indices (in the simulation output fits files) corresponding
    # to each of the wavelengths in the specified list (expressed in micron). The function searches the simulation's
    # wavelength grid for the wavelength nearest to the requested value. It raises an error if the simulation
    # wavelengths are not available.
    def frameindices(self, wavelengths):
        # get the wavelength grid
        grid = np.array(self.wavelengths())
        # loop over the specified wavelengths
        return [ np.argmin(np.abs(grid-wave)) for wave in wavelengths ]

    # -----------------------------------------------------------------

    ## This function returns an appropriate axis label for the flux described in the simulation output sed files,
    # including a description of the physical quantity and the corresponding units. If there are no sed output files,
    # or if the units are not recognized, the function returns the string "Flux".
    def fluxlabel(self):
        # get the paths of the sed output files
        sedpaths = self.seddatpaths()
        if len(sedpaths)>0:
            # get the second line of the file, which contains the description of the flux column
            sedfile = arch.opentext(sedpaths[0])
            sedfile.readline()
            fluxdescription = sedfile.readline()
            sedfile.close()
            # select the appropriate label based on the units given in the description
            if "(W/m2)" in fluxdescription: return r"$\lambda\,F_\lambda\,(\mathrm{W}\,\mathrm{m}^{-2})$"
            if "(W/m3)" in fluxdescription: return r"$F_\lambda\,(\mathrm{W}\,\mathrm{m}^{-3})$"
            if "(W/m2/micron)" in fluxdescription: return r"$F_\lambda\,(\mathrm{W}\,\mathrm{m}^{-2}\,\mu \mathrm{m}^{-1})$"
            if "(W/m2/Hz)" in fluxdescription: return r"$F_\nu\,(\mathrm{W}\,\mathrm{m}^{-2}\,\mathrm{Hz}^{-1})$"
            if "(Jy)" in fluxdescription: return r"$F_\nu\,(\mathrm{Jy})$"
        # failed
        return "Flux"

    # -----------------------------------------------------------------

    ## This function retrieves a field from the specified simulation text output file, and returns its value
    # converted to the specified units. The field is located by a trigger (a text string that must occur on a line
    # before the one containing the field) and a header (a text string that must occur on the line containing
    # the field). The last text segment on the line represents the units of the value in the file, and the segment
    # before the units represents the value itself. The value is converted from the units in the file to the
    # requested units. If the function can't locate the field, it returns zero.
    def getfieldfromfile(self, filesuffix, trigger, header, units):
        filepath = self.outfilepath(filesuffix)
        if arch.isfile(filepath):
            triggered = False
            for line in arch.opentext(filepath):
                if trigger in line: triggered = True
                if triggered and header in line:
                    segments = line.split()
                    if len(segments)>2:
                        return self.units().convert(segments[-2], from_unit=segments[-1], to_unit=units)
        return 0;

    ## This function returns the total dust mass in the simulation's configuration space, in solar masses.
    # The function retrieves the 'expected' dust mass value listed in the convergence check data file.
    # It raises an error if the convergence check data file is not available or if the dust mass is zero.
    def dustmass(self):
        result = self.getfieldfromfile("ds_convergence.dat", "total dust mass", "expected value", "Msun")
        if result <= 0: raise ValueError("Can't determine dust mass")
        return result

    ## This function returns the total dust mass in the simulation's dust grid, in solar masses.
    # The function retrieves the 'actual' dust mass value listed in the convergence check data file.
    # It raises an error if the convergence check data file is not available or if the dust grid mass is zero.
    def dustgridmass(self):
        result = self.getfieldfromfile("ds_convergence.dat", "total dust mass", "actual value", "Msun")
        if result <= 0: raise ValueError("Can't determine dust grid mass")
        return result

    ## This function returns the total mass of the cold gass represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry
    # written by the SPH dust distribution. It raises an error if this entry is not found or if the mass is zero.
    def coldgasmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH gas", "Total gas mass", "Msun")
        if result <= 0: raise ValueError("Can't determine SPH cold gas mass")
        return result

    ## This function returns the total mass of the metallic gas represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry
    # written by the SPH dust distribution. It raises an error if this entry is not found or if the mass is zero.
    def metallicgasmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH gas", "Total metal mass", "Msun")
        if result <= 0: raise ValueError("Can't determine SPH metallic gas mass")
        return result

    ## This function returns the total initial stellar mass (i.e. the mass at the time of birth) represented by
    # the set of SPH particles imported for the simulation, in solar masses. The function retrieves this information
    # from the log file entry written by the SPH stellar component. It raises an error if this entry is not found or
    # if the mass is zero.
    def initialstellarmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH star", "Total stellar mass", "Msun")
        if result <= 0: raise ValueError("Can't determine SPH initial stellar mass")
        return result

    ## This function returns the total mass in hii regions represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry written
    # by the SPH starburst component. It raises an error if this entry is not found or if the mass is zero.
    def hiiregionmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH HII region", "Total stellar mass", "Msun")
        if result <= 0: raise ValueError("Can't determine SPH HII region mass")
        return result

    ## This function returns the total luminosity represented by the set of SPH stellar particles imported for
    # the simulation, in solar bolometric luminosity units. The function retrieves this information from the
    # log file entry written by the SPH stellar component. It raises an error if this entry is not found or if
    # the total stellar luminosity is zero.
    def stellarluminosity(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH star", "Total luminosity", "Lsun")
        if result <= 0: raise ValueError("Can't determine SPH stellar luminosity")
        return result

    ## This function returns the total luminosity represented by the set of SPH region II particles imported for
    # the simulation, in solar bolometric luminosity units. The function retrieves this information from the
    # log file entry written by the SPH starburst component. It raises an error if this entry is not found or if
    # the total stellar luminosity is zero.
    def hiiregionluminosity(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH HII region", "HII luminosity", "Lsun")
        if result <= 0: raise ValueError("Can't determine SPH HII region luminosity")
        return result

    # -----------------------------------------------------------------

    ## This function returns a tuple with three relevant properties of the dust grid used in the simulation:
    # the total number of dust cells, the largest optical depth for any one cell, and the optical depth
    # at the 90% percentile point (i.e. 90% of the cells that actually contain dust have an optical depth below
    # this value). The function retrieves this information from the dust cell properties data file optionally
    # written by the dust system. It raises an error if this file is not found.
    def dustcellstats(self):
        # load the optical depths from the file
        filepath = self.outfilepath("ds_cellprops.dat")
        depths = np.loadtxt(arch.opentext(filepath), usecols=(3,))
        # calculate and return the statistics
        nonzerodepths = depths[depths>0]
        if len(nonzerodepths) > 0:
            return ( len(depths), np.amax(nonzerodepths), np.percentile(nonzerodepths, 90) )
        else:
            return ( len(depths), 0, 0 )

    # -----------------------------------------------------------------

    # ============ re-evaluate below this line

    ## This function returns the monochromatic stellar flux \f$\lambda F_{*,\lambda}\f$ at the wavelength
    # \f$\lambda\f$ (specified in micron) and at the distance of the simulation's instruments, in W/m2.
    # The stellar flux is defined as the flux generated by the set of SPH particles imported for the simulation
    # without taking into account the presence of dust (i.e. the "transparent flux").
    # The implementation requires the presence of the output files "wavelengths.dat" and "luminosities.dat",
    # and it assumes that the wavelengths in these files are specified in micron and luminosities in \f$L_\odot\f$.
    def starlambdaflux(self, wavelength):
        # get the wavelengths and the bin widths for the simulation's wavelength grid
        datapath = self.outfilepath("wavelengths.dat")
        if not arch.isfile(datapath): raise ValueError("Can't locate wavelengths file")
        wavelengths, binwidths = np.loadtxt(arch.opentext(datapath), usecols=(0,1), unpack=True)

        # get the stellar luminosities, which should be sampled on the same wavelength grid
        datapath = self.outfilepath("luminosities.dat")
        if not arch.isfile(datapath): raise ValueError("Can't locate stellar luminosities file")
        luminosities = np.loadtxt(arch.opentext(datapath), usecols=(1,))
        if len(luminosities) != len(wavelengths): raise ValueError("Luminosities not sampled on wavelength grid")

        # convert luminosity (in Lsun) to monochromatic flux (in W/m2)
        distance = pc * self.instrumentdistance()
        fluxes = (Lsun / (4.*np.pi*distance*distance)) * luminosities * wavelengths / binwidths
        return self.interpolate(wavelength, wavelengths, fluxes)

    ## This function returns the monochromatic flux \f$\lambda F_{\lambda}\f$ at the wavelength
    # \f$\lambda\f$ (specified in micron) received by the specified instrument, in W/m2.
    # The implementation requires the presence of the output file "prefix_name_sed.dat",
    # and it assumes that the wavelengths in this file are specified in micron and the fluxes in W/m2.
    def lambdaflux(self, name, wavelength):
        wavelength = wavelengthforband(wavelength)
        # get the wavelengths and fluxes, and interpolate
        datapath = self.outfilepath(name + "_sed.dat")
        if not arch.isfile(datapath): raise ValueError("Can't locate sed file: " + datapath)
        wavelengths, fluxes = np.loadtxt(arch.opentext(datapath), usecols=(0,1), unpack=True)
        return self.interpolate(wavelength, wavelengths, fluxes)

    ## This function returns the integrated flux \f[\int_{\lambda_1}^{\lambda_2} F_{\lambda}\,\text{d}\lambda\f]
    # over a wavelength range \f$[\lambda_1,\lambda_2]\f$ (specified in micron) received by the specified instrument,
    # in W/m2. The function simply sums the fluxes for all wavelength bins touched by the specified wavelength range,
    # so the result is accurate only when the wavelength range spans many bins in the simulation's wavelength grid.
    # The implementation requires the presence of the output files "wavelengths.dat" and "prefix_name_sed.dat",
    # and it assumes that the wavelengths in these files are specified in micron and the fluxes in W/m2.
    def integratedflux(self, name, wavelength1, wavelength2):
        # get the wavelengths and the bin widths for the simulation's wavelength grid
        datapath = self.outfilepath("wavelengths.dat")
        if not arch.isfile(datapath): raise ValueError("Can't locate wavelengths file")
        wavelengths, binwidths = np.loadtxt(arch.opentext(datapath), usecols=(0,1), unpack=True)

        # get the fluxes, which should be sampled on the same wavelength grid
        datapath = self.outfilepath(name + "_sed.dat")
        if not arch.isfile(datapath): raise ValueError("Can't locate sed file: " + datapath)
        fluxes = np.loadtxt(arch.opentext(datapath), usecols=(1,))
        if len(fluxes) != len(wavelengths): raise ValueError("Fluxes not sampled on wavelength grid")

        # locate the first and last wavelength bin touching the specified range
        leftbounds = wavelengths-0.5*binwidths
        index1 = max(0, np.searchsorted(leftbounds, wavelength1)-1 )
        index2 = max(0, np.searchsorted(leftbounds, wavelength2)-1 )

        # perform the integration
        integrated = binwidths*fluxes/wavelengths
        return integrated[index1:index2+1].sum()

    ## This function returns the luminosity in \f$L_\odot\f$ corresponding to the specified flux in W/m2,
    # assuming the source is located at the distance of the simulation's instruments.
    def luminosityforflux(self, flux):
        distance = pc * self.instrumentdistance()
        return flux * 4.*np.pi*distance*distance / Lsun

    # -----------------------------------------------------------------

    ## This function returns the flux density \f$F_{\lambda}\f$ at the wavelength
    # \f$\lambda\f$ (specified in micron) received by the specified instrument, in W/m2/m.
    # The implementation requires the presence of the output file "prefix_name_sed.dat",
    # and it assumes that the wavelengths in this file are specified in micron and the fluxes in W/m2.
    def fluxdensityWave(self, name, wavelength):
        wavelength = wavelengthforband(wavelength)
        return self.lambdaflux(name, wavelength) / (wavelength*1e-6)

    ## This function returns the flux density \f$F_{\nu}\f$ at the wavelength
    # \f$\lambda\f$ (specified in micron) received by the specified instrument, in Yansky. We have
    # \f$\lambda F_\lambda=\nu F_\nu\f$ and \f$\lambda\nu=c\f$ so that \f$F_\nu=(\lambda F_\lambda)\times\lambda/c\f$.
    # And finally, \f$1\,\textrm{Jy}=10^{-26}\,\textrm{W}\textrm{m}^{-2}\textrm{Hz}^{-1}\f$.
    # The implementation requires the presence of the output file "prefix_name_sed.dat",
    # and it assumes that the wavelengths in this file are specified in micron and the fluxes in W/m2.
    def fluxdensityFreq(self, name, wavelength):
        wavelength = wavelengthforband(wavelength)
        lambdaflux = self.lambdaflux(name, wavelength)
        return 1e26*lambdaflux*wavelength*1e-6/c

    ## This function returns the weighed flux density in W/m2/m received by the specified instrument, calculated as
    # \f[ \frac{\sum_i w(\lambda_i) F_{\lambda}(\lambda_i)} {\sum_i w(\lambda_i)} \f] where the list of
    # wavelengths \f$\lambda_i\f$ (in micron) and the list of weights \f$w(\lambda_i)\f$ are specified as arguments.
    # The implementation requires the presence of the output file "prefix_name_sed.dat",
    # and it assumes that the wavelengths in this file are specified in micron and the fluxes in W/m2.
    def weighedfluxdensityWave(self, name, wavelengths, weights):
        result = 0.
        for wavelength, weight in zip(wavelengths, weights):
            result +=  weight * self.fluxdensityWave(name, wavelength)
        return result / weights.sum()

    ## This function returns the luminosity density, in W/Hz, at the wavelength
    # \f$\lambda\f$ (specified in micron) corresponding to the flux received by the specified instrument
    # assuming the source is located at the distance of the simulation's instruments.
    # The implementation requires the presence of the output file "prefix_name_sed.dat",
    # and it assumes that the wavelengths in this file are specified in micron and the fluxes in W/m2.
    def luminositydensityFreq(self, name, wavelength):
        flux = 1e-26*self.fluxdensityFreq(name, wavelength)
        distance = pc * self.instrumentdistance()
        return flux * 4.*np.pi*distance*distance

    ## This function returns the absolute AB magnitude of the modeled source as seen by the specified instrument
    # at the wavelength \f$\lambda\f$ (specified in micron).
    # Given a flux density \f$F_\nu\f$, measured in ergs per second per square cm per Hz, the corresponding
    # AB magnitude is defined as \f$\text{AB}=-2.5\log_{10} F_\nu -48.60\f$. Furthermore, we have
    # \f$\lambda F_\lambda=\nu F_\nu\f$ and \f$\lambda\nu=c\f$ so that \f$F_\nu=(\lambda F_\lambda)\times\lambda/c\f$,
    # and also \f$F_\nu^{(cgs)}=10^3 F_\nu^{(SI)}\f$. Finally, the resulting apparent magnitude is converted to the
    # absolute magnitude using the standard formula \f$M=m-5\log_{10}d^{(pc)}+5\f$.
    def absmagnitude(self, name, wavelength):
        wavelength = wavelengthforband(wavelength)
        lambdaflux = self.lambdaflux(name, wavelength)
        apparent = -2.5*np.log10(1e3*lambdaflux*wavelength*1e-6/c) - 48.60
        distance = self.instrumentdistance()
        absolute = apparent - 5*np.log10(distance) + 5
        return absolute

    ## This function returns the specified color index of the modeled source as seen by the specified instrument.
    # for the specified wavelength \f$\lambda\f$ (specified in micron). The desired color index is given as a string
    # with two segments separated by a hyphen. Each segment specifies a standard band (such as r, V or NUV) or a
    # wavelength in micron. Spaces are not allowed. Examples: "g-r", "U-V", "NUV-r", "250-350".
    # The color index is defined as the difference between the AB magnitudes for the two bands/wavelengths.
    def color(self, name, colorspec):
        segments = colorspec.split("-")
        if len(segments) != 2: raise ValueError("Invalid color spec: " + colorspec)
        return self.absmagnitude(name, segments[0]) - self.absmagnitude(name, segments[1])

# -----------------------------------------------------------------

c = 2.99792458e8     # light speed in m/s
pc = 3.08568025e16   # parsec in m
Lsun = 3.845e26      # solar bolometric luminosity in W
Msun = 1.98892e30    # solar mass in kg

bands = {            # central wavelength in micron for standard bands
          'FUV': 0.152, 'NUV': 0.231,                                           # GALEX
          'U': 0.365,   'B': 0.445,   'V': 0.551,   'R': 0.658,   'I': 0.806,   # Johnson
          'Z': 0.882,   'Y': 1.031,   'J': 1.248,   'H': 1.631,   'K': 2.201,   # UKIDSS
          'u': 0.3551,  'g': 0.4686,  'r': 0.6165,  'i': 0.7481,  'z': 0.8931,  # SDSS ugriz
        }

## This function returns the central wavelength in micron for the usual GALEX, Johnson, UKIDSS and SDSS ugriz bands,
# specified as a string such as "V", "FUV" or "r". If the argument does not specify a known band, it is converted
# to a floating point number and returned (i.e. it is interpreted as the wavelength in micron).
def wavelengthforband(bandspec):
    if bandspec in bands: return bands[bandspec]
    return float(bandspec)

# -----------------------------------------------------------------
