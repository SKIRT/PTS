#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.simulation.simulation Representing the files related to a SKIRT simulation
#
# An instance of the SkirtSimulation class in this module represents all input and output files
# related to a single performed SKIRT simulation.

# -----------------------------------------------------------------

# Import standard modules
import os
import os.path
import types
import numpy as np
import importlib
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools import serialization
from ..tools import filesystem as fs
from .skifile import SkiFile
from .logfile import LogFile
from ..tools import archive as arch
from ..launch.options import AnalysisOptions

# -----------------------------------------------------------------

## This function returns a list of SkirtSimulation objects depending on the type and value of the \em source argument:
# - missing: all simulations for which there is a log file in the current directory
# - empty string: all simulations for which there is a log file in the current directory
# - string containing a slash: all simulations for which there is a log file in the specified directory
# - nonempty string without a slash: the simulation with the specified prefix in the current directory
# - simulation object: the simulation represented by the specified object
# - list of strings and/or simulation objects: all simulations in the listed objects, as defined above
#
def createsimulations(source="", single=False):
    simulations = []
    sourcelist = source if isinstance(source, (types.TupleType,types.ListType)) else [ source ]
    for source in sourcelist:
        if isinstance(source, types.StringTypes):
            if source == "" or "/" in source:
                dirpath = os.path.realpath(os.path.expanduser(source))
                logfiles = arch.listdir(dirpath, "_log.txt")
                for logfile in logfiles:
                    prefix = logfile[:-8]
                    ski_path = fs.join(dirpath, prefix + ".ski")
                    ski_path = ski_path if fs.is_file(ski_path) else None
                    simulations.append(SkirtSimulation(prefix=logfile[:-8], outpath=dirpath, ski_path=ski_path))
            else:
                if os.path.exists(source + "_log.txt"):
                    ski_path = source + ".ski"
                    ski_path = ski_path if fs.is_file(ski_path) else None
                    simulations.append(SkirtSimulation(prefix=source, ski_path=ski_path))
        elif isinstance(source, SkirtSimulation):
            simulations.append(source)
        else:
            raise ValueError("Unsupported source type for simulation")

    # If a single simulation is expected
    if single:

        if not len(simulations) == 1: raise ValueError("Multiple simulations were found")
        else: return simulations[0]

    # If multiple simulations are expected
    else: return simulations

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
class SkirtSimulation(object):

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
    def __init__(self, prefix="", inpath="", outpath="", ski_path=None, parameters=None, name=None):

        # Set the full path to the input directory or the paths to the input files
        if inpath is not None:
            if isinstance(inpath, list): self._inpath = inpath
            else: self._inpath = os.path.realpath(os.path.expanduser(inpath))
        else: self._inpath = None
        self._outpath = os.path.realpath(os.path.expanduser(outpath if outpath is not None else ""))
        self._prefix = prefix
        if self._prefix.endswith(".ski"):
            self._prefix = self._prefix[0:len(self._prefix)-len(".ski")]
        if self._prefix == "":
            logfiles = arch.listdir(self._outpath, "_log.txt")
            if len(logfiles) == 0: raise ValueError("No log file in path: " + self._outpath)
            if len(logfiles) > 1: raise ValueError("Multiple log files in path: " + self._outpath)
            self._prefix = logfiles[0][0:len(logfiles[0])-len("_log.txt")]

        self.ski_path = ski_path
        self.input_path = self._inpath
        self.output_path = self._outpath

        if self.ski_path is None: self.ski_path = self.outfilepath("parameters.xml")

        # The base path = the directory where the ski file is located
        self.base_path = fs.directory_of(self.ski_path)

        # provide placeholders for caching frequently-used objects
        self._parameters = parameters
        self._units = None
        self._processes = None
        self._threads = None

        # A name given to the simulation
        self._name = name

        # The options for analysing the simulation output
        self.analysis = AnalysisOptions()

    ## This function returns whether the simulation requires input
    @property
    def has_input(self):
        return self.input_path is not None

    ## This function returns the simulation name, used as a prefix for output filenames
    def prefix(self):
        return self._prefix

    ## This function returns the simulation name, which is the prefix is a name was not set
    @property
    def name(self):
        return self._name if self._name is not None else self.prefix()

    ## This functions sets the name
    @name.setter
    def name(self, value):
        self._name = value

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

    ## This function returns a LogFile object created from the simulation's log file
    @property
    def log_file(self):
        return LogFile(self.outfilepath("log.txt"))

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

    ## This function returns the number of processes used for this simulation
    def processes(self):
        if self._processes is None:
            with open(self.logfilepath()) as logfile:
                for line in logfile:
                    if "Starting simulation" in line:
                        if "with" in line:
                            self._processes = int(line.split(' with ')[1].split()[0])
                        else: self._processes = 1
            if self._processes is None: raise ValueError("Cannot determine the number of processes from the log file")
        return self._processes

    ## This function returns the number of threads used for this simulation
    def threads(self):
        if self._threads is None:
            with open(self.logfilepath()) as logfile:
                triggered = False
                max_thread_number = 0
                for line in logfile:
                    if "Initializing random number generator" in line:
                        triggered = True
                        max_thread_number = int(line.split("thread number ")[1].split(" with seed")[0])
                    elif triggered: self._threads = max_thread_number+1
            if self._threads is None: raise ValueError("Cannot determine the number of threads from the log file")
        return self._threads

    ## This function returns a SkiFile object representing the parameter file for this simulation.
    def parameters(self):
        if self._parameters is None:
            self._parameters = SkiFile(self.ski_path)
        return self._parameters

    @property
    def ski_file(self):
        return self.parameters()

    ## This function returns a SkirtUnits object representing the default SKIRT units for this simulation.
    def units(self):
        if self._units is None:
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

    ## This function returns a list of absolute filepaths for all log files produced by the simulation, including
    # the master log file and any log files produced by parallel (MPI) processes. The list includes only paths for
    # log files that actually exist, and the paths are listed in order of process rank.
    def logfilepaths(self):
        logname = self._prefix + "_log"
        logfiles = sorted(filter(lambda fn: fn.startswith(logname), arch.listdir(self._outpath,".txt")))
        return [ os.path.join(self._outpath, logfile) for logfile in logfiles ]

    ## This function returns a list of LogFile objects, created from the log files produced by the simulation.
    def logfiles(self):
        return [LogFile(logfilepath) for logfilepath in self.logfilepaths()]

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
    # requested units. If the function can't locate the field, it returns -1.
    def getfieldfromfile(self, filesuffix, trigger, header, units=None):
        filepath = self.outfilepath(filesuffix)
        if arch.isfile(filepath):
            triggered = False
            for line in arch.opentext(filepath):
                if trigger in line: triggered = True
                if triggered and header in line:
                    segments = line.split()
                    if len(segments)>2:
                        if units!=None:
                            return self.units().convert(segments[-2], from_unit=segments[-1], to_unit=units)
                        else:
                            return float(segments[-1])
        return -1

    ## This function returns the total dust mass in the simulation's configuration space, in solar masses.
    # The function retrieves the 'expected' dust mass value listed in the convergence check data file.
    # It raises an error if the convergence check data file is not available or if the dust mass is zero.
    def dustmass(self):
        result = self.getfieldfromfile("ds_convergence.dat", "total dust mass", "expected value", "Msun")
        if result < 0: raise ValueError("Can't determine dust mass")
        return result

    ## This function returns the total dust mass in the simulation's dust grid, in solar masses.
    # The function retrieves the 'actual' dust mass value listed in the convergence check data file.
    # It raises an error if the convergence check data file is not available or if the dust grid mass is zero.
    def dustgridmass(self):
        result = self.getfieldfromfile("ds_convergence.dat", "total dust mass", "actual value", "Msun")
        if result < 0: raise ValueError("Can't determine dust grid mass")
        return result

    ## This function returns the total mass of the cold gass represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry
    # written by the SPH dust distribution. It raises an error if this entry is not found.
    def coldgasmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH gas", "Total gas mass", "Msun")
        if result < 0: raise ValueError("Can't determine SPH cold gas mass")
        return result

    ## This function returns the total mass of the metallic gas represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry
    # written by the SPH dust distribution. It raises an error if this entry is not found.
    def metallicgasmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH gas", "Total metal mass", "Msun")
        if result < 0: raise ValueError("Can't determine SPH metallic gas mass")
        return result

    ## This function returns the number of "cold" SPH gas particles imported for the simulation, i.e. the gas
    # particles that actually contain dust. The function retrieves this information from the log file entry
    # written by the SPH dust distribution. It raises an error if this entry is not found.
    def coldgasparticles(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH gas", "gas particles containing dust")
        if result < 0: raise ValueError("Can't determine SPH cold gas particles")
        return int(result)

    ## This function returns the total initial stellar mass (i.e. the mass at the time of birth) represented by
    # the set of SPH particles imported for the simulation, in solar masses. The function retrieves this information
    # from the log file entry written by the SPH stellar component. It raises an error if this entry is not found or
    # if the mass is zero.
    def initialstellarmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH star", " mass:", "Msun")
        if result < 0: raise ValueError("Can't determine SPH initial stellar mass")
        return result

    ## This function returns the total mass in hii regions represented by the set of SPH particles imported for
    # the simulation, in solar masses. The function retrieves this information from the log file entry written
    # by the SPH starburst component. It raises an error if this entry is not found or if the mass is zero.
    def hiiregionmass(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH HII region", " mass:", "Msun")
        if result < 0: raise ValueError("Can't determine SPH HII region mass")
        return result

    ## This function returns the total luminosity represented by the set of SPH stellar particles imported for
    # the simulation, in solar bolometric luminosity units. The function retrieves this information from the
    # log file entry written by the SPH stellar component. It raises an error if this entry is not found or if
    # the total stellar luminosity is zero.
    def stellarluminosity(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH star", " luminosity:", "Lsun")
        if result < 0: raise ValueError("Can't determine SPH stellar luminosity")
        return result

    ## This function returns the total luminosity represented by the set of SPH HII region particles imported for
    # the simulation, in solar bolometric luminosity units. The function retrieves this information from the
    # log file entry written by the SPH starburst component. It raises an error if this entry is not found or if
    # the total stellar luminosity is zero.
    def hiiregionluminosity(self):
        result = self.getfieldfromfile("log.txt", "Reading SPH HII region", " luminosity:", "Lsun")
        if result < 0: raise ValueError("Can't determine SPH HII region luminosity")
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

    ## This function returns a numpy array with the wavelengths used by the simulation, if available.
    # The wavelengths are given in micron, and are sorted in increasing order.
    # For an oligochromatic simulation, the wavelengths are obtained from the ski file.
    # For a panchromatic simulation, the wavelengths are read from the "wavelength.dat" file optionally
    # written by the WavelengthGrid class, or from one of the "sed.dat" files written by instruments.
    # If none of these files is present, the function raises an error.
    def wavelengths(self):
        # first try the ski file (for oligochromatic simulations)
        result = self.parameters().wavelengths()
        if len(result) > 0: return np.sort(result)

        # if that fails, try an SED data file or the wavelengths data file
        sedpaths = self.seddatpaths()
        if len(sedpaths) > 0:
            filepath = sedpaths[0]
        else:
            filepath = self.outfilepath("wavelengths.dat")
        if arch.isfile(filepath):
            result = np.loadtxt(arch.opentext(filepath), usecols=(0,))
            if len(result) > 0:
                return self.units().convert(result, to_unit='micron', quantity='wavelength')

        # if everything fails, raise an error
        raise ValueError("Can't determine wavelengths for simulation")

    ## This function returns a list of the frame indices (in the simulation output fits files) corresponding
    # to each of the wavelengths in the specified list (expressed in micron). The function searches the simulation's
    # wavelength grid for the wavelength nearest to the requested value. It raises an error if the simulation
    # wavelengths are not available.
    def frameindices(self, wavelengths):
        # get the wavelength grid
        grid = self.wavelengths()
        # loop over the specified wavelengths
        return [ np.argmin(np.abs(grid-wave)) for wave in wavelengths ]

    ## This function returns a numpy array representing the wavelength grid used by the simulation, including both
    # the wavelength bin centers and the corresponding bin widths. The returned array has a shape of (2,N) where N
    # is the number of wavelength bins, so that it can be easily unpacked in two separate arrays. The returned values
    # are expressed in the specified units ('micron' by default). The function requires the presence of the
    # "wavelength.dat" file. If this file is not found, the function raises an error.
    def wavelengthbins(self, unit='micron'):
        filepath = self.outfilepath("wavelengths.dat")
        result = np.loadtxt(arch.opentext(filepath), usecols=(0,1), unpack=True)
        return self.units().convert(result, to_unit=unit, quantity='wavelength')

    ## This function returns a numpy array containing the luminosities of the SPH stellar particles imported for
    # the simulation, specified for each wavelength bin. Thus the returned array is expected to have the same length
    # as the wavelength grid. The returned values are expressed in the specified units ('Lsun' by default).
    # The function requires the presence of the "luminosities.dat" file. If this file is not found, the function
    # raises an error.
    def stellarluminosities(self, unit='Lsun'):
        filepath = self.outfilepath("star_luminosities.dat")
        result = np.loadtxt(arch.opentext(filepath), usecols=(1,))
        return self.units().convert(result, to_unit=unit, quantity='luminosity')

    ## This function returns a numpy array containing the luminosities of the SPH HII region particles imported for
    # the simulation, specified for each wavelength bin. Thus the returned array is expected to have the same length
    # as the wavelength grid. The returned values are expressed in the specified units ('Lsun' by default).
    # The function requires the presence of the "luminosities.dat" file. If this file is not found, the function
    # raises an error.
    def hiiregionluminosities(self, unit='Lsun'):
        filepath = self.outfilepath("hii_luminosities.dat")
        result = np.loadtxt(arch.opentext(filepath), usecols=(1,))
        return self.units().convert(result, to_unit=unit, quantity='luminosity')

    ## This function returns a numpy array containing the flux densities received by the specified instrument, at
    # each wavelength bin center. Thus the returned array is expected to have the same length as the wavelength grid.
    # The returned values are expressed in the specified units ('Jy' by default), with automatic conversion between
    # flux styles (neutral, wavelength, frequency) as needed. The function requires the presence of the
    # "prefix_name_sed.dat" file. If this file is not found, the function raises an error.
    def fluxdensities(self, name, unit='Jy'):
        filepath = self.outfilepath(name + "_sed.dat")
        wavelengths, fluxes = np.loadtxt(arch.opentext(filepath), usecols=(0,1), unpack=True)
        wavelengths = self.units().convert(wavelengths, to_unit='micron', quantity='wavelength')
        return self.units().convert(fluxes, to_unit=unit, quantity='fluxdensity', wavelength=wavelengths)

    ## This function allows setting the analysis options from a dictionary (or an actual AnalysisOptions object)
    def set_analysis_options(self, options):

        # If the options is an actual AnalysisOptions object, set the analysis attribute directly
        if isinstance(options, AnalysisOptions): self.analysis = options

        # Load the options into the AnalysisOptions object
        else: self.analysis.set_options(options)

    ## This function checks the analysis options, adapts them if necessary and also adapts the logging options if possible
    def check_analysis_options(self, logging_options=None):
        self.analysis.check(logging_options)

    @property
    def from_batch(self):
        return self.analysis.timing_table_path is not None or self.analysis.memory_table_path is not None

    @property
    def from_scaling_test(self):
        return self.analysis.scaling_run_name is not None

    @property
    def from_modeling(self):
        return self.analysis.modeling_path is not None

# -----------------------------------------------------------------

class RemoteSimulation(SkirtSimulation):

    """
    This class ...
    """

    def __init__(self, ski_path, input_path, output_path):

        """
        The constructor ...
        :return:
        """

        # Determine the simulation prefix
        prefix = fs.strip_extension(fs.name(ski_path))

        # Call the constructor of the base class
        super(RemoteSimulation, self).__init__(prefix, input_path, output_path, ski_path)

        # -- Attributes --

        # Properties of the remote host on which the simulation was run
        self.host_id = None
        self.cluster_name = None

        # The simulation file path
        self.path = None

        # Basic properties
        self.id = None
        self.remote_ski_path = None
        self.remote_simulation_path = None
        self.remote_input_path = None
        self.remote_output_path = None
        self.submitted_at = None

        # Options for retrieval
        self.retrieve_types = None

        # The parallelization properties
        self.parallelization = None

        # Options for removing remote or local input and output
        self.remove_remote_input = True                 # After retrieval
        self.remove_remote_output = True                # After retrieval
        self.remove_remote_simulation_directory = True  # After retrieval
        self.remove_local_output = False                # After analysis

        # The execution handle
        self.handle = None

        # Flag indicating whether this simulation has been retrieved or not
        self.retrieved = False

        # Flag indicating whether this simulation has been analysed or not
        self.analysed = False

        # The paths to the extra simulation analysers
        self.analyser_paths = []

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load the simulation object from file
        simulation = serialization.load(path)

        # Set the path of the simulation file
        simulation.path = path

        # Return the simulation object
        return simulation

    # -----------------------------------------------------------------

    def add_analyser(self, clspath):

        """
        This function ...
        :param clspath:
        :return:
        """

        self.analyser_paths.append(clspath)

    # -----------------------------------------------------------------

    @property
    def analyser_classes(self):

        """
        This function ...
        :return:
        """

        # The list of classes
        classes = []

        # Loop over the class paths
        for class_path in self.analyser_paths:
            module_path, class_name = class_path.rsplit('.', 1)

            # Get the class of the configurable of which an instance has to be created
            module = importlib.import_module(module_path)
            cls = getattr(module, class_name)

            # Add the class to the list of classes
            classes.append(cls)

        # Return the list of classes
        return classes

    # -----------------------------------------------------------------

    def to_file(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Serialize and dump the simulation object
        serialization.dump(self, path, method="pickle")

        # Set the simulation file path
        self.path = path

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Check whether a path is defined for the simulation file
        if self.path is None: raise RuntimeError("The simulation file does not exist yet")

        # Save to the original path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Set the new path
        self.path = path

        # Set the _parameters to None to avoid an error when trying to pickle the SkiFile instance
        self._parameters = None

        # Serialize and dump the simulation object
        serialization.dump(self, self.path, method="pickle")

# -----------------------------------------------------------------
