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
import copy
import os.path
import types as systypes
import numpy as np
import importlib
import warnings

# Import the relevant PTS classes and modules
from ..tools import serialization
from ..tools import filesystem as fs
from .skifile import SkiFile
from .logfile import LogFile
from ..tools import archive as arch
from ..launch.options import AnalysisOptions
from .status import LogSimulationStatus
from .input import SimulationInput
from .output import SimulationOutput, ExtractionOutput, PlottingOutput, MiscOutput
from ..tools import types
from ..tools.stringify import tostr
from ..tools import strings
from ..tools.utils import create_lazified_class

# -----------------------------------------------------------------

## This function returns a list of SkirtSimulation objects depending on the type and value of the \em source argument:
# - missing: all simulations for which there is a log file in the current directory
# - empty string: all simulations for which there is a log file in the current directory
# - string containing a slash: all simulations for which there is a log file in the specified directory
# - nonempty string without a slash: the simulation with the specified prefix in the current directory
# - simulation object: the simulation represented by the specified object
# - list of strings and/or simulation objects: all simulations in the listed objects, as defined above
#
def createsimulations(source="", single=False, name=None, cls=None):

    # Initialize list for simulations
    simulations = []
    sourcelist = source if isinstance(source, (systypes.TupleType,systypes.ListType)) else [ source ]

    # Set class
    if cls is None: cls = SkirtSimulation
    elif not issubclass(cls, SkirtSimulation): raise ValueError("Passed class must be subclass of SkirtSimulation")

    # Loop over the sources
    for source in sourcelist:

        if isinstance(source, systypes.StringTypes):

            if source == "" or "/" in source:

                dirpath = os.path.realpath(os.path.expanduser(source))
                logfiles = arch.listdir(dirpath, "_log.txt")

                for logfile in logfiles:

                    prefix = logfile[:-8]
                    ski_path = fs.join(dirpath, prefix + ".ski")
                    ski_path = ski_path if fs.is_file(ski_path) else None

                    # Create the simulation object
                    sim = cls(prefix=logfile[:-8], outpath=dirpath, ski_path=ski_path)

                    # Add the simulation
                    simulations.append(sim)

            else:

                if os.path.exists(source + "_log.txt"):

                    ski_path = source + ".ski"
                    ski_path = ski_path if fs.is_file(ski_path) else None

                    # Create the simulation object
                    sim = cls(prefix=source, ski_path=ski_path)

                    # Add the simulation
                    simulations.append(sim)

        elif isinstance(source, SkirtSimulation): simulations.append(source)
        else: raise ValueError("Unsupported source type for simulation")

    # If a single simulation is expected
    if single:

        if len(simulations) == 0: raise ValueError("No simulations were found matching the source '" + str(source) + "'")
        elif len(simulations) > 1: raise ValueError("Multiple simulations were found for source '" + str(source) + "': " + ", ".join([simulation.prefix() + " in " + simulation.output_path for simulation in simulations]))
        else:

            simulation = simulations[0]
            if name is not None: simulation.name = name
            return simulation

    # If multiple simulations are expected
    elif name is not None: raise ValueError("Cannot specify name if single is False")
    else: return simulations

# -----------------------------------------------------------------
#  SkirtSimulation class
# -----------------------------------------------------------------

# Dictionary of default attribute values
# Add new attributes here, will be set in constructor,
# AND also in from_file when missing (when simulation file is outdated w.r.t. the additions in this class)
default_attributes = dict()

# The simulation file path
default_attributes["path"] = None

# The parallelization properties
default_attributes["parallelization"] = None

# The options for analysing the simulation output
default_attributes["analysis"] = AnalysisOptions()

# The paths to the extra simulation analysers
#default_attributes["analyser_paths"] = []

# Flag indicating whether this simulation has been analysed or not
default_attributes["analysed"] = False

# Performed extraction steps
default_attributes["analysed_extraction"] = []

# Performed plotting steps
default_attributes["analysed_plotting"] = []

# Performed misc steps
default_attributes["analysed_misc"] = []

# Performed batch analysis
default_attributes["analysed_batch"] = False

# Performed scaling analysis
default_attributes["analysed_scaling"] = False

# Performed extra analysis
default_attributes["analysed_extra"] = [] # class names

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
            if types.is_sequence(inpath): self._inpath = inpath
            elif types.is_string_type(inpath): self._inpath = os.path.realpath(os.path.expanduser(inpath))
            elif isinstance(inpath, SimulationInput): self._inpath = inpath
            elif types.is_dictionary(inpath): self._inpath = SimulationInput.unchecked(**inpath)
            else: raise ValueError("Invalid value for 'inpath': " + tostr(inpath) + " (" + str(type(inpath)) + ")")
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
        # BETTER AS PROPERTIES SO THAT SETTING THEM WILL NOT LEAVE INCONSISTENT STATE!
        #self.input_path = self._inpath
        #self.output_path = self._outpath

        # Try to obtain ski path from file with name prefix.ski in current working directory
        if self.ski_path is None and fs.is_file(fs.join(fs.cwd(), prefix + ".ski")): self.ski_path = fs.join(fs.cwd(), prefix + ".ski")

        # The base path = the directory where the ski file is located
        if self.ski_path is not None: self.base_path = fs.directory_of(self.ski_path)
        else:
            self.base_path = None # if ski path is not specified, base path is unknown
            self.ski_path = self.outfilepath("parameters.xml")
            if not fs.is_file(self.ski_path):
                warnings.warn("No parameters file can be found for this simulation")
                self.ski_path = None

        # Set parameters, if passed
        self._parameters = parameters

        # Provide placeholders for caching frequently-used objects
        self._units = None
        self._processes = None
        self._threads = None

        # A name given to the simulation
        self._name = name

        # Set default values for other attributes
        for attr_name in default_attributes:

            # Create a copy
            value = copy.deepcopy(default_attributes[attr_name])

            # Set the attribute
            setattr(self, attr_name, value)

    @property
    def input_path(self):
        return self._inpath

    @input_path.setter
    def input_path(self, value):
        self._inpath = value

    @property
    def output_path(self):
        return self._outpath

    @output_path.setter
    def output_path(self, value):
        self._outpath = value

    @property
    def analyser_paths(self):

        """
        This function ...
        :return:
        """

        return self.analysis.analyser_paths

    ## Set analysed flags
    def set_analysed(self, value=True, all=False):

        changed = False
        old_analysed = self.analysed
        self.analysed = value
        if old_analysed != value: changed = True

        if not value and all: changed |= self.unset_analysed_extraction()
        if not value and all: changed |= self.unset_analysed_plotting()
        if not value and all: changed |= self.unset_analysed_misc()
        if not value and all: changed |= self.unset_analysed_batch()
        if not value and all: changed |= self.unset_analysed_scaling()
        if not value and all: changed |= self.unset_analysed_extra()

        return changed

    ## Unset analysed extraction
    def unset_analysed_extraction(self):

        from ..tools import sequences

        # Extraction
        if not sequences.is_empty(self.analysed_extraction): changed = True
        else: changed = False
        self.analysed_extraction = []
        return changed

    @property
    def analysed_any_extraction(self):
        from ..tools import sequences
        return not sequences.is_empty(self.analysed_extraction)

    @property
    def extraction_features(self):
        features = []
        # extraction_names = [progress_name, timeline_name, memory_name]
        from ..launch.options import progress_name, timeline_name, memory_name
        if self.analysis.extraction.progress: features.append(progress_name)
        if self.analysis.extraction.timeline: features.append(timeline_name)
        if self.analysis.extraction.memory: features.append(memory_name)
        return features

    @property
    def analysed_all_extraction(self):
        from ..tools import sequences
        #from ..launch.options import extraction_names
        return sequences.contains_all(self.analysed_extraction, self.extraction_features)

    ## Unset analysed plotting
    def unset_analysed_plotting(self):

        from ..tools import sequences

        # Plotting
        if not sequences.is_empty(self.analysed_plotting): changed = True
        else: changed = False
        self.analysed_plotting = []
        return changed

    @property
    def analysed_any_plotting(self):
        from ..tools import sequences
        return not sequences.is_empty(self.analysed_plotting)

    @property
    def plotting_features(self):
        features = []
        # plotting_names = [progress_name, timeline_name, memory_name, seds_name, grids_name]
        from ..launch.options import progress_name, timeline_name, memory_name, seds_name, grids_name
        if self.analysis.plotting.progress: features.append(progress_name)
        if self.analysis.plotting.timeline: features.append(timeline_name)
        if self.analysis.plotting.memory: features.append(memory_name)
        if self.analysis.plotting.seds: features.append(seds_name)
        if self.analysis.plotting.grids: features.append(grids_name)
        return features

    @property
    def analysed_all_plotting(self):
        from ..tools import sequences
        #from ..launch.options import plotting_names
        #return sequences.contains_all(self.analysed_plotting, plotting_names)
        return sequences.contains_all(self.analysed_plotting, self.plotting_features)

    ## Unset analysed misc
    def unset_analysed_misc(self):

        from ..tools import sequences

        # Misc
        if not sequences.is_empty(self.analysed_misc): changed = True
        else: changed = False
        self.analysed_misc = []
        return changed

    @property
    def analysed_any_misc(self):
        from ..tools import sequences
        return not sequences.is_empty(self.analysed_misc)

    @property
    def misc_features(self):
        features = []
        # misc_names = [rgb_name, animations_name, fluxes_name, fluxes_from_images_name, images_name]
        from ..launch.options import rgb_name, animations_name, fluxes_name, fluxes_from_images_name, images_name
        if self.analysis.misc.rgb: features.append(rgb_name)
        if self.analysis.misc.animations: features.append(animations_name)
        if self.analysis.misc.fluxes: features.append(fluxes_name)
        if self.analysis.misc.fluxes_from_images: features.append(fluxes_from_images_name)
        if self.analysis.misc.images: features.append(images_name)
        return features

    @property
    def analysed_all_misc(self):
        from ..tools import sequences
       # from ..launch.options import misc_names
        return sequences.contains_all(self.analysed_misc, self.misc_features)

    ## Unset analysed batch
    def unset_analysed_batch(self):

        # Batch
        if self.analysed_batch: changed = True
        else: changed = False
        self.analysed_batch = False
        return changed

    ## Unset analysed scaling
    def unset_analysed_scaling(self):

        # Scaling
        if self.analysed_scaling: changed = True
        else: changed = False
        self.analysed_scaling = False
        return changed

    ## Unset analysed extra
    def unset_analysed_extra(self):

        from ..tools import sequences

        # Extra
        if not sequences.is_empty(self.analysed_extra): changed = True
        else: changed = False
        self.analysed_extra = []
        return changed

    @property
    def analysed_any_extra(self):
        from ..tools import sequences
        return not sequences.is_empty(self.analysed_extra)

    @property
    def analysed_all_extra(self):
        for analyser_class in self.analyser_classes:
            class_name = analyser_class.__name__
            if class_name not in self.analysed_extra: return False
        return True

    @property
    def analysed_any(self):
        return self.analysed_any_extraction or self.analysed_any_plotting or self.analysed_any_misc or self.analysed_batch or self.analysed_scaling or self.analysed_any_extra

    ## This property returns a SingleSimulationDefinition object
    @property
    def definition(self):
        from .definition import SingleSimulationDefinition
        return SingleSimulationDefinition(self.ski_path, self.output_path, self.input_path, name=self.name)

    ## This function returns a SkirtArguments object
    def get_arguments(self, logging_options=None, parallelization=None):
        from .arguments import SkirtArguments
        return SkirtArguments.from_definition(self.definition, logging_options=logging_options, parallelization=parallelization)

    ## This property returns a SimulationInput object
    @property
    def input(self):
        if not self.has_input: return None
        elif isinstance(self.input_path, SimulationInput): return self.input_path
        else: return SimulationInput.from_any(self.input_path)

    ## This property returns whether the simulation input is defined in terms of a single directory path
    @property
    def has_input_directory(self):
        return self.input_path is not None and types.is_string_type(self.input_path)

    ## This function returns whether a ski or parameters file is found for this simulation
    @property
    def has_ski(self):
        return self.ski_path is not None

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
        if isinstance(self._inpath, basestring): return os.path.join(self._inpath, name)
        elif isinstance(self._inpath, SimulationInput): return self._inpath[name]
        else: raise ValueError("The input path is of invalid type")

    ## This function returns the absolute path for a simulation output file, given the file's partial name
    # (the partial name does not include the prefix and the subsequent underscore).
    def outfilepath(self, partialname):
        return os.path.join(self._outpath, self._prefix + "_" + partialname)

    ## This function returns the absolute path for the simulation log file
    def logfilepath(self):
        return self.outfilepath("log.txt")

    ## This property returns whether the log file is present (yet)
    @property
    def has_logfile(self):
        return fs.is_file(self.logfilepath())

    ## This function returns a LogFile object created from the simulation's log file
    @property
    def log_file(self):
        path = self.logfilepath()
        if not fs.is_file(path): raise IOError("The log file is not present at '" + path + "'")
        return LogFile(path)

    # -----------------------------------------------------------------

    @property
    def extraction_path(self):
        return self.analysis.extraction.path

    @extraction_path.setter
    def extraction_path(self, value):
        self.analysis.extraction.path = value

    @property
    def plotting_path(self):
        return self.analysis.plotting.path

    @plotting_path.setter
    def plotting_path(self, value):
        self.analysis.plotting.path = value

    @property
    def misc_path(self):
        return self.analysis.misc.path

    @misc_path.setter
    def misc_path(self, value):
        self.analysis.misc.path = value

    # -----------------------------------------------------------------

    @property
    def output(self):
        return SimulationOutput.from_directory(self.outpath(), prefix=self.prefix())

    @property
    def extraction_output(self):
        if self.has_extraction_output is not None: return ExtractionOutput.from_directory(self.extraction_path)
        else: return None

    @property
    def plotting_output(self):
        if self.has_plotting_output is not None: return PlottingOutput.from_directory(self.plotting_path)
        else: return None

    @property
    def misc_output(self):
        if self.has_misc_output is not None: return MiscOutput.from_directory(self.misc_path)
        else: return None

    # -----------------------------------------------------------------

    ## This function returns a SimulationStatus object, that can be refreshed when desired
    def get_status(self):
        logpath = self.logfilepath()
        return LogSimulationStatus(logpath)

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
        depths = np.loadtxt(arch.opentext(filepath), usecols=(8,))
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

    ## This function updates the analysis options with extra options that have been added after this simulation object was created
    def update_analysis_options(self):
        self.analysis = AnalysisOptions(**self.analysis.to_dict())

    @property
    def from_batch(self):
        return self.analysis.timing_table_path is not None or self.analysis.memory_table_path is not None

    @property
    def from_scaling_test(self):
        return self.analysis.scaling_run_name is not None

    @property
    def from_modeling(self):
        return self.analysis.modeling_path is not None

    @property
    def has_extraction_output(self):
        return self.extraction_path is not None and fs.is_directory(self.extraction_path) and not fs.is_empty(self.extraction_path)

    @property
    def has_plotting_output(self):
        return self.plotting_path is not None and fs.is_directory(self.plotting_path) and not fs.is_empty(self.plotting_path)

    @property
    def has_misc_output(self):
        return self.misc_path is not None and fs.is_directory(self.misc_path) and not fs.is_empty(self.misc_path)

    ## This function adds an analyser class to the simulation
    def add_analyser(self, clspath):
        self.analyser_paths.append(clspath)

    @property
    def analyser_class_names(self):

        """
        Thisn function ...
        :return:
        """

        # The list of class names
        names = []

        # Loop over the class paths
        for class_path in self.analyser_paths:
            module_path, class_name = class_path.rsplit('.', 1)
            names.append(class_name)

        # Return the list of class names
        return names

    @property
    def analyser_classes(self):

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

    @property
    def log_file_path(self):

        """
        This function ...
        :return:
        """

        return self.outfilepath("log.txt")

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

        # Check the ID
        filename = fs.strip_extension(fs.name(path))
        if strings.is_integer(filename):
            simulation_id = int(filename)
            if simulation_id != simulation.id:
                warnings.warn("The ID of the simulation object doesn't match the ID in the filename, fixing ...")
                simulation.id = simulation_id
                simulation.save()

        # Loop over the attribute names, check if defined
        for attr_name in default_attributes:

            # This attribute is present: OK
            if hasattr(simulation, attr_name): continue

            # Get the default value
            value = copy.copy(default_attributes[attr_name])

            # Set the attribute
            setattr(simulation, attr_name, value)

        # THIS IS A HACK, FOR WHEN A FILTER OBJECT COULD STILL CONTAIN AN LXML.ETREE ELEMENT (NOW SOLVED)
        if hasattr(simulation, "analysis_plotting_ignore_filter_names"):
            from ..filter.filter import parse_filter
            # Check
            assert simulation.analysis.plotting.ignore_filters == []
            # Set the filters to the analysis options
            simulation.analysis.plotting.ignore_filters = [parse_filter(name) for name in simulation.analysis_plotting_ignore_filter_names]
            # Remove the 'hack' attribute, make the simulation object 'normal' again
            delattr(simulation, "analysis_plotting_ignore_filter_names")

        # THIS IS A FIX FOR SIMULATIONS THAT STILL HAVE THE ANALYSER_PATHS AS AN ATTRIBUTE, INSTEAD OF BEING DEFINED
        # IN THE ANALYSIS OPTIONS
        if "analyser_paths" in simulation.__dict__:

            # Get the paths
            analyser_paths = simulation.__dict__.pop("analyser_paths")

            # Check that the analyser paths in the analysis options are not defined
            if "analyser_paths" in simulation.analysis: raise RuntimeError("Something is wrong")

            # Update the analysis options
            simulation.update_analysis_options()

            # Set the analyser paths in the analysis options
            simulation.analysis.analyser_paths = analyser_paths

        # Return the simulation object
        return simulation

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
        if self.path is None:
            warnings.warn("Not saving this local simulation object")
            return

        # Save to the original path
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Set the new path
        if update_path: self.path = path

        # Set the _parameters to None to avoid an error when trying to pickle the SkiFile instance
        parameters = self._parameters
        self._parameters = None

        # THIS IS A HACK, FOR WHEN A FILTER OBJECT COULD STILL CONTAIN AN LXML.ETREE ELEMENT (NOW SOLVED)
        if len(self.analysis.plotting.ignore_filters) > 0:
            filters = self.analysis.plotting.ignore_filters
            self.analysis.plotting.ignore_filters = []
            filter_names = [str(fltr) for fltr in filters]
            self.analysis_plotting_ignore_filter_names = filter_names

        # Serialize and dump the simulation object
        serialization.dump(self, self.path, method="pickle")

        # Set the parameters
        self._parameters = parameters

# -----------------------------------------------------------------

# Dictionary of default attribute values of remote simulation
default_remote_attributes = dict()

# Properties of the remote host on which the simulation was run
default_remote_attributes["host_id"] = None
default_remote_attributes["cluster_name"] = None

# Basic properties
default_remote_attributes["id"] = None
default_remote_attributes["remote_ski_path"] = None
default_remote_attributes["remote_simulation_path"] = None
default_remote_attributes["remote_input_path"] = None
default_remote_attributes["remote_output_path"] = None
default_remote_attributes["submitted_at"] = None

# Options for retrieval
default_remote_attributes["retrieve_types"] = None

# Options for removing remote or local input and output
default_remote_attributes["remove_remote_input"] = True                 # After retrieval
default_remote_attributes["remove_remote_output"] = True                # After retrieval
default_remote_attributes["remove_remote_simulation_directory"] = True  # After retrieval
default_remote_attributes["remove_local_output"] = False                # After analysis

# The execution handle
default_remote_attributes["handle"] = None

# Flag indicating whether this simulation has finished or not
default_remote_attributes["finished"] = False

# Flag indicating whether this simulation has been retrieved or not
default_remote_attributes["retrieved"] = False

# -----------------------------------------------------------------

class RemoteSimulation(SkirtSimulation):

    """
    This class ...
    """

    def __init__(self, ski_path, input_path, output_path, **kwargs):

        """
        The constructor ...
        :param ski_path:
        :param input_path:
        :param output_path:
        :param kwargs:
        :return:
        """

        # Determine the simulation prefix
        prefix = fs.strip_extension(fs.name(ski_path))

        # Call the constructor of the base class
        super(RemoteSimulation, self).__init__(prefix, input_path, output_path, ski_path, name=kwargs.pop("name", None))

        # -- Attributes --

        # Set default values for other attributes
        for attr_name in default_remote_attributes:

            # Create a copy
            if kwargs.get(attr_name, None) is not None: value = kwargs.pop(attr_name)
            else: value = copy.deepcopy(default_remote_attributes[attr_name])

            # Set the attribute
            setattr(self, attr_name, value)

        # Remote
        self._remote = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Load from file
        simulation = super(RemoteSimulation, cls).from_file(path)

        # Loop over the attribute names, check if defined
        for attr_name in default_remote_attributes:

            # This attribute is present: OK
            if hasattr(simulation, attr_name): continue

            # Get the default value
            value = copy.copy(default_remote_attributes[attr_name])

            # Set the attribute
            setattr(simulation, attr_name, value)

        # Return the remote simulation
        return simulation

    # -----------------------------------------------------------------

    def set_finished(self, value=True):

        """
        This function ...
        :param value:
        :return:
        """

        old_finished = self.finished
        self.finished = value
        return old_finished != self.finished

    # -----------------------------------------------------------------

    def set_retrieved(self, value=True):

        """
        This function ...
        :param value:
        :return:
        """

        old_retrieved = self.retrieved
        self.retrieved = value
        return old_retrieved != self.retrieved

    # -----------------------------------------------------------------

    @property
    def host(self):

        """
        This function ...
        :return:
        """

        from ..remote.host import load_host
        if self._remote is not None: return self._remote.host
        elif self.host_id is not None: return load_host(self.host_id, clustername=self.cluster_name)
        else: return None

    # -----------------------------------------------------------------

    @property
    def remote(self):

        """
        This function ...
        :return:
        """

        from ..remote.remote import Remote
        if self._remote is None: self._remote = Remote(host_id=self.host_id)
        return self._remote

    # -----------------------------------------------------------------

    @remote.setter
    def remote(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        from ..remote.remote import load_remote

        # Load remote
        value = load_remote(value)

        # Check host ID
        if self.host_id is not None:
            if value.host_id != self.host_id: raise ValueError("This is the wrong remote instance: host ID should be '" + self.host_id + "'")
        else: self.host_id = value.host_id

        # Set the remote
        self._remote = value

    # -----------------------------------------------------------------

    def get_remote_output(self, remote=None):

        """
        This function ...
        :param remote:
        :return:
        """

        # Get remote
        if remote is None: remote = self.remote
        else: self.remote = remote

        # Return the simulation output
        return SimulationOutput.from_remote_directory(self.remote_output_path, remote, prefix=self.prefix())

    # -----------------------------------------------------------------

    @property
    def remote_output(self):

        """
        Thisf unction ...
        :return:
        """

        return self.get_remote_output()

    # -----------------------------------------------------------------

    def get_remote_input(self, remote=None):

        """
        This function ...
        :param remote:
        :return:
        """

        # Get remote
        if remote is None: remote = self.remote
        else: self.remote = remote

        # Return the simulation input
        if not self.has_input: return None
        else: return SimulationInput.from_remote_directory(self.remote_input_path, remote, prefix=self.prefix())

    # -----------------------------------------------------------------

    @property
    def remote_input(self):

        """
        This function ...
        :return:
        """

        return self.get_remote_input()

    # -----------------------------------------------------------------

    def remote_input_file_path(self, name):

        """
        This function returns the absolute path for a simulation input file, given the file's name
        """

        return fs.join(self.remote_input_path, name)

    # -----------------------------------------------------------------

    def remote_output_file_path(self, partialname):

        """
        This function returns the absolute path for a simulation output file, given the file's partial name
        (the partial name does not include the prefix and the subsequent underscore).
        :param partialname:
        :return:
        """

        return fs.join(self.remote_output_path, self._prefix + "_" + partialname)

    # -----------------------------------------------------------------

    @property
    def remote_log_file_path(self):

        """
        This function ...
        :return:
        """

        return self.remote_output_file_path("log.txt")

    # -----------------------------------------------------------------

    def get_remote_log_file_paths(self, remote):

        """
        This function returns a list of absolute filepaths for all log files produced by the simulation, including
        the master log file and any log files produced by parallel (MPI) processes. The list includes only paths for
        log files that actually exist, and the paths are listed in order of process rank.
        :return:
        """

        logname = self._prefix + "_log"
        filepaths = remote.files_in_path(self.remote_output_path, startswith=logname, extension="txt")
        nfiles = len(filepaths)
        if nfiles == 0: return []

        had_logzero = False
        if self.remote_log_file_path in filepaths:
            had_logzero = True
            filepaths.remove(self.remote_log_file_path)

        # Sort on process rank
        filepaths = list(sorted(filepaths, key=lambda path: int(path.split("_logP")[1].split(".")[0])))
        if had_logzero: filepaths = [self.remote_log_file_path] + filepaths
        return filepaths

    # -----------------------------------------------------------------

    def remove_from_remote(self, remote, full=False):

        """
        This function ...
        :param remote:
        :param full:
        :return:
        """

        # Remove the remote input, if present, if requested
        if (self.remove_remote_input or full) and self.has_input and remote.is_directory(self.remote_input_path): remote.remove_directory(self.remote_input_path)

        # Remove the remote output, if requested
        if (self.remove_remote_output or full) and remote.is_directory(self.remote_output_path): remote.remove_directory(self.remote_output_path)

        # If both the input and output directories have to be removed, the remote simulation directory
        # can be removed too
        if (self.remove_remote_simulation_directory or full) and remote.is_directory(self.remote_simulation_path): remote.remove_directory(self.remote_simulation_path)

    # -----------------------------------------------------------------

    def saveto(self, path, update_path=True):

        """
        This function ...
        :param path:
        :param update_path:
        :return:
        """

        # Set the remote to None
        remote = self._remote
        self._remote = None

        # Call the saveto function of the base class
        super(RemoteSimulation, self).saveto(path, update_path=update_path)

        # Set the remote back
        self._remote = remote

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

StaticSkirtSimulation = create_lazified_class(SkirtSimulation, "StaticSkirtSimulation")

# -----------------------------------------------------------------
