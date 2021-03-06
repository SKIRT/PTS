#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.skirt Contains the SKIRTTestSuite class, performing a suite of SKIRT test cases.
#
# An instance of the SkirtTestSuite class in this module represents a suite of SKIRT test cases, stored as
# a nested structure of files and directories according to a specific layout, and provides facilities to
# perform the tests, verify the results, and prepare a summary test report.

# -----------------------------------------------------------------

# Import standard modules
import filecmp
import sys
import os
import os.path
import re
import time
import datetime
import numpy as np
import multiprocessing
try: import pyfits
except ImportError: import astropy.io.fits as pyfits

# Import the relevant PTS classes and modules
from ..simulation.execute import SkirtExec
from ..basics.log import log

# -----------------------------------------------------------------
#  SkirtTestSuite class
# -----------------------------------------------------------------

## An instance of the SkirtTestSuite class represents a suite of <tt>SKIRT</tt> test cases, stored as
# a nested structure of files and directories according to a specific layout, and provides facilities to
# perform the tests, verify the results, and prepare a summary test report.
#
# A test suite consists of a set of independent test cases (i.e. test cases can be executed in arbitrary order)
#
# Each test case in a test suite is defined by a collection of files and directories as follows:
#  - a directory with arbitrary name containing all test case files and directories, called the "case directory"
#  - immediately inside the case directory there is:
#    - exactly one \em ski file with an arbitrary name (with the \c .ski filename extension) specifying the simulation
#      to be performed for the test case
#    - a directory named \c in containing the input files for the simulation, if any
#    - a directory named \c ref containing the reference files for the test, i.e. a copy of the output files
#      generated by a correct simulation run
#    - a directory named \c out to receive the actual output files when the test is performed; this directory
#      and its contents are automatically removed and created when running the test case
#    - everything else is ignored, as long as there are no additional files with a \c .ski filename extension
#
# A test suite is defined by a collection of files and directories as follows:
#  - a directory directly or indirectly containing all test cases, called the "suite directory";
#    a test suite is named after this directory
#  - each ski file directly or indirectly contained in the suite directory defines a test case that
#    must adhere to the description above (no other ski files in the same directory, special directories
#    next to the \em ski file, etc.)
#
# For example, a test suite may be structured with nested sub-suites as follows (where each \c CaseN directory
# contains a ski file plus \c ref, \c in, and \c out directories):
# \verbatim
# SKIRT Tests
#   SPH simulations
#       Case1
#       Case2
#   Geometries
#     Radial
#         Case1
#         Case2
#     Cilindrical
#         Case1
#         Case2
#         Case3
#     Full 3D
#         Case1
#         Case2
#   Instruments
# \endverbatim
#
# It is also allowed to nest test cases inside another test case, but this is not recommended.
#
class SkirtTestSuite(object):

    ## The constructor accepts three arguments:
    #
    #   - suitepath: the path of the directory containing the complete test suite
    #   - subsuite: the name of a certain subsuite (the name of one particular test case or an overarching name)
    #   - parallel: this flag indicates whether the test suite is to be executed in parallel mode or not
    #   - skirtpath: this optional argument specifies the path to the skirt executable. If none is given,
    #     the skirt version used will be the one found in the standard system path.
    #
    #  Paths may be absolute, relative to a user's home folder, or relative to the
    #  current working directory.
    #
    def __init__(self, suitepath, subsuite=None, parallel=False, skirtpath=None):

        # Set the basic characteristics of this test run
        self._suitepath = os.path.realpath(os.path.expanduser(suitepath))
        self._suitename = os.path.basename(self._suitepath)
        self._subsuitepath = findsubdirectory(self._suitepath, subsuite)
        self._parallel = parallel

        # Create a SKIRT execution context
        self._skirt = SkirtExec(skirtpath)

        # Initialize some data structures for the test run
        self._statistics = dict()
        self._simulations = []
        self._modes = []
        self._modenames = []
        self._skipatterns = []

        # Find out the number of CPU cores on this system and determine the number of threads in singleprocessing mode,
        # the number of processes and the number of threads per process in multiprocessing mode (based on whether
        # parallel mode is enabled).
        cores = multiprocessing.cpu_count()
        if self._parallel:
            threads = cores
            processes = 2
            threadspp = cores/processes
        else:
            threads = 1
            #processes = cores  # If we would choose to build the 'reference test case output' seperately on each computer
            processes = 4  # If we want the 'reference test case output' to be portable from one computer to another, this has to be fixed
            threadspp = 1

        # Check in which modes the test suite should be executed (singleprocessing and/or multiprocessing mode)
        # and create the appropriate ski file pattern(s).
        if self._subsuitepath == self._suitepath:
            # Add configurations for all three modes
            self._modes = [(cores, threads, 1, False), (1, threadspp, processes, False), (1, threadspp, processes, True)]
            self._modenames = ["in single-processing mode ", "in multi-processing mode ", "in data-parallelization mode "]
            self._skipatterns = [os.path.join(self._suitepath, "Singleprocessing", "*.ski"),
                                 os.path.join(self._suitepath, "Multiprocessing", "*.ski"),
                                 os.path.join(self._suitepath, "Dataparallelization", "*.ski")]

        elif "Singleprocessing" in self._subsuitepath:
            self._modes = [(cores, threads, 1, False)]
            self._modenames = [""]
            self._skipatterns = [os.path.join(self._subsuitepath, "*.ski")]

        elif "Multiprocessing" in self._subsuitepath:
            self._modes = [(1, threadspp, processes, False)]
            self._modenames = [""]
            self._skipatterns = [os.path.join(self._subsuitepath, "*.ski")]

        elif "Dataparallelization" in self._subsuitepath:
            self._modes = [(1, threadspp, processes, True)]
            self._modenames = [""]
            self._skipatterns = [os.path.join(self._subsuitepath, "*.ski")]

    ## This function performs all tests in the test suite, verifies the results, and prepares a summary test report.
    #  As an argument, it can take the time in seconds to sleep before checking for simulation completion again.
    #  The default value is 60 seconds.
    #
    def perform(self, sleepsecs="60"):

        # Define a name identifying this test run
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d--%H-%M-%S")
        self._testname = os.path.basename(self._subsuitepath) + "_" + timestamp

        # Inform the user of the fact that the test suite has been initiated
        log.info("Starting report for test suite " + self._subsuitepath)
        log.info("Using " + self._skirt.version() + " in " + self._skirt.root_directory)

        # Create a report file to contain a detailed report of the test run
        self._createreportfile()

        # Clean the "out" directories
        self._clean()

        # Set the number of finished simulations to zero
        self._finished = 0

        # Perform singleprocessing and multiprocessing mode sequentially
        self._numsimulations = 0
        for mode, modename, skipattern in zip(self._modes, self._modenames, self._skipatterns):

            # Start performing the simulations
            simulations = self._skirt.execute(skipattern, recursive=True, inpath="in", outpath="out", skirel=True,
                            parallel=mode[0], threads=mode[1], processes=mode[2], dataparallel=mode[3], wait=False)
            numsimulations = len(simulations)

            # Inform the user on the number of test cases (in this mode)
            log.info("Number of test cases " + modename + ": " + str(numsimulations))
            self._report.write("Number of test cases " + modename + ": " + str(numsimulations) + "<br>\n")

            # Add the new simulations to the list
            self._simulations += simulations
            self._numsimulations += numsimulations

            # Verify the results for each test case
            self._verify(sleepsecs)

            # Wait for the skirt execution context to finish
            self._skirt.wait()

        # Write statistics about the number of successful test cases
        self._writestatistics()

        # Close the report file
        self._report.close()

    ## This function creates a HTML file that contains a detailed report of the test run
    def _createreportfile(self):

        csscommands = """<style type="text/css">
                             [id^="togList"],                        /* HIDE CHECKBOX */
                             [id^="togList"] ~ .list,                /* HIDE LIST */
                             [id^="togList"] + label  span + span,   /* HIDE "Collapse" */
                             [id^="togList"]:checked + label span{   /* HIDE "Expand" (IF CHECKED) */
                               display:none;
                             }
                             [id^="togList"]:checked + label span + span{
                               display:inline-block;                 /* SHOW "Collapse" (IF CHECKED) */
                             }
                             [id^="togList"]:checked ~ .list{
                               display:block;                        /* SHOW LIST (IF CHECKED) */
                             }
                             </style>"""

        # Open the report file
        filepath = os.path.join(self._suitepath, "report_" + self._testname + ".html")
        self._report = open(filepath, 'w')

        # Write some general info to the report file
        self._report.write("<html>\n<head>\n</head>\n<body>\n")
        self._report.write(csscommands + "\n")
        self._report.write("Report file for test suite " + self._subsuitepath + "<br>\n")
        self._report.write("Using " + self._skirt.version() + " in " + self._skirt.root_directory + "<br>\n")

    ## This function cleans up the contents of all "out" directories that reside next to a ski file
    def _clean(self):

        # Find every directory in the suite directory that has the name "out" and subsequently delete its content
        for dirpath, dirs, files in os.walk(self._subsuitepath):
            if dirpath.endswith("/out"):
                for skifile in filter(lambda fn: fn.endswith(".ski"), os.listdir(dirpath[:-4])):
                    prefix = skifile[:-4] + "_"
                    for name in filter(lambda fn: fn.startswith(prefix), files):
                        os.remove(os.path.join(dirpath, name))

    ## This function verifies the results for each test case
    def _verify(self, sleepsecs):

        # Verify the results of each simulation, once it finishes (processed items are removed from the list)
        while True:
            for simulation in self._simulations[:]:
                if simulation.status().endswith("shed"): # "Finished" or "Crashed"
                    self._reportsimulation(simulation)
                    self._simulations.remove(simulation)
            if len(self._simulations) == 0 or not self._skirt.isrunning(): break
            time.sleep(sleepsecs)

        # Report potentially unstarted simulations (which are still in the list of simulations and therefore not yet processed)
        for simulation in self._simulations:
            self._reportsimulation(simulation)
            self._simulations.remove(simulation)

    ## This function writes statistics about the number of successful test cases
    def _writestatistics(self):

        log.info("Summary for total of: "+ str(self._numsimulations))
        self._report.write("Summary for total of: " + str(self._numsimulations) + "<br>\n")

        for key,value in self._statistics.iteritems():

            log.info("  " + key + ": " + str(value))
            self._report.write("  " + key + ": " + str(value) + "<br>\n")

        log.info("Finished report for test suite " + self._subsuitepath)

    ## This function verifies and reports on the test result of the given simulation.
    # It writes a line to the console and it updates the statistics.
    def _reportsimulation(self, simulation):

        # Get the full path of the simulation directory and the name of this directory
        casedirpath = os.path.dirname(simulation.outpath())
        casename = os.path.relpath(casedirpath, self._suitepath) if self._parallel else os.path.basename(casedirpath)

        # Determine the most relevant string to identify this simulation within the current test suite
        residual = casedirpath
        while (not self._parallel):
            if casename == os.path.basename(self._subsuitepath): break
            if (os.path.dirname(residual) == self._subsuitepath):
                break
            else:
                residual = os.path.dirname(residual)
                casename = os.path.join(os.path.basename(residual), casename)

        # Report the status of this simulation
        status = simulation.status()
        if status == "Finished":

            # Increment the number of finished simulations
            self._finished += 1

            extra, missing, differ = self._finddifference(casedirpath)

            if len(extra) == 0 and len(missing) == 0 and len(differ) == 0:

                status = "Succeeded"
                log.success("Test case " + casename + ": succeeded")

                # Write to the report file
                self._report.write("<span style='color:LightGreen'>Test case " + casename + ": succeeded</span><br>\n")

            else:

                status = "Failed"
                log.error("Test case " + casename + ": failed")

                # Write to the report file
                self._report.write("<div class='row'>\n")
                self._report.write("  <input id='togList" + str(self._finished) + "' type='checkbox'>\n")
                self._report.write("  <label for='togList" + str(self._finished) + "'>\n")
                self._report.write("    <span style='color:red'>Test case " + casename + ": failed [click to expand] </span>\n")
                self._report.write("    <span style='color:red'>Test case " + casename + ": failed [click to collapse]</span>\n")
                self._report.write("  </label>\n")
                self._report.write("  <div class='list'>\n")
                self._report.write("    <ul>\n")

                if len(extra) > 0:

                    if len(extra) == 1: self._report.write("      <li>The output contains an extra file:<br>")
                    else: self._report.write("      <li>The output contains extra files:<br>")

                    for filename in extra:

                        self._report.write("  - " + filename + "<br>")

                    self._report.write("</li>\n")

                if len(missing) > 0:

                    if len(missing) == 1: self._report.write("      <li>The output misses a file:<br>")
                    else: self._report.write("      <li>The output misses files:<br>")

                    for filename in missing:

                        self._report.write("  - " + filename + "<br>")

                    self._report.write("</li>\n")

                if len(differ) > 0:

                    if len(differ) == 1: self._report.write("      <li>The following file differs:<br>")
                    else: self._report.write("      <li>The following files differ:<br>")

                    for filename in differ:

                        self._report.write("  - " + filename + "<br>")

                    self._report.write("</li>\n")

                self._report.write("    </ul>\n  </div>\n</div>\n")

        self._statistics[status] = self._statistics.get(status,0) + 1

    ## This function looks for relevant differences between the contents of the output and reference directories
    # of a test case. The test case is specified as an absolute case directory path.
    # The function returns a brief message describing the first relevant difference found, or the empty string
    # if there are no relevant differences.
    def _finddifference(self, casedirpath):

        # Get the full output and reference paths
        outpath = os.path.join(casedirpath, "out")
        refpath = os.path.join(casedirpath, "ref")

        # Check for the presence of the reference directory
        if not os.path.isdir(refpath): return "Test case has no reference directory"

        # Initialize lists to contain the extra files, missing files and differing files
        extra = []
        missing = []
        differ = []

        # Verify list of filenames
        if len(filter(lambda fn: os.path.isdir(fn), os.listdir(outpath))) > 0: return "Output contains a directory"
        dircomp = filecmp.dircmp(outpath, refpath, ignore=['.DS_Store'])

        # Create a list of files that were not present in the reference directory
        for file in dircomp.left_only: extra.append(file)

        # Create a list of files that are missing from the output directory
        for file in dircomp.right_only: missing.append(file)

        # Compare files, focusing on those that aren't trivially equal.
        matches, mismatches, errors = filecmp.cmpfiles(outpath, refpath, dircomp.common, shallow=False)

        for filename in mismatches + errors:

            reffile = os.path.join(refpath, filename)
            outfile = os.path.join(outpath, filename)

            # For a soft comparison between files, check if both files are similar enough
            if self._parallel:
                if not similarfiles(reffile, outfile): differ.append(filename)

            # In the other case, check if both files are equal (except for potential timestamps)
            else:
                if not equalfiles(reffile, outfile): differ.append(filename)

        # Return the list of extra files, the list of missing files and the list of differing files
        return extra, missing, differ

# -----------------------------------------------------------------

# This functions searches for directories within a certain parent directories that have a specific name
def findsubdirectory(parent, name):

    # If the name of the subdirectory is empty, simply return the parent directory.
    if not name: return parent

    # Recursively find all subdirectories of the parent directory.
    dirlist = os.walk(parent)

    # Find those subdirectories whose name correponds to the argument "name".
    # The path of the first subdirectory found with that name is returned.
    for dirpath, dirnames, filenames in dirlist:
        for dirname in dirnames:
            # Compare only the name of the directory "dirname" with the "name" string.
            if (dirname.lower() == name.lower()): return os.path.join(dirpath, dirname)
            # Compare the "dirname", prefixed with the name of its parent directory, to the "name" string. (e.g. "Instruments/Simple")
            if (os.path.join(os.path.basename(dirpath), dirname).lower() == name.lower()): return os.path.join(dirpath, dirname)

    # If no match is found, show an error message and quit
    log.error(name + " not found in " + parent)
    sys.exit()

# -----------------------------------------------------------------

## This function returns True if the specified files are similar (meaning that there relative difference is lower than
#  a certain threshold.
def similarfiles(filepath1, filepath2, threshold=0.1):

    # Don't compare log files because it is too complicated (time & duration differences, changes to messages, ...)
    if filepath1.endswith("_log.txt"): return True

    # Supported file types
    if filepath1.endswith(".fits"): return similarfitsfiles(filepath1, filepath2, threshold)
    if filepath1.endswith("_parameters.xml"): return equaltextfiles(filepath1, filepath2, 1)
    if filepath1.endswith("_parameters.tex"): return equaltextfiles(filepath1, filepath2, 2)
    if filepath1.endswith("_sed.dat"): return similarseds(filepath1, filepath2, threshold)
    if filepath1.endswith("_convergence.dat"): return similarconvergence(filepath1, filepath1, threshold)
    if filepath1.endswith("_isrf.dat"): return True
    if filepath1.endswith("_abs.dat"): return True
    if filepath1.endswith("_gridxz.dat") or filepath1.endswith("_gridyz.dat") or filepath1.endswith("_gridxy.dat") or filepath1.endswith("_gridxyz.dat"):
        return True

    # Unsupported file type
    return False

## This function returns True if the specified files are equal except for irrelevant differences (such as
# time stamps in a file header or prolog), and False otherwise. Since the structure of the allowed differences
# varies with the type of file, this function dispatches the comparison to various other functions depending on
# the last portion of the first filename (i.e a generalized filename extension).
def equalfiles(filepath1, filepath2):

    # Don't compare log files because it is too complicated (time & duration differences, changes to messages, ...)
    if filepath1.endswith("_log.txt"): return True

    # Supported file types
    if filepath1.endswith(".fits"): return equalfitsfiles(filepath1, filepath2)
    if filepath1.endswith("_parameters.xml"): return equaltextfiles(filepath1, filepath2, 1)
    if filepath1.endswith("_parameters.tex"): return equaltextfiles(filepath1, filepath2, 2)

    # Unsupported file type
    return False

## This function returns True if the specified text files are equal except for time stamp information in
#  a number of lines not larger than \em allowedDiffs, and False otherwise.
def equaltextfiles(filepath1, filepath2, allowedDiffs):
    return equallists(readlines(filepath1), readlines(filepath2), allowedDiffs)

## This function returns True if the specified fits files are equal except for time stamp information in a single
#  header record, and False otherwise.
def equalfitsfiles(filepath1, filepath2):
    return equallists(readblocks(filepath1), readblocks(filepath2), 1)

## This function returns True if the specified FITS files are similar, meaning that their relative difference
#  is lower than the specified threshold.
def similarfitsfiles(filepath1, filepath2, threshold):

    # Obtain the datacubes for both FITS files
    data1 = pyfits.getdata(filepath1)
    data2 = pyfits.getdata(filepath2)

    # Compare the dimension of the datacubes
    if data1.shape != data2.shape: return False

    # Determine a cutoff for the lowest fluxes
    cutoff = min(np.amax(data1)/100.0,np.amax(data2)/100.0)

    # Compare the pixel values
    return np.isclose(data1, data2, rtol=threshold, atol=cutoff).any()

## This function returns True if the specified SED data files are similar, meaning that their relative difference
#  is lower than the specified threshold.
def similarseds(filepath1, filepath2, threshold):

    # Obtain the fluxes of both SED files
    fluxes1 = np.loadtxt(filepath1, usecols=[1], ndmin=1)
    fluxes2 = np.loadtxt(filepath2, usecols=[1], ndmin=1)

    # Compare the number of flux points, if
    if len(fluxes1) != len(fluxes2): return False

    # Compare the fluxes
    return np.isclose(fluxes1, fluxes2, rtol=threshold, atol=.0).any()

## This function returns \c True if the specified convergence files contain similar surface densities and dust masses
def similarconvergence(filepath1, filepath2, threshold):

    # Get the lines from the first file
    with open(filepath1) as file1: lines1 = file1.readlines()

    # Get the lines from the second file
    with open(filepath2) as file2: lines2 = file2.readlines()

    # Compare the number of lines in both files
    if len(lines1) != len(lines2): return False

    # Compare the values in the convergence files, by looping over all lines
    for i in range(len(lines1)):

        line1 = lines1[i]
        line2 = lines2[i]

        # Check if this line states a value
        if "value = " in line1:

            # Get the values in the corresponding line for both files
            value1 = float(line1.split("value = ")[1].split()[0])
            value2 = float(line2.split("value = ")[1].split()[0])

            # If the values are not similar, return False
            if not np.isclose(value1, value2, rtol=threshold): return False

    # If no dissimilarities were found, return True
    return True

## This function returns True if the specified lists are equal except for possible time information in
#  a number of items not larger than \em allowedDiffs, and False otherwise.
def equallists(list1, list2, allowedDiffs):

    # The lists must have the same length, which must be at least 2 (to avoid everything being read into 1 line)
    length = len(list1)
    if length < 2 or length != len(list2): return False

    # compare the lists item by item
    diffs = 0
    for index in range(length):
        if list1[index] != list2[index]:
            # verify against allowed number of differences
            diffs += 1
            if diffs > allowedDiffs: return False
            # verify that the differing items are identical up to certain numerics
            pattern = re.compile(r"([0-9a-f]{7,7}(-dirty){0,1})|\d{1,4}")
            item1 = re.sub(pattern, "*", list1[index])
            item2 = re.sub(pattern, "*", list2[index])
            if item1 != item2: return False

    # no relevant differences
    return True

# -----------------------------------------------------------------

## This function reads the lines of the specified text file into a list of strings, and returns the list.
def readlines(filepath):
    with open(filepath) as f: result = f.readlines()
    return result

## This function reads blocks of 80 bytes from the specified binary file into a list of strings, and returns the list.
def readblocks(filepath):
    result = []
    with open(filepath) as f:
        while (True):
            block = f.read(80)
            if len(block) == 0: break
            result.append(block)
    return result

# -----------------------------------------------------------------
