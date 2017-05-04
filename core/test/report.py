#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.report Contains the TestReport class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# -----------------------------------------------------------------

class TestReport(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        :param config:
        """

        # Call the constructor of the base class
        super(TestReport, self).__init__()

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

    # -----------------------------------------------------------------

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

    ## This function writes statistics about the number of successful test cases
    def _writestatistics(self):

        log.info("Summary for total of: "+ str(self._numsimulations))
        self._report.write("Summary for total of: " + str(self._numsimulations) + "<br>\n")

        for key,value in self._statistics.iteritems():

            log.info("  " + key + ": " + str(value))
            self._report.write("  " + key + ": " + str(value) + "<br>\n")

        log.info("Finished report for test suite " + self._subsuitepath)

# -----------------------------------------------------------------
