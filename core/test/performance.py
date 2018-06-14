#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.test.performance Contains the PeformanceTest class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ..simulation.execute import SkirtExec
from ..basics.log import log
from ..basics.configurable import Configurable
from ..launch.batch import BatchLauncher

# -----------------------------------------------------------------

class PerformanceTest(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(PerformanceTest, self).__init__(*args, **kwargs)

        # The local SKIRT execution environment
        self.skirt = SkirtExec()

        # The SKIRT batch launcher
        self.launcher = BatchLauncher()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        super(PerformanceTest, self).setup(**kwargs)

        # Set launcher options
        self.launcher.config.remotes = [self.config.remote]

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

# -----------------------------------------------------------------
