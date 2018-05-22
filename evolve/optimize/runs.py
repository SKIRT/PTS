#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.evolve.optimize.runs Contains the RunsOptimizer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from .stepwise import StepWiseOptimizer

# -----------------------------------------------------------------

class RunsOptimizer(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(RunsOptimizer, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 3. Perform the runs
        self.perform()

        # 4. Show
        if self.config.show: self.show()

        # 5. Write
        if self.config.write: self.write()

        # 6. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(RunsOptimizer, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def perform(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Performing the runs ...")

        # Loop over the different runs
        for run in range(self.config.runs):

            # Determine run ID
            run_id = "run" + str(run)

            # Create the optimizer
            optimizer = StepWiseOptimizer()

            # Set options
            optimizer.config.run_id = run_id

            # Inform the user
            log.info("Starting run '" + run_id + "' ...")

            # Run the optimizer
            optimizer.run()

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# -----------------------------------------------------------------
