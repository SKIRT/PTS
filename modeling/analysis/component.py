#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.component Contains the AnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ..core.component import ModelingComponent
from ...core.tools import filesystem as fs
from ...core.launch.timing import TimingTable
from ...core.launch.memory import MemoryTable

# -----------------------------------------------------------------

class AnalysisComponent(ModelingComponent):
    
    """
    This class...
    """

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisComponent, self).__init__(config)

        # -- Attributes --

        # The path to the timing table
        self.timing_table_path = None

        # The path to the memory table
        self.memory_table_path = None

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisComponent, self).setup()

        ## TIMING TABLE

        # Set the path to the timing table
        self.timing_table_path = fs.join(self.analysis_path, "timing.dat")

        # Initialize the timing table if necessary
        if not fs.is_file(self.timing_table_path):

            # Create the table and save it
            timing_table = TimingTable()
            timing_table.saveto(self.timing_table_path)

        ## MEMORY TABLE

        # Set the path to the memory table
        self.memory_table_path = fs.join(self.analysis_path, "memory.dat")

        # Initialize the memory table if necessary
        if not fs.is_file(self.memory_table_path):

            # Create the table and save it
            memory_table = MemoryTable()
            memory_table.saveto(self.memory_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def timing_table(self):

        """
        This function ...
        :return:
        """

        return TimingTable.from_file(self.timing_table_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def memory_table(self):

        """
        This function ...
        :return:
        """

        return MemoryTable.from_file(self.memory_table_path)

# -----------------------------------------------------------------
