#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.cache Contains the AnalysisRunCacher class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import AnalysisComponent
from ...core.basics.log import log
from ...core.remote.remote import Remote
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.basics.configuration import prompt_proceed

# -----------------------------------------------------------------

class AnalysisRunCacher(AnalysisComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(AnalysisRunCacher, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The remote host
        self.remote = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Cache
        self.cache()

        # 3. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(AnalysisRunCacher, self).setup(**kwargs)

        # Connect to the remote host
        self.remote = Remote()
        if not self.remote.setup(self.config.remote):
            raise RuntimeError("Could not connect to the remote host")

    # -----------------------------------------------------------------

    @property
    def run_name(self):

        """
        This function ...
        :return:
        """

        return self.config.run

    # -----------------------------------------------------------------

    @lazyproperty
    def run_path(self):

        """
        This function ...
        :return:
        """

        return self.analysis_runs.get_path(self.config.run)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_name(self):

        """
        Thisf unction ...
        :return:
        """

        return self.galaxy_name + "_analysis"

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.remote.home_directory, self.cache_directory_name)
        if not self.remote.is_directory(path): self.remote.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_path_run(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.cache_directory_path, self.run_name)
        if not self.remote.is_directory(path): self.remote.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def cache(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Caching ...")

        # Upload
        self.remote.upload_directory_to(self.run_path, self.cache_directory_path)

        # Prompt for proceed
        if not prompt_proceed("was the upload succesfull?"): exit()

        # Remove the local directory
        fs.remove_directory(self.run_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the cached runs table
        self.write_table()

    # -----------------------------------------------------------------

    def write_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the cached runs table ...")

        # Add entry to the cached table
        self.cached_table.add_entry(self.config.run, self.config.remote)

        # Save the table
        self.cached_table.save()

# -----------------------------------------------------------------
