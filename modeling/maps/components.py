#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.components Contains the ComponentMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapsComponent
from ...core.basics.configuration import prompt_string_list
from ...core.tools.utils import lazyproperty
from ...magic.core.frame import Frame

# -----------------------------------------------------------------

class ComponentMapsMaker(MapsComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ComponentMapsMaker, self).__init__(*args, **kwargs)

        # The selections
        self.old_selection = None
        self.young_selection = None
        self.ionizing_selection = None
        self.dust_selection = None

        # The maps
        self.old_maps = dict()
        self.young_maps = dict()
        self.ionizing_maps = dict()
        self.dust_maps = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Prompt
        self.prompt()

        # 3. Load the maps
        self.load_maps()

        # 4. Process the maps
        self.process_maps()

        # Deproject the maps
        self.deproject()

        # Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsMaker, self).setup(**kwargs)

        # Set selections
        if self.config.old is not None: self.old_selection = self.config.old
        if self.config.young is not None: self.young_selection = self.config.young
        if self.config.ionizing is not None: self.ionizing_selection = self.config.ionizing
        if self.config.dust is not None: self.dust_selection = self.config.dust

        # All maps?
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # All
        if self.config.all_dust: self.old_selection = self.dust_map_names
        if self.config.all_young: self.young_selection = self.young_map_names
        if self.config.all_ionizing: self.ionizing_selection = self.ionizing_map_names
        if self.config.all_dust: self.dust_selection = self.dust_map_names

    # -----------------------------------------------------------------

    @property
    def has_old_selection(self):

        """
        This function ...
        :return:
        """

        return self.old_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_young_selection(self):

        """
        This function ...
        :return:
        """

        return self.young_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_ionizing_selection(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_selection is not None

    # -----------------------------------------------------------------

    @property
    def has_dust_selection(self):

        """
        Thisf unctino ...
        :return:
        """

        return self.dust_selection is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def old_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_old_stellar_disk_map_paths()

    # -----------------------------------------------------------------

    @property
    def old_map_names(self):

        """
        This function ...
        :return:
        """

        return self.old_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def young_map_paths(self):

        """
        This function ...
        :return:
        """

        return self.maps_collection.get_young_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def young_map_names(self):

        """
        This function ...
        :return:
        """

        return self.young_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_map_paths(self):

        """
        Thisn function ...
        :return:
        """

        # Get map paths
        return self.maps_collection.get_ionizing_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_names(self):

        """
        This function ...
        :return:
        """

        return self.ionizing_map_paths.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map_paths(self):

        """
        Thisn function ...
        :return:
        """

        # Get paths
        # map_paths = self.maps_collection.get_dust_map_paths(flatten=True)
        return self.maps_collection.get_not_hot_dust_map_paths(flatten=True)

    # -----------------------------------------------------------------

    @property
    def dust_map_names(self):

        """
        This function ...
        :return:
        """

        return self.dust_map_paths.keys()

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("")

        # Old
        if not self.has_old_selection: self.prompt_old()

        # Young
        if not self.has_young_selection: self.prompt_young()

        # Ionizing
        if not self.has_ionizing_selection: self.prompt_ionizing()

        # Dust
        if not self.has_dust_selection: self.prompt_dust()

    # -----------------------------------------------------------------

    def prompt_old(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the old stellar maps ...")

        # Ask for the old stellar maps
        self.old_selection = prompt_string_list("old_maps", "selection of old stellar disk maps", choices=self.old_map_names)

    # -----------------------------------------------------------------

    def prompt_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the young stellar maps ...")

        # Ask for the young stellar maps
        self.young_selection = prompt_string_list("young_maps", "selection of young stellar disk maps", choices=self.young_map_names)

    # -----------------------------------------------------------------

    def prompt_ionizing(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the ionizing stellar maps ...")

        # Ask for the ionizing stellar map
        self.ionizing_selection = prompt_string_list("ionizing_maps", "selection of ionizing stellar disk maps", choices=self.ionizing_map_names)

    # -----------------------------------------------------------------

    def prompt_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the dust maps ...")

        # Ask for the dust map to use
        self.dust_selection = prompt_string_list("dust_maps", "selection of dust disk maps", choices=self.dust_map_names)

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected maps ...")

        # Load the old stellar maps
        self.load_old()

        # Load the young stellar maps
        self.load_young()

        # Load the ionizing stellar maps
        self.load_ionizing()

        # Load the dust maps
        self.load_dust()

    # -----------------------------------------------------------------

    def load_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected old stellar maps ...")

        # Load
        for name in self.old_selection: self.old_maps[name] = Frame.from_file(self.old_map_paths[name])

    # -----------------------------------------------------------------

    def load_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected young stellar maps ...")

        # Load
        for name in self.young_selection: self.young_maps[name] = Frame.from_file(self.young_map_paths[name])

    # -----------------------------------------------------------------

    def load_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected ionizing stellar maps ...")

        # Load
        for name in self.ionizing_selection: self.ionizing_maps[name] = Frame.from_file(self.ionizing_map_paths[name])

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected dust maps ...")

        # Load
        for name in self.dust_selection: self.dust_maps[name] = Frame.from_file(self.dust_map_paths[name])

    # -----------------------------------------------------------------

    def process_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the maps ...")

        # Truncate
        self.truncate_maps()

        #self.
        #self.soften_edges()

    # -----------------------------------------------------------------

    def truncate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the maps ...")



    # -----------------------------------------------------------------

    def deproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the maps ...")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        self.write_maps()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Old
        self.write_old()

        self.write_young()

        self.write_ionizing()

        self.write_dust()

    # -----------------------------------------------------------------

    def write_old(self):

        """
        This function ...
        :return:
        """

    # -----------------------------------------------------------------

    def write_young(self):

    # -----------------------------------------------------------------

    def write_ionizing(self):

    # -----------------------------------------------------------------

    def write_dust(self):

        """
        This function ...
        :return:
        """

# -----------------------------------------------------------------