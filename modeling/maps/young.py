#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.youngstars Contains the YoungStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .component import MapMakingComponent
from ...magic.maps.youngstars.young import YoungStellarMapsMaker
from ...core.tools import filesystem as fs
from ...core.tools.stringify import tostr
from ...core.basics.containers import create_subdict
from .component import select_maps
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

methods = None

# -----------------------------------------------------------------

class YoungStellarMapMaker(MapMakingComponent):

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
        super(YoungStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The input FUV and FUV error maps
        self.fuv = None
        self.fuv_errors = None

        # The map of the old stellar disk
        self.old = None

        # The maps of FUV attenuation
        self.fuv_attenuations = None

        # The origins
        self.old_origin = None
        self.fuv_attenuations_origins = None

        # Methods
        self.old_method = None
        self.fuv_attenuations_methods = None

        # Nans
        self.fuv_attenuations_nans = None

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return self.maps_young_path

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the necessary input maps
        self.load_input()

        # 3. Make the map of young stars
        self.make_maps()

        # 4. Writing
        self.write()

        # 5. Plotting
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(YoungStellarMapMaker, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

    # -----------------------------------------------------------------

    def load_input(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the necessary input ...")

        # Load the GALEX FUV image and error map
        self.load_fuv()

        # Load FUV attenuation map
        self.load_fuv_attenuation_maps()

        # Load old stellar map
        self.load_old_stellar_map()

    # -----------------------------------------------------------------

    def load_fuv(self):

        """
        This function ...
        :return:
        """

        # Get FUV frame and error map
        self.fuv = self.dataset.get_frame("GALEX FUV") # in original MJy/sr units
        self.fuv_errors = self.dataset.get_errormap("GALEX FUV") # in original MJy/sr units

    # -----------------------------------------------------------------

    def load_fuv_attenuation_maps(self):

        """
        THis function ...
        :return:
        """

        # Inform the user
        log.info("Loading the maps of the FUV attenuation ...")

        # Get the FUV attenuation maps
        fuv_attenuations, fuv_attenuations_origins, fuv_attenuations_methods, fuv_attenuations_nans = self.get_fuv_attenuation_maps_origins_methods_and_nans(flatten=True, cortese=self.config.use_cortese, buat=self.config.use_buat)

        # Debugging
        log.debug("Making selection ...")

        # Set things
        names = self.config.attenuation_maps
        prompt = self.config.select_attenuation
        title = "FUV attenuation maps to correct FUV emission"

        # Get only certain FUV attenuation maps
        if names is not None:
            fuv_attenuations = create_subdict(fuv_attenuations, names)
            fuv_attenuations_origins = create_subdict(fuv_attenuations_origins, names)
            fuv_attenuations_methods = create_subdict(fuv_attenuations_methods, names)
            fuv_attenuations_nans = create_subdict(fuv_attenuations_nans, names)

        # Select interactively
        if prompt:
            fuv_attenuations, attenuation_names = select_maps(fuv_attenuations, title, return_names=True)
            fuv_attenuations_origins = create_subdict(fuv_attenuations_origins, attenuation_names)
            fuv_attenuations_methods = create_subdict(fuv_attenuations_methods, attenuation_names)
            fuv_attenuations_nans = create_subdict(fuv_attenuations_nans, attenuation_names)

        # Set
        self.fuv_attenuations = fuv_attenuations
        self.fuv_attenuations_origins = fuv_attenuations_origins
        self.fuv_attenuations_methods = fuv_attenuations_methods
        self.fuv_attenuations_nans = fuv_attenuations_nans

    # -----------------------------------------------------------------

    @property
    def use_old_bulge(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.old_component == "bulge"

    # -----------------------------------------------------------------

    @property
    def use_old_total(self):

        """
        This function ...
        :return:
        """

        return self.config.old_component == "total"

    # -----------------------------------------------------------------

    @property
    def use_old_disk(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.old_component == "disk"

    # -----------------------------------------------------------------

    def load_old_stellar_map(self):

        """
        This function ...
        """

        # Inform the user
        log.info("Loading the map of old stars ...")

        # Load map
        if self.use_old_total: self.old = self.get_old_stellar_total_map(self.config.old)
        elif self.use_old_disk: self.old = self.get_old_stellar_disk_map(self.config.old)
        elif self.use_old_bulge: self.old = self.get_old_stellar_bulge_map(self.config.old)
        else: raise ValueError("Invalid option for 'old_component'")

        # Set origin
        self.old_origin = self.config.old

        # Set the old method
        self.old_method = self.config.old_component

    # -----------------------------------------------------------------

    @lazyproperty
    def factors(self):

        """
        This function ...
        :return:
        """

        return self.config.factor_range.linear(self.config.factor_nvalues, as_list=True)

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return: 
        """

        # Inform the user
        log.info("Making the maps ...")

        # Clear?
        if self.config.clear: self.clear_current_all()

        # Get the current maps
        if self.config.remake: current = dict()
        else: current = self.get_current_maps(factors=self.factors)

        # Create the map maker
        maker = YoungStellarMapsMaker()

        # Debugging
        log.debug("Using the following factors for subtracting diffuse emission: " + tostr(self.factors, delimiter=", "))

        # Run the map maker
        maker.run(fuv=self.fuv, fuv_errors=self.fuv_errors, old=self.old, fuv_attenuations=self.fuv_attenuations,
                  factors=self.factors, old_origin=self.old_origin, fuv_attenuations_origins=self.fuv_attenuations_origins,
                  old_method=self.old_method, fuv_attenuations_methods=self.fuv_attenuations_methods, maps=current,
                  region_of_interest=self.truncation_ellipse, fuv_attenuations_nans=self.fuv_attenuations_nans)

        # Set the maps
        self.maps = maker.maps

        # Set the origins
        self.origins = maker.origins

        # Set the methods
        self.methods = maker.methods

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the maps
        self.write_maps()

        # 2. Write origins
        self.write_origins()

        # 3. Write the methods
        self.write_methods()

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot
        self.plot_young()

# -----------------------------------------------------------------
