#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.maps.dust.attenuation Contains the AttenuationDustMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from copy import copy

# Import the relevant PTS classes and modules
from ....magic.core.frame import Frame
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ...calibrations.cortese import CorteseAttenuationCalibration
from .tir_to_uv import make_tir_to_uv
from ....core.filter.filter import parse_filter
from ....core.tools import sequences
from ...core.list import NamedFrameList
from ...tools import plotting
from ....magic.core.image import Image
from ....magic.core.mask import union

# -----------------------------------------------------------------

def make_map(fuv, tir, ssfr, ssfr_colour, return_tir_to_fuv=False):

    """
    This function ...
    :param fuv:
    :param tir:
    :param ssfr:
    :param ssfr_colour:
    :param return_tir_to_fuv:
    :return: 
    """

    # Create the attenuation map maker
    maker = CorteseAttenuationMapsMaker()

    # Set input
    tirs = {"standard": tir}
    ssfrs = {ssfr_colour: ssfr}

    # Run
    maker.run(fuv=fuv, tirs=tirs, ssfrs=ssfrs)

    # Get the map
    if return_tir_to_fuv: return maker.single_map, maker.tirtofuvs["standard"]
    else: return maker.single_map

# -----------------------------------------------------------------

class CorteseAttenuationMapsMaker(Configurable):

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
        super(CorteseAttenuationMapsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The FUV map
        self.fuv = None

        # The TIR maps
        self.tirs = None

        # The ssfr maps
        self.ssfrs = None

        # Origins
        self.tirs_origins = None
        self.ssfrs_origins = None

        # Methods
        self.tirs_methods = None
        self.ssfrs_methods = None

        # NaNs
        self.tirs_nans = None
        self.ssfrs_nans = None

        # The table describing the calibration parameters from Cortese et. al 2008
        # Title of table: Relations to convert the TIR/FUV ratio in A(FUV) for different values of tau and
        # FUV − NIR/optical colours.
        self.cortese = None

        # The SSFR maps (the FUV/optical-NIR colour maps)
        self.ssfrs = dict()

        # The attenuation maps (for different FUV/optical-NIR colours)
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

        # Method name
        self.method_name = None

        # The TIR to FUV maps
        self.tirtofuvs = dict()

        # The region of interest
        self.region_of_interest = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make the dust map
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(CorteseAttenuationMapsMaker, self).setup(**kwargs)

        # Get input
        self.fuv = kwargs.pop("fuv")
        self.tirs = kwargs.pop("tirs")
        self.ssfrs = kwargs.pop("ssfrs")

        # Get origins
        self.tirs_origins = kwargs.pop("tirs_origins", None)
        self.ssfrs_origins = kwargs.pop("ssfrs_origins", None)

        # Get methods
        self.tirs_methods = kwargs.pop("tirs_methods", None)
        self.ssfrs_methods = kwargs.pop("ssfrs_methods", None)

        # Get NaN maps
        self.tirs_nans = kwargs.pop("tirs_nans", None)
        self.ssfrs_nans = kwargs.pop("ssfrs_nans", None)

        # Get method name
        self.method_name = kwargs.pop("method_name", None)
        if self.has_methods and self.method_name is None: raise ValueError("Method name should be specified when methods are given")

        # Get already calculated maps
        self.maps = kwargs.pop("maps", dict())

        # Get already calculated TIR to FUV maps
        self.tirtofuvs = kwargs.pop("tir_to_fuvs", dict())

        # The region of interest
        self.region_of_interest = kwargs.pop("region_of_interest", None)

        # Create the Cortese instance
        self.cortese = CorteseAttenuationCalibration()

    # -----------------------------------------------------------------

    @property
    def has_origins(self):

        """
        This function ...
        :return: 
        """

        return self.tirs_origins is not None and self.ssfrs_origins is not None

    # -----------------------------------------------------------------

    @property
    def has_methods(self):

        """
        This function ...
        :return:
        """

        return self.tirs_methods is not None and self.ssfrs_methods is not None

    # -----------------------------------------------------------------

    @property
    def has_tirs_nans(self):

        """
        Thisf unction ...
        :return:
        """

        return self.tirs_nans is not None

    # -----------------------------------------------------------------

    @property
    def has_ssfrs_nans(self):

        """
        Thisf unction ...
        :return:
        """

        return self.ssfrs_nans is not None

    # -----------------------------------------------------------------

    @property
    def has_nans(self):

        """
        This function ...
        :return:
        """

        return self.tirs_nans is not None and self.ssfrs_nans is not None

    # -----------------------------------------------------------------

    def has_nans_for_tir(self, tir_name):

        """
        This function ...
        :param tir_name:
        :return:
        """

        return self.has_tirs_nans and tir_name in self.tirs_nans and self.tirs_nans[tir_name] is not None

    # -----------------------------------------------------------------

    def has_nans_for_ssfr(self, ssfr_name):

        """
        Thisnfunction ...
        :param ssfr_name:
        :return:
        """

        return self.has_ssfrs_nans and ssfr_name in self.ssfrs_nans and self.ssfrs_nans[ssfr_name] is not None

    # -----------------------------------------------------------------

    @property
    def has_all_maps(self):

        """
        This function ...
        :return:
        """

        # Loop over the different TIR maps
        for name in self.tirs:

            # Loop over the different colour options
            for ssfr_colour in self.ssfrs:

                # Determine name
                key = name + "__" + ssfr_colour

                # If map is not yet present, return False
                if key not in self.maps: return False

        # ALl maps are already present
        return True

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the attenuation maps ...")

        # Dust = FUV attenuation = function of (ratio of TIR and FUV luminosity)

        # CHECK if ALL MAPS ARE ALREADY PRESENT: IN THAT CASE, WE DON'T HAVE TO CREATE THE TIR_TO_FUV MAP!
        if self.has_all_maps:
            log.debug("All maps are already present. Not creating the TIR to FUV map and the attenuation maps but setting origins and methods.")
            need_any = False
        else: need_any = True

        # Loop over the different TIR maps
        for name in self.tirs:

            # Debugging
            log.debug("Creating attenuation maps with the '" + name + "' TIR map ...")

            # Already calcualted TIR to FUV map
            if name in self.tirtofuvs:

                log.success("The '" + name + "' TIR to FUV map is already created: not creating it again")
                tir_to_fuv = self.tirtofuvs[name]
                if isinstance(tir_to_fuv, Image): tir_to_fuv = tir_to_fuv.primary

            # Not yet created
            else:

                # Debugging
                log.debug("Creating the '" + name + "' TIR to FUV map ...")

                # Make the TIR to FUV map
                if need_any: tir_to_fuv = make_tir_to_uv(self.tirs[name], self.fuv)
                #log_tir_to_fuv = Frame(np.log10(tir_to_fuv), wcs=tir_to_fuv.wcs) # unit is lost: cannot do rebinning because 'frame.unit.is_per_pixelsize' is not accessible ...
                else: tir_to_fuv = None

                # Replace NaNs and add the TIR to FUV map to the dictionary
                if tir_to_fuv is not None:

                    # Interpolate NaNs in TIR to FUV
                    tirfuv_nans = tir_to_fuv.interpolate_nans_if_below(min_max_in=self.region_of_interest)

                    # Create image with mask
                    tir_to_fuv_image = Image()
                    tir_to_fuv_image.add_frame(tir_to_fuv, "tir_to_fuv")
                    if tirfuv_nans is not None: tir_to_fuv_image.add_mask(tirfuv_nans, "nans")

                    # Add the image
                    self.tirtofuvs[name] = tir_to_fuv_image

            # Loop over the different colour options
            for ssfr_colour in self.ssfrs:

                # Determine name
                key = name + "__" + ssfr_colour

                # Set origins
                if self.has_origins:

                    origins = copy(self.tirs_origins[name])
                    origins_ssfr = copy(self.ssfrs_origins[ssfr_colour])
                    sequences.extend_unique(origins, origins_ssfr)
                    sequences.append_unique(origins, parse_filter("FUV"))
                    self.origins[key] = origins

                # Set methods
                if self.has_methods:

                    methods = copy(self.tirs_methods[name])
                    methods_ssfr = copy(self.ssfrs_methods[ssfr_colour])
                    sequences.extend_unique(methods, methods_ssfr)
                    methods.append(self.method_name)
                    self.methods[key] = methods

                # Check whether a map is already present
                if key in self.maps:
                    log.success("The '" + key + "' attenuation map is already created: not creating it again")
                    continue

                # Debugging
                log.debug("Creating the '" + key + "' attenuation map ...")

                # Get the ssfr map
                ssfr = self.ssfrs[ssfr_colour]
                if isinstance(ssfr, Image): ssfr_frame = ssfr.primary
                elif isinstance(ssfr, Frame): ssfr_frame = ssfr
                else: raise ValueError("Something went wrong")

                # Rebin and convolve the TIR-to-FUV, FUV and the sSFR maps
                frames = NamedFrameList(fuv=self.fuv, ssfr=ssfr_frame, tirtofuv=tir_to_fuv)
                frames.convolve_and_rebin()
                
                # Create log of TIRtoFUV frame
                log_tir_to_fuv = Frame(np.log10(frames["tirtofuv"].data), wcs=frames["tirtofuv"].wcs)

                # Plot log TIR to FUV
                if self.config.plot: plotting.plot_box(log_tir_to_fuv, "log10(TIR/FUV)")

                # Create the FUV attenuation map according to the calibration in Cortese et. al 2008
                fuv_attenuation = make_fuv_attenuation_map(self.cortese, ssfr_colour, log_tir_to_fuv, frames["ssfr"], plot=self.config.plot)

                # Set properties
                fuv_attenuation.unit = None # no unit for attenuation
                fuv_attenuation.filter = None # no filter for attenuation
                fuv_attenuation.wcs = frames.wcs
                fuv_attenuation.distance = frames.distance
                fuv_attenuation.pixelscale = frames.pixelscale
                fuv_attenuation.psf_filter = frames.psf_filter
                fuv_attenuation.fwhm = frames.fwhm

                # Interpolate
                nans = fuv_attenuation.interpolate_nans_if_below(min_max_in=self.region_of_interest)
                image = Image()
                image.add_frame(fuv_attenuation, "fuv_attenuation")

                nan_masks = []
                if nans is not None: nan_masks.append(nans)
                if self.has_nans_for_ssfr(ssfr_colour): nan_masks.append(self.ssfrs_nans[ssfr_colour])
                if self.has_nans_for_tir(name): nan_masks.append(self.tirs_nans[name])
                if len(nan_masks) > 0:
                    nans = union(*nan_masks, rebin=True)
                    image.add_mask(nans, "nans")

                # Set attenuation to zero where the original FUV map is smaller than zero
                fuv_attenuation[frames["fuv"] < 0.0] = 0.0

                # Make positive: replace NaNs and negative pixels by zeros
                # Set negatives and NaNs to zero
                #fuv_attenuation.replace_nans(0.0)
                fuv_attenuation.replace_negatives(0.0)

                # Add the attenuation map to the dictionary
                #self.maps[key] = fuv_attenuation
                self.maps[key] = image

    # -----------------------------------------------------------------

    @property
    def single_map(self):

        """
        This function ...
        :return: 
        """

        if len(self.maps) != 1: raise ValueError("Not a single map")
        return self.maps[self.maps.keys()[0]]

# -----------------------------------------------------------------

def make_fuv_attenuation_map(cortese, ssfr_colour, log_tir_to_fuv, ssfr, plot=False):

    """
    This function ...
    :param cortese:
    :param ssfr_colour:
    :param log_tir_to_fuv:
    :param ssfr:
    :param plot:
    :return:
    """

    # Inform the user
    log.info("Creating the FUV attenuation map according to the relation to the TIR/FUV ratio as described in Cortese et. al 2008 ...")

    # Calculate powers of log(tir_to_fuv)
    log_tir_to_fuv2 = np.power(log_tir_to_fuv.data, 2.0)
    log_tir_to_fuv3 = np.power(log_tir_to_fuv.data, 3.0)
    log_tir_to_fuv4 = np.power(log_tir_to_fuv.data, 4.0)

    # Create an empty image
    a_fuv_cortese = Frame.zeros_like(log_tir_to_fuv)

    # Use the parameters for the lowest tau (highest sSFR colour) where the sSFR colour value exceeds this maximum
    tau_min, colour_range, parameters = cortese.minimum_tau_range_and_parameters(ssfr_colour)

    # Debugging
    log.debug("Setting FUV attenuation values for tau < " + str(tau_min) + " ...")
    log.debug("This corresponds to sSFR colour values above " + str(colour_range.max))

    # Where?
    where_above = ssfr > colour_range.max

    # Show the number of pixels with this value
    log.debug("There are " + str(np.sum(where_above)) + " pixels of the sSFR map above this value")

    # Plot which pixels
    if plot: plotting.plot_mask(where_above, title="Pixels where sSFR > " + str(colour_range.max) + " (tau < " + str(tau_min) + ")")
    a_fuv_cortese[where_above] = parameters[0] + parameters[1] * log_tir_to_fuv[where_above] + parameters[2] * log_tir_to_fuv2[where_above] + \
                                                   parameters[3] * log_tir_to_fuv3[where_above] + parameters[4] * log_tir_to_fuv4[where_above]

    # Create the FUV attenuation map
    for tau, colour_range, parameters in cortese.taus_ranges_and_parameters(ssfr_colour):

        # Debugging
        log.debug("Setting FUV attenuation values for tau = " + str(tau) + " ...")

        # Debugging
        if colour_range.min == float("-inf"): log.debug("This corresponds to sSFR colour values below " + str(colour_range.max))
        else: log.debug("This corresponds to a sSFR colour range between " + str(colour_range.min) + " and " + str(colour_range.max))

        # Set mask
        where = (ssfr >= colour_range.min) * (ssfr < colour_range.max)

        # Show the number of pixels with this value
        log.debug("There are " + str(np.sum(where)) + " pixels of the sSFR map within this range")

        # Plot which pixels
        if plot:
            if colour_range.min == float("-inf"): title = "Pixels where sSFR < " + str(colour_range.max) + " (tau >= " + str(tau) + ")"
            else: title = "Pixels where " + str(colour_range.min) + " <= sSFR < " + str(colour_range.max) + " (tau = " + str(tau) + ")"
            plotting.plot_mask(where, title=title)

        # Set the appropriate pixels
        a_fuv_cortese[where] = parameters[0] + parameters[1] * log_tir_to_fuv[where] + parameters[2] * log_tir_to_fuv2[where] + \
                               parameters[3] * log_tir_to_fuv3[where] + parameters[4] * log_tir_to_fuv4[where]

    # Set attenuation to zero where tir_to_fuv is NaN
    a_fuv_cortese[np.isnan(log_tir_to_fuv)] = 0.0

    # Set attenuation to zero where sSFR colour is smaller than zero
    a_fuv_cortese[ssfr < 0.0] = 0.0

    # Return the A(FUV) map
    return a_fuv_cortese

# -----------------------------------------------------------------
