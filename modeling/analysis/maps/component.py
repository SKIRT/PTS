#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.component Contains the MapsAnalysisComponent class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta
import numpy as np

# Import the relevant PTS classes and modules
from ..component import AnalysisRunComponent
from ....core.basics.log import log
from ...maps.component import MapMakerBase
from ...maps.collection import MapsCollection, StaticMapsCollection
from ....core.filter.filter import parse_filter
from ....magic.core.list import NamedFrameList, uniformize
from ....core.basics.distribution import Distribution
from ....core.tools import filesystem as fs
from ....core.tools import types

# -----------------------------------------------------------------

class MapsAnalysisComponent(AnalysisRunComponent, MapMakerBase):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        AnalysisRunComponent.__init__(self, no_config=True)
        MapMakerBase.__init__(self, *args, **kwargs)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        #super(MapsAnalysisComponent, self).setup(**kwargs)
        AnalysisRunComponent.setup(self, **kwargs)
        MapMakerBase.setup(self, **kwargs)

    # -----------------------------------------------------------------

    @property
    def simulated_dataset(self):
        return self.analysis_run.simulated_dataset

    # -----------------------------------------------------------------

    def get_frame_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.analysis_run.get_simulated_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_frames_for_filters(self, *fltrs, **kwargs):

        """
        This function ...
        :param fltrs:
        :param kwargs:
        :return:
        """

        # Get the frames
        frames = [self.get_frame_for_filter(fltr) for fltr in fltrs]

        # Uniformize
        if kwargs.pop("uniformize", False): frames = uniformize(*frames)

        # Return
        return frames

    # -----------------------------------------------------------------

    def get_frame(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Parse filter
        fltr = parse_filter(fltr)

        # Get frame
        return self.get_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    @property
    def frame_list(self):
        return self.analysis_run.simulated_frame_list

    #  -----------------------------------------------------------------

    @property
    def named_frame_list(self):
        return self.analysis_run.simulated_named_frame_list

    # -----------------------------------------------------------------

    @property
    def colours_path(self):
        return self.analysis_run.colour_maps_path

    # -----------------------------------------------------------------

    @property
    def ssfr_path(self):
        return self.analysis_run.ssfr_maps_path

    # -----------------------------------------------------------------

    @property
    def tir_path(self):
        return self.analysis_run.tir_maps_path

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):
        return self.analysis_run.attenuation_maps_path

    # -----------------------------------------------------------------

    @property
    def old_path(self):
        return self.analysis_run.old_maps_path

    # -----------------------------------------------------------------

    @property
    def young_path(self):
        return self.analysis_run.young_maps_path

    # -----------------------------------------------------------------

    @property
    def ionizing_path(self):
        return self.analysis_run.ionizing_maps_path

    # -----------------------------------------------------------------

    @property
    def dust_path(self):
        return self.analysis_run.dust_maps_path

    # -----------------------------------------------------------------

    @property
    def rgb_path(self):
        return self.analysis_run.rgb_maps_path

    # -----------------------------------------------------------------

    def load_collection(self):
        return MapsCollection.from_analysis_run(self.analysis_run)

    # -----------------------------------------------------------------

    def load_static_collection(self):
        return StaticMapsCollection.from_analysis_run(self.analysis_run)

    # -----------------------------------------------------------------

    def create_residuals_for_map(self, original_map, map_name, method_name=None):

        """
        This function ...
        :param original_map:
        :param map_name:
        :param method_name:
        :return:
        """

        # Get the new map
        new_map = self.maps[method_name][map_name] if method_name is not None else self.maps[map_name]

        # Rebin and normalize both maps
        maps = NamedFrameList(original=original_map, new=new_map)
        maps.rebin_to_highest_pixelscale()
        maps.normalize()

        # Create residual map
        residuals = maps["new"] / maps["original"] - 1.

        # Replace infs and nans
        residuals.replace_infs(0.0)
        residuals.replace_nans(0.0)

        # Return the residuals frame
        return residuals

    # -----------------------------------------------------------------

    def create_residuals_distribution(self, residuals, min_value=-20, max_value=20, nbins=20):

        """
        This function ...
        :param residuals:
        :param min_value:
        :param max_value:
        :param nbins:
        :return:
        """

        # Get the values within the truncation ellipse
        values = residuals.values_in(self.truncation_ellipse)

        # REMOVE EXACT ZEROES
        indices = np.argwhere(values == 0)
        values = np.delete(values, indices)

        # REMOVE TOO LOW OR TOO HIGH (PROBABLY NOISE)
        indices = np.argwhere(values < min_value)
        values = np.delete(values, indices)
        indices = np.argwhere(values > max_value)
        values = np.delete(values, indices)

        # Create distribution
        distribution = Distribution.from_values("Residual", values, nbins=nbins)

        # Return the distribution
        return distribution

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # Loop over the methods
        for method in self.maps:

            # Depending on whether subdictionaries
            if types.is_dictionary(self.maps[method]):

                # Loop over the maps
                for name in self.maps[method]:

                    # Determine path
                    plot_path = self.get_path_for_map(name, method, extension="png")

                    # If map already exists and we don't have to remake
                    if fs.is_file(plot_path) and not self.config.remake: continue

                    # Save
                    self.make_rgba_plot(name, self.maps[method][name], plot_path)

            # No different methods
            else:

                # Determine path
                plot_path = self.get_path_for_map(method, extension="png")

                # If map already exists and we don't have to remake
                if fs.is_file(plot_path) and not self.config.remake: continue

                # Save
                self.make_rgba_plot(method, self.maps[method], plot_path)

    # -----------------------------------------------------------------

    def make_rgba_plot(self, name, frame, filepath, cropping_factor=1.3):

        """
        This function ...
        :param name:
        :param frame:
        :param filepath:
        :param cropping_factor:
        :return:
        """

        # Debugging
        log.debug("Making an RGBA plot from the '" + name + "' map at '" + filepath + "' ...")

        # Crop the frame
        frame = frame.cropped_to(self.truncation_box, factor=cropping_factor)

        # Get the truncation mask and mask out the pixel beyond the truncation limit
        wcs, xsize, ysize = frame.wcs, frame.xsize, frame.ysize
        ellipse = self.truncation_ellipse.to_pixel(wcs)
        mask = ellipse.to_mask(xsize, ysize).inverse()

        frame[mask] = 0.0

        # Make RGBA image
        rgba = frame.to_rgba(scale=self.config.scale, colours=self.config.colours)
        #rgba.soften_edges(self.softening_ellipse.to_pixel(wcs), self.softening_range)

        # Save
        rgba.saveto(filepath)

# -----------------------------------------------------------------
