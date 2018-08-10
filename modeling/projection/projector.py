#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.projection.projector Contains the Projector class and derived classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from .data import DataProjections
from ...core.basics.log import log
from ..core.data import Data3D, show_data_properties
from ...magic.core.image import Image
from ...core.tools.utils import lazyproperty
from ...magic.tools.plotting import plot_map

# -----------------------------------------------------------------

faceon_name = "faceon"
edgeon_name = "edgeon"
orientations = [faceon_name, edgeon_name]

# -----------------------------------------------------------------

class Projector(Configurable):
    pass

# -----------------------------------------------------------------

stddev_name = "stddev"
ncells_name = "ncells"

# -----------------------------------------------------------------

class DataProjector(Projector):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(DataProjector, self).__init__(*args, **kwargs)

        # The 3D data
        self.data = None

        # The projections instance
        self.projections = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Show
        self.show()

        # Project
        self.project()

        # Interpolate
        if self.config.interpolate: self.interpolate()

        # Write
        self.write()

        # Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DataProjector, self).setup(**kwargs)

        # Load the data
        self.load_data()

    # -----------------------------------------------------------------

    def load_data(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the 3D data from file ...")

        # Load
        self.data = Data3D.from_file(self.config.filename)

    # -----------------------------------------------------------------

    @property
    def do_faceon(self):
        return faceon_name in self.config.orientations

    # -----------------------------------------------------------------

    @property
    def do_edgeon(self):
        return edgeon_name in self.config.orientations

    # -----------------------------------------------------------------

    @property
    def do_any(self):
        return self.do_faceon or self.do_edgeon

    # -----------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing data properties ...")

        # Show properties of the data
        show_data_properties(self.data)

    # -----------------------------------------------------------------

    def project(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the data ...")

        # Create projections
        self.projections = DataProjections(self.data, faceon=self.do_faceon, edgeon=self.do_edgeon,
                                           faceon_height=self.config.height, edgeon_width=self.config.width,
                                           faceon_spacing=self.config.spacing, edgeon_spacing=self.config.spacing,
                                           faceon_spacing_factor=self.config.spacing_factor,
                                           edgeon_spacing_factor=self.config.spacing_factor)

    # -----------------------------------------------------------------

    @property
    def faceon_projection(self):
        return self.projections.projection_faceon

    # -----------------------------------------------------------------

    @property
    def edgeon_projection(self):
        return self.projections.projection_edgeon

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon(self):
        return Image.from_frames(self.projections.faceon, stddev=self.projections.faceon_stddev, ncells=self.projections.faceon_ncells)

    # -----------------------------------------------------------------

    @property
    def faceon_frame(self):
        return self.faceon.primary

    # -----------------------------------------------------------------

    @property
    def faceon_ncells(self):
        return self.faceon.frames.ncells

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon(self):
        return Image.from_frames(self.projections.edgeon, stddev=self.projections.edgeon_stddev, ncells=self.projections.edgeon_ncells)

    # -----------------------------------------------------------------

    @property
    def edgeon_frame(self):
        return self.edgeon.primary

    # -----------------------------------------------------------------

    @property
    def edgeon_ncells(self):
        return self.edgeon.frames.ncells

    # -----------------------------------------------------------------

    def interpolate_map(self, frame, ncells):

        """
        This function ...
        :param frame:
        :param ncells:
        :return:
        """

        # Copy
        #interpolated = frame.copy()

        # Get outside nans
        outside_nans = frame.nans.largest()
        # plotting.plot_mask(outside_nans, title="outside nans")
        #outside_nans.saveto(fs.join(self.cell_heating_path, "outside_nans.fits"))
        not_nans = outside_nans.inverse()
        not_nans.disk_dilate(radius=self.config.not_nans_dilation_radius)
        # not_nans.fill_holes()
        #not_nans.saveto(fs.join(self.cell_heating_path, "not_nans.fits"))
        do_nans = not_nans.largest().inverse()

        # Get mask
        where = ncells.where_smaller_than(self.config.min_ncells)

        # plotting.plot_mask(where, title="where smaller than " + str(self.config.min_ncells))
        # plotting.plot_mask(self.map_interpolated.nans, title="nans")

        # Put pixels to NaN
        frame.replace_by_nans(where)

        # plotting.plot_mask(self.map_interpolated.nans, title="nans")

        # Interpolate nans
        frame.interpolate_nans(sigma=2.)
        frame.replace_by_nans(do_nans)

        # Return the interpolated frame
        #return interpolated

    # -----------------------------------------------------------------

    def interpolate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the maps ...")

        # Faceon
        if self.do_faceon: self.interpolate_faceon()

        # Edgeon
        if self.do_edgeon: self.interpolate_edgeon()

    # -----------------------------------------------------------------

    def interpolate_faceon(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Interpolating the face-on map ...")

        # Interpolate
        self.interpolate_map(self.faceon_frame, self.faceon_ncells)

    # -----------------------------------------------------------------

    def interpolate_edgeon(self):

        """
        This functino ...
        :return:
        """

        # Debugging
        log.debug("Interpolating the edge-on map ...")

        # Interpolate
        self.interpolate_map(self.edgeon_frame, self.edgeon_ncells)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the projections?
        self.write_projections()

        # Write the projected maps
        self.write_maps()

    # -----------------------------------------------------------------

    def write_projections(self):

        """
        This function ...
        :return:
        """

        # Faceon
        if self.do_faceon: self.write_faceon_projection()

        # Edgeon
        if self.do_edgeon: self.write_edgeon_projection()

    # -----------------------------------------------------------------

    def write_faceon_projection(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("faceon.proj")

        # Save
        self.faceon_projection.saveto(path)

    # -----------------------------------------------------------------

    def write_edgeon_projection(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("edgeon.proj")

        # Save
        self.edgeon_projection.saveto(path)

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Faceon
        if self.do_faceon: self.write_faceon_map()

        # Edgeon
        if self.do_edgeon: self.write_edgeon_map()

    # -----------------------------------------------------------------

    def write_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("faceon.fits")

        # Save
        self.faceon.saveto(path)

    # -----------------------------------------------------------------

    def write_edgeon_map(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("edgeon.fits")

        # Save
        self.edgeon.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the maps
        self.plot_maps()

    # -----------------------------------------------------------------

    def plot_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the maps ...")

        # Faceon
        if self.do_faceon: self.plot_faceon_map()

        # Edgeon
        if self.do_edgeon: self.plot_edgeon_map()

    # -----------------------------------------------------------------

    @property
    def plotting_interval(self):

        """
        This function ...
        :return:
        """

        if self.config.plotting.minmax is not None: return self.config.plotting.minmax
        else: return self.config.plotting.interval

    # -----------------------------------------------------------------

    def plot_faceon_map(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("faceon.pdf")

        # Plot
        plot_map(self.faceon_frame, path=path, interval=self.plotting_interval, contours=self.config.plotting.contours,
                 ncontours=self.config.plotting.ncontours, contours_color=self.config.plotting.contours_color)

    # -----------------------------------------------------------------

    def plot_edgeon_map(self):

        """
        This function ...
        :return:
        """

        # Determine path
        path = self.output_path_file("edgeon.pdf")

        # Plot
        plot_map(self.edgeon_frame, path=path, interval=self.plotting_interval, contours=self.config.plotting.contours,
                 ncontours=self.config.plotting.ncontours, contours_color=self.config.plotting.contours_color)

# -----------------------------------------------------------------
