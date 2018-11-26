#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.maps.tir Contains the TIRMapsAnalyser class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....magic.maps.tir.single import SingleBandTIRMapMaker
from ....magic.maps.tir.multi import MultiBandTIRMapMaker
from ....magic.core.list import FrameList
from ....core.tools.utils import lazyproperty
from ...maps.tir import singleband_filter_names, multiband_filter_names
from ....core.filter.filter import parse_filter
from .component import MapsAnalysisComponent
from ....core.tools import time
from ....core.units.parsing import parse_quantity as q

# -----------------------------------------------------------------

total_name = "total"
young_name = "young"
ionizing_name = "ionizing"
unevolved_name = "unevolved"

# -----------------------------------------------------------------

class TIRMapsAnalyser(MapsAnalysisComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(TIRMapsAnalyser, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        Thisf ucntion ...
        :param kwargs:
        :return:
        """

        # Make the maps
        self.make_maps()

        # Write
        self.write()

        # Plot
        self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(TIRMapsAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    @lazyproperty
    def available_filters_singleband(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the colours
        for filter_name in singleband_filter_names:

            # Parse fltr
            fltr = parse_filter(filter_name)

            # If no image is avilalbe for this filters, skip
            if not self.has_frame_for_filter(fltr): continue

            # otherwise, add to the list of filters
            filters.append(fltr)

        # Return the available filters
        return filters

    # -----------------------------------------------------------------

    @lazyproperty
    def available_filters_multiband(self):

        """
        This function ...
        :return:
        """

        filters = []

        # Loop over the colours
        for filter_name in multiband_filter_names:

            # Parse filter
            fltr = parse_filter(filter_name)

            # If no image is avilalbe for this filters, skip
            if not self.has_frame_for_filter(fltr): continue

            # otherwise, add to the list of filters
            filters.append(fltr)

        # Return the available filters
        return filters

    # -----------------------------------------------------------------

    def load_frames_singleband(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        frames = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_singleband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frame.distance = self.galaxy_distance
            frames.append(frame, fltr)

        # Return the frames
        return frames

    # -----------------------------------------------------------------

    def load_frames_multiband(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the data ...")

        # Frames and error maps
        frames = FrameList()

        # Loop over the filters
        for fltr in self.available_filters_multiband:

            # Debugging
            log.debug("Loading the '" + str(fltr) + "' frame ...")

            # Load the frame
            frame = self.get_frame_for_filter(fltr)
            frame.distance = self.galaxy_distance
            frames.append(frame, fltr)

        # Return the frames
        return frames

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making TIR map ...")

        # Make maps based on a single band
        self.make_tir_maps_single()

        # Make maps based on multiple bands
        self.make_tir_maps_multi()

        # Make maps based on integrating the full dust SED
        self.make_tir_maps_integration()

    # -----------------------------------------------------------------

    def make_tir_maps_single(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on a single band ...")

        # Set the method name
        method_name = "single"

        # Create
        maker = SingleBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames = self.load_frames_singleband()
        print(frames)
        maker.run(frames=frames, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    def make_tir_maps_multi(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on multiple bands ...")

        # Set method name
        method_name = "multi"

        # Create
        maker = MultiBandTIRMapMaker()

        # Get current
        if self.config.remake: current = dict()
        else: current = self.get_current_maps_method(method_name)

        # Run
        frames = self.load_frames_multiband()
        print(frames)
        maker.run(frames=frames, maps=current, method_name=method_name)

        # Set the maps
        self.maps[method_name] = maker.maps

        # Set the origins
        self.origins[method_name] = maker.origins

        # Set the methods
        self.methods[method_name] = maker.methods

    # -----------------------------------------------------------------

    @property
    def total_cube(self):
        return self.model.total_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def young_cube(self):
        return self.model.young_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def ionizing_cube(self):
        return self.model.sfr_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @property
    def unevolved_cube(self):
        return self.model.unevolved_bolometric_luminosity_cube_earth

    # -----------------------------------------------------------------

    @lazyproperty
    def tir_min_wavelength(self):
        return q("8 micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def tir_max_wavelength(self):
        return q("1000 micron")

    # -----------------------------------------------------------------

    def make_tir_maps_integration(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making maps based on integrating the dust emission ...")

        # Set method name
        method_name = "integration"

        # The TIR maps
        maps = dict()

        # Create
        maps[total_name] = integrate_datacube(self.total_cube, min_wavelength=self.tir_min_wavelength, max_wavelength=self.tir_max_wavelength)
        maps[young_name] = integrate_datacube(self.young_cube, min_wavelength=self.tir_min_wavelength, max_wavelength=self.tir_max_wavelength)
        maps[ionizing_name] = integrate_datacube(self.ionizing_cube, min_wavelength=self.tir_min_wavelength, max_wavelength=self.tir_max_wavelength)
        maps[unevolved_name] = integrate_datacube(self.unevolved_cube, min_wavelength=self.tir_min_wavelength, max_wavelength=self.tir_max_wavelength)

        # Set the maps
        self.maps[method_name] = maps

        # Set the origins
        #self.origins[method_name] = origins

        # Set the methods
        #self.methods[method_name] = methods

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):
        return self.tir_path

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the colour maps
        self.write_maps()

        # Write origins
        self.write_origins()

        # Write the methods
        self.write_methods()

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

def integrate_datacube(cube, min_wavelength=None, max_wavelength=None, show_time=False):

    """
    This function ...
    :param cube:
    :param min_wavelength:
    :param max_wavelength:
    :param show_time:
    :return:
    """

    # Splice the datacube
    cube = cube.splice(min_wavelength=min_wavelength, max_wavelength=max_wavelength)
    wavelength_unit = cube.wavelength_unit

    # Get wavelength grid
    wavelength_grid = cube.wavelength_grid
    #wavelengths = wavelength_grid.wavelengths(unit=wavelength_unit, asarray=True)
    deltas = wavelength_grid.deltas(unit=wavelength_unit, asarray=True)

    # Get converted to wavelength density
    cube = cube.converted_to_corresponding_wavelength_density_unit(wavelength_unit=wavelength_unit)

    # Get 3D array
    array = cube.asarray(axis=2) # wavelength along second axis

    # Perform the integration
    with time.elapsed_timer() as elapsed:
        integrated2 = np.trapz(y=array, dx=deltas)
        integrated = np.sum(array * deltas, axis=2)
        if show_time: print("Integration performed in " + str(elapsed()) + " seconds")

    # Return
    return integrated

# -----------------------------------------------------------------
