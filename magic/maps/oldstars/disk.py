#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.oldstars.disk Contains the DiskOldStellarMapMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ....core.basics.configurable import Configurable
from ...core.list import FrameList, NamedFrameList
from ....core.tools.stringify import tostr

# -----------------------------------------------------------------

def make_map(frame, bulge):

    """
    This function ...
    :param frame:
    :param bulge:
    :return: 
    """

    # Create the maker
    maker = DiskOldStellarMapMaker()

    # Set input
    frames = FrameList(frame)
    bulges = FrameList(bulge)

    # Run the maker
    maker.run(frames=frames, bulges=bulges)

    # Return the map
    return maker.single_map

# -----------------------------------------------------------------

class DiskOldStellarMapMaker(Configurable):

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
        super(DiskOldStellarMapMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The input
        self.frames = None
        self.bulges = None

        # The method name
        self.method_name = None

        # The maps
        self.maps = dict()

        # The origins
        self.origins = dict()

        # The methods
        self.methods = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 4. Make the map of old stars
        self.make_maps()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(DiskOldStellarMapMaker, self).setup(**kwargs)

        # Get input
        self.frames = kwargs.pop("frames")
        self.bulges = kwargs.pop("bulges")

        # Get already created maps
        self.maps = kwargs.pop("maps", dict())

        # Get method name
        self.method_name = kwargs.pop("method_name", None)

    # -----------------------------------------------------------------

    @property
    def filters(self):

        """
        This function ...
        :return: 
        """

        return self.frames.filters

    # -----------------------------------------------------------------

    @property
    def has_method_name(self):

        """
        This function ...
        :return:
        """

        return self.method_name is not None

    # -----------------------------------------------------------------

    def make_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the maps of old stars ...")

        # Loop over the frames
        #for name in self.frames:
        for fltr in self.filters:

            # Set name
            name = tostr(fltr, delimiter="_")

            # Set origin
            self.origins[name] = [fltr]

            # Set methods
            if self.has_method_name: self.methods[name] = [self.method_name]

            # Check if already present
            if name in self.maps:
                log.warning("The " + name + " old stellar disk map is already created: not creating it again")
                continue

            # Old stars = IRAC3.6 - bulge
            # From the IRAC 3.6 micron map, we must subtract the bulge component to only retain the disk emission

            # The relative contribution of the bulge to the 3.6mu emission
            #bulge_rel_contribution = self.parameters.bulge.f

            # Total flux of the IRAC 3.6mu image
            #total_flux = np.sum(self.images["3.6mu"].frames.primary)

            # Calculate factor
            #factor = bulge_rel_contribution * total_flux / np.sum(self.bulge)

            # Create the old stars map
            #old_stars = self.images["3.6mu"].frames.primary - factor * self.bulge

            #assert str(self.masked_bulge_frame.unit) == "Jy"

            frame = self.frames[fltr]
            bulge = self.bulges[fltr]

            # REBIN TO THE SAME PIXELSCALE (AND CONVOLVE?)
            frames = NamedFrameList(observation=frame, bulge=bulge)
            frames.convolve_and_rebin()

            # Subtract bulge from the IRAC I1 image
            minus_bulge = frames["observation"] - frames["bulge"]

            #bulge_residual = self.images["3.6mu"].frames.primary - self.disk
            #bulge_residual_path = fs.join(self.maps_intermediate_path, "bulge_residual.fits")
            #bulge_residual.save(bulge_residual_path)

            # Set the old stars map zero for pixels with low signal-to-noise in the 3.6 micron image
            #old_stars[self.irac < self.config.old_stars.irac_snr_level*self.irac_errors] = 0.0

            # Create copy
            #map = self.i1_jy_minus_bulge.copy()

            # Make sure all pixel values are larger than or equal to zero
            minus_bulge[minus_bulge < 0.0] = 0.0

            # Normalize the old stellar map
            minus_bulge.normalize()

            # Add
            self.maps[name] = minus_bulge

            # Mask pixels outside of the low signal-to-noise contour
            #old_stars[self.mask] = 0.0

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
