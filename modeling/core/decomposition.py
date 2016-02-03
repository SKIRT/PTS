#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.decomposition Contains the GalaxyDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from .component import ModelingComponent

# -----------------------------------------------------------------

class GalaxyDecomposer(ModelingComponent):
    
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
        super(GalaxyDecomposer, self).__init__(config)

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Do the decomposition
        self.decompose()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyDecomposer, self).setup()

    # -----------------------------------------------------------------

    def decompose(self):

        """
        This function ...
        :return:
        """

        # TODO: do the bulge/disk fitting here
        # Run the decomposition
        #decomposer.run()

        # Set path of bulge and disk images
        bulge_path = os.path.join(self.prep_path, "Bulge", "bulge.fits")
        disk_path = os.path.join(self.prep_path, "Disk", "disk.fits")

        # Create list
        paths = {"Bulge": bulge_path, "Disk": disk_path}

        # For bulge and disk ...
        for name, path in paths.items():

            # Open the frame
            frame = Frame.from_file(path)

            # Convolve the frame to the PACS 160 resolution
            kernels_dir = os.path.expanduser("~/Kernels")
            kernel_path = os.path.join(kernels_dir, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")
            kernel = Frame.from_file(kernel_path)
            frame.convolve(kernel)

            # Rebin the convolved frame to the PACS 160 frame (just as we did with the other images)
            reference_path = os.path.join(self.data_path, "PACS160.fits")
            reference = Frame.from_file(reference_path)
            frame.rebin(reference)

            # Finally, crop the image
            frame.frames.primary.crop(350, 725, 300, 825)

            # Save the convolved, rebinned and cropped bulge or disk frame
            component_path = os.path.join(self.prep_path, name, 'final.fits')
            frame.save(component_path)

# -----------------------------------------------------------------
