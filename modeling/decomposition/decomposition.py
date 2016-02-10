#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.decomposition.decomposition Contains the GalaxyDecomposer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import os

# Import the relevant PTS classes and modules
from ..core import ModelingComponent
from ...core.basics.map import Map
from ...core.tools import inspection

# -----------------------------------------------------------------

# Bulge properties from S4G
bulge = Map()
bulge.type = "Sersic"
bulge.xc = 1455.27
bulge.yc = 1905.39
bulge.mag = 7.0877
bulge.re = 102.9578
bulge.n = 3.5566
bulge.ar = 0.6540
bulge.pa = -34.9801

# Disk properties from S4G
disk = Map()
disk.type = "ExpDisk"
disk.mag = 6.9121
disk.rs = 204.2269
disk.ar = 0.5426
disk.pa = -23.6950

# -----------------------------------------------------------------

s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

local_table_path = os.path.join(inspection.pts_dat_dir("modeling"), "s4G", "s4g_p4_table8.dat")

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

        # 2. Get the parameters
        self.get_parameters()

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

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        #The table columns are:
        # (1) the running number (1-2352),
        # (2) the galaxy name,
        # (3) the type of final decomposition model,
        # (4) the number N of components in the final model , and
        # (5) the quality flag Q.

        # 6- 8 for unresolved central component ('psf'; phys, frel,  mag),
        # 9-15 for inner sersic-component ('sersic1'; phys,  frel,  mag,    re,    ar,      pa,    n),
        # 16-22 for inner disk-component ('expo1'; phys,  frel,  mag,    hr,    ar,      pa,   mu0),
        # 23-28 for inner ferrers-component ('ferrers1'; phys,  frel,  mu0,   rout,   ar,     pa),
        # 29-34 for inner edge-on disk component ('edgedisk1';  phys, frel, mu0,  rs,    hs,      pa ).
        # 35-41 for outer sersic-component ('sersic2'; phys,  frel,  mag,    re,    ar,      pa,    n),
        # 42-49 for outer disk-component ('expo2'; phys,  frel,  mag,    hr,    ar,      pa,   mu0),
        # 50-55 for outer ferrers-component ('ferrers2'; phys,  frel,  mu0,   rout,   ar,     pa),
        # 56-61 for outer edge-on disk component ('edgedisk2';  phys, frel, mu0,  rs,    hs,      pa ).

        sersic_1_index = 8

        #For each function:

        #the first entry stands for the physical intepretation of the component:
        #'N' for a central source (or unresolved bulge), 'B' for a bulge (or elliptical), 'D' for a disk, 'BAR' for a bar, and 'Z' for an edge-on disk.

        # 'rel    =  the relative contribution of the component to the total model flux,
        #  mag    =  the component's total 3.6 micron AB magnitude,
        #  mu0    =  the central  surface brightness (mag/arcsec^2;  de-projected central surface brightness for expdisk and edgedisk, and
        #                                                          sky-plane central surface brightness for ferrer)
        #  ar     =  axial ratio
        #  pa     =  position angle (degrees ccw from North)
        #  n      =  sersic index
        #  hr     =  exponential scale lenght (arcsec)
        #  rs     =  radial scale lenght (arcsec)
        #  hs     =  vertical scale height (arcsec)
        #  rout   =  bar outer truncation radius (arcsec)

        with open(local_table_path, 'r') as s4g_table:

            for line in s4g_table:

                splitted = line.split()

                name = splitted[1]

                if name != "NGC3031": continue

                model_type = splitted[2]
                number_of_components = splitted[3]
                quality = splitted[4]

                sersic_1_physical_interpretation = splitted[sersic_1_index]
                sersic_1_rel = splitted[sersic_1_index+1]
                sersic_1_mag = splitted[sersic_1_index+2]
                sersic_1_mu0 = splitted[sersic_1_index+3]
                sersic_1_ar = splitted[sersic_1_index+4]
                sersic_1_pa = splitted[sersic_1_index+5]
                sersic_1_n = splitted[sersic_1_index+6]
                sersic_1_hr = splitted[sersic_1_index+7]
                sersic_1_rs = splitted[sersic_1_index+8]
                sersic_1_hs = splitted[sersic_1_index+9]
                sersic_1_rout = splitted[sersic_1_index+10]

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
