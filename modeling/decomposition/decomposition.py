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

# Import astronomical modules
from astroquery.vizier import Vizier
from astropy.units import Unit
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from ..core import ModelingComponent
from ...core.basics.map import Map
from ...core.tools import inspection, filesystem
from ...core.tools.logging import log
from ...core.simulation.skifile import SkiFile

# Import the relevant AstroMagic classes and modules
from ...magic.tools import catalogs
from ...magic.basics import Extent
from ...magic.basics.skygeometry import SkyEllipse, SkyCoord
from ...magic.basics import SkyRegion

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

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = os.path.join(inspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

template_path = os.path.join(inspection.pts_dat_dir("modeling"), "ski")

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

        # The NGC name of the galaxy
        self.ngc_id = None
        self.ngc_id_nospaces = None

        # The path to the disk and bulge directories
        self.bulge_path = None
        self.disk_path = None

        # The bulge and disk parameters
        self.parameters = Map()

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new GalaxyDecomposer instance
        decomposer = cls(arguments.config)

        # Set the input and output path
        decomposer.config.path = arguments.path
        decomposer.config.input_path = os.path.join(arguments.path, "prep")
        decomposer.config.output_path = os.path.join(arguments.path, "components")

        # Return the new instance
        return decomposer

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

        # 2. Simulate the bulge and disk images
        self.simulate()

        # 3. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyDecomposer, self).setup()

        # Get the NGC name of the galaxy
        self.ngc_id = catalogs.get_ngc_name(self.galaxy_name)
        self.ngc_id_nospaces = self.ngc_id.replace(" ", "")

        # Determine the path to the bulge and disk directories
        self.bulge_path = self.full_output_path("bulge")
        self.disk_path = self.full_output_path("disk")

    # -----------------------------------------------------------------

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        # Query the S4G catalog using Vizier for general parameters
        self.get_general_parameters()

        # Parse the S4G table 8 to get the decomposition parameters
        self.get_component_parameters()

    # -----------------------------------------------------------------

    def get_general_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Querying the S4G catalog ...")

        vizier = Vizier(keywords=["galaxies"])

        # Get parameters from S4G catalog
        result = vizier.query_object(self.galaxy_name, catalog=["J/PASP/122/1397/s4g"])
        table = result[0]

        # Galaxy name for S4G catalog
        self.parameters.galaxy_name = table["Name"][0]

        # Galaxy center from decomposition (?)
        ra_center = table["_RAJ2000"][0]
        dec_center = table["_DEJ2000"][0]
        center = SkyCoord(ra=ra_center, dec=dec_center, unit=(Unit("deg"), Unit("deg")), frame='fk5')
        self.parameters.center = center

        # Major axis, ellipticity, position angle
        self.parameters.major = table["amaj"][0] * Unit("arcsec")
        self.parameters.ellipticity = table["ell"][0]
        self.parameters.position_angle = Angle(table["PA"][0] - 90.0, Unit("deg"))

        # Distance
        self.parameters.distance = table["Dmean"][0] * Unit("Mpc")
        self.parameters.distance_error = table["e_Dmean"][0] * Unit("Mpc")

        # Magnitudes
        asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        self.parameters.i1_mag = asymptotic_ab_magnitude_i1
        self.parameters.i1_mag_error = asymptotic_ab_magnitude_i1_error
        self.parameters.i2_mag = asymptotic_ab_magnitude_i2
        self.parameters.i2_mag_error = asymptotic_ab_magnitude_i2_error

        # Other ...
        #absolute_magnitude_i1 = table["M3.6"][0]
        #absolute_magnitude_i2 = table["M4.5"][0]
        #stellar_mass = 10.0**table["logM_"][0] * u.Unit("Msun")

    # -----------------------------------------------------------------

    def get_component_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Parsing S4G table 8 to get the decomposition parameters ...")

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
        disk_1_index = 15

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

        # SERSIC: phys,  frel,  mag,    re,    ar,      pa,    n
        # DISK: phys,  frel,  mag,    hr,    ar,      pa,   mu0

        with open(local_table_path, 'r') as s4g_table:

            for line in s4g_table:

                splitted = line.split()

                if len(splitted) < 2: continue

                name = splitted[1]

                # Only look at the line corresponding to the galaxy
                if name != self.ngc_id_nospaces: continue

                self.parameters.model_type = splitted[2]
                self.parameters.number_of_components = splitted[3]
                self.parameters.quality = splitted[4]

                # BULGE
                self.parameters.bulge = Map()
                self.parameters.bulge.interpretation = splitted[sersic_1_index].split("|")[1] # COMMON with disk
                self.parameters.bulge.rel = splitted[sersic_1_index + 1] # COMMON with disk
                self.parameters.bulge.mag = splitted[sersic_1_index + 2] # COMMON with disk
                self.parameters.bulge.re = splitted[sersic_1_index + 3] # re ??? effective radius?
                self.parameters.bulge.ar = splitted[sersic_1_index + 4] # COMMON with disk
                self.parameters.bulge.pa = splitted[sersic_1_index + 5] # COMMON with disk
                self.parameters.bulge.n = splitted[sersic_1_index + 6]

                # DISK
                self.parameters.disk = Map()
                self.parameters.disk.interpretation = splitted[disk_1_index].split("|")[1] # COMMON
                self.parameters.disk.rel = splitted[disk_1_index + 1] # COMMON
                self.parameters.disk.mag = splitted[disk_1_index + 2] # COMMON
                self.parameters.disk.hr = splitted[disk_1_index + 3]
                self.parameters.disk.ar = splitted[disk_1_index + 4] # COMMON
                self.parameters.disk.pa = splitted[disk_1_index + 5] # COMMON
                self.parameters.disk.mu0 = splitted[disk_1_index + 6]

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Simulate the stellar bulge
        self.simulate_bulge()

        # Simulate the stellar disk
        self.simulate_disk()

    # -----------------------------------------------------------------

    def simulate_bulge(self):

        """
        This function ...
        :return:
        """

        # Load the bulge ski file template
        bulge_template_path = filesystem.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

    # -----------------------------------------------------------------

    def simulate_disk(self):

        """
        This function ...
        :return:
        """

        # Load the disk ski file template
        disk_template_path = filesystem.join(template_path, "disk.ski")
        ski = SkiFile(disk_template_path)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Write out the parameters in a data file
        self.write_parameters()

        # Write out the disk ellipse
        self.write_disk_ellipse()

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Determine the full path to the parameters file
        path = self.full_output_path("parameters.dat")

        # Create the parameters file
        with open(path, 'w') as parameter_file:

            # Add general info
            print("Galaxy name:", self.parameters.galaxy_name, file=parameter_file)
            print("Center:", str(self.parameters.center), file=parameter_file)
            print("Major axis length:", str(self.parameters.major), file=parameter_file)
            print("Ellipticity:", self.parameters.ellipticity, file=parameter_file)
            print("Position angle:", self.parameters.position_angle, file=parameter_file)
            print("Distance:", str(self.parameters.distance) + " +/- " + str(self.parameters.distance_error), file=parameter_file)
            print("I1 magnitude:", self.parameters.i1_mag, "+/-", self.parameters.i1_mag_error, file=parameter_file)
            print("I2 magnitude:", self.parameters.i2_mag, "+/-", self.parameters.i2_mag_error, file=parameter_file)

            print("Model type:", self.parameters.model_type, file=parameter_file)
            print("Number of components:", self.parameters.number_of_components, file=parameter_file)
            print("Quality:", self.parameters.quality, file=parameter_file)

            # Add components parameters
            for component in ["bulge", "disk"]:

                print(component.title()+":", file=parameter_file)
                print("    Interpretation:", self.parameters[component].interpretation, file=parameter_file)
                print("    Relative contribution of the component to the total model flux:", self.parameters[component].rel, file=parameter_file)
                print("    Total 3.6 micron AB magnitude:", self.parameters[component].mag, file=parameter_file)
                print("    Axial ratio:", self.parameters[component].ar, file=parameter_file)
                print("    Position angle (degrees ccw from North):", self.parameters[component].pa, file=parameter_file)

                if component == "bulge":

                    print("    Effective radius ???:", self.parameters[component].re, file=parameter_file)
                    print("    Sersic index:", self.parameters[component].n, file=parameter_file)

                elif component == "disk":

                    print("    Central surface brightness (mag/arcsec2):", self.parameters[component].mu0, file=parameter_file)
                    print("    Exponential scale lenght (arcsec):", self.parameters[component].hr, file=parameter_file)

    # -----------------------------------------------------------------

    def write_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        minor = (1.0 - self.parameters.ellipticity) * self.parameters.major

        # Ellipse radius
        radius = Extent(self.parameters.major, minor)

        # Create sky ellipse
        sky_ellipse = SkyEllipse(self.parameters.center, radius, self.parameters.position_angle)

        # Create region
        region = SkyRegion()
        region.append(sky_ellipse)
        region_path = self.full_output_path("disk.reg")
        region.save(region_path)

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
