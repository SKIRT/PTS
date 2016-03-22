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

# Import astronomical modules
from astroquery.vizier import Vizier
from astropy.units import Unit, dimensionless_angles
from astropy.coordinates import Angle

# Import the relevant PTS classes and modules
from .component import DecompositionComponent
from ...core.basics.map import Map
from ...core.tools import inspection, filesystem
from ...core.tools.logging import log
from ...core.simulation.skifile import SkiFile
from ..preparation import unitconversion
from ...core.basics.errorbar import ErrorBar
from ...core.simulation.arguments import SkirtArguments
from ...core.simulation.execute import SkirtExec
from ...magic.tools import catalogs
from ...magic.basics.vector import Extent, Position
from ...magic.basics.skygeometry import SkyEllipse, SkyCoordinate
from ...magic.basics.skyregion import SkyRegion
from ...magic.core.frame import Frame
from ...magic.basics.coordinatesystem import CoordinateSystem

# -----------------------------------------------------------------

# Online table url
s4g_decomposition_table_link = "http://www.oulu.fi/astronomy/S4G_PIPELINE4/s4g_p4_table8.dat"

# Local table path
local_table_path = filesystem.join(inspection.pts_dat_dir("modeling"), "s4g", "s4g_p4_table8.dat")

# -----------------------------------------------------------------

# The path to the template ski files directory
template_path = filesystem.join(inspection.pts_dat_dir("modeling"), "ski")

# -----------------------------------------------------------------

# Reference image
reference_image = "Pacs red"

# -----------------------------------------------------------------

class GalaxyDecomposer(DecompositionComponent):
    
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
        self.bulge_directory = None
        self.disk_directory = None

        # The bulge and disk parameters
        self.parameters = Map()

        # The SKIRT execution context
        self.skirt = SkirtExec()

        # The bulge and disk image
        self.bulge = None
        self.disk = None

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
        self.bulge_directory = self.full_output_path("bulge")
        self.disk_directory = self.full_output_path("disk")

        # Create the bulge and disk directories
        filesystem.create_directories([self.bulge_directory, self.disk_directory])

    # -----------------------------------------------------------------

    def get_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the structural galaxy parameters from the S4G catalog ...")

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
        center = SkyCoordinate(ra=ra_center, dec=dec_center, unit="deg", frame='fk5')
        self.parameters.center = center

        # Distance
        self.parameters.distance = table["Dmean"][0] * Unit("Mpc")
        self.parameters.distance_error = table["e_Dmean"][0] * Unit("Mpc")

        # Major axis, ellipticity, position angle
        self.parameters.major_arcsec = table["amaj"][0] * Unit("arcsec")
        self.parameters.major = (self.parameters.distance * self.parameters.major_arcsec).to("pc", equivalencies=dimensionless_angles())

        # Ellipticity
        self.parameters.ellipticity = table["ell"][0]
        self.parameters.position_angle = Angle(table["PA"][0] + 90.0, Unit("deg"))

        # Magnitudes
        asymptotic_ab_magnitude_i1 = table["__3.6_"][0]
        asymptotic_ab_magnitude_i2 = table["__4.5_"][0]
        asymptotic_ab_magnitude_i1_error = table["e__3.6_"][0]
        asymptotic_ab_magnitude_i2_error = table["e__4.5_"][0]

        self.parameters.i1_mag = asymptotic_ab_magnitude_i1
        self.parameters.i1_mag_error = asymptotic_ab_magnitude_i1_error
        self.parameters.i2_mag = asymptotic_ab_magnitude_i2
        self.parameters.i2_mag_error = asymptotic_ab_magnitude_i2_error

        self.parameters.i1_fluxdensity = unitconversion.ab_to_jansky(self.parameters.i1_mag)
        i1_fluxdensity_lower = unitconversion.ab_to_jansky(self.parameters.i1_mag + self.parameters.i1_mag_error)
        i1_fluxdensity_upper = unitconversion.ab_to_jansky(self.parameters.i1_mag - self.parameters.i1_mag_error)
        i1_error = ErrorBar(i1_fluxdensity_lower, i1_fluxdensity_upper, at=self.parameters.i1_fluxdensity)
        self.parameters.i1_error = i1_error.average

        self.parameters.i2_fluxdensity = unitconversion.ab_to_jansky(self.parameters.i2_mag)
        i2_fluxdensity_lower = unitconversion.ab_to_jansky(self.parameters.i2_mag + self.parameters.i2_mag_error)
        i2_fluxdensity_upper = unitconversion.ab_to_jansky(self.parameters.i2_mag - self.parameters.i2_mag_error)
        i2_error = ErrorBar(i2_fluxdensity_lower, i2_fluxdensity_upper, at=self.parameters.i2_fluxdensity)
        self.parameters.i2_error = i2_error.average

        # Other ...
        #absolute_magnitude_i1 = table["M3.6"][0]
        #absolute_magnitude_i2 = table["M4.5"][0]
        #stellar_mass = 10.0**table["logM_"][0] * u.Unit("Msun")

        # Inform the user
        log.info("Querying the catalog of radial profiles for 161 face-on spirals ...")

        # Radial profiles for 161 face-on spirals (Munoz-Mateos+, 2007)
        radial_profiles_result = vizier.query_object(self.galaxy_name, catalog="J/ApJ/658/1006")

        distance = float(radial_profiles_result[0][0]["Dist"])
        inclination = Angle(float(radial_profiles_result[0][0]["i"]), "deg")

        # Set the inclination
        self.parameters.inclination = inclination

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
                self.parameters.bulge.interpretation = splitted[sersic_1_index].split("|")[1]
                self.parameters.bulge.rel = float(splitted[sersic_1_index + 1])
                self.parameters.bulge.mag = float(splitted[sersic_1_index + 2])
                self.parameters.bulge.fluxdensity = unitconversion.ab_to_jansky(self.parameters.bulge.mag) * Unit("Jy")

                # Effective radius in pc
                re_arcsec = float(splitted[sersic_1_index + 3]) * Unit("arcsec")
                self.parameters.bulge.re = (self.parameters.distance * re_arcsec).to("pc", equivalencies=dimensionless_angles())

                self.parameters.bulge.ar = float(splitted[sersic_1_index + 4])
                self.parameters.bulge.pa = float(splitted[sersic_1_index + 5])
                self.parameters.bulge.n = float(splitted[sersic_1_index + 6])

                # DISK
                self.parameters.disk = Map()
                self.parameters.disk.interpretation = splitted[disk_1_index].split("|")[1]
                self.parameters.disk.rel = float(splitted[disk_1_index + 1])
                self.parameters.disk.mag = float(splitted[disk_1_index + 2])
                self.parameters.disk.fluxdensity = unitconversion.ab_to_jansky(self.parameters.disk.mag) * Unit("Jy")

                # Scale length in pc
                hr_arcsec = float(splitted[disk_1_index + 3]) * Unit("arcsec")
                self.parameters.disk.hr = (self.parameters.distance * hr_arcsec).to("pc", equivalencies=dimensionless_angles())

                self.parameters.disk.ar = float(splitted[disk_1_index + 4]) # axial ratio
                self.parameters.disk.pa = float(splitted[disk_1_index + 5])
                self.parameters.disk.mu0 = float(splitted[disk_1_index + 6])

    # -----------------------------------------------------------------

    def simulate(self):

        """
        This function ...
        :return:
        """

        # Get the parameters describing the pixel grid of the prepared images
        reference_path = filesystem.join(self.prep_path, reference_image, "result.fits")
        reference_wcs = CoordinateSystem.from_file(reference_path)

        # Simulate the stellar bulge
        self.simulate_bulge(reference_wcs)

        # Simulate the stellar disk
        self.simulate_disk(reference_wcs)

    # -----------------------------------------------------------------

    def simulate_bulge(self, reference_wcs):

        """
        This function ...
        :param reference_wcs:
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the bulge image ...")

        # Load the bulge ski file template
        bulge_template_path = filesystem.join(template_path, "bulge.ski")
        ski = SkiFile(bulge_template_path)

        # Change the ski file parameters
        # component_id, index, radius, y_flattening=1, z_flattening=1
        ski.set_stellar_component_sersic_geometry(0, self.parameters.bulge.n, self.parameters.bulge.re, z_flattening=self.parameters.bulge.ar)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add a new SimpleInstrument
        distance = self.parameters.distance
        inclination = self.parameters.inclination
        azimuth = 0.0
        position_angle = self.parameters.position_angle
        pixels_x = reference_wcs.xsize
        pixels_y = reference_wcs.ysize
        #center_x = reference_wcs.center_pixel.x
        #center_y = reference_wcs.center_pixel.y
        pixel_center = self.parameters.center.to_pixel(reference_wcs)
        center = Position(0.5*pixels_x - pixel_center.x, 0.5*pixels_y - pixel_center.y)
        field_x_angular = reference_wcs.pixelscale.x.to("deg/pix") * pixels_x * Unit("pix")
        field_y_angular = reference_wcs.pixelscale.y.to("deg/pix") * pixels_y * Unit("pix")
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())
        ski.add_simple_instrument("earth", distance, inclination, azimuth, position_angle, field_x_physical, field_y_physical, pixels_x, pixels_y, center.x, center.y)

        # Determine the path to the ski file
        ski_path = filesystem.join(self.bulge_directory, "bulge.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = filesystem.join(self.bulge_directory, "out")

        # Create the output directory
        filesystem.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the bulge simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug else True)

        # Determine the path to the output FITS file
        bulge_image_path = filesystem.join(out_path, "bulge_earth_total.fits")

        # Check if the output contains the "bulge_earth_total.fits" file
        if not filesystem.is_file(bulge_image_path):
            raise RuntimeError("Something went wrong with the bulge simulation: output FITS file missing")

        # Open the bulge image
        self.bulge = Frame.from_file(bulge_image_path)

        # Convolve the frame to the PACS 160 resolution
        #kernels_dir = os.path.expanduser("~/Kernels")
        #kernel_path = os.path.join(kernels_dir, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")
        #kernel = Frame.from_file(kernel_path)
        #frame.convolve(kernel)

        # Rebin the convolved frame to the PACS 160 frame (just as we did with the other images)
        #reference_path = os.path.join(self.data_path, "PACS160.fits")
        #reference = Frame.from_file(reference_path)
        #frame.rebin(reference)

        # Finally, crop the image
        #frame.frames.primary.crop(350, 725, 300, 825)

    # -----------------------------------------------------------------

    def simulate_disk(self, reference_wcs):

        """
        This function ...
        :param reference_wcs:
        :return:
        """

        # Inform the user
        log.info("Creating ski file to simulate the disk image ...")

        # Load the disk ski file template
        disk_template_path = filesystem.join(template_path, "disk.ski")
        ski = SkiFile(disk_template_path)

        # Change the ski file parameters ...
        radial_scale = self.parameters.disk.hr
        axial_scale = self.parameters.disk.ar * radial_scale
        ski.set_stellar_component_expdisk_geometry(0, radial_scale, axial_scale, radial_truncation=0, axial_truncation=0, inner_radius=0)

        # Remove all existing instruments
        ski.remove_all_instruments()

        # Add a new SimpleInstrument
        distance = self.parameters.distance
        inclination = self.parameters.inclination
        azimuth = 0.0
        position_angle = self.parameters.position_angle
        pixels_x = reference_wcs.xsize
        pixels_y = reference_wcs.ysize
        #center_x = reference_wcs.center_pixel.x
        #center_y = reference_wcs.center_pixel.y
        pixel_center = self.parameters.center.to_pixel(reference_wcs)
        center = Position(0.5*pixels_x - pixel_center.x, 0.5*pixels_y - pixel_center.y)
        field_x_angular = reference_wcs.pixelscale.x.to("deg/pix") * pixels_x * Unit("pix")
        field_y_angular = reference_wcs.pixelscale.y.to("deg/pix") * pixels_y * Unit("pix")
        field_x_physical = (field_x_angular * distance).to("pc", equivalencies=dimensionless_angles())
        field_y_physical = (field_y_angular * distance).to("pc", equivalencies=dimensionless_angles())
        ski.add_simple_instrument("earth", distance, inclination, azimuth, position_angle, field_x_physical, field_y_physical, pixels_x, pixels_y, center.x, center.y)

        # Determine the path to the ski file
        ski_path = filesystem.join(self.disk_directory, "disk.ski")

        # Save the ski file to the new path
        ski.saveto(ski_path)

        # Determine the path to the simulation output directory
        out_path = filesystem.join(self.disk_directory, "out")

        # Create the output directory
        filesystem.create_directory(out_path)

        # Create a SkirtArguments object
        arguments = SkirtArguments()

        # Adjust the parameters
        arguments.ski_pattern = ski_path
        arguments.output_path = out_path
        arguments.single = True   # we expect a single simulation from the ski pattern

        # Inform the user
        log.info("Running the disk simulation ...")

        # Run the simulation
        simulation = self.skirt.run(arguments, silent=False if log.is_debug else True)

        # Determine the path to the output FITS file
        disk_image_path = filesystem.join(out_path, "disk_earth_total.fits")

        # Check if the output contains the "disk_earth_total.fits" file
        if not filesystem.is_file(disk_image_path):
            raise RuntimeError("Something went wrong with the disk simulation: output FITS file missing")

        # Open the disk image
        self.disk = Frame.from_file(disk_image_path)

        # Convolve the frame to the PACS 160 resolution
        #kernels_dir = os.path.expanduser("~/Kernels")
        #kernel_path = os.path.join(kernels_dir, "Kernel_HiRes_Moffet_00.5_to_PACS_160.fits")
        #kernel = Frame.from_file(kernel_path)
        #frame.convolve(kernel)

        # Rebin the convolved frame to the PACS 160 frame (just as we did with the other images)
        #reference_path = os.path.join(self.data_path, "PACS160.fits")
        #reference = Frame.from_file(reference_path)
        #frame.rebin(reference)

        # Finally, crop the image
        #frame.frames.primary.crop(350, 725, 300, 825)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write out the final bulge and disk images
        self.write_images()

        # Write out the parameters in a data file
        self.write_parameters()

        # Write out the disk ellipse
        self.write_disk_ellipse()

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the bulge image and save it
        final_bulge_path = filesystem.join(self.components_path, "bulge.fits")
        self.bulge.save(final_bulge_path)

        # Determine the path to the disk image and save it
        final_disk_path = filesystem.join(self.components_path, "disk.fits")
        self.disk.save(final_disk_path)

    # -----------------------------------------------------------------

    def write_parameters(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing data file with parameters ...")

        # Determine the full path to the parameters file
        path = self.full_output_path("parameters.dat")

        # Write the parameters to the specified location
        write_parameters(self.parameters, path)

    # -----------------------------------------------------------------

    def write_disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing regions file with disk ellipse ...")

        minor = (1.0 - self.parameters.ellipticity) * self.parameters.major_arcsec

        # Ellipse radius
        radius = Extent(self.parameters.major_arcsec, minor)

        # Create sky ellipse
        sky_ellipse = SkyEllipse(self.parameters.center, radius, self.parameters.position_angle)

        # Create region
        region = SkyRegion()
        region.append(sky_ellipse)
        region_path = self.full_output_path("disk.reg")
        region.save(region_path)

# -----------------------------------------------------------------

def write_parameters(parameters, path):

    """
    This function ...
    :param parameters:
    :param path:
    :return:
    """

    # Create the parameters file
    with open(path, 'w') as parameter_file:

        # Add general info
        print("Name:", parameters.galaxy_name, file=parameter_file)
        print("Center RA:", str(parameters.center.ra.to("deg").value) + " deg", file=parameter_file)
        print("Center DEC:", str(parameters.center.dec.to("deg").value) + " deg", file=parameter_file)
        print("Major axis length:", str(parameters.major), file=parameter_file)
        print("Ellipticity:", parameters.ellipticity, file=parameter_file)
        print("Position angle:", str(parameters.position_angle.to("deg").value) + " deg", file=parameter_file)
        print("Distance:", str(parameters.distance), file=parameter_file)
        print("Distance error:", str(parameters.distance_error), file=parameter_file)
        print("Inclination:", str(parameters.inclination.to("deg").value) + " deg", file=parameter_file)
        print("IRAC 3.6um flux density:", parameters.i1_fluxdensity, file=parameter_file)
        print("IRAC 3.6um flux density error:", parameters.i1_error, file=parameter_file)
        print("IRAC 4.5um flux density:", parameters.i2_fluxdensity, file=parameter_file)
        print("IRAC 4.5um flux density error:", parameters.i2_error, file=parameter_file)

        #print("Model type:", parameters.model_type, file=parameter_file)
        #print("Number of components:", parameters.number_of_components, file=parameter_file)
        #print("Quality:", parameters.quality, file=parameter_file)

        # Add components parameters
        for component in ["bulge", "disk"]:

            #print(component.title() + ": Interpretation:", parameters[component].interpretation, file=parameter_file)
            print(component.title() + ": Relative contribution:", parameters[component].rel, file=parameter_file)
            #print(component.title() + ": Total IRAC 3.6um AB magnitude:", parameters[component].mag, file=parameter_file)
            print(component.title() + ": IRAC 3.6um flux density:", parameters[component].fluxdensity, file=parameter_file)
            print(component.title() + ": Axial ratio:", parameters[component].ar, file=parameter_file)
            print(component.title() + ": Position angle:", parameters[component].pa, file=parameter_file) # (degrees ccw from North)

            if component == "bulge":

                print(component.title() + ": Effective radius:", str(parameters[component].re), file=parameter_file)
                print(component.title() + ": Sersic index:", parameters[component].n, file=parameter_file)

            elif component == "disk":

                print(component.title() + ": Central surface brightness:", parameters[component].mu0, file=parameter_file) # (mag/arcsec2)
                print(component.title() + ": Exponential scale length:", str(parameters[component].hr), file=parameter_file)

# -----------------------------------------------------------------

def load_parameters(path):

    """
    This function ...
    :param path:
    :return:
    """

    # Create parameters structure
    parameters = Map()
    parameters.bulge = Map()
    parameters.disk = Map()

    ra = None
    dec = None

    # Read the parameter file
    with open(path, 'r') as parameter_file:

        # Loop over all lines in the file
        for line in parameter_file:

            splitted = line.split(": ")

            # Bulge parameters
            if splitted[0] == "Bulge":

                if splitted[1] == "Relative contribution": parameters.bulge.rel = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.bulge.fluxdensity = get_quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.bulge.ar = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.bulge.pa = get_angle(splitted[2])
                elif splitted[1] == "Effective radius": parameters.bulge.re = get_quantity(splitted[2])
                elif splitted[1] == "Sersic index": parameters.bulge.n = float(splitted[2])

            # Disk parameters
            elif splitted[0] == "Disk":

                if splitted[1] == "Relative contribution": parameters.disk.rel = float(splitted[2])
                elif splitted[1] == "IRAC 3.6um flux density": parameters.disk.fluxdensity = get_quantity(splitted[2])
                elif splitted[1] == "Axial ratio": parameters.disk.ar = float(splitted[2])
                elif splitted[1] == "Position angle": parameters.disk.pa = get_angle(splitted[2])
                elif splitted[1] == "Central surface brightness": parameters.disk.mu0 = float(splitted[2])
                elif splitted[1] == "Exponential scale length": parameters.disk.hr = get_quantity(splitted[2])

            # Other parameters
            elif len(splitted) == 2:

                if splitted[0] == "Name": parameters.galaxy_name = splitted[1]
                elif splitted[0] == "Center RA": ra = get_quantity(splitted[1])
                elif splitted[0] == "Center DEC": dec = get_quantity(splitted[1])
                elif splitted[0] == "Major axis length": parameters.major = get_quantity(splitted[1])
                elif splitted[0] == "Ellipticity": parameters.ellipticity = float(splitted[1])
                elif splitted[0] == "Position angle": parameters.position_angle = get_angle(splitted[1])
                elif splitted[0] == "Distance": parameters.distance = get_quantity(splitted[1])
                elif splitted[0] == "Distance error": parameters.distance_error = get_quantity(splitted[1])
                elif splitted[0] == "Inclination": parameters.inclination = get_angle(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density": parameters.i1_fluxdensity = float(splitted[1])
                elif splitted[0] == "IRAC 3.6um flux density error": parameters.i1_error = float(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density": parameters.i2_fluxdensity = float(splitted[1])
                elif splitted[0] == "IRAC 4.5um flux density error": parameters.i2_error = float(splitted[1])

    # Add the center coordinate
    parameters.center = SkyCoordinate(ra=ra, dec=dec)

    # Return the parameters
    return parameters

# -----------------------------------------------------------------

def get_quantity(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create a quantity object and return it
    if unit is not None: value = value * Unit(unit)
    return value

# -----------------------------------------------------------------

def get_angle(entry, default_unit=None):

    """
    This function ...
    :param entry:
    :param default_unit:
    :return:
    """

    splitted = entry.split()
    value = float(splitted[0])
    try: unit = splitted[1]
    except IndexError: unit = default_unit

    # Create an Angle object and return it
    if unit is not None: value = Angle(value, unit)
    return value

# -----------------------------------------------------------------
