#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization Contains the InputInitializer

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import astronomical modules
from astropy.units import Unit, dimensionless_angles

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.tools import inspection, tables, filesystem
from ...core.simulation.skifile import SkiFile
from ...core.basics.filter import Filter
from ..basics.models import SersicModel, DeprojectionModel
from ...magic.basics.vector import Position
from ..decomposition.decomposition import load_parameters

# -----------------------------------------------------------------

template_ski_path = filesystem.join(inspection.pts_dat_dir("modeling"), "ski", "template.ski")

# -----------------------------------------------------------------

class InputInitializer(FittingComponent):
    
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
        super(InputInitializer, self).__init__(config)

        # -- Attributes --

        # The ski file
        self.ski = None

        # The wavelength grid
        self.wavelength_grid = None

        # The structural parameters
        self.parameters = None

        # The geometric bulge model
        self.bulge = None

        # The deprojection model
        self.deprojection = None

        # The fluxes table
        self.fluxes = None

        # Filters
        self.i1 = None
        self.fuv = None

    # -----------------------------------------------------------------

    @classmethod
    def from_arguments(cls, arguments):

        """
        This function ...
        :param arguments:
        :return:
        """

        # Create a new InputInitializer instance
        initializer = cls(arguments.config)

        # Set the output path
        initializer.config.path = arguments.path

        # Set minimum and maximum wavelength of the total grid
        if arguments.lambda_minmax is not None:
            initializer.config.wavelengths.min = arguments.lambda_minmax[0]
            initializer.config.wavelengths.max = arguments.lambda_minmax[1]

        # Set minimum and maximum wavelength of the zoomed-in grid
        if arguments.lambda_minmax_zoom is not None:
            initializer.config.wavelengths.min_zoom = arguments.lambda_minmax_zoom[0]
            initializer.config.wavelengths.max_zoom = arguments.lambda_minmax_zoom[1]

        # Set the number of wavelength points
        if arguments.nlambda is not None:
            # Based on npoints = 1.1 * npoints_zoom
            initializer.config.wavelengths.npoints = 1.1 * arguments.nlambda / 2.1
            initializer.config.wavelengths.npoints_zoom = arguments.nlambda / 2.1

        # Set the number of photon packages per wavelength
        if arguments.packages is not None: initializer.config.packages = arguments.packages

        # Set selfabsorption
        initializer.config.selfabsorption = arguments.selfabsorption

        # Return the new instance
        return initializer

    # -----------------------------------------------------------------

    def run(self):

        """
        This function ...
        :return:
        """

        # 1. Call the setup function
        self.setup()

        # 2. Load the template ski file
        self.load_template()

        # 3. Load the structural parameters for the galaxy
        self.load_parameters()

        # Load the fluxes
        self.load_fluxes()

        # 4. Create the wavelength grid
        self.create_wavelength_grid()

        # Create the bulge model
        self.create_bulge_model()

        # 5. Adjust the ski file
        self.adjust_ski()

        # 4. Place the input
        self.place_input()

        # 5. Place ski file
        self.place_ski_file()

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(InputInitializer, self).setup()

        # Create filters
        self.i1 = Filter.from_string("I1")
        self.fuv = Filter.from_string("FUV")

    # -----------------------------------------------------------------

    def load_template(self):

        """
        This function ...
        :return:
        """

        # Open the template ski file
        self.ski = SkiFile(template_ski_path)

    # -----------------------------------------------------------------

    def load_parameters(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the parameters file
        path = filesystem.join(self.components_path, "parameters.dat")

        # Load the parameters
        self.parameters = load_parameters(path)

    # -----------------------------------------------------------------

    def load_fluxes(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the fluxes table
        fluxes_path = filesystem.join(self.phot_path, "fluxes.dat")

        # Load the fluxes table
        self.fluxes = tables.from_file(fluxes_path)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Verify the grid parameters
        if self.config.wavelengths.npoints < 2: raise ValueError("the number of points in the low-resolution grid should be at least 2")
        if self.config.wavelengths.npoints_zoom < 2: raise ValueError("the number of points in the high-resolution subgrid should be at least 2")
        if self.config.wavelengths.min <= 0: raise ValueError("the shortest wavelength should be positive")
        if (self.config.wavelengths.min_zoom <= self.config.wavelengths.min
            or self.config.wavelengths.max_zoom <= self.config.wavelengths.min_zoom
            or self.config.wavelengths.max <= self.config.wavelengths.max_zoom):
                raise ValueError("the high-resolution subgrid should be properly nested in the low-resolution grid")

        # Build the high- and low-resolution grids independently
        base_grid = np.logspace(float(self.config.wavelengths.min), float(self.config.wavelengths.max), num=self.config.wavelengts.npoints, endpoint=True, base=10.)
        zoom_grid = np.logspace(float(self.config.wavelengths.min_zoom), float(self.config.wavelengths.max_zoom), num=self.config.wavelengths.npoints_zoom, endpoint=True, base=10.)

        # Merge the two grids
        total_grid = []

        # Add the wavelengths of the low-resolution grid before the first wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength < self.config.wavelengths.min_zoom: total_grid.append(wavelength)

        # Add the wavelengths of the high-resolution grid
        for wavelength in zoom_grid: total_grid.append(wavelength)

        # Add the wavelengths of the low-resolution grid after the last wavelength of the high-resolution grid
        for wavelength in base_grid:
            if wavelength > self.config.wavelengths.max_zoom: total_grid.append(wavelength)

        # -- Add the wavelengths of the bands of interest --

        # Open the fluxes.dat table to get the filters that are used for the SED
        fluxes_table_path = filesystem.join(self.phot_path, "fluxes.dat")
        fluxes_table = tables.from_file(fluxes_table_path)

        # Loop over the entries in the fluxes table, get the filter
        for entry in fluxes_table:

            # Get the filter
            filter_id = entry["Instrument"] + "." + entry["Band"]
            filter = Filter.from_string(filter_id)

            # Get the wavelength in micron
            wavelength = filter.pivotwavelength()

            # Insert the wavelength at the appropriate place
            for i in range(len(total_grid)):
                if total_grid[i] > wavelength:
                    total_grid.insert(i, wavelength)
                    break
            # If no break is encountered, no value in the grid was greater than our filter wavelength,
            # so add the filter wavelength at the end
            else: total_grid.append(wavelength)

        # Create table for the wavelength grid
        self.wavelength_grid = tables.new([total_grid], names=["Wavelength"])

    # -----------------------------------------------------------------

    def create_bulge_model(self):

        """
        :return:
        """

        # Create a Sersic model for the bulge
        self.bulge = SersicModel.from_galfit(self.parameters.bulge, self.parameters.inclination, self.parameters.disk.PA)

    # -----------------------------------------------------------------

    def create_deprojection_model(self):

        """
        This function ...
        :return:
        """

        filename = None
        hz = None

        distance = self.parameters.distance
        inclination = self.parameters.inclination
        #azimuth = 0.0
        pa = self.parameters.disk.PA
        pixels_x = self.reference_wcs.xsize
        pixels_y = self.reference_wcs.ysize
        pixel_center = self.parameters.center.to_pixel(self.reference_wcs)

        center = Position(0.5*pixels_x - pixel_center.x - 1, 0.5*pixels_y - pixel_center.y - 1)
        center_x = center.x * Unit("pix")
        center_y = center.y * Unit("pix")
        center_x = (center_x * self.reference_wcs.pixelscale.x.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())
        center_y = (center_y * self.reference_wcs.pixelscale.y.to("deg/pix") * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the pixelscale in physical units
        pixelscale_angular = self.reference_wcs.xy_average_pixelscale * Unit("pix") # in deg
        pixelscale = (pixelscale_angular * distance).to("pc", equivalencies=dimensionless_angles())

        # Get the x and y size
        x_size = self.reference_wcs.xsize
        y_size = self.reference_wcs.ysize

        # Create the deprojection model
        self.deprojection = DeprojectionModel(filename, pixelscale, pa, inclination, x_size, y_size, xc, yc, hz)

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add a new SEDInstrument
        # name, distance, inclination, azimuth, position_angle
        self.ski.add_sed_instrument("earth", self.parameters.distance, self.parameters.inclination, 0.0, self.parameters.disk.PA)

        # Set the number of photon packages
        self.ski.setpackages(self.config.packages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # Set the stellar and dust components
        self.set_components()

        # Set transient dust emissivity
        self.ski.set_transient_dust_emissivity()

        min_x = -15. * Unit("kpc")
        max_x = 15. * Unit("kpc")
        min_y = -15. * Unit("kpc")
        max_y = 15. * Unit("kpc")
        min_z = -3. * Unit("kpc")
        max_z = 3. * Unit("kpc")
        # Set the dust grid
        self.ski.set_binary_tree_dust_grid(min_x, max_x, min_y, max_y, min_z, max_z)

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Dust self-absorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

        # Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Set the evolved stellar bulge component
        self.set_bulge_component()

        # Set the evolved stellar disk component
        self.set_old_stellar_component()

        # Set the young stellar component
        self.set_young_stellar_component()

        # Set the ionizing stellar component
        self.set_ionizing_stellar_component()

        # The dust component
        self.set_dust_component()

    # -----------------------------------------------------------------

    def set_bulge_component(self):

        """
        This function ...
        :return:
        """

        bulge_template = "BruzualCharlot"
        bulge_age = 12
        bulge_metallicity = 0.02

        fluxdensity =

        # Set the parameters of the bulge
        self.ski.set_stellar_component_geometry("Evolved stellar bulge", self.bulge)
        self.ski.set_stellar_component_sed("Evolved stellar bulge", bulge_template, bulge_age, bulge_metallicity) # SED
        self.ski.set_stellar_component_luminosity("Evolved stellar bulge", luminosity, self.i1) # normalization

    # -----------------------------------------------------------------

    def set_old_stellar_component(self):

        """
        This function ...
        :return:
        """

        disk_template = "BruzualCharlot"
        disk_age = 8
        disk_metallicity = 0.02

        scale_height = 538 * Unit("pc")  # Andromeda

        # Set the parameters of the evolved stellar component
        self.deprojection.file_name = "old_stars.fits"
        self.deprojection.scale_height = scale_height
        self.ski.set_stellar_component_geometry("Evolved stellar disk", self.deprojection)
        self.ski.set_stellar_component_sed("Evolved stellar disk", disk_template, disk_age, disk_metallicity) # SED
        self.ski.set_stellar_component_luminosity("Evolved stellar disk", luminosity, self.i1) # normalization

    # -----------------------------------------------------------------

    def set_young_stellar_component(self):

        """
        This function ...
        :return:
        """

        young_template = "BruzualCharlot"
        young_age = 0.1
        young_metallicity = 0.02

        scale_height = 190 * Unit("pc") # Andromeda

        # Set the parameters of the young stellar component
        self.deprojection.file_name = "young_stars.fits"
        self.deprojection.scale_height = scale_height
        self.ski.set_stellar_component_geometry("Young stars", self.deprojection)
        self.ski.set_stellar_component_sed("Young stars", young_template, young_age, young_metallicity) # SED
        self.ski.set_stellar_component_luminosity("Young stars", luminosity, self.fuv) # normalization

    # -----------------------------------------------------------------

    def set_ionizing_stellar_component(self):

        """
        This function ...
        :return:
        """

        ionizing_metallicity = 0.02
        ionizing_compactness = 6
        ionizing_pressure = 1e12 * Unit("K/m3")
        ionizing_covering_factor = 0.2

        scale_height = 190 * Unit("pc") # Andromeda

        # Set the parameters of the ionizing stellar component
        self.deprojection.file_name = "ionizing_stars.fits"
        self.deprojection.scale_height = scale_height
        self.ski.set_stellar_component_geometry("Ionizing stars", self.deprojection)
        self.ski.set_stellar_component_mappingssed("Ionizing stars", ionizing_metallicity, ionizing_compactness, ionizing_pressure, ionizing_covering_factor) # SED
        self.ski.set_stellar_component_luminosity("Ionizing stars", luminosity, self.fuv) # normalization

    # -----------------------------------------------------------------

    def set_dust_component(self):

        """
        This function ...
        :return:
        """

        scale_height = 238 * Unit("pc") # Andromeda
        dust_mass = 4.3e7 * Unit("Msun") # Andromeda

        # Set the parameters of the dust component
        self.deprojection.file_name = "dust.fits"
        self.deprojection.scale_height = scale_height
        self.ski.set_dust_component_geometry(0, self.deprojection)
        self.ski.set_dust_component_themis_mix(0, hydrocarbon_pops, enstatite_pops, forsterite_pops) # dust mix
        self.ski.set_dust_component_mass(0, dust_mass) # dust mass

    # -----------------------------------------------------------------

    def place_input(self):

        """
        This function ...
        :return:
        """

        # -- The wavelength grid --

        # Determine the path to the wavelength grid file
        grid_path = filesystem.join(self.fit_in_path, "wavelengths.txt")

        # Write the wavelength grid
        tables.write(self.wavelength_grid, grid_path)

        # -- The old stars map --

        # Determine the path to the old stars map
        old_stars_path = filesystem.join(self.maps_path, "old_stars.fits")

        # Copy the map
        filesystem.copy_file(old_stars_path, self.fit_in_path)

        # -- The young stars map --

        # Determine the path to the young stars map
        young_stars_path = filesystem.join(self.maps_path, "young_stars.fits")

        # Copy the map
        filesystem.copy_file(young_stars_path, self.fit_in_path)

        # -- The ionizing stars map --

        # Determine the path to the ionizing stars map
        ionizing_stars_path = filesystem.join(self.maps_path, "ionizing_stars.fits")

        # Copy the map
        filesystem.copy_file(ionizing_stars_path, self.fit_in_path)

        # -- The dust map --

        # Determine the path to the dust map
        dust_path = filesystem.join(self.maps_path, "dust.fits")

        # Copy the map
        filesystem.copy_file(dust_path, self.fit_in_path)

    # -----------------------------------------------------------------

    def place_ski_file(self):

        """
        This function ...
        :return:
        """

        # Save the ski file to the specified location
        self.ski.saveto(self.fit_ski_path)

# -----------------------------------------------------------------
