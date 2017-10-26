#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.io.fits import Header
from ipyvolume.embed import embed_html

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.tools import introspection
from pts.modeling.basics.models import load_3d_model
from pts.modeling.basics.properties import GalaxyProperties
from pts.magic.core.frame import Frame
from pts.magic.basics.coordinatesystem import CoordinateSystem
from pts.modeling.basics.models import DeprojectionModel3D
from pts.modeling.basics.models import load_2d_model
from pts.core.units.quantity import PhotometricQuantity
from pts.core.units.parsing import parse_quantity
from pts.core.units.parsing import parse_unit as u
from pts.modeling.config.parameters import default_units
from pts.modeling.plotting.model import plot_galaxy_components, generate_html
from pts.core.tools.utils import lazyproperty

# -----------------------------------------------------------------

# Determine the path to the dropbox path and the path of the directory with the data for M81
m81_data_path = fs.join(introspection.get_dropbox_tests_pts_path_for_subproject("modeling"), "M81")

# -----------------------------------------------------------------

models_path = fs.join(m81_data_path, "models")
disk2d_path = fs.join(models_path, "disk.mod")
bulge2d_path = fs.join(models_path, "bulge.mod")

# -----------------------------------------------------------------

# Determine paths to the components
path = fs.join(m81_data_path, "components")
bulge_path = fs.join(path, "bulge.mod")
disk_path = fs.join(path, "disk.mod")

# -----------------------------------------------------------------

instrument_name = "earth"

# -----------------------------------------------------------------

titles = dict()
titles["bulge"] = "Evolved stellar bulge"
titles["old"] = "Evolved stellar disk"
titles["young"] = "Young stars"
titles["ionizing"] = "Ionizing stars"
titles["dust"] = "Dust disk"

# -----------------------------------------------------------------

# Bulge
bulge_template = "BruzualCharlot"
bulge_age = 10
bulge_metallicity = 0.03

# -----------------------------------------------------------------

# Old stellar disk
disk_template = "BruzualCharlot"
disk_age = 8
# disk_metallicity = 0.02
disk_metallicity = 0.03

# -----------------------------------------------------------------

# Young stellar disk
young_template = "BruzualCharlot"
young_age = 0.1
# young_metallicity = 0.02
young_metallicity = 0.03

# -----------------------------------------------------------------

# Ionizing stellar disk
ionizing_metallicity = 0.03  # XU KONG et al. 2000
ionizing_compactness = 6.
ionizing_pressure = 1e12 * u("K/m3")
ionizing_covering_factor = 0.2

# Convert the SFR into a FUV luminosity
sfr = 0.8  # The star formation rate # see Perez-Gonzalez 2006 (mentions Devereux et al 1995)

# -----------------------------------------------------------------

# fuv_young:6.0068695608165e+36 W/micron
# fuv_ionizing:2.4590756925069244e+33 W/micron]
# 15450820.890962543 Msun

fuv_young = PhotometricQuantity(1e36, "W/micron")
fuv_ionizing = PhotometricQuantity(2.5e33, "W/micron")
dust_mass = parse_quantity("1.5e7 Msun")

#dust_mass = 1.5e7 * u("Msun")

hydrocarbon_pops = 25
silicate_pops = 25

# -----------------------------------------------------------------

old_filename = "old_stars.fits"
young_filename = "young_stars.fits"
ionizing_filename = "ionizing_stars.fits"
dust_filename = "dust.fits"

# -----------------------------------------------------------------

# Define the possible free parameters
possible_free_parameters = ["dust_mass", "fuv_young", "fuv_ionizing", "distance"]
default_free_parameters = ["dust_mass", "fuv_young", "fuv_ionizing"]

# Define the free parameter types
free_parameter_types = dict()

# Parameter types
free_parameter_types["dust_mass"] = "mass"
free_parameter_types["fuv_young"] = "spectral luminosity density"
free_parameter_types["fuv_ionizing"] = "spectral luminosity density"
free_parameter_types["distance"] = "length"

# Define free parameter units
free_parameter_units = dict()
for label in free_parameter_types:
    parameter_type = free_parameter_types[label]
    free_parameter_units[label] = default_units[parameter_type]

# Define the number of digits
parameter_ndigits = dict()
for label in free_parameter_types:
    parameter_ndigits[label] = 3

# Absolute ski file parameters
free_parameters_absolute_paths = dict()

# Stellar component parameters
free_parameters_relative_stellar_component_paths = dict()
free_parameters_relative_stellar_component_paths["fuv_young"] = ("normalization/SpectralLuminosityStellarCompNormalization/luminosity", titles["young"])
free_parameters_relative_stellar_component_paths["fuv_ionizing"] = ("normalization/SpectralLuminosityStellarCompNormalization/luminosity", titles["ionizing"])

# Dust component parameters
free_parameters_relative_dust_component_paths = dict()
free_parameters_relative_dust_component_paths["dust_mass"] = ("normalization/DustMassDustCompNormalization/dustMass", titles["dust"])

# Instrument parameters
free_parameters_relative_instruments_paths = dict()
free_parameters_relative_instruments_paths["distance"] = ("distance", instrument_name)

# Free parameter descriptions
free_parameter_descriptions = dict()
free_parameter_descriptions["dust_mass"] = "total dust mass"
free_parameter_descriptions["fuv_young"] = "FUV spectral luminosity of the young stellar component"
free_parameter_descriptions["fuv_ionizing"] = "FUV spectral luminosity of the ionizing stellar component"
free_parameter_descriptions["distance"] = "galaxy distance"

# -----------------------------------------------------------------

seds_path = fs.join(m81_data_path, "seds")
dustpedia_sed_path = fs.join(seds_path, "DustPedia.dat")

# -----------------------------------------------------------------

reference_wavelength_grid_filename = "wavelengths.txt"

# -----------------------------------------------------------------

# Determine the path
properties_path = fs.join(m81_data_path, "properties.dat")

# -----------------------------------------------------------------

# Determine path to maps directory
maps_path = fs.join(m81_data_path, "maps")
old_map_path = fs.join(maps_path, old_filename)
young_map_path = fs.join(maps_path, young_filename)
ionizing_map_path = fs.join(maps_path, ionizing_filename)
dust_map_path = fs.join(maps_path, dust_filename)

# -----------------------------------------------------------------

header_path = fs.join(maps_path, "header.txt")

# -----------------------------------------------------------------

class M81TestData(object):

    """
    This class ...
    """

    def __init__(self):

        """
        This function ...
        """

    # -----------------------------------------------------------------

    @lazyproperty
    def properties(self):

        """
        This function ...
        :return:
        """

        # Load
        return GalaxyProperties.from_file(properties_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk2d_model(self):

        """
        This function ....
        :return:
        """

        return load_2d_model(disk2d_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge2d_model(self):

        """
        This function ...
        :return:
        """

        return load_2d_model(bulge2d_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_fluxdensity(self):

        """
        This function ...
        :return:
        """

        # Get the flux density of the bulge
        return self.bulge2d_model.fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def old_scale_height(self):

        """
        This function ...
        :return:
        """

        # Get the scale height
        return self.disk2d_model.scalelength / 8.26  # De Geyter et al. 2014

    # -----------------------------------------------------------------

    @lazyproperty
    def old_fluxdensity(self):

        """
        This function ...
        :return:
        """

        # Get the 3.6 micron flux density with the bulge subtracted
        total_i1_fluxdensity = PhotometricQuantity(10.6552814592, "Jy")
        return total_i1_fluxdensity - self.bulge_fluxdensity

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_scale_height(self):

        """
        This function ...
        :return:
        """

        # scale_height = 260.5 * Unit("pc") # first models
        # dust_scale_height = 200. * u("pc")  # M51
        return 0.25 * self.old_scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def young_scale_height(self):

        """
        This function ...
        :return:
        """

        # Get the scale height
        # scale_height = 150 * Unit("pc") # first models
        # young_scale_height = 100. * u("pc")  # M51
        return 0.5 * self.old_scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_scale_height(self):

        """
        This function ...
        :return:
        """

        # Get the scale height
        # scale_height = 150 * Unit("pc") # first models
        # ionizing_scale_height = 100. * u("pc")  # M51
        return 0.25 * self.old_scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge(self):

        """
        This function ...
        :return:
        """

        bulge = load_3d_model(bulge_path)
        # No y flattening: this is a mistake in the file
        bulge.y_flattening = 1.
        return bulge

    # -----------------------------------------------------------------

    @lazyproperty
    def disk(self):

        """
        This function ...
        :return:
        """

        return load_3d_model(disk_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def header(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the header file
        return Header.fromtextfile(header_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def wcs(self):

        """
        Thisfunction ...
        :return:
        """

        return CoordinateSystem(header=self.header)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_map(self):

        """
        This function ...
        :return:
        """

        # Old stars

        old_map = Frame.from_file(old_map_path)
        old_map.wcs = self.wcs
        return old_map

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_map(self):

        """
        This function ...
        :return:
        """

        # young stars

        young_map = Frame.from_file(young_map_path)
        young_map.wcs = self.wcs
        return young_map

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_map(self):

        """
        This function ....
        :return:
        """

        # Ionizing stars

        ionizing_map = Frame.from_file(ionizing_map_path)
        ionizing_map.wcs = self.wcs
        return ionizing_map

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):

        """
        This function ...
        :return:
        """

        # Dust

        dust_map = Frame.from_file(dust_map_path)
        dust_map.wcs = self.wcs
        return dust_map

    # -----------------------------------------------------------------

    @lazyproperty
    def old_deprojection(self):

        """
        This function ...
        :return:
        """

        return DeprojectionModel3D.from_wcs(self.wcs, self.properties.center, self.properties.distance,
                                                        self.properties.position_angle, self.properties.inclination,
                                                        old_map_path,
                                                        self.old_scale_height)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_deprojection(self):

        """
        This function
        :return:
        """

        return DeprojectionModel3D.from_wcs(self.wcs, self.properties.center, self.properties.distance,
                                                          self.properties.position_angle, self.properties.inclination,
                                                          young_map_path, self.young_scale_height)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_deprojection(self):

        """
        This function ...
        :return:
        """

        return DeprojectionModel3D.from_wcs(self.wcs, self.properties.center, self.properties.distance,
                                                             self.properties.position_angle, self.properties.inclination,
                                                             ionizing_map_path, self.ionizing_scale_height)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_deprojection(self):

        """
        This function ...
        :return:
        """

        return DeprojectionModel3D.from_wcs(self.wcs, self.properties.center, self.properties.distance,
                                                         self.properties.position_angle, self.properties.inclination,
                                                         dust_map_path, self.dust_scale_height)

    # -----------------------------------------------------------------

    @lazyproperty
    def components(self):

        """
        This function ...
        :return:
        """

        # RETURN
        components = {"disk": self.disk, "bulge": self.bulge, "old": self.old_deprojection, "ionizing": self.ionizing_deprojection,
                      "young": self.young_deprojection, "dust": self.dust_deprojection}
        return components

    # -----------------------------------------------------------------

    def show_components(self):

        """
        This function ...
        :return:
        """

        old_components = {"disk": self.components["disk"], "bulge": self.components["bulge"]}

        # Determine the filepath for rendering
        filename = "render.html"
        filepath = fs.join(introspection.pts_temp_dir, filename)

        # Plot, create HTML
        box = plot_galaxy_components(old_components, draw=True, show=False)
        embed_html(filepath, box)

        # Open HTML
        fs.open_file(filepath)

    # -----------------------------------------------------------------

    def render_components_html(self, components, output_path, plot_kwargs=None, render_kwargs=None):

        """
        This function ...
        :param components:
        :param output_path:
        :param plot_kwargs:
        :param render_kwargs:
        :return:
        """

        if plot_kwargs is None: plot_kwargs = {}
        if render_kwargs is None: render_kwargs = {}

        title = "Bulge and disk"
        box = plot_galaxy_components(components, draw=True, show=False, **plot_kwargs)
        return generate_html(box, title, output_path, **render_kwargs)

# -----------------------------------------------------------------
