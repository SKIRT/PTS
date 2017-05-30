#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.component.galaxy Contains the GalaxyModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import math
import gc
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.data.sed import ObservedSED
from ...magic.core.dataset import DataSet
from ...magic.core.frame import Frame
from ...magic.region.list import SkyRegionList
from ..basics.models import load_3d_model, load_2d_model
from ..basics.properties import GalaxyProperties
from ...core.tools.logging import log
from ...magic.prepare.statistics import PreparationStatistics
from .component import ModelingComponent
from ...magic.region.ellipse import SkyEllipseRegion
from ...magic.basics.stretch import SkyStretch
from ...core.tools import types
from ...core.filter.filter import parse_filter

# -----------------------------------------------------------------

class GalaxyModelingComponent(ModelingComponent):
    
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
        super(GalaxyModelingComponent, self).__init__(*args, **kwargs)

        # Attributes
        self.galaxy_name = None

        # The path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = None

        # The path to the galaxy properties file
        self.galaxy_properties_path = None

        # The path to the components/models directory
        self.components_models_path = None

        # The path to the components/projections directory
        self.components_projections_path = None

        # The paths to the bulge and disk models
        self.bulge_model_path = None
        self.disk_model_path = None

        # The path to the components/images directory
        self.components_images_path = None

        # The paths to the final bulge, disk and model images
        self.bulge_image_path = None
        self.disk_image_path = None
        self.model_image_path = None

        # The path to the truncation/masks directory
        #self.truncation_masks_path = None

        # The path to the data/seds directory
        self.data_seds_path = None

        # The path to the DustPedia observed SED for the galaxy
        self.dustpedia_sed_path = None

        # The path to the data/images directory
        self.data_images_path = None

        # The path to the disk region file
        self.disk_region_path = None

        # The path to the initial data set file
        self.initial_dataset_path = None

        # The path to the prepared data set file
        self.prepared_dataset_path = None

        # The path to the preparation statistics file
        self.preparation_statistics_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyModelingComponent, self).setup(**kwargs)

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

        # Set the path to the galaxy properties file
        self.galaxy_properties_path = fs.join(self.data_path, "properties.dat")

        # Set the path to the components/models directory
        self.components_models_path = fs.create_directory_in(self.components_path, "models")

        # Set the path to the components/projections directory
        self.components_projections_path = fs.create_directory_in(self.components_path, "projections")

        # Set the paths to the bulge and disk models
        self.bulge_model_path = fs.join(self.components_models_path, "bulge.mod")
        self.disk_model_path = fs.join(self.components_models_path, "disk.mod")

        # Set the path to the components/images directory
        self.components_images_path = fs.create_directory_in(self.components_path, "images")

        # Set the path to the final bulge, disk and model images
        self.bulge_image_path = fs.join(self.components_images_path, "bulge.fits")
        self.disk_image_path = fs.join(self.components_images_path, "disk.fits")
        self.model_image_path = fs.join(self.components_images_path, "model.fits")

        # Set the path to the truncation/masks directory
        #self.truncation_masks_path = fs.create_directory_in(self.truncation_path, "masks")

        # Set ...
        self.data_seds_path = fs.create_directory_in(self.data_path, "SEDs")

        # The DustPedia SED path
        self.dustpedia_sed_path = fs.join(self.data_seds_path, "DustPedia.dat")

        # Set ...
        self.data_images_path = fs.create_directory_in(self.data_path, "images")

        # Set the path to the disk region file
        self.disk_region_path = fs.join(self.components_path, "disk.reg")

        # Set the path to the initial dataset file
        self.initial_dataset_path = fs.join(self.prep_path, "initial_dataset.dat")

        # Set the path to the prepared dataset file
        self.prepared_dataset_path = fs.join(self.prep_path, "dataset.dat")

        # Set the path to the preparation statistics file
        self.preparation_statistics_path = fs.join(self.prep_path, "statistics.dat")

    # -----------------------------------------------------------------

    @property
    def data_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.data_path

    # -----------------------------------------------------------------

    @property
    def prep_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.prep_path

    # -----------------------------------------------------------------

    @property
    def inspect_path(self):

        """
        This function ...
        :return: 
        """

        return self.environment.inspect_path

    # -----------------------------------------------------------------

    @property
    def truncation_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.truncation_path

    # -----------------------------------------------------------------

    @property
    def phot_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.phot_path

    # -----------------------------------------------------------------

    @property
    def maps_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.maps_path

    # -----------------------------------------------------------------

    @property
    def components_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.components_path

    # -----------------------------------------------------------------

    @property
    def deprojection_path(self):

        """
        This function ...
        :return:
        """

        return self.environment.deprojection_path

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_name(self):

        """
        This function ...
        :return:
        """

        # Get the NGC name from the configuration
        return self.modeling_configuration.ngc_name

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_name_nospaces(self):

        """
        This function ...
        :return:
        """

        return self.ngc_name.replace(" ", "")

    # -----------------------------------------------------------------

    @lazyproperty
    def hyperleda_name(self):

        """
        This function ...
        :return:
        """

        # Get the HYPERLEDA name from the configuration
        return self.modeling_configuration.hyperleda_name

    # -----------------------------------------------------------------

    @lazyproperty
    def initial_dataset(self):

        """
        This function ...
        :return:
        """

        return DataSet.from_file(self.initial_dataset_path, check=False) # don't check whether the file are actually present (caching on remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_statistics(self):

        """
        This function ...
        :return:
        """

        return PreparationStatistics.from_file(self.preparation_statistics_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def wcs_reference_filter(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return self.preparation_statistics.rebinning_filter

    # -----------------------------------------------------------------

    #@lazyproperty
    #def wcs_reference_image_name(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return self.dataset.get_name_for_filter(self.wcs_reference_filter)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def fwhm_reference_filter(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return self.preparation_statistics.convolution_filter

    # -----------------------------------------------------------------

    #@lazyproperty
    #def fwhm_reference_image_name(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return self.dataset.get_name_for_filter(self.fwhm_reference_filter)

    # -----------------------------------------------------------------

    @lazyproperty
    def dataset(self):

        """
        This funtion ...
        :return: 
        """

        return DataSet.from_file(self.prepared_dataset_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def dataset(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Initialize the dataset
        #dataset = DataSet()

        # Loop over all directories in the preparation directory
        #for path, name in fs.directories_in_path(self.prep_path, not_contains="Halpha", returns=["path", "name"]):

            # Check whether the 'result' file exists
            #result_path = fs.join(path, "result.fits")
            #if not fs.is_file(result_path): raise RuntimeError("The " + name + " result image does not exist")

            # Add the image path to the dataset
            #dataset.add_path(name, result_path)

            # Check whether a truncation mask is available
            ##mask_path = self.truncation_mask_path(name)

            # Add the mask path
            ##if mask_path is not None: dataset.add_mask_path(name, mask_path)

        # Return the dataset
        #return dataset

    # -----------------------------------------------------------------

    #@lazyproperty
    #def masked_halpha_frame(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Get the frame
        #frame = self.halpha_frame.copy()

        # Check whether the reference truncation mask is present
        #if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the image
        #mask = Mask.from_file(self.reference_mask_path)
        #frame[mask] = 0.0

        # Return the frame
        #return frame

    # -----------------------------------------------------------------

    #@lazyproperty
    #def halpha_frame(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Determine the path to the image
        #path = fs.join(self.prep_path, "Mosaic Halpha", "result.fits")

        # Check whether the Halpha image is present
        #if not fs.is_file(path): raise IOError("The prepared H-alpha image is missing")

        # Load and return the frame
        #frame = Frame.from_file(path)

        # Return the frame
        #return frame

    # -----------------------------------------------------------------

    #@lazyproperty
    #def halpha_errors(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Determine the path to the image
        #path = fs.join(self.prep_path, "Mosaic Halpha", "result.fits")

        # Check whether the Halpha image is present
        #if not fs.is_file(path): raise IOError("The prepared H-alpha image is missing")

        # Load the errors frame
        #errors = Frame.from_file(path, plane="errors")

        # Return the errors frame
        #return errors

    # -----------------------------------------------------------------

    #@lazyproperty
    #def halpha_relative_errors(self):

        #"""
        #This function ...
        #:return:
        #"""

        #frame = self.halpha_frame
        #errors = self.halpha_errors

        # Calculate the relative errors frame and return it
        #return errors / frame

    # -----------------------------------------------------------------

    #@lazyproperty
    #def halpha_significance(self):

        #"""
        #This function ...
        #:return:
        #"""

        #frame = self.halpha_frame
        #errors = self.halpha_errors

        # Calculate the significance map and return it
        #return frame / errors

    # -----------------------------------------------------------------

    #def get_halpha_significance_levels(self, levels, below_levels_value=float("nan")):

        #"""
        #This function ...
        #:param levels:
        #:param below_levels_value:
        #:return:
        #"""

        # Get the significance map
        #significance = self.halpha_significance

        # Create a frame full of nans
        #significance_levels = Frame.filled_like(significance, below_levels_value)

        # Loop over the levels
        #for level in levels: significance_levels[significance > level] = level

        # Return the significance levels map
        #return significance_levels

    # -----------------------------------------------------------------

    #def get_halpha_significance_mask(self, level):

        #"""
        #This function ...
        #:param level:
        #:return:
        #"""

        #return self.halpha_significance > level

    # -----------------------------------------------------------------

    #@lazyproperty
    #def masked_disk_frame(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Get the disk frame
        #frame = self.disk_frame.copy()

        # Check whether the reference truncation mask is present
        #if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the disk frame
        #mask = Mask.from_file(self.reference_mask_path)
        #frame[mask] = 0.0

        # Return the frame
        #return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_frame(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.disk_image_path): raise IOError("The disk image is not present. Run the 'decompose' step to create it")

        # Load the frame
        frame = Frame.from_file(self.disk_image_path)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    #@lazyproperty
    #def masked_bulge_frame(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Get the bulge frame
        #frame = self.bulge_frame.copy()

        # Check whether the reference truncation mask is present
        #if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the bulge frame
        #mask = Mask.from_file(self.reference_mask_path)
        #frame[mask] = 0.0

        # Return the frame
        #return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_frame(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.bulge_image_path): raise IOError("The bulge image is not present. Run the 'decompose' step to create it")

        # Load the frame
        frame = Frame.from_file(self.bulge_image_path)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_ellipse(self):

        """
        This function ...
        :return:
        """

        # Open the region
        region = SkyRegionList.from_file(self.disk_region_path)

        # Return the first and only shape
        return region[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.disk_ellipse.angle

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse(self):

        """
        This function ...
        :return:
        """

        # Load the ellipse
        path = fs.join(self.truncation_path, "ellipse.reg")
        region = SkyRegionList.from_file(path)
        ellipse = region[0]

        # Return the (sky) ellipse
        return ellipse

    # -----------------------------------------------------------------

    @lazyproperty
    def physical_truncation_ellipse(self):

        """
        This fucntion ...
        :return: 
        """

        raise NotImplementedError("Not implemented yet")
        #return self.truncation_area

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_area(self):

        """
        This function ...
        :return:
        """

        # Convert the semi minor and semi major axis lengths from angular to physical sizes
        semimajor = (self.truncation_ellipse.semimajor * self.galaxy_distance).to("kpc", equivalencies=dimensionless_angles())
        semiminor = (self.truncation_ellipse.semiminor * self.galaxy_distance).to("kpc", equivalencies=dimensionless_angles())

        # Calculate the area in kpc^2
        # A = pi * a * b
        area = math.pi * semimajor * semiminor

        # Return the area
        return area

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_box(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse.bounding_box

    # -----------------------------------------------------------------

    #@lazyproperty
    #def reference_wcs_path(self):

        #"""
        #This function ...
        #:return:
        #"""

        #reference_path = fs.join(self.prep_path, self.wcs_reference_image_name, "result.fits")
        #return reference_path

    # -----------------------------------------------------------------

    #@lazyproperty
    #def reference_wcs(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return CoordinateSystem.from_file(self.reference_wcs_path)

    # -----------------------------------------------------------------

    def highest_resolution_wcs(self):

        """
        This function ...
        :return:
        """

        return self.dataset.min_pixelscale_wcs

    # -----------------------------------------------------------------

    def lowest_resolution_wcs(self):

        """
        This function ...
        :return:
        """

        return self.dataset.max_pixelscale_wcs

    # -----------------------------------------------------------------

    def wcs_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.dataset.get_wcs_for_filter(fltr)

    # -----------------------------------------------------------------

    def sky_annulus_region(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        path = fs.join(self.prep_path, image_name, "sky", "annulus.reg")
        return SkyRegionList.from_file(path)

    # -----------------------------------------------------------------

    def sky_annulus_outer(self, image_name):

        """
        This function ...
        :return:
        """

        # Get the sky annulus region for this image
        region = self.sky_annulus_region(image_name)

        # Return the ellipse with the largest radius
        #return max(region, key=lambda ellipse: ellipse.radius)

        return region[0].base

    # -----------------------------------------------------------------

    def sky_annulus_inner(self, image_name):

        """
        This function ...
        :return:
        """

        # Get the sky annulus region for this image
        region = self.sky_annulus_region(image_name)

        # Return the ellipse with the smallest radius
        #return min(region, key=lambda ellipse: ellipse.radius)

        return region[0].exclude

    # -----------------------------------------------------------------

    #def truncation_mask_path(self, image_name):

        #"""
        #This function ...
        #:param image_name:
        #:return:
        #"""

        # Check whether mask is present with image name, or else use the reference mask file
        #path = fs.join(self.truncation_masks_path, image_name + ".fits")
        #if not fs.is_file(path): path = self.reference_mask_path

        # Return None if truncation has not been performed yet
        #if not fs.is_file(path): return None
        #else: return path

    # -----------------------------------------------------------------

    #def truncation_mask(self, image_name):

        #"""
        #This function ...
        #:param image_name:
        #:return:
        #"""

        # Get the path to the truncation mask
        #path = self.truncation_mask_path(image_name)

        # Return None if no mask is present
        #if path is None: return None

        # Else, return the mask
        #return Mask.from_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.galaxy_properties_path): raise IOError("The galaxy properties file is not present. Perform 'fetch_properties' to create this file'")

        # Load the properties
        properties = GalaxyProperties.from_file(self.galaxy_properties_path)

        # Return the property map
        return properties

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_ellipse(self): # from properties

        """
        This function ...
        :return: 
        """

        # Get properties
        center = self.galaxy_properties.center
        major = self.galaxy_properties.major_arcsec
        position_angle = self.galaxy_properties.position_angle
        ellipticity = self.galaxy_properties.ellipticity

        # 1 / axial_ratio = 1 - ellipticity
        axial_ratio = 1. / (1. - ellipticity)

        # Set radius
        minor = major * axial_ratio
        radius = SkyStretch(major, minor)

        # Create and return the region
        region = SkyEllipseRegion(center, radius, position_angle)
        return region

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_distance(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_center(self):

        """
        This function ...
        :return: 
        """

        return self.galaxy_properties.center

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_inclination(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.inclination

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_position_angle(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_redshift(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.redshift

    # -----------------------------------------------------------------

    #@lazyproperty
    #def earth_projection(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.earth_projection_path): raise IOError("The earth projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        #projection = GalaxyProjection.from_file(self.earth_projection_path)

        # Return the projection
        #return projection

    # -----------------------------------------------------------------

    #@lazyproperty
    #def edgeon_projection(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.edgeon_projection_path): raise IOError("The edgeon projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        #projection = GalaxyProjection.from_file(self.edgeon_projection_path)

        # Return the projection
        #return projection

    # -----------------------------------------------------------------

    #@lazyproperty
    #def faceon_projection(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.faceon_projection_path): raise IOError("The faceon projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        #projection = GalaxyProjection.from_file(self.faceon_projection_path)

        # Return the projection
        #return projection

    # -----------------------------------------------------------------

    #def create_instrument(self, instrument_type, projection):

        #"""
        #This function ...
        #:param instrument_type: "frame", "SED", "simple", or "full"
        #:param projection: "earth", "faceon", or "edgeon"
        #:return:
        #"""

        # Determine the instrument class
        #if instrument_type == "SED": instrument_class = SEDInstrument
        #elif instrument_type == "frame": instrument_class = FrameInstrument
        #elif instrument_type == "simple": instrument_class = SimpleInstrument
        #elif instrument_type == "full": instrument_class = FullInstrument
        #else: raise ValueError("Invalid instrument type: " + instrument_type)

        # Create the instrument and return it
        #if projection == "earth": return instrument_class.from_projection(self.earth_projection)
        #elif projection == "faceon": return instrument_class.from_projection(self.faceon_projection)
        #elif projection == "edgeon": return instrument_class.from_projection(self.edgeon_projection)
        #else: raise ValueError("Invalid projection: " + projection)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_model(self):

        """
        This function returns the bulge model
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.bulge_model_path): raise IOError("The bulge model file is not present. Run the 'decompose' step to create it")

        # Load the model
        return load_3d_model(self.bulge_model_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_model(self):

        """
        This function returns the disk model
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.disk_model_path): raise IOError("The disk model file is not present. Run the 'decompose' step to create it")

        # Load the model
        return load_3d_model(self.disk_model_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge2d_model(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.components_path, "2D", "S4G", "bulge.mod")
        return load_2d_model(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk2d_model(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.components_path, "2D", "S4G", "disk.mod")
        return load_2d_model(path)

    # -----------------------------------------------------------------

    #@property
    #def input_map_paths(self):

        #"""
        #This function returns the paths to the input maps of stellar and dust distribution
        #:return:
        #"""

        # Check whether the FITS files exist
        #if not fs.is_file(self.old_stellar_map_path): raise RuntimeError("The map of old stars is not present. Run make_old_map first.")
        #if not fs.is_file(self.young_stellar_map_path): raise RuntimeError("The map of young stars is not present. Run make_young_map first.")
        #if not fs.is_file(self.ionizing_stellar_map_path): raise RuntimeError("The map of ionizing stars is not present. Run make_ionizing_map first.")
        #if not fs.is_file(self.dust_map_path): raise RuntimeError("The dust map is not present. Run make_dust_map first.")

        # Return the paths to the maps of stars and dust
        #return [self.old_stellar_map_path, self.young_stellar_map_path, self.ionizing_stellar_map_path, self.dust_map_path]

    # -----------------------------------------------------------------

    #@property
    #def old_stellar_map_filename(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return fs.name(self.old_stellar_map_path)

    # -----------------------------------------------------------------

    #@property
    #def young_stellar_map_filename(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return fs.name(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    #@property
    #def ionizing_stellar_map_filename(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return fs.name(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    #@property
    #def dust_map_filename(self):

        #"""
        #This function ...
        #:return:
        #"""

        #return fs.name(self.dust_map_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def old_stars_map(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.old_stellar_map_path): raise IOError("The map of old stars is not present. Run 'make_old_map' to create it")

        # Open and return the frame of the old stellar distribution
        #return Frame.from_file(self.old_stellar_map_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def young_stars_map(self):

        #""""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.young_stellar_map_path): raise IOError("The map of young stars is not present. Run 'make_young_map' to create it")

        # Open and return the frame of the young stellar distribution
        #return Frame.from_file(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def ionizing_stars_map(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.ionizing_stellar_map_path): raise IOError("The map of ionizing stars is not present. Run 'make_ionizing_map' to create it")

        # Open and return the frame of the ionizing stellar distribution
        #return Frame.from_file(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    #@lazyproperty
    #def dust_map(self):

        #"""
        #This function ...
        #:return:
        #"""

        # Check whether the file is present
        #if not fs.is_file(self.dust_map_path): raise IOError("The dust map is not present. Run 'make_dust_map' to create it")

        # Open and return the frame of the dust distribution
        #return Frame.from_file(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the reference SEDs ...")

        # Loop over the SEDs in the data/SEDs directory
        for path, name in fs.files_in_path(self.data_seds_path, extension="dat", returns=["path", "name"], not_contains="Lines"):

            # Open the observed SED
            sed = ObservedSED.from_file(path)

            # Add the SED to the dictionary
            self.reference_seds[name] = sed

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_sed_labels(self):

        """
        This function ...
        :return:
        """

        return self.reference_seds.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def dustpedia_sed(self):

        """
        This function ...
        :return:
        """

        # Open the SED
        sed = ObservedSED.from_file(self.dustpedia_sed_path)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    @staticmethod
    def needs_poisson_errors(fltr):

        """
        This function ...
        :param fltr: 
        :return: 
        """

        if types.is_string_type(fltr): fltr = parse_filter(fltr)
        filter_string = str(fltr)
        return "GALEX" in filter_string or "SDSS" in filter_string

    # -----------------------------------------------------------------

    def get_data_image_and_error_paths(self):

        """
        This function ...
        :return: 
        """

        paths = dict()
        error_paths = dict()

        # Loop over the images
        for image_path, image_name in fs.files_in_path(self.data_images_path, extension="fits", not_contains="poisson",
                                                       returns=["path", "name"], recursive=True, recursion_level=1):

            # Determine directory path
            path = fs.directory_of(image_path)

            # Load the primary image frame
            frame = load_image_frame(image_path)

            # Determine name
            name = frame.filter_name

            # Add the image path
            paths[name] = image_path

            # Determine path to poisson error map
            poisson_path = fs.join(path, image_name + "_poisson.fits")

            # Set the path to the poisson error map
            if fs.is_file(poisson_path):

                # Debugging
                log.debug("Poisson error frame found for '" + name + "' image ...")
                error_paths[name] = poisson_path

            # Poisson frame not present
            elif self.needs_poisson_errors(frame.filter): raise RuntimeError("Poisson error frame not found for the " + name + " image. Run the appropriate command to create the mosaics and poisson frames.")

            # Free memory
            gc.collect()

        # Return the paths and error paths
        return paths, error_paths

    # -----------------------------------------------------------------

    def get_data_image_paths(self):

        """
        This function ...
        :return: 
        """

        paths = dict()

        # Loop over the images
        for image_path, image_name in fs.files_in_path(self.data_images_path, extension="fits", not_contains="poisson", returns=["path", "name"], recursive=True):

            # Load the primary image frame
            frame = load_image_frame(image_path)

            # Determine name
            name = frame.filter_name

            # Add the image path
            paths[name] = image_path

            # Free memory
            gc.collect()

        # Return the image paths
        return paths

# -----------------------------------------------------------------

def load_preparation_statistics(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Determine path
    prep_path = fs.join(modeling_path, "prep")
    path = fs.join(prep_path, "statistics.dat")

    # Load and return the statistics
    if fs.is_file(path): return PreparationStatistics.from_file(path)
    else: return None

# -----------------------------------------------------------------

def get_observed_sed_file_path(modeling_path):

    """
    This function ...
    :return:
    """

    return fs.join(modeling_path, "phot", "fluxes.dat")

# -----------------------------------------------------------------

def get_observed_sed(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return ObservedSED.from_file(get_observed_sed_file_path(modeling_path))

# -----------------------------------------------------------------

def get_galaxy_properties_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    data_path = fs.join(modeling_path, "data")
    return fs.join(data_path, "properties.dat")

# -----------------------------------------------------------------

def get_data_seds_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    data_path = fs.join(modeling_path, "data")
    return fs.join(data_path, "SEDs")

# -----------------------------------------------------------------

def get_data_images_path(modeling_path):

    """
    This function ....
    :param modeling_path:
    :return:
    """

    data_path = fs.join(modeling_path, "data")
    return fs.join(data_path, "images")

# -----------------------------------------------------------------

def get_reference_seds(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    seds = dict()

    # Loop over the SEDs in the data/SEDs directory
    for path, name in fs.files_in_path(get_data_seds_path(modeling_path), extension="dat", returns=["path", "name"], not_contains="Lines"):

        # Open the observed SED
        sed = ObservedSED.from_file(path)

        # Add the SED to the dictionary
        seds[name] = sed

    # Return the SEDs
    return seds

# -----------------------------------------------------------------

def get_dustpedia_sed_path(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    return fs.join(get_data_seds_path(modeling_path), "DustPedia.dat")

# -----------------------------------------------------------------

def get_dustpedia_sed(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    # Open the SED and return it
    return ObservedSED.from_file(get_dustpedia_sed_path(modeling_path))

# -----------------------------------------------------------------

def load_image_frame(path):

    """
    This function ...
    :param path: 
    :return: 
    """

    name = fs.strip_extension(fs.name(path))

    frame = None
    # Try opening
    try:
        # Open the image frame
        frame = Frame.from_file(path)
    except IOError:
        log.warning("The file '" + path + "' is probably damaged. Removing the file and exitting. Run the command again.")
        fs.remove_file(path)
        exit()

    # Determine the preparation name
    #if frame.filter is not None: prep_name = str(frame.filter)
    #else: prep_name = image_name
    if frame.filter is None:

        #log.warning("Did not recognize the filter of the '" + image_name + "' image: skipping")
        #continue
        raise RuntimeError("Did not recognize the filter for the '" + name + "' image")

    # Return the frame
    return frame

# -----------------------------------------------------------------
