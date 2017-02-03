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
from ..basics.projection import GalaxyProjection
from ..basics.properties import GalaxyProperties
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.mask import Mask
from ...core.tools.logging import log
from ..basics.instruments import SEDInstrument, FrameInstrument, SimpleInstrument, FullInstrument
from ...magic.prepare.batch import PreparationStatistics
from .component import ModelingComponent

# -----------------------------------------------------------------

class GalaxyModelingComponent(ModelingComponent):
    
    """
    This class...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, config=None):

        """
        The constructor ...
        :param config:
        :return:
        """

        # Call the constructor of the base class
        super(GalaxyModelingComponent, self).__init__(config)

        # Attributes
        self.galaxy_name = None

        # Modeling directories
        self.data_path = None
        self.prep_path = None
        self.truncation_path = None
        self.phot_path = None
        self.maps_path = None
        self.components_path = None
        self.deprojection_path = None

        # The path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = None

        # The path to the maps
        self.old_stellar_map_path = None
        self.young_stellar_map_path = None
        self.ionizing_stellar_map_path = None
        self.dust_map_path = None

        # The paths to the significance maps
        self.old_stellar_significance_path = None
        self.young_stellar_significance_path = None
        self.ionizing_stellar_significance_path = None
        self.dust_significance_path = None

        # The paths to the cutoff masks
        self.old_stellar_cutoff_path = None
        self.young_stellar_cutoff_path = None
        self.ionizing_stellar_cutoff_path = None
        self.dust_cutoff_path = None

        # The path to the galaxy properties file
        self.galaxy_properties_path = None

        # The path to the components/models directory
        self.components_models_path = None

        # The path to the components/projections directory
        self.components_projections_path = None

        # The paths to the bulge and disk models
        self.bulge_model_path = None
        self.disk_model_path = None

        # The paths to the different galaxy projections
        self.earth_projection_path = None
        self.edgeon_projection_path = None
        self.faceon_projection_path = None

        # The path to the components/images directory
        self.components_images_path = None

        # The paths to the final bulge, disk and model images
        self.bulge_image_path = None
        self.disk_image_path = None
        self.model_image_path = None

        # The path to the truncation/masks directory
        self.truncation_masks_path = None

        # The path to the truncation mask of the reference image (and rebinned images in the dataset)
        self.reference_mask_path = None

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
        super(GalaxyModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.data_path = fs.create_directory_in(self.config.path, "data")
        self.prep_path = fs.create_directory_in(self.config.path, "prep")
        self.truncation_path = fs.create_directory_in(self.config.path, "truncated")
        self.phot_path = fs.create_directory_in(self.config.path, "phot")
        self.maps_path = fs.create_directory_in(self.config.path, "maps")
        self.components_path = fs.create_directory_in(self.config.path, "components")
        self.deprojection_path = fs.create_directory_in(self.config.path, "deprojection")

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

        # Set the paths to the input maps
        self.old_stellar_map_path = fs.join(self.maps_path, "old_stars.fits")
        self.young_stellar_map_path = fs.join(self.maps_path, "young_stars.fits")
        self.ionizing_stellar_map_path = fs.join(self.maps_path, "ionizing_stars.fits")
        self.dust_map_path = fs.join(self.maps_path, "dust.fits")

        # The paths to the significance masks
        self.old_stellar_significance_path = fs.join(self.maps_path, "old_stars_significance.fits")
        self.young_stellar_significance_path = fs.join(self.maps_path, "young_stars_significance.fits")
        self.ionizing_stellar_significance_path = fs.join(self.maps_path, "ionizing_stars_significance.fits")
        self.dust_significance_path = fs.join(self.maps_path, "dust_significance.fits")

        # The paths to the significance masks
        self.old_stellar_cutoff_path = fs.join(self.maps_path, "old_stars_cutoff.fits")
        self.young_stellar_cutoff_path = fs.join(self.maps_path, "young_stars_cutoff.fits")
        self.ionizing_stellar_cutoff_path = fs.join(self.maps_path, "ionizing_stars_cutoff.fits")
        self.dust_cutoff_path = fs.join(self.maps_path, "dust_cutoff.fits")

        # Set the path to the galaxy properties file
        self.galaxy_properties_path = fs.join(self.data_path, "properties.dat")

        # Set the path to the components/models directory
        self.components_models_path = fs.create_directory_in(self.components_path, "models")

        # Set the path to the components/projections directory
        self.components_projections_path = fs.create_directory_in(self.components_path, "projections")

        # Set the paths to the bulge and disk models
        self.bulge_model_path = fs.join(self.components_models_path, "bulge.mod")
        self.disk_model_path = fs.join(self.components_models_path, "disk.mod")

        # The paths to the different galaxy projections
        self.earth_projection_path = fs.join(self.components_projections_path, "earth.proj")
        self.edgeon_projection_path = fs.join(self.components_projections_path, "edgeon.proj")
        self.faceon_projection_path = fs.join(self.components_projections_path, "faceon.proj")

        # Set the path to the components/images directory
        self.components_images_path = fs.create_directory_in(self.components_path, "images")

        # Set the path to the final bulge, disk and model images
        self.bulge_image_path = fs.join(self.components_images_path, "bulge.fits")
        self.disk_image_path = fs.join(self.components_images_path, "disk.fits")
        self.model_image_path = fs.join(self.components_images_path, "model.fits")

        # Set the path to the truncation/masks directory
        self.truncation_masks_path = fs.create_directory_in(self.truncation_path, "masks")

        # The path to the truncation mask of the reference image (and rebinned images in the dataset)
        self.reference_mask_path = fs.join(self.truncation_masks_path, "reference.fits")

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

        return DataSet.from_file(self.initial_dataset_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def preparation_statistics(self):

        """
        This function ...
        :return:
        """

        return PreparationStatistics.from_file(self.preparation_statistics_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def wcs_reference_filter(self):

        """
        This function ...
        :return:
        """

        return self.preparation_statistics.rebinning_filter

    # -----------------------------------------------------------------

    @lazyproperty
    def wcs_reference_image_name(self):

        """
        This function ...
        :return:
        """

        return self.dataset.get_name_for_filter(self.wcs_reference_filter)

    # -----------------------------------------------------------------

    @lazyproperty
    def fwhm_reference_filter(self):

        """
        This function ...
        :return:
        """

        return self.preparation_statistics.convolution_filter

    # -----------------------------------------------------------------

    @lazyproperty
    def fwhm_reference_image_name(self):

        """
        This function ...
        :return:
        """


        return self.dataset.get_name_for_filter(self.fwhm_reference_filter)

    # -----------------------------------------------------------------

    @lazyproperty
    def dataset(self):

        """
        This function ...
        :return:
        """

        # Initialize the dataset
        dataset = DataSet()

        # Loop over all directories in the preparation directory
        for path, name in fs.directories_in_path(self.prep_path, not_contains="Halpha", returns=["path", "name"]):

            # Check whether the 'result' file exists
            result_path = fs.join(path, "result.fits")
            if not fs.is_file(result_path): raise RuntimeError("The " + name + " result image does not exist")

            # Add the image path to the dataset
            dataset.add_path(name, result_path)

            # Check whether a truncation mask is available
            mask_path = self.truncation_mask_path(name)

            # Add the mask path
            if mask_path is not None: dataset.add_mask_path(name, mask_path)

        # Return the dataset
        return dataset

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_halpha_frame(self):

        """
        This function ...
        :return:
        """

        # Get the frame
        frame = self.halpha_frame.copy()

        # Check whether the reference truncation mask is present
        if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the image
        mask = Mask.from_file(self.reference_mask_path)
        frame[mask] = 0.0

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def halpha_frame(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the image
        path = fs.join(self.prep_path, "Mosaic Halpha", "result.fits")

        # Check whether the Halpha image is present
        if not fs.is_file(path): raise IOError("The prepared H-alpha image is missing")

        # Load and return the frame
        frame = Frame.from_file(path)

        # Return the frame
        return frame

    # -----------------------------------------------------------------

    @lazyproperty
    def halpha_errors(self):

        """
        This function ...
        :return:
        """

        # Determine the path to the image
        path = fs.join(self.prep_path, "Mosaic Halpha", "result.fits")

        # Check whether the Halpha image is present
        if not fs.is_file(path): raise IOError("The prepared H-alpha image is missing")

        # Load the errors frame
        errors = Frame.from_file(path, plane="errors")

        # Return the errors frame
        return errors

    # -----------------------------------------------------------------

    @lazyproperty
    def halpha_relative_errors(self):

        """
        This function ...
        :return:
        """

        frame = self.halpha_frame
        errors = self.halpha_errors

        # Calculate the relative errors frame and return it
        return errors / frame

    # -----------------------------------------------------------------

    @lazyproperty
    def halpha_significance(self):

        """
        This function ...
        :return:
        """

        frame = self.halpha_frame
        errors = self.halpha_errors

        # Calculate the significance map and return it
        return frame / errors

    # -----------------------------------------------------------------

    def get_halpha_significance_levels(self, levels, below_levels_value=float("nan")):

        """
        This function ...
        :param levels:
        :param below_levels_value:
        :return:
        """

        # Get the significance map
        significance = self.halpha_significance

        # Create a frame full of nans
        significance_levels = Frame.filled_like(significance, below_levels_value)

        # Loop over the levels
        for level in levels: significance_levels[significance > level] = level

        # Return the significance levels map
        return significance_levels

    # -----------------------------------------------------------------

    def get_halpha_significance_mask(self, level):

        """
        This function ...
        :param level:
        :return:
        """

        return self.halpha_significance > level

    # -----------------------------------------------------------------

    @lazyproperty
    def masked_disk_frame(self):

        """
        This function ...
        :return:
        """

        # Get the disk frame
        frame = self.disk_frame.copy()

        # Check whether the reference truncation mask is present
        if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the disk frame
        mask = Mask.from_file(self.reference_mask_path)
        frame[mask] = 0.0

        # Return the frame
        return frame

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

    @lazyproperty
    def masked_bulge_frame(self):

        """
        This function ...
        :return:
        """

        # Get the bulge frame
        frame = self.bulge_frame.copy()

        # Check whether the reference truncation mask is present
        if not fs.is_file(self.reference_mask_path): raise IOError("The truncation mask has not been created")

        # Mask the bulge frame
        mask = Mask.from_file(self.reference_mask_path)
        frame[mask] = 0.0

        # Return the frame
        return frame

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

    @lazyproperty
    def reference_wcs_path(self):

        """
        This function ...
        :return:
        """

        reference_path = fs.join(self.prep_path, self.wcs_reference_image_name, "result.fits")
        return reference_path

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):

        """
        This function ...
        :return:
        """

        #return self.dataset.get_wcs(self.reference_image)
        return CoordinateSystem.from_file(self.reference_wcs_path)

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

    def truncation_mask_path(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        # Check whether mask is present with image name, or else use the reference mask file
        path = fs.join(self.truncation_masks_path, image_name + ".fits")
        if not fs.is_file(path): path = self.reference_mask_path

        # Return None if truncation has not been performed yet
        if not fs.is_file(path): return None
        else: return path

    # -----------------------------------------------------------------

    def truncation_mask(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        # Get the path to the truncation mask
        path = self.truncation_mask_path(image_name)

        # Return None if no mask is present
        if path is None: return None

        # Else, return the mask
        return Mask.from_file(path)

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
    def galaxy_distance(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_redshift(self):

        """
        This function ...
        :return:
        """

        return self.galaxy_properties.redshift

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.earth_projection_path): raise IOError("The earth projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        projection = GalaxyProjection.from_file(self.earth_projection_path)

        # Return the projection
        return projection

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.edgeon_projection_path): raise IOError("The edgeon projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        projection = GalaxyProjection.from_file(self.edgeon_projection_path)

        # Return the projection
        return projection

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.faceon_projection_path): raise IOError("The faceon projection file is not present. Run the 'decompose' step to create it")

        # Load the projection
        projection = GalaxyProjection.from_file(self.faceon_projection_path)

        # Return the projection
        return projection

    # -----------------------------------------------------------------

    def create_instrument(self, instrument_type, projection):

        """
        This function ...
        :param instrument_type: "frame", "SED", "simple", or "full"
        :param projection: "earth", "faceon", or "edgeon"
        :return:
        """

        # Determine the instrument class
        if instrument_type == "SED": instrument_class = SEDInstrument
        elif instrument_type == "frame": instrument_class = FrameInstrument
        elif instrument_type == "simple": instrument_class = SimpleInstrument
        elif instrument_type == "full": instrument_class = FullInstrument
        else: raise ValueError("Invalid instrument type: " + instrument_type)

        # Create the instrument and return it
        if projection == "earth": return instrument_class.from_projection(self.earth_projection)
        elif projection == "faceon": return instrument_class.from_projection(self.faceon_projection)
        elif projection == "edgeon": return instrument_class.from_projection(self.edgeon_projection)
        else: raise ValueError("Invalid projection: " + projection)

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

    @property
    def input_map_paths(self):

        """
        This function returns the paths to the input maps of stellar and dust distribution
        :return:
        """

        # Check whether the FITS files exist
        if not fs.is_file(self.old_stellar_map_path): raise RuntimeError("The map of old stars is not present. Run make_old_map first.")
        if not fs.is_file(self.young_stellar_map_path): raise RuntimeError("The map of young stars is not present. Run make_young_map first.")
        if not fs.is_file(self.ionizing_stellar_map_path): raise RuntimeError("The map of ionizing stars is not present. Run make_ionizing_map first.")
        if not fs.is_file(self.dust_map_path): raise RuntimeError("The dust map is not present. Run make_dust_map first.")

        # Return the paths to the maps of stars and dust
        return [self.old_stellar_map_path, self.young_stellar_map_path, self.ionizing_stellar_map_path, self.dust_map_path]

    # -----------------------------------------------------------------

    @property
    def old_stellar_map_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.old_stellar_map_path)

    # -----------------------------------------------------------------

    @property
    def young_stellar_map_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_stellar_map_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    @property
    def dust_map_filename(self):

        """
        This function ...
        :return:
        """

        return fs.name(self.dust_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_stars_map(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.old_stellar_map_path): raise IOError("The map of old stars is not present. Run 'make_old_map' to create it")

        # Open and return the frame of the old stellar distribution
        return Frame.from_file(self.old_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_stars_map(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.young_stellar_map_path): raise IOError("The map of young stars is not present. Run 'make_young_map' to create it")

        # Open and return the frame of the young stellar distribution
        return Frame.from_file(self.young_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_stars_map(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.ionizing_stellar_map_path): raise IOError("The map of ionizing stars is not present. Run 'make_ionizing_map' to create it")

        # Open and return the frame of the ionizing stellar distribution
        return Frame.from_file(self.ionizing_stellar_map_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_map(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.dust_map_path): raise IOError("The dust map is not present. Run 'make_dust_map' to create it")

        # Open and return the frame of the dust distribution
        return Frame.from_file(self.dust_map_path)

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
