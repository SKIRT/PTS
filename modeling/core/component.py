#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.core.component Contains the ModelingComponent class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from abc import ABCMeta

# Import astronomical modules
from astropy.utils import lazyproperty

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.tools import introspection
from ...core.tools import filesystem as fs
from ..core.sed import ObservedSED
from ...core.basics.filter import Filter
from ...magic.core.dataset import DataSet
from ...magic.core.frame import Frame
from ...magic.basics.skyregion import SkyRegion
from ...magic.basics.region import Region
from ..basics.models import load_3d_model, load_2d_model
from ..basics.projection import GalaxyProjection
from ..basics.properties import GalaxyProperties
from ...magic.tools import catalogs
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.mask import Mask
from ...core.tools.logging import log
from ...core.basics.configuration import Configuration

# -----------------------------------------------------------------

class ModelingComponent(Configurable):
    
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
        super(ModelingComponent, self).__init__(config)

        # Attributes
        self.galaxy_name = None

        # The modeling meta file
        self.meta_file_path = None

        # Modeling directories
        self.data_path = None
        self.prep_path = None
        self.truncation_path = None
        self.phot_path = None
        self.maps_path = None
        self.components_path = None
        self.fit_path = None
        self.analysis_path = None
        self.reports_path = None
        self.visualisation_path = None
        self.plot_path = None
        self.log_path = None
        self.config_path = None
        self.show_path = None

        # PTS directories
        self.kernels_path = None

        # Reference image
        self.reference_image = "Pacs red"

        # The path to the observed SEDs
        self.observed_sed_path = None
        self.observed_sed_dustpedia_path = None

        # The path to the fitting configuration file
        self.fitting_configuration_path = None

        # The path to the maps
        self.old_stellar_map_path = None
        self.young_stellar_map_path = None
        self.ionizing_stellar_map_path = None
        self.dust_map_path = None

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

    # -----------------------------------------------------------------

    def setup(self):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ModelingComponent, self).setup()

        # -- Attributes --

        # Get the name of the galaxy (the name of the base directory)
        self.galaxy_name = fs.name(self.config.path)

        # Determine the path to the modeling meta file
        self.meta_file_path = fs.join(self.config.path, "modeling.meta")

        # Check for the presence of the meta file
        if not fs.is_file(self.meta_file_path): raise ValueError("The current working directory is not a radiative transfer modeling directory (the meta file is missing)")

        # Get the full paths to the necessary subdirectories and CREATE THEM
        self.data_path = fs.create_directory_in(self.config.path, "data")
        self.prep_path = fs.create_directory_in(self.config.path, "prep")
        self.truncation_path = fs.create_directory_in(self.config.path, "truncated")
        self.phot_path = fs.create_directory_in(self.config.path, "phot")
        self.maps_path = fs.create_directory_in(self.config.path, "maps")
        self.components_path = fs.create_directory_in(self.config.path, "components")
        self.fit_path = fs.create_directory_in(self.config.path, "fit")
        self.analysis_path = fs.create_directory_in(self.config.path, "analysis")
        self.reports_path = fs.create_directory_in(self.config.path, "reports")
        self.visualisation_path = fs.create_directory_in(self.config.path, "visualisation")
        self.plot_path = fs.create_directory_in(self.config.path, "plot")
        self.log_path = fs.create_directory_in(self.config.path, "log")
        self.config_path = fs.create_directory_in(self.config.path, "config")
        self.show_path = fs.create_directory_in(self.config.path, "show")

        # Determine the path to the kernels user directory
        self.kernels_path = fs.join(introspection.pts_user_dir, "kernels")

        # Set the path to the observed SED
        self.observed_sed_path = fs.join(self.phot_path, "fluxes.dat")

        # Set the path to the DustPedia observed SED
        self.observed_sed_dustpedia_path = fs.join(self.data_path, "fluxes.dat")

        # Set the path to the fitting configuration file
        self.fitting_configuration_path = fs.join(self.fit_path, "configuration.cfg")

        # Set the paths to the input maps
        self.old_stellar_map_path = fs.join(self.maps_path, "old_stars.fits")
        self.young_stellar_map_path = fs.join(self.maps_path, "young_stars.fits")
        self.ionizing_stellar_map_path = fs.join(self.maps_path, "ionizing_stars.fits")
        self.dust_map_path = fs.join(self.maps_path, "dust.fits")

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

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_id(self):

        """
        This function ...
        :return:
        """

        # If the galaxy properties file has not been created yet
        if not fs.is_file(self.galaxy_properties_path):

            # Determine the NGC id of the galaxy
            ngc_id = catalogs.get_ngc_name(self.galaxy_name)

            # Return the ID
            return ngc_id

        # If the galaxy properties file has been created
        else: return self.galaxy_properties.ngc_id

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_id_nospaces(self):

        """
        This function ...
        :return:
        """

        return self.ngc_id.replace(" ", "")

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_sed(self):

        """
        This function ...
        :return:
        """

        return ObservedSED.from_file(self.observed_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filters(self):

        """
        This function ...
        :return:
        """

        return self.observed_sed.filters()

    # -----------------------------------------------------------------

    @lazyproperty
    def observed_filter_names(self):

        """
        This function ...
        :return:
        """

        return [str(fltr) for fltr in self.observed_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filters(self):

        """
        This function ...
        :return:
        """

        return map(Filter.from_string, self.fitting_filter_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_configuration(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.fitting_configuration_path) if fs.is_file(self.fitting_configuration_path) else None

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.filters if self.fitting_configuration is not None else None

    # -----------------------------------------------------------------

    @lazyproperty
    def free_parameter_labels(self):

        """
        This function ...
        :return:
        """

        return self.fitting_configuration.free_parameters if self.fitting_configuration is not None else None

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
        region = SkyRegion.from_file(self.disk_region_path)

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
        region = SkyRegion.from_file(path)
        ellipse = region[0]

        # Return the (sky) ellipse
        return ellipse

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_filter(self):

        """
        This function ...
        :return:
        """

        return Filter.from_string(self.reference_image)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):

        """
        This function ...
        :return:
        """

        #return self.dataset.get_wcs(self.reference_image)
        reference_path = fs.join(self.prep_path, self.reference_image, "result.fits")
        return CoordinateSystem.from_file(reference_path)

    # -----------------------------------------------------------------

    def sky_annulus_region(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        path = fs.join(self.prep_path, image_name, "sky", "annulus.reg")
        return Region.from_file(path)

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
