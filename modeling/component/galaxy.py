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
import gc
from abc import ABCMeta

# Import astronomical modules
from astropy.io.fits import getheader

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.data.sed import ObservedSED
from ...magic.core.dataset import DataSet
from ...magic.core.frame import Frame
from ...magic.region.list import SkyRegionList
from ..basics.models import load_3d_model, load_2d_model
from ..basics.properties import GalaxyProperties
from ...core.basics.log import log
from ...magic.prepare.statistics import PreparationStatistics
from .component import ModelingComponent
from ...core.tools import types
from ...core.filter.filter import parse_filter
from ...core.remote.remote import Remote
from ..core.steps import cached_directory_path_for_single_command
from ..core.environment import GalaxyModelingEnvironment
from ...magic.core.remote import get_filter_name
from ...magic.tools import headers
from ...core.tools.utils import lazyproperty
from ..basics.projection import get_npixels, get_field, get_center
from ...magic.core.mask import Mask
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

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
        self.bulge2d_image_path = None
        self.bulge_image_path = None
        self.disk_image_path = None
        self.model_image_path = None

        # The path to the DustPedia observed SED for the galaxy
        self.dustpedia_sed_path = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(GalaxyModelingComponent, self).setup(**kwargs)

        # -- Attributes --

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
        self.bulge2d_image_path = fs.join(self.components_images_path, "bulge2d.fits")
        self.bulge_image_path = fs.join(self.components_images_path, "bulge.fits")
        self.disk_image_path = fs.join(self.components_images_path, "disk.fits")
        self.model_image_path = fs.join(self.components_images_path, "model.fits")

        # The DustPedia SED path
        self.dustpedia_sed_path = fs.join(self.data_seds_path, "DustPedia.dat")

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):
        return self.environment.galaxy_name

    # -----------------------------------------------------------------

    @property
    def data_path(self):
        return self.environment.data_path

    # -----------------------------------------------------------------

    @property
    def truncated_sed_path(self):
        return self.environment.truncated_sed_path

    # -----------------------------------------------------------------

    @property
    def asymptotic_sed_path(self):
        return self.environment.asymptotic_sed_path

    # -----------------------------------------------------------------

    @property
    def phot_images_path(self):
        return self.environment.phot_images_path

    # -----------------------------------------------------------------

    @property
    def photometry_dataset(self):
        return self.environment.photometry_dataset

    # -----------------------------------------------------------------

    @property
    def static_photometry_dataset(self):
        return self.environment.static_photometry_dataset

    # -----------------------------------------------------------------

    @property
    def photometry_image_names(self):
        return self.environment.photometry_image_names

    # -----------------------------------------------------------------

    @property
    def photometry_image_paths(self):
        return self.environment.photometry_image_paths

    # -----------------------------------------------------------------

    @property
    def photometry_image_paths_for_filters(self):
        return self.environment.photometry_image_paths_for_filters

    # -----------------------------------------------------------------

    @property
    def photometry_image_paths_for_filter_names(self):
        return self.environment.photometry_image_paths_for_filter_names

    # -----------------------------------------------------------------

    def get_photometry_image_path(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return self.environment.get_photometry_image_path(name)

    # -----------------------------------------------------------------

    def get_photometry_image_path_for_filter(self, fltr):

        """
        Thisfunction ...
        :param fltr:
        :return:
        """

        return self.environment.get_photometry_image_path_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_photometry_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.environment.get_photometry_image(name)

    # -----------------------------------------------------------------

    def get_photometry_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return self.environment.get_photometry_image_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_photometry_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.environment.get_photometry_frame(name)

    # -----------------------------------------------------------------

    def get_photometry_frame_for_filter(self, fltr):
        return self.environment.get_photometry_frame_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_photometry_background(self, name):
        return self.environment.get_photometry_background(name)

    # -----------------------------------------------------------------

    def get_photometry_background_for_filter(self, fltr):
        return self.environment.get_photometry_background_for_filter(fltr)

    # -----------------------------------------------------------------

    def get_photometry_mask(self, name):
        return self.environment.get_photometry_mask(name)

    # -----------------------------------------------------------------

    def get_photometry_mask_for_filter(self, fltr):
        return self.environment.get_photometry_mask_for_filter(fltr)

    # -----------------------------------------------------------------

    @property
    def observed_sed_dustpedia_path(self):
        return self.environment.observed_sed_dustpedia_path

    # -----------------------------------------------------------------

    @property
    def galaxy_properties_path(self):
        return self.environment.galaxy_properties_path

    # -----------------------------------------------------------------

    @property
    def galaxy_info_path(self):
        return self.environment.galaxy_info_path

    # -----------------------------------------------------------------

    @property
    def data_seds_path(self):
        return self.environment.data_seds_path

    # -----------------------------------------------------------------

    @property
    def data_images_path(self):
        return self.environment.data_images_path

    # -----------------------------------------------------------------

    @property
    def prep_path(self):
        return self.environment.prep_path

    # -----------------------------------------------------------------

    @property
    def inspect_path(self):
        return self.environment.inspect_path

    # -----------------------------------------------------------------

    @property
    def truncation_path(self):
        return self.environment.truncation_path

    # -----------------------------------------------------------------

    @property
    def phot_path(self):
        return self.environment.phot_path

    # -----------------------------------------------------------------

    @property
    def maps_path(self):
        return self.environment.maps_path

    # -----------------------------------------------------------------

    @property
    def maps_raw_path(self):
        return self.environment.maps_raw_path

    # -----------------------------------------------------------------

    @property
    def maps_components_path(self):
        return self.environment.maps_components_path

    # -----------------------------------------------------------------

    @property
    def components_path(self):
        return self.environment.components_path

    # -----------------------------------------------------------------

    @property
    def playground_path(self):
        return self.environment.playground_path

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_name(self):
        # Get the NGC name from the configuration
        return self.modeling_configuration.ngc_name

    # -----------------------------------------------------------------

    @lazyproperty
    def ngc_name_nospaces(self):
        return self.ngc_name.replace(" ", "")

    # -----------------------------------------------------------------

    @lazyproperty
    def hyperleda_name(self):
        # Get the HYPERLEDA name from the configuration
        return self.modeling_configuration.hyperleda_name

    # -----------------------------------------------------------------

    @property
    def initial_dataset_path(self):
        return self.environment.initial_dataset_path

    # -----------------------------------------------------------------

    @property
    def initial_dataset(self):
        return self.environment.initial_dataset

    # -----------------------------------------------------------------

    @property
    def preparation_names(self):
        return self.environment.preparation_names

    # -----------------------------------------------------------------

    @property
    def prep_names(self):
        return self.preparation_names

    # -----------------------------------------------------------------

    @property
    def preparation_paths(self):
        return self.environment.preparation_paths

    # -----------------------------------------------------------------

    @property
    def dataset_path(self):
        return self.environment.prepared_dataset_path

    # -----------------------------------------------------------------

    @property
    def dataset(self):
        return self.environment.dataset

    # -----------------------------------------------------------------

    @property
    def frame_list(self):
        return self.environment.frame_list

    # -----------------------------------------------------------------

    def get_names_for_filters(self, filters):

        """
        Thisf unction ...
        :param filters:
        :return:
        """

        return self.environment.get_names_for_filters(filters)

    # -----------------------------------------------------------------

    def get_frames_for_filters(self, filters, framelist=False, named=True):

        """
        Thisf unction ...
        :param filters:
        :param framelist:
        :param named:
        :return:
        """

        return self.environment.get_frames_for_filters(filters, framelist=framelist, named=named)

    # -----------------------------------------------------------------

    @property
    def named_frame_list(self):
        return self.environment.named_frame_list

    # -----------------------------------------------------------------

    @property
    def errormap_list(self):
        return self.environment.errormap_list

    # -----------------------------------------------------------------

    @property
    def named_errormap_list(self):
        return self.environment.named_errormap_list

    # -----------------------------------------------------------------

    @property
    def frame_path_list(self):
        return self.environment.frame_path_list

    # -----------------------------------------------------------------

    @property
    def named_frame_path_list(self):
        return self.environment.named_frame_path_list

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_header(self):

        """
        This function ...
        :return:
        """

        # Check wether the file is present
        if not fs.is_file(self.disk_image_path): raise IOError("The disk image is not present. Run the 'decompose' step to create it")

        # Get the header
        return headers.get_header(self.disk_image_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def bulge_header(self):

        """
        This function ...
        :return:
        """

        # Check whether the file is present
        if not fs.is_file(self.bulge_image_path): raise IOError("The bulge image is not present. Run the 'decompose' step to create it")

        # Get the header
        return headers.get_header(self.bulge_image_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def decomposition_filter(self):

        """
        Thisf unction ...
        :return:
        """

        # Obtain the PSF filter of this image
        psf_filter = headers.get_psf_filter(self.disk_header)

        # Return the filter
        return psf_filter

    # -----------------------------------------------------------------

    @lazyproperty
    def decomposition_filter_name(self):
        return str(self.decomposition_filter)

    # -----------------------------------------------------------------

    @lazyproperty
    def decomposition_filter_no_spaces(self):
        return tostr(self.decomposition_filter, delimiter="_")

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

    @property
    def disk_region_path(self):
        return self.environment.disk_region_path

    # -----------------------------------------------------------------

    @property
    def disk_ellipse(self):
        return self.environment.disk_ellipse

    # -----------------------------------------------------------------

    @property
    def disk_position_angle(self):
        return self.environment.disk_position_angle

    # -----------------------------------------------------------------

    @property
    def disk_axial_ratio(self):
        return self.environment.disk_axial_ratio

    # -----------------------------------------------------------------

    @property
    def disk_ellipticity(self):
        return self.environment.disk_ellipticity

    # -----------------------------------------------------------------

    @property
    def disk_inclination(self):
        return self.environment.disk_inclination

    # -----------------------------------------------------------------

    @property
    def significance_levels_path(self):
        return self.environment.significance_levels_path

    # -----------------------------------------------------------------

    @property
    def truncation_ellipse_path(self):
        return self.environment.truncation_ellipse_path

    # -----------------------------------------------------------------

    @property
    def has_truncation_ellipse(self):
        return self.environment.has_truncation_ellipse

    # -----------------------------------------------------------------

    @property
    def truncation_ellipse(self):
        return self.environment.truncation_ellipse

    # -----------------------------------------------------------------

    @property
    def truncation_radius(self):
        return self.environment.truncation_radius

    # -----------------------------------------------------------------

    @property
    def physical_truncation_radius(self):
        return self.environment.physical_truncation_radius

    # -----------------------------------------------------------------

    @property
    def truncation_factor(self):
        return self.environment.truncation_factor

    # -----------------------------------------------------------------

    @property
    def truncation_box_axial_ratio(self):
        return self.environment.truncation_box_axial_ratio

    # -----------------------------------------------------------------

    def get_pixel_truncation_ellipse(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        return self.truncation_ellipse.to_pixel(wcs)

    # -----------------------------------------------------------------

    def get_truncation_mask(self, wcs, invert=True):

        """
        This function ...
        :param wcs:
        :param invert:
        :return:
        """

        ellipse = self.get_pixel_truncation_ellipse(wcs)
        return ellipse.to_mask(wcs.xsize, wcs.ysize, invert=invert)

    # -----------------------------------------------------------------

    @property
    def physical_truncation_ellipse(self):
        return self.environment.physical_truncation_ellipse

    # -----------------------------------------------------------------

    @property
    def truncation_area(self):
        return self.environment.truncation_area

    # -----------------------------------------------------------------

    @property
    def truncation_box(self):
        return self.environment.truncation_box

    # -----------------------------------------------------------------

    @property
    def has_significance_levels(self):
        return self.environment.has_significance_levels

    # -----------------------------------------------------------------

    @property
    def significance_levels(self):
        return self.environment.significance_levels

    # -----------------------------------------------------------------

    def get_significance_level(self, image_name):

        """
        This function ...
        :param image_name:
        :return:
        """

        return self.significance_levels[image_name]

    # -----------------------------------------------------------------

    def get_significance_level_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return self.get_significance_level(str(fltr))

    # -----------------------------------------------------------------

    def get_significance_map(self, observed, errors):

        """
        This function ...
        :param observed:
        :param errors:
        :return:
        """

        return observed / errors

    # -----------------------------------------------------------------

    def get_significance_mask(self, observed, errors, invert=True, min_npixels=1, connectivity=4):

        """
        This function ...
        :param observed:
        :param errors:
        :param invert:
        :param min_npixels:
        :param connectivity:
        :return:
        """

        # Get the level
        level = self.get_significance_level(observed.filter_name)

        # Create the mask
        significance = self.get_significance_map(observed, errors)
        mask = Mask.above(significance, level, wcs=observed.wcs)

        # Only keep largest patch
        mask = mask.largest(npixels=min_npixels, connectivity=connectivity)

        # Fill holes
        mask.fill_holes()

        # Return the mask
        if invert: mask.invert()
        return mask

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

    @property
    def galaxy_properties(self):
        return self.environment.galaxy_properties

    # -----------------------------------------------------------------

    @property
    def galaxy_info(self):
        return self.environment.galaxy_info

    # -----------------------------------------------------------------

    @property
    def hubble_type(self):
        return self.environment.hubble_type

    # -----------------------------------------------------------------

    @property
    def hubble_stage(self):
        return self.environment.hubble_stage

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_stage_type(self):

        """
        This function returns the Hubble type, determined from the stage number
        :return:
        """

        from ..core.stellar_mass import hubble_stage_to_type
        return hubble_stage_to_type(self.hubble_stage)

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_stage_subtype(self):

        """
        This function returns the Hubble subtype, determined from the stage number
        :return:
        """

        from ..core.stellar_mass import hubble_stage_to_type
        return hubble_stage_to_type(self.hubble_stage, add_subtype=True)[1]

    # -----------------------------------------------------------------

    @property
    def galaxy_ellipse(self): # from properties
        return self.environment.galaxy_ellipse

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):
        return self.environment.galaxy_distance

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_center(self):
        return self.environment.galaxy_center

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_inclination(self):
        return self.environment.galaxy_inclination

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_position_angle(self):
        return self.environment.galaxy_position_angle

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_redshift(self):
        return self.environment.galaxy_redshift

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
        path = fs.join(self.components_path, "2D", "S4G", "bulge.mod")
        return load_2d_model(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk2d_model(self):
        path = fs.join(self.components_path, "2D", "S4G", "disk.mod")
        return load_2d_model(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the reference SEDs ...")

        # The seds
        seds = dict()

        # Loop over the SEDs in the data/SEDs directory
        for path, name in fs.files_in_path(self.data_seds_path, extension="dat", returns=["path", "name"], not_contains="Lines"):

            # Open the observed SED
            sed = ObservedSED.from_file(path)

            # Add the SED to the dictionary
            seds[name] = sed

        # Return the SEDs
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_sed_labels(self):
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

    @property
    def has_cache_host(self):
        return self.cache_host_id is not None

    # -----------------------------------------------------------------

    @property
    def cache_host_id(self):
        return self.environment.cache_host_id

    # -----------------------------------------------------------------

    @property
    def cache_host(self):
        return self.environment.cache_host

    # -----------------------------------------------------------------

    @property
    def cache_remote(self):
        return self.environment.cache_remote

    # -----------------------------------------------------------------

    def get_data_image_paths(self):

        """
        This function ...
        :return:
        """

        return get_data_image_paths(self.config.path)

    # -----------------------------------------------------------------

    def get_data_image_paths_with_cached(self, lazy=False, origins=None, not_origins=None, attached=False):

        """
        This function ...
        :param lazy:
        :parma origins:
        :parma not_origins:
        :param attached:
        :return:
        """

        return get_data_image_paths_with_cached(self.config.path, self.cache_host_id, lazy=lazy, origins=origins,
                                                not_origins=not_origins, attached=attached)

    # -----------------------------------------------------------------

    def get_data_image_paths_with_cached_for_origin(self, origin, lazy=False, attached=False):

        """
        This function ...
        :param origin:
        :param lazy:
        :param attached:
        :return:
        """

        return get_data_image_paths_with_cached_for_origin(self.config.path, origin, self.cache_host_id, lazy=lazy, attached=attached)

    # -----------------------------------------------------------------

    def get_data_image_and_error_paths(self):

        """
        This function ...
        :return: 
        """

        return get_data_image_and_error_paths(self.config.path)

    # -----------------------------------------------------------------

    def get_data_image_and_error_paths_with_cached(self, lazy=False):

        """
        This function ...
        :param lazy:
        :return:
        """

        return get_data_image_and_error_paths_with_cached(self.config.path, self.cache_host_id, lazy=lazy)

    # -----------------------------------------------------------------

    @property
    def analysis_context(self):
        return self.environment.analysis_context

    # -----------------------------------------------------------------

    @property
    def analysis_runs(self):
        return self.environment.analysis_runs

    # -----------------------------------------------------------------

    @property
    def cached_analysis_runs(self):
        return self.environment.cached_analysis_runs

    # -----------------------------------------------------------------

    def deproject_maps(self, maps, scale_height, root_path, method="pts", edgeon=False, write=True,
                       return_deprojections=False, downsample_factor=2.):

        """
        This function ...
        :param maps:
        :param scale_height:
        :param root_path:
        :param method:
        :param edgeon:
        :param write:
        :param return_deprojections:
        :param downsample_factor:
        :return:
        """

        from ..misc.deprojector import Deprojector

        # Create the deprojector
        deprojector = Deprojector(cwd=root_path)

        # Check edgeon setting
        if edgeon and method != "skirt": raise ValueError("Edgeon is not possible when method is not 'skirt'")

        # Set modeling path
        deprojector.config.path = self.config.path

        # Set settings
        deprojector.config.method = method
        deprojector.config.writing.deprojections = write
        deprojector.config.writing.maps = False
        deprojector.config.downsample_factor = downsample_factor

        # Run the deprojector
        deprojector.run(maps=maps, scale_height=scale_height, root_path=root_path)

        # Return the deprojected maps
        if edgeon:
            if return_deprojections: return deprojector.deprojections, deprojector.deprojected, deprojector.edgeon
            else: return deprojector.deprojected, deprojector.edgeon
        else:
            if return_deprojections: return deprojector.deprojections, deprojector.deprojected
            else: return deprojector.deprojected

    # -----------------------------------------------------------------

    def deproject_models(self, deprojections, root_path, method="pts", edgeon=False, downsample_factor=2.):

        """
        Thisf unction ...
        :param deprojections:
        :param root_path:
        :param method:
        :param edgeon:
        :param downsample_factor:
        :return:
        """

        from ..misc.deprojector import Deprojector

        # Create the deprojector
        deprojector = Deprojector(cwd=root_path)

        # Check edgeon setting
        if edgeon and method != "skirt": raise ValueError("Edgeon is not possible when method is not 'skirt'")

        # Set the modeling path
        deprojector.config.path = self.config.path

        # Set settings
        deprojector.config.method = method
        deprojector.config.writing.deprojections = False
        deprojector.config.writing.maps = False
        deprojector.config.downsample_factor = downsample_factor

        # Run the deprojector
        deprojector.run(deprojections=deprojections, root_path=root_path)

        # Return the deprojected maps
        if edgeon: return deprojector.deprojected, deprojector.edgeon
        else: return deprojector.deprojected

    # -----------------------------------------------------------------

    def project_models_faceon(self, models, npixels, pixelscale, root_path, write=True, npackages=int(1e7), map_aliases=None):

        """
        This function projects models to face-on views
        :param models:
        :param root_path:
        :param npixels:
        :param pixelscale:
        :param write:
        :param npackages:
        :param map_aliases:
        :return:
        """

        from ..misc.projector import Projector

        # Create the projector
        projector = Projector(cwd=root_path)

        # Get number of pixels
        npixels = get_npixels(npixels)

        # Get field of view
        field = get_field(pixelscale, npixels, self.galaxy_distance)

        # Get the center pixel
        center = get_center(npixels)

        # Set the modeling path
        projector.config.path = self.config.path

        # Set settings
        projector.config.writing.projections = write
        projector.config.faceon = True
        projector.config.edgeon = False
        projector.config.distance = self.galaxy_distance
        projector.config.center = center
        projector.config.npixels = npixels
        projector.config.field = field
        projector.config.npackages = npackages

        # Run the projector
        projector.run(models=models, root_path=root_path, map_aliases=map_aliases)

        # Return the projected maps
        return projector.faceon

    # -----------------------------------------------------------------

    def project_models_edgeon(self, models, npixels, pixelscale, root_path, write=True, npackages=int(1e7), map_aliases=None):

        """
        This function projects models to edge-on views
        :param models:
        :param root_path:
        :param npixels:
        :param pixelscale:
        :param write:
        :param npackages:
        :param map_aliases:
        :return:
        """

        from ..misc.projector import Projector

        # Create the projector
        projector = Projector(cwd=root_path)

        # Get number of pixels
        npixels = get_npixels(npixels)

        # Get field of view
        field = get_field(pixelscale, npixels, self.galaxy_distance)

        # Get the center pixel
        center = get_center(npixels)

        # Set the modeling path
        projector.config.path = self.config.path

        # Set settings
        projector.config.writing.projections = write
        projector.config.faceon = False
        projector.config.edgeon = True
        projector.config.distance = self.galaxy_distance
        projector.config.center = center
        projector.config.npixels = npixels
        projector.config.field = field
        projector.config.npackages = npackages

        # Run the projector
        projector.run(models=models, root_path=root_path, map_aliases=map_aliases)

        # Return the projected maps
        return projector.edgeon

    # -----------------------------------------------------------------

    def project_models_to_wcs(self, models, wcs, root_path, write=True, npackages=int(1e7)):

        """
        This function projects models to an arbritrary view
        :param models:
        :param wcs:
        :param root_path:
        :param write:
        :param npackages:
        :return:
        """

        from ..misc.projector import Projector

        # Set azimuth
        from ...core.units.parsing import parse_angle
        azimuth = parse_angle("0 deg")

        # Create the projector
        projector = Projector(cwd=root_path)

        # Set the modeling path
        projector.config.path = self.config.path

        # Set settings
        projector.config.writing.projections = write
        projector.config.faceon = False
        projector.config.edgeon = False
        projector.config.distance = self.galaxy_distance
        projector.config.center = self.galaxy_center.to_sky(wcs)
        #projector.config.npixels = npixels # only for creating edgeon and faceon
        #projector.config.field = field # only for creating edgeon and faceon
        projector.config.inclination = self.disk_inclination
        projector.config.azimuth = azimuth
        projector.config.position_angle = self.disk_position_angle
        projector.config.npackages = npackages

        # Name for the projection
        name = "single"

        # Run the projector
        projector.run(models=models, name=name, wcs=wcs, root_path=root_path)

        # Return the projected maps for the only projection
        return projector.projected[name]

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

def get_galaxy_properties(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    # get path
    path = get_galaxy_properties_path(modeling_path)

    # Check whether the file is present
    if not fs.is_file(path): raise IOError("The galaxy properties file is not present. Perform 'fetch_properties' to create this file'")

    # Load the properties
    properties = GalaxyProperties.from_file(path)

    # Return the property map
    return properties

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

def get_disk_position_angle(modeling_path):

    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    # Determine the path to the regions file
    components_path = fs.join(modeling_path, "components")
    disk_region_path = fs.join(components_path, "disk.reg")

    # Open the region
    region = SkyRegionList.from_file(disk_region_path)

    # Return the first and only shape
    disk_ellipse = region[0]

    # Return the orientation angle
    return disk_ellipse.angle

# -----------------------------------------------------------------

def get_initial_dataset_path(modeling_path):
    
    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    return fs.join(modeling_path, "prep", "initial_dataset.dat")

# -----------------------------------------------------------------

def get_initial_dataset(modeling_path, check=True):

    """
    This function ...
    :param modeling_path:
    :param check:
    :return:
    """

    return DataSet.from_file(get_initial_dataset_path(modeling_path), check=check)

# -----------------------------------------------------------------

def get_prepared_dataset_path(modeling_path):
    
    """
    This function ...
    :param modeling_path: 
    :return: 
    """

    return fs.join(modeling_path, "prep", "dataset.dat")

# -----------------------------------------------------------------

def get_prepared_dataset(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    return DataSet.from_file(get_prepared_dataset_path(modeling_path))

# -----------------------------------------------------------------

def get_data_image_paths(modeling_path, origins=None, not_origins=None):

    """
    This function ...
    :param modeling_path:
    :param origins:
    :param not_origins:
    :return:
    """

    data_images_path = get_data_images_path(modeling_path)

    paths = dict()

    # Loop over the images
    for image_path, image_name in fs.files_in_path(data_images_path, extension="fits", not_contains="poisson",
                                                   returns=["path", "name"], recursive=True, recursion_level=1):

        # Determine the origin
        origin = fs.name(fs.directory_of(image_path))

        # Check
        if origins is not None and origin not in origins: continue
        if not_origins is not None and origin in not_origins: continue

        # Load the primary image frame
        #frame = load_image_frame(image_path)

        # Get filter name
        header = getheader(image_path)
        fltr = headers.get_filter(image_name, header)

        # Determine name
        name = str(fltr)

        # Add the image path
        paths[name] = image_path

        # Free memory
        #gc.collect()

    # Return the paths
    return paths

# -----------------------------------------------------------------

def get_data_image_paths_for_origin(modeling_path, origin):

    """
    This function ...
    :param modeling_path:
    :param origin:
    :return:
    """

    # Set path
    data_images_path = get_data_images_path(modeling_path)
    origin_path = fs.join(data_images_path, origin)

    # Check path
    if not fs.is_directory(origin_path): return dict()

    paths = dict()

    # Loop over the images
    for image_path, image_name in fs.files_in_path(origin_path, extension="fits", not_contains="poisson", returns=["path", "name"]):

        # Get filter name
        header = getheader(image_path)
        fltr = headers.get_filter(image_name, header)

        # Determine name
        name = str(fltr)

        # Add the image path
        paths[name] = image_path

    # Return the paths
    return paths

# -----------------------------------------------------------------

def get_data_image_paths_with_cached(modeling_path, host_id, lazy=False, origins=None, not_origins=None, attached=False):

    """
    This function ...
    :param modeling_path:
    :param host_id:
    :param lazy:
    :param origins:
    :param not_origins:
    :return:
    """

    paths = get_data_image_paths(modeling_path, origins=origins, not_origins=not_origins)
    paths.update(**get_cached_data_image_paths(modeling_path, host_id, origins=origins, not_origins=not_origins, lazy=lazy, attached=attached))
    return paths

# -----------------------------------------------------------------

def get_data_image_paths_with_cached_for_origin(modeling_path, origin, host_id, lazy=False, attached=False):

    """
    This function ...
    :param modeling_path:
    :param origin:
    :param host_id:
    :param lazy:
    :param attached:
    :return:
    """

    paths = get_data_image_paths_for_origin(modeling_path, origin)
    paths.update(**get_cached_data_image_paths_for_origin(modeling_path, origin, host_id, lazy=lazy, attached=attached))
    return paths

# -----------------------------------------------------------------

def get_cached_data_image_paths(modeling_path, host_id, origins=None, not_origins=None, lazy=False, attached=False):

    """
    This function ...
    :param modeling_path:
    :param host_id:
    :param origins:
    :param not_origins:
    :param lazy:
    :param attached:
    :return:
    """

    # Create the remote and start (detached) python session
    remote = Remote(host_id=host_id)
    if not lazy: session = remote.start_python_session(output_path=remote.pts_temp_path, attached=attached)
    else: session = None

    # Load the environment
    environment = GalaxyModelingEnvironment(modeling_path)

    paths = dict()

    command_name = "initialize_preparation"

    # Get the remote path
    remote_data_path = cached_directory_path_for_single_command(environment, command_name, remote)

    # Loop over the images
    for image_path, image_name in remote.files_in_path(remote_data_path, extension="fits", not_contains="poisson",
                                                   returns=["path", "name"], recursive=True): #, recursion_level=1):

        # Determine the origin
        origin = fs.name(fs.directory_of(image_path))

        # Check
        if origins is not None and origin not in origins: continue
        if not_origins is not None and origin in not_origins: continue

        # Get filter name
        if lazy:
            #name = str(parse_filter(image_name))
            fltr = headers.get_filter(image_name)
            if fltr is None:
                #raise RuntimeError("Could not determine the filter for the '" + image_name + "' image")
                #log.warning("Could not determine the filter for the '" + image_name + "' image: skipping ...")
                #continue
                if session is None: session = remote.start_python_session(output_path=remote.pts_temp_path, attached=attached)
                name = get_filter_name(image_path, session)
            else: name = str(fltr)
        else: name = get_filter_name(image_path, session)

        if name is None: raise RuntimeError("Could not determine the filter name for the '" + image_name + "' image")

        # Add the image path
        paths[name] = image_path

    # Return the paths
    return paths

# -----------------------------------------------------------------

def get_cached_data_image_paths_for_origin(modeling_path, origin, host_id, lazy=False, attached=False):

    """
    This function ...
    :param modeling_path:
    :param origin:
    :param host_id:
    :param lazy:
    :param attached:
    :return:
    """

    # Create the remote and start (detached) python session
    remote = Remote(host_id=host_id)
    if not lazy: session = remote.start_python_session(output_path=remote.pts_temp_path, attached=attached)
    else: session = None

    # Load the environment
    environment = GalaxyModelingEnvironment(modeling_path)

    paths = dict()

    command_name = "initialize_preparation"

    # Get the remote path
    remote_data_path = cached_directory_path_for_single_command(environment, command_name, remote)
    remote_origin_path = fs.join(remote_data_path, origin)

    # Check path
    if not remote.is_directory(remote_origin_path): return dict() #raise ValueError("The origin directory does not exist")

    # Loop over the images
    for image_path, image_name in remote.files_in_path(remote_origin_path, extension="fits", not_contains="poisson", returns=["path", "name"]):

        # Get filter name
        if lazy:
            # name = str(parse_filter(image_name))
            fltr = headers.get_filter(image_name)
            if fltr is None:
                # raise RuntimeError("Could not determine the filter for the '" + image_name + "' image")
                #log.warning("Could not determine the filter for the '" + image_name + "' image: skipping ...")
                #continue
                if session is None: session = remote.start_python_session(output_path=remote.pts_temp_path, attached=attached)
                name = get_filter_name(image_path, session)
            name = str(fltr)
        else: name = get_filter_name(image_path, session)

        if name is None: raise RuntimeError("Could not determine the filter name for the '" + image_name + "' image")

        # Add the image path
        paths[name] = image_path

    # Return the paths
    return paths

# -----------------------------------------------------------------

def get_data_image_and_error_paths(modeling_path):

    """
    This function ...
    :return: 
    """

    data_images_path = get_data_images_path(modeling_path)

    paths = dict()
    error_paths = dict()

    # Loop over the images
    for image_path, image_name in fs.files_in_path(data_images_path, extension="fits", not_contains="poisson",
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
        elif needs_poisson_errors(frame.filter): raise RuntimeError("Poisson error frame not found for the " + name + " image. Run the appropriate command to create the mosaics and poisson frames.")

        # Free memory
        gc.collect()

    # Return the paths and error paths
    return paths, error_paths

# -----------------------------------------------------------------

def get_data_image_and_error_paths_with_cached(modeling_path, host_id, lazy=False):

    """
    This function ...
    :param modeling_path:
    :param host_id:
    :param lazy:
    :return:
    """

    paths, error_paths = get_data_image_and_error_paths(modeling_path)
    cached_paths, cached_error_paths = get_cached_data_image_and_error_paths(modeling_path, host_id, lazy=lazy)
    paths.update(**cached_paths)
    error_paths.update(**cached_error_paths)
    return paths, error_paths

# -----------------------------------------------------------------

def get_cached_data_image_and_error_paths(modeling_path, host_id, lazy=False):

    """
    This function ...
    :param modeling_path:
    :param host_id:
    :param lazy:
    :return:
    """

    # Create the remote and start (detached) python session
    if isinstance(host_id, Remote): remote = host_id
    else: remote = Remote(host_id=host_id)
    if not lazy: session = remote.start_python_session(output_path=remote.pts_temp_path)
    else: session = None

    # Load the environment
    environment = GalaxyModelingEnvironment(modeling_path)

    paths = dict()
    error_paths = dict()

    command_name = "initialize_preparation"

    # Get the remote path
    remote_data_path = cached_directory_path_for_single_command(environment, command_name, remote)

    # Loop over the images
    for image_path, image_name in remote.files_in_path(remote_data_path, extension="fits", not_contains="poisson",
                                                   returns=["path", "name"], recursive=True): #, recursion_level=1):

        # Determine directory path
        path = fs.directory_of(image_path)

        # Get filter name
        if lazy:
            fltr = headers.get_filter(image_name)
            if fltr is None:
                log.warning("Could not determine the filter of the '" + image_name + "' image: skipping ...")
                continue
            else: name = str(fltr)
        else: name = get_filter_name(image_path, session)

        if name is None: raise RuntimeError("Could not determine the filter name for the '" + image_name + "' image")

        # Add the image path
        paths[name] = image_path

        # Determine path to poisson error map
        poisson_path = fs.join(path, image_name + "_poisson.fits")

        # Set the path to the poisson error map
        if remote.is_file(poisson_path):

            # Debugging
            log.debug("Poisson error frame found for '" + name + "' image ...")
            error_paths[name] = poisson_path

        # Poisson frame not present
        elif needs_poisson_errors(name): raise RuntimeError("Poisson error frame not found for the " + name + " image. Run the appropriate command to create the mosaics and poisson frames.")

    # Return the paths
    return paths, error_paths

# -----------------------------------------------------------------

def get_bulge2d_model_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    components_path = fs.join(modeling_path, "components")
    return fs.join(components_path, "2D", "S4G", "bulge.mod")

# -----------------------------------------------------------------

def has_bulge2d_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_bulge2d_model_path(modeling_path)
    return fs.is_file(path)

# -----------------------------------------------------------------

def get_bulge2d_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_bulge2d_model_path(modeling_path)
    return load_2d_model(path)

# -----------------------------------------------------------------

def get_disk2d_model_path(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    components_path = fs.join(modeling_path, "components")
    return fs.join(components_path, "2D", "S4G", "disk.mod")

# -----------------------------------------------------------------

def has_disk2d_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_disk2d_model_path(modeling_path)
    return fs.is_file(path)

# -----------------------------------------------------------------

def get_disk2d_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_disk2d_model_path(modeling_path)
    return load_2d_model(path)

# -----------------------------------------------------------------

def get_bulge_model_path(modeling_path):

    """
    Thisf unction ...
    :param modeling_path:
    :return:
    """

    components_path = fs.join(modeling_path, "components")
    return fs.join(components_path, "models", "bulge.mod")

# -----------------------------------------------------------------

def has_bulge_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_bulge_model_path(modeling_path)
    return fs.is_file(path)

# -----------------------------------------------------------------

def get_bulge_model(modeling_path):

    """
    Thisj function ...
    :param modeling_path:
    :return:
    """

    path = get_bulge_model_path(modeling_path)
    return load_3d_model(path)

# -----------------------------------------------------------------

def get_disk_model_path(modeling_path):

    """
    Thisj function ...
    :param modeling_path:
    :return:
    """

    components_path = fs.join(modeling_path, "components")
    return fs.join(components_path, "models", "disk.mod")

# -----------------------------------------------------------------

def has_disk_model(modeling_path):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    path = get_disk_model_path(modeling_path)
    return fs.is_file(path)

# -----------------------------------------------------------------

def get_disk_model(modeling_path):

    """
    Thisf unction ...
    :param modeling_path:
    :return:
    """

    path = get_disk_model_path(modeling_path)
    return load_3d_model(path)

# -----------------------------------------------------------------
