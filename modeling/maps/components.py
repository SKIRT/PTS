#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.components Contains the ComponentMapsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import astronomical modules
from astropy.units import dimensionless_angles

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .selectioncomponent import MapsSelectionComponent
from ...core.basics.configuration import prompt_real
from ...core.tools.utils import lazyproperty
from ...magic.core.frame import Frame
from ...core.tools import filesystem as fs
from ...core.tools import numbers
from ...core.basics.range import RealRange
from ...core.tools import sequences
from ...core.remote.remote import Remote
from ...magic.core.mask import Mask, intersection
from ...magic.core.alpha import AlphaMask, load_mask_or_alpha_mask
from ...magic.core.detection import Detection
from ...core.tools.stringify import tostr
from ...core.basics.range import QuantityRange
from ...core.units.parsing import parse_unit as u
from ...magic.basics.vector import PixelShape
from ...core.filter.filter import parse_filter
from ...core.tools.parsing import real
from ...magic.core.image import Image
from ...magic.tools import colours
from ...magic.basics.mask import Mask as oldMask
from ...core.basics.containers import ordered_by_key
from ...magic.core.fits import get_plane_names
from ...core.tools import formatting as fmt
from ...core.tools.serialization import write_dict, load_dict
from ...magic.basics.mask import MaskBase
from ...magic.region.region import PixelRegion
from ...magic.tools import plotting
from .collection import MapsCollection, StaticMapsCollection

# -----------------------------------------------------------------

negatives_plane_name = "negatives"
nans_plane_name = "nans"

# -----------------------------------------------------------------

steps_name = "steps"
masks_name = "masks"
deprojected_name = "deprojected"
deprojected_skirt_name = "deprojected_skirt"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

correct_step = "corrected"
crop_step = "cropped"
interpolate_negatives_step = "interpolated_negatives"
interpolate_step = "interpolated"
truncate_step = "truncated"
clip_step = "clipped"
softened_step = "softened"
steps = [correct_step, crop_step, interpolate_negatives_step, interpolate_step, truncate_step, clip_step, softened_step]

# -----------------------------------------------------------------

clip_suffix = "clip"
softening_suffix = "softening"

# -----------------------------------------------------------------

def steps_before(step):

    """
    This function ...
    :param step:
    :return:
    """

    before = []
    for stepi in steps:
        if stepi == step: break
        before.append(stepi)
    return before

# -----------------------------------------------------------------

def steps_before_and_including(step):

    """
    This function ...
    :param step:
    :return:
    """

    return steps_before(step) + [step]

# -----------------------------------------------------------------

def steps_after(step):

    """
    This function ...
    :param step:
    :return:
    """

    after = []
    for stepi in reversed(steps):
        if stepi == step: break
        after.append(stepi)
    return list(reversed(after))

# -----------------------------------------------------------------

def steps_after_and_including(step):

    """
    This function ...
    :param step:
    :return:
    """

    return [step] + steps_after(step)

# -----------------------------------------------------------------

cropped_negatives_filename = "cropped_negatives.fits"
central_ellipse_filename = "central_ellipse.reg"
negatives_mask_filename = "negatives.fits"
negatives_dilated_filename = "negatives_dilated.fits"
negatives_filled_filename = "negatives_filled.fits"

interpolation_mask_filename = "interpolation_mask.fits"
interpolation_negatives_mask_filename = "interpolation_negatives_mask.fits"
interpolation_background_filename = "interpolation_background.fits"
interpolation_negatives_background_filename = "interpolation_negatives_background.fits"

# -----------------------------------------------------------------

interpolation_ellipse_filename = "interpolation_ellipse.reg"

# -----------------------------------------------------------------

class ComponentMapsMaker(MapsSelectionComponent):

    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ComponentMapsMaker, self).__init__(*args, **kwargs)

        # Path to save intermediate results
        self.old_steps_path = None
        self.young_steps_path = None
        self.ionizing_steps_path = None
        self.dust_steps_path = None

        # Mask directory paths
        self.old_masks_path = None
        self.young_masks_path = None
        self.ionizing_masks_path = None
        self.dust_masks_path = None

        # Deprojection directory paths
        self.old_deprojection_path = None
        self.young_deprojection_path = None
        self.ionizing_deprojection_path = None
        self.dust_deprojection_path = None

        # Deprojection (SKIRT) directory paths
        self.old_deprojection_skirt_path = None
        self.young_deprojection_skirt_path = None
        self.ionizing_deprojection_skirt_path = None
        self.dust_deprojection_skirt_path = None

        # Edgeon directory paths
        self.old_edgeon_path = None
        self.young_edgeon_path = None
        self.ionizing_edgeon_path = None
        self.dust_edgeon_path = None

        # Negative masks
        self.old_negative_masks = dict()
        self.young_negative_masks = dict()
        self.ionizing_negative_masks = dict()
        self.dust_negative_masks = dict()

        # The clip masks
        self.old_clip_masks = dict()
        self.young_clip_masks = dict()
        self.ionizing_clip_masks = dict()
        self.dust_clip_masks = dict()

        # The softening masks
        self.old_softening_masks = dict()
        self.young_softening_masks = dict()
        self.ionizing_softening_masks = dict()
        self.dust_softening_masks = dict()

        # The deprojections
        self.old_deprojections = None
        self.young_deprojections = None
        self.ionizing_deprojections = None
        self.dust_deprojections = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Load the selection
        self.load_selection()

        # 3. Set the S/N levels
        #self.set_levels()

        # 4. Set rerun
        if self.rerun or self.rerun_all: self.set_rerun()

        # 5. Set redeproject
        if self.redeproject: self.set_redeproject()

        # 6. Set redeproject SKIRT
        if self.redeproject_skirt: self.set_redeproject_skirt()

        # 7. Set reproject
        if self.reproject: self.set_reproject()

        # 8. Remove other
        if self.remove: self.remove_other()

        # 9. Load the maps
        self.load_maps()

        # 10. Load the masks
        self.load_masks()

        # 11. Process the maps
        self.process_maps()

        # 12. Deproject the maps
        self.deproject()

        # 13. Project the maps to the edge-on view
        self.project()

        # 14. Writing
        self.write()

        # 15. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def load_collection(self):
        return MapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    def load_static_collection(self):
        return StaticMapsCollection(self.maps_raw_path)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsMaker, self).setup(**kwargs)

        # Set the number of allowed open file handles
        fs.set_nallowed_open_files(self.config.nopen_files)

        # Clear all?
        if self.config.clear_all: self.config.clear_old = self.config.clear_young = self.config.clear_ionizing = self.config.clear_dust = True

        # Clear?
        if self.config.clear_old: fs.clear_directory(self.old_component_maps_path)
        if self.config.clear_young: fs.clear_directory(self.young_component_maps_path)
        if self.config.clear_ionizing: fs.clear_directory(self.ionizing_component_maps_path)
        if self.config.clear_dust: fs.clear_directory(self.dust_component_maps_path)

        # Set the steps paths
        self.old_steps_path = fs.join(self.old_component_maps_path, steps_name)
        self.young_steps_path = fs.join(self.young_component_maps_path, steps_name)
        self.ionizing_steps_path = fs.join(self.ionizing_component_maps_path, steps_name)
        self.dust_steps_path = fs.join(self.dust_component_maps_path, steps_name)

        # Create
        if self.config.steps and not fs.is_directory(self.old_steps_path): fs.create_directory(self.old_steps_path)
        if self.config.steps and not fs.is_directory(self.young_steps_path): fs.create_directory(self.young_steps_path)
        if self.config.steps and not fs.is_directory(self.ionizing_steps_path): fs.create_directory(self.ionizing_steps_path)
        if self.config.steps and not fs.is_directory(self.dust_steps_path): fs.create_directory(self.dust_steps_path)

        # Masks directories
        self.old_masks_path = fs.create_directory_in(self.old_component_maps_path, masks_name)
        self.young_masks_path = fs.create_directory_in(self.young_component_maps_path, masks_name)
        self.ionizing_masks_path = fs.create_directory_in(self.ionizing_component_maps_path, masks_name)
        self.dust_masks_path = fs.create_directory_in(self.dust_component_maps_path, masks_name)

        # Deprojected directories
        self.old_deprojection_path = fs.create_directory_in(self.old_component_maps_path, deprojected_name)
        self.young_deprojection_path = fs.create_directory_in(self.young_component_maps_path, deprojected_name)
        self.ionizing_deprojection_path = fs.create_directory_in(self.ionizing_component_maps_path, deprojected_name)
        self.dust_deprojection_path = fs.create_directory_in(self.dust_component_maps_path, deprojected_name)

        # Deprojected with SKIRT directories
        self.old_deprojection_skirt_path = fs.create_directory_in(self.old_component_maps_path, deprojected_skirt_name)
        self.young_deprojection_skirt_path = fs.create_directory_in(self.young_component_maps_path, deprojected_skirt_name)
        self.ionizing_deprojection_skirt_path = fs.create_directory_in(self.ionizing_component_maps_path, deprojected_skirt_name)
        self.dust_deprojection_skirt_path = fs.create_directory_in(self.dust_component_maps_path, deprojected_skirt_name)

        # Edgeon directories
        self.old_edgeon_path = fs.create_directory_in(self.old_component_maps_path, edgeon_name)
        self.young_edgeon_path = fs.create_directory_in(self.young_component_maps_path, edgeon_name)
        self.ionizing_edgeon_path = fs.create_directory_in(self.ionizing_component_maps_path, edgeon_name)
        self.dust_edgeon_path = fs.create_directory_in(self.dust_component_maps_path, edgeon_name)

        # Clear results?
        if self.config.clear_results: self.config.clear_results_old = self.config.clear_results_young = self.config.clear_results_ionizing = self.config.clear_results_dust = True

        # Clear results?
        if self.config.clear_results_old: self.clear_results_old()
        if self.config.clear_results_young: self.clear_results_young()
        if self.config.clear_results_ionizing: self.clear_results_ionizing()
        if self.config.clear_results_dust: self.clear_results_dust()

        # Levels
        if self.from_previous_levels: self.load_previous_levels()
        elif self.config.levels is not None: self.levels = self.config.levels

        # Set rerun options
        if self.config.rerun is not None:
            self.config.rerun_old = self.config.rerun
            self.config.rerun_young = self.config.rerun
            self.config.rerun_ionizing = self.config.rerun
            self.config.rerun_dust = self.config.rerun

        # Set rerun all options
        if self.config.rerun_all: self.config.rerun_all_old = self.config.rerun_all_young = self.config.rerun_all_ionizing = self.config.rerun_all_dust = True

        # Set redeproject options
        if self.config.redeproject:
            self.config.redeproject_old = True
            self.config.redeproject_young = True
            self.config.redeproject_ionizing = True
            self.config.redeproject_dust = True

        # Set redeproject SKIRT options
        if self.config.redeproject_skirt:
            self.config.redeproject_skirt_old = True
            self.config.redeproject_skirt_young = True
            self.config.redeproject_skirt_ionizing = True
            self.config.redeproject_skirt_dust = True

        # Set reproject options
        if self.config.reproject:
            self.config.reproject_old = True
            self.config.reproject_young = True
            self.config.reproject_ionizing = True
            self.config.reproject_dust = True

    # -----------------------------------------------------------------

    def clear_results_old(self):

        """
        This fucntion ...
        :return:
        """

        # Inform the user
        log.info("Clearing results from old stellar component maps ...")

        # Remove all FITS files in the directory
        fs.remove_files_in_path(self.old_component_maps_path, extension="fits")

        # Clear masks
        fs.clear_directory(self.old_masks_path)

        # Clear deprojected
        fs.clear_directory(self.old_deprojection_path)

        # Clear deprojected with SKIRT
        fs.clear_directory(self.old_deprojection_skirt_path)

        # Clear edgeon
        fs.clear_directory(self.old_edgeon_path)

    # -----------------------------------------------------------------

    def clear_results_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing results from young stellar component maps ...")

        # Remove all FITS files in the directory
        fs.remove_files_in_path(self.young_component_maps_path, extension="fits")

        # Clear masks
        fs.clear_directory(self.young_masks_path)

        # Clear deprojected
        fs.clear_directory(self.young_deprojection_path)

        # Clear deprojected with SKIRT
        fs.clear_directory(self.young_deprojection_skirt_path)

        # Clear edgeon
        fs.clear_directory(self.young_edgeon_path)

    # -----------------------------------------------------------------

    def clear_results_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing results from ionizing stellar component maps ...")

        # Remove all FITS files in the directory
        fs.remove_files_in_path(self.ionizing_component_maps_path, extension="fits")

        # Clear masks
        fs.clear_directory(self.ionizing_masks_path)

        # Clear deprojected
        fs.clear_directory(self.ionizing_deprojection_path)

        # Clear deprojected with SKIRT
        fs.clear_directory(self.young_deprojection_skirt_path)

        # Clear edgeon
        fs.clear_directory(self.ionizing_edgeon_path)

    # -----------------------------------------------------------------

    def clear_results_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing results from dust component maps ...")

        # Remove all FITS files in the directory
        fs.remove_files_in_path(self.dust_component_maps_path, extension="fits")

        # Clear masks
        fs.clear_directory(self.dust_masks_path)

        # Clear deprojected
        fs.clear_directory(self.dust_deprojection_path)

        # Clear deprojected with SKIRT
        fs.clear_directory(self.dust_deprojection_skirt_path)

        # Clear edgeon
        fs.clear_directory(self.dust_edgeon_path)

    # -----------------------------------------------------------------

    def load_selection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.debug("Loading maps selections ...")

        # Determine the path
        path = fs.join(self.maps_components_path, "selection_" + str(self.config.selection) + ".dat")

        # Load the selection file
        selection = load_dict(path)

        # Set the selections
        self.old_selection = selection["old"]
        self.young_selection = selection["young"]
        self.ionizing_selection = selection["ionizing"]
        self.dust_selection = selection["dust"]

    # -----------------------------------------------------------------

    def set_levels(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Setting levels ...")

        # Prompt for levels
        if not self.has_levels: self.prompt_levels()

        # Check given levels
        else: self.check_levels()

        # Write levels
        if not self.from_previous_levels: self.write_levels()

    # -----------------------------------------------------------------

    def load_previous_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading previous sigma levels ...")

        # Determine the path
        path = fs.join(self.maps_components_path, "levels_" + str(self.config.previous_levels) + ".dat")

        # Load the levels file
        self.levels = load_dict(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def rerun_steps_old(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_old is None: return []
        else: return steps_after_and_including(self.config.rerun_old)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_rerun_steps_old(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.rerun_steps_old) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def not_rerun_steps_old(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_old is None: return steps
        else: return steps_before(self.config.rerun_old)

    # -----------------------------------------------------------------

    @lazyproperty
    def rerun_steps_young(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_young is None: return []
        else: return steps_after_and_including(self.config.rerun_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_rerun_steps_young(self):

        """
        This function ...
        :return:
        """

        return len(self.rerun_steps_young) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def not_rerun_steps_young(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_young is None: return steps
        else: return steps_before(self.config.rerun_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def rerun_steps_ionizing(self):

        """
        Thisn function ...
        :return:
        """

        if self.config.rerun_ionizing is None: return []
        else: return steps_after_and_including(self.config.rerun_ionizing)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_rerun_steps_ionizing(self):

        """
        Thisf unction ...
        :return:
        """

        return len(self.rerun_steps_ionizing) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def not_rerun_steps_ionizing(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_ionizing is None: return steps
        else: return steps_before(self.config.rerun_ionizing)

    # -----------------------------------------------------------------

    @lazyproperty
    def rerun_steps_dust(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_dust is None: return []
        else: return steps_after_and_including(self.config.rerun_dust)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_rerun_steps_dust(self):

        """
        This function ...
        :return:
        """

        return len(self.rerun_steps_dust) > 0

    # -----------------------------------------------------------------

    @lazyproperty
    def not_rerun_steps_dust(self):

        """
        This function ...
        :return:
        """

        if self.config.rerun_dust is None: return steps
        else: return steps_before(self.config.rerun_dust)

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        if self.config.remote is None: return None
        else: return Remote(host_id=self.config.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.old_selection: origins.update(self.old_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.young_selection: origins.update(self.young_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.ionizing_selection: origins.update(self.ionizing_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        for name in self.dust_selection: origins.update(self.dust_map_origins[name])
        return list(origins)

    # -----------------------------------------------------------------

    @lazyproperty
    def selection_origins(self):

        """
        This function ...
        :return:
        """

        origins = set()
        origins.update(self.old_selection_origins)
        origins.update(self.young_selection_origins)
        origins.update(self.ionizing_selection_origins)
        origins.update(self.dust_selection_origins)
        return list(origins)

    # -----------------------------------------------------------------

    def auto_select(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate maps ...")

        # Make selection
        old_name, young_name, ionizing_name, dust_name = self.auto_select_maps()

        # Show selections
        log.info("Selected the following maps automaticaly:")
        log.info("")
        log.info(" - " + fmt.bold + "old stellar disk: " + fmt.reset + old_name)
        log.info(" - " + fmt.bold + "young stellar disk: " + fmt.reset + young_name)
        log.info(" - " + fmt.bold + "ionizing stellar disk: " + fmt.reset + ionizing_name)
        log.info(" - " + fmt.bold + "dust disk: " + fmt.reset + dust_name)
        log.info("")

        # Set selections
        self.old_selection = [old_name]
        self.young_selection = [young_name]
        self.ionizing_selection = [ionizing_name]
        self.dust_selection = [dust_name]

    # -----------------------------------------------------------------

    def auto_select_maps(self):

        """
        This function ...
        :return:
        """

        # sSFR
        ssfr_method, ssfr_name = self.auto_select_ssfr_map()

        # TIR
        tir_method, tir_name = self.auto_select_tir_map()

        # Attenuation
        attenuation_method, attenuation_name = self.auto_select_attenuation_map(tir_method, tir_name, ssfr_method, ssfr_name)

        # Old
        old_method, old_name = self.auto_select_old_map()

        # Dust
        dust_method, dust_name = self.auto_select_dust_map(attenuation_method, attenuation_name)

        # Young
        young_name = self.auto_select_young_map(attenuation_method, attenuation_name, old_name)

        # Hot dust
        hot_dust_method, hot_dust_name = self.auto_select_hot_dust_map(old_name)

        # Ionizing
        ionizing_name = self.auto_select_ionizing_map(hot_dust_name)

        # Make full names
        #full_old_name = old_method + "_" + old_name
        full_old_name = old_name # because we only consider disk maps
        full_young_name = young_name
        full_ionizing_name = ionizing_name
        full_dust_name = dust_method + "_" + dust_name

        # Return the full names
        return full_old_name, full_young_name, full_ionizing_name, full_dust_name

    # -----------------------------------------------------------------

    def auto_select_ssfr_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate sSFR map ...")

        preferred_method = "colours"
        preferred_colours = ["FUV-r", "FUV-H", "FUV-i", "FUV-g"]

        # Select sSFR
        if not self.has_ssfr_maps: raise IOError("No sSFR maps are present")
        if not self.ssfr_has_methods: raise IOError("Place the contents of the sSFR maps inside a 'colours' directory")

        # Get sSFR method names
        ssfr_methods = self.ssfr_map_methods
        if preferred_method not in ssfr_methods: raise RuntimeError("Cannot make automatic choice when sSFR maps are not based on colours")

        # Get the sSFR colour map names
        map_names = self.ssfr_map_names_for_method("colours")

        # Select the preferred sSFR colour map name
        ssfr_map_name = None

        # Loop over the preferred colours
        for colour in preferred_colours:

            map_name = find_matching_colour_map_name(colour, map_names)
            if map_name is not None:
                ssfr_map_name = map_name
                break

        # Check
        if ssfr_map_name is None: raise RuntimeError("Cannot make automatic choice: none of the expected sSFR colour maps are present")

        # Return method and map name
        return preferred_method, ssfr_map_name

    # -----------------------------------------------------------------

    def auto_select_tir_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate TIR map ...")

        # Preferred methods
        preferred_methods = ["multi", "single"]

        # Check
        if not self.has_tir_maps: raise IOError("No TIR maps are present")
        if not self.tir_has_methods: raise IOError("No methods for the TIR maps")

        # Get TIR method names
        tir_methods = self.tir_map_methods

        # The best combination
        best = None
        best_nfilters = None
        best_has_spire = None

        # Loop over the preferred methods
        for method in preferred_methods:

            # No maps for this method
            if method not in tir_methods: continue

            # Get the map names for this method
            map_names = self.tir_map_names_for_method(method)

            # Loop over the maps
            for name in map_names:

                # Get the filters (origins) for this map
                filters = self.tir_origins[method][name]
                nfilters = len(filters)
                has_spire = sequences.contains_any(filters, self.spire_filters)

                if best is None:

                    best = (method, name)
                    best_nfilters = nfilters
                    best_has_spire = has_spire

                elif nfilters > best_nfilters:

                    best = (method, name)
                    best_nfilters = nfilters
                    best_has_spire = has_spire

                elif nfilters == best_nfilters and best_has_spire and not has_spire:

                    best = (method, name)
                    best_nfilters = nfilters
                    best_has_spire = has_spire

        # Return the best map
        return best

    # -----------------------------------------------------------------

    def auto_select_attenuation_map(self, tir_method, tir_name, ssfr_method, ssfr_name):

        """
        This function ...
        :param tir_method:
        :param tir_name:
        :param ssfr_method:
        :param ssfr_name:
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate FUV attenuation map ...")

        # Check
        if not self.has_attenuation_maps: raise IOError("No atttenuation maps are present")
        if not self.attenuation_has_methods: raise IOError("No methods for the attenuation maps")

        # Get attenuation method names
        attenuation_methods = self.attenuation_map_methods

        # Prefer Cortese, Buat otherwise
        if "cortese" in attenuation_methods: return self.auto_select_cortese_attenuation_map(tir_method, tir_name, ssfr_method, ssfr_name)
        elif "buat" in attenuation_methods: return self.auto_select_buat_attenuation_map(tir_method, tir_name)
        else: raise ValueError("Cannot find a proper attenuation map method")

    # -----------------------------------------------------------------

    def auto_select_cortese_attenuation_map(self, tir_method, tir_name, ssfr_method, ssfr_name):

        """
        This function ...
        :param tir_method:
        :param tir_name:
        :param ssfr_method:
        :param ssfr_name:
        :return:
        """

        # Inform the user
        log.debug("Automatically selecting appropriate Cortese attenuation map ...")

        # Loop over the map names
        map_name = None
        for name in self.attenuation_map_names_for_method("cortese"):

            if tir_name in name and ssfr_name in name:
                map_name = name
                break

        # Not found?
        if map_name is None: raise RuntimeError("Something went wrong finding the required map")

        # Return
        return "cortese", map_name

    # -----------------------------------------------------------------

    def auto_select_buat_attenuation_map(self, tir_method, tir_name):

        """
        This function ...
        :param tir_method:
        :param tir_name:
        :return:
        """

        # Inform the user
        log.debug("Automatically selecting appropriate Buat attenuation map ...")

        preferred_filters = ["FUV", "NUV"]

        map_name = None

        # Loop over the preferred filters
        for uv_filter_name in preferred_filters:

            # Loop over the maps for this UV filter
            for name in self.attenuation_map_names_for_method("buat"):
                if not name.startswith(uv_filter_name): continue

                if tir_name in name:
                    map_name = name
                    break

        # Not found?
        if map_name is None: raise RuntimeError("Something went wrong finding the required map")

        # Return
        return "buat", map_name

    # -----------------------------------------------------------------

    def auto_select_old_map(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate old stellar map ...")

        preferred_filters = ["IRAC I1", "IRAC I2", "2MASS H", "2MASS K", "WISE W1", "WISE W2"]

        # Check
        if not self.has_old_maps: raise IOError("No old stellar maps are present")
        if not self.old_has_methods: raise IOError("No methods for the old stellar maps")

        method_name = "disk"
        if method_name not in self.old_map_methods: raise ValueError("'disk' method is not present among the old stellar maps")

        # Loop over the preferred filters
        map_name = None
        for filter_name in preferred_filters:

            fltr = parse_filter(filter_name)

            # Loop over the maps
            for name in self.old_map_names_for_method(method_name):

                map_filter = parse_filter(name)
                if map_filter == fltr:
                    map_name = name
                    break

        # No map?
        if map_name is None: raise ValueError("No appropriate old stellar disk map was found")

        # Return
        return method_name, map_name

    # -----------------------------------------------------------------

    def auto_select_dust_map(self, attenuation_method, attenuation_name):

        """
        Thisfunction ...
        :param attenuation_method:
        :param attenuation_name:
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate dust map ...")

        preferred_method = "attenuation"

        # Check
        if not self.has_dust_maps: raise IOError("No dust maps are present")
        if not self.dust_has_methods: raise IOError("No methods for the dust maps")

        if preferred_method not in self.dust_map_methods: raise ValueError("'attenuation' method is not present among the dust maps")

        # Loop over the attenuation dust maps
        map_name = None
        for name in self.dust_map_names_for_method(preferred_method):

            # Check whether the right one
            if attenuation_method in name and attenuation_name in name:
                map_name = name
                break

        # Check if found
        if map_name is None: raise RuntimeError("Appropriate dust map not found")

        # Return
        return preferred_method, map_name

    # -----------------------------------------------------------------

    def auto_select_young_map(self, attenuation_method, attenuation_name, old_name):

        """
        This function ...
        :param attenuation_method:
        :param attenuation_name:
        :param old_name:
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate young stellar map ...")

        # Check
        if not self.has_young_maps: raise IOError("No young stellar maps are present")
        if self.young_has_methods: raise IOError("Didn't expect different methods for the young stellar maps")

        # Get the paths to the young stellar maps
        paths = self.get_young_map_paths()

        # Initialize dictionary to store the number of negatives for each factor
        nnegatives_dict = dict()
        names_dict = dict()

        # Loop over the young stellar maps made with the specific attenuation map
        # AND WITH THE SPECIFIC OLD STELLAR FILTER
        for name in self.young_map_names_no_methods:
            if not (attenuation_method in name and attenuation_name in name): continue
            if not old_name in name: continue

            # Get the factor
            factor = real(name.split("__")[-1])

            # Open the map to get the number of negatives in a central ellipse
            map_path = paths[name]
            image = Image.from_file(map_path)
            if "negatives" not in image.masks:
                log.warning("Negatives mask not present in the '" + name + "' young stellar map image: skipping ...")
                continue

            # Count the number of negatives
            mask = image.masks["negatives"]
            if isinstance(mask, oldMask): mask = Mask(mask, wcs=image.wcs) # Fix type
            nnegatives = mask.relative_nmasked_in(self.central_ellipse)

            # Add to dictionary
            nnegatives_dict[factor] = nnegatives
            names_dict[factor] = name

        # Check if anything is found
        if len(nnegatives_dict) == 0: raise RuntimeError("No appropriate young stellar maps were found")

        # If only one is found
        if len(nnegatives_dict) == 1:
            log.warning("Map was found for only one factor")
            factor = names_dict.keys()[0]
            return names_dict[factor]

        # Sort
        nnegatives_dict = ordered_by_key(nnegatives_dict)

        # Find the factor
        factor = find_factor_max_nnegatives(nnegatives_dict, self.config.young_max_nnegatives)

        # Get the corresponding map name
        map_name = names_dict[factor]

        # Return the map name
        return map_name

    # -----------------------------------------------------------------

    def auto_select_hot_dust_map(self, old_name):

        """
        This function ...
        :param old_name:
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate hot dust map ...")

        method_name = "hot"

        # Check
        if not self.has_dust_maps: raise IOError("No dust maps are present")
        if not self.dust_has_methods: raise IOError("No methods for the dust maps")

        # No hot dust maps
        if method_name not in self.dust_map_methods:
            log.warning("No hot dust maps are present")
            return None, None

        # Get the paths to the hot dust maps
        paths = self.get_dust_map_paths(flatten=False, method=method_name)

        # Iinitialize dictionary to store the number of negatives for each factor
        nnegatives_dict = dict()
        names_dict = dict()

        # Loop over the hot dust maps correct with the same old stellar filter as the old stellar disk map
        for name in self.dust_map_names_for_method(method_name):
            if not name.startswith(old_name): continue

            # Get the factor
            factor = real(name.split("__")[1])

            # Open the map to get the number of negatives in a central ellipse
            map_path = paths[name]
            image = Image.from_file(map_path)
            if "negatives" not in image.masks:
                log.warning("Negatives mask not present in the '" + name + "' hot dust map image: skipping ...")
                continue

            # Count the number of negatives
            mask = image.masks["negatives"]
            if isinstance(mask, oldMask): mask = Mask(mask, wcs=image.wcs) # Fix type
            nnegatives = mask.relative_nmasked_in(self.central_ellipse)

            # Add to dictionary
            nnegatives_dict[factor] = nnegatives
            names_dict[factor] = name

        # Check
        if len(nnegatives_dict) == 0: raise RuntimeError("No appropriate hot dust maps were found")

        # Only one is found
        if len(nnegatives_dict) == 1:
            log.warning("Maps was found for only one factor")
            factor = names_dict.keys()[0]
            return method_name, names_dict[factor]

        # Sort
        nnegatives_dict = ordered_by_key(nnegatives_dict)

        # Find the factor
        factor = find_factor_max_nnegatives(nnegatives_dict, self.config.hot_dust_max_nnegatives)

        # Get the corresponding map name
        map_name = names_dict[factor]

        # Return
        return method_name, map_name

    # -----------------------------------------------------------------

    def auto_select_ionizing_map(self, hot_dust_name):

        """
        This function ...
        :param hot_dust_name:
        :return:
        """

        # Inform the user
        log.info("Automatically selecting appropriate ionizing stellar map ...")

        # Check
        if not self.has_ionizing_maps: raise IOError("No ionizing stellar maps are present")
        if self.ionizing_has_methods: raise IOError("Didn't expect different methods for the ionizing stellar maps")

        # No hot dust maps could be made
        if hot_dust_name is None:

            name = "halpha"
            if name not in self.ionizing_map_names_no_methods: raise RuntimeError("Could not find appropriate ionizing stellar map: no hot dust and no halpha map")
            else: return name

        # Hot dust map was found
        else:

            # Loop over the ionizing stellar map names
            map_name = None
            for name in self.ionizing_map_names_no_methods:

                if hot_dust_name in name:
                    map_name = name
                    break

            # Check
            if map_name is None: raise RuntimeError("Could not find the appropriate ionizing stellar map")

            # Return
            return map_name

    # -----------------------------------------------------------------

    def prompt_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for the significance levels ...")

        # Initialize the levels dictionary
        self.levels = dict()

        # Loop over the images (origins)
        for fltr in self.selection_origins:

            # Set filter name
            filter_name = str(fltr)
            id = filter_name.lower().replace(" ", "_")

            # Get the default level
            default_level = self.get_significance_level_for_filter(fltr)

            # Prompt for the level
            if not self.config.all_default_levels: level = prompt_real(id + "_level", "sigma level for the " + filter_name + " image", default=default_level)
            else: level = default_level

            # Set the level
            self.levels[fltr] = level

    # -----------------------------------------------------------------

    def check_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Checking the specified significance levels ...")

        # Check whether all filters are present
        for fltr in self.selection_origins:

            # Check
            if fltr not in self.levels: raise ValueError("The level for the '" + str(fltr) + "' filter is not specified")

            # Convert to float to make sure
            self.levels[fltr] = float(self.levels[fltr])

    # -----------------------------------------------------------------

    def has_old_steps(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        return fs.is_directory(map_path) and not fs.is_empty(map_path)

    # -----------------------------------------------------------------

    def has_young_steps(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        return fs.is_directory(map_path) and not fs.is_empty(map_path)

    # -----------------------------------------------------------------

    def has_ionizing_steps(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        return fs.is_directory(map_path) and not fs.is_empty(map_path)

    # -----------------------------------------------------------------

    def has_dust_steps(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        return fs.is_directory(map_path) and not fs.is_empty(map_path)

    # -----------------------------------------------------------------

    def has_old_step(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        filepath = fs.join(map_path, step + ".fits")
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_young_step(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        filepath = fs.join(map_path, step + ".fits")
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_ionizing_step(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        filepath = fs.join(map_path, step + ".fits")
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def has_dust_step(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        filepath = fs.join(map_path, step + ".fits")
        return fs.is_file(filepath)

    # -----------------------------------------------------------------

    def get_last_old_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        for step in reversed(steps):
            filepath = fs.join(self.old_steps_path_for_map(name), step + ".fits")
            if fs.is_file(filepath): return step, filepath

        #raise RuntimeError("We shouldn't get here")
        return None, None

    # -----------------------------------------------------------------

    def get_last_young_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        for step in reversed(steps):
            filepath = fs.join(self.young_steps_path_for_map(name), step + ".fits")
            if fs.is_file(filepath): return step, filepath

        #raise RuntimeError("We shouldn't get here")
        return None, None

    # -----------------------------------------------------------------

    def get_last_ionizing_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        for step in reversed(steps):
            filepath = fs.join(self.ionizing_steps_path_for_map(name), step + ".fits")
            if fs.is_file(filepath): return step, filepath

        #raise RuntimeError("We shouldn't get here")
        return None, None

    # -----------------------------------------------------------------

    def get_last_dust_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        for step in reversed(steps):
            filepath = fs.join(self.dust_steps_path_for_map(name), step + ".fits")
            if fs.is_file(filepath): return step, filepath

        #raise RuntimeError("We shouldn't get here")
        return None, None

    # -----------------------------------------------------------------

    @property
    def from_previous_levels(self):

        """
        This function ...
        :return:
        """

        return self.config.previous_levels is not None

    # -----------------------------------------------------------------

    def write_levels(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the levels ...")

        # Determine path for the levels file
        current_indices = fs.files_in_path(self.maps_components_path, extension="dat", returns="name", startswith="levels", convert=int, convert_split_pattern="_", convert_split_index=1)
        index = numbers.lowest_missing_integer(current_indices)
        levels_path = fs.join(self.maps_components_path, "levels_" + str(index) + ".dat")

        # Write levels dictionary
        write_dict(self.levels, levels_path)

    # -----------------------------------------------------------------

    @property
    def rerun(self):

        """
        This function ...
        :return:
        """

        return self.rerun_old or self.rerun_young or self.rerun_ionizing or self.rerun_dust

    # -----------------------------------------------------------------

    @property
    def rerun_all(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.rerun_all_old or self.config.rerun_all_young or self.config.rerun_all_ionizing or self.config.rerun_all_dust

    # -----------------------------------------------------------------

    @property
    def rerun_old(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_old is not None or self.config.rerun_all_old

    # -----------------------------------------------------------------

    @property
    def rerun_young(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_young is not None or self.config.rerun_all_young

    # -----------------------------------------------------------------

    @property
    def rerun_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_ionizing is not None or self.config.rerun_all_ionizing

    # -----------------------------------------------------------------

    @property
    def rerun_dust(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_dust is not None or self.config.rerun_all_dust

    # -----------------------------------------------------------------

    def set_rerun(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for rerunning certain processing steps ...")

        # Old
        if self.rerun_old: self.set_rerun_old()

        # Young
        if self.rerun_young: self.set_rerun_young()

        # Ionizing
        if self.rerun_ionizing: self.set_rerun_ionizing()

        # Dust
        if self.rerun_dust: self.set_rerun_dust()

    # -----------------------------------------------------------------

    def set_rerun_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting rerun for the old stellar maps ...")

        # Debugging
        if self.has_rerun_steps_old: log.debug("Rerun steps: " + tostr(self.rerun_steps_old))

        # Loop over the old stellar maps
        for name in self.old_selection:

            # Loop over the rerun steps
            for step in self.rerun_steps_old:

                # Debugging
                log.debug("Removing the intermediate output from the '" + step + "' step for the '" + name + "' old stellar map ...")

                # Determine the path for the intermediate map and potential mask
                map_path = self.old_step_path_for_map(name, step)
                mask_path = self.old_step_path_for_mask(name, step)

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(mask_path): fs.remove_file(mask_path)

            # Remove the end result
            if self.has_rerun_steps_old:

                # Debugging
                log.debug("Removing the end results for the '" + name + "' old stellar map ...")

                # Set the paths
                map_path = fs.join(self.old_component_maps_path, name + ".fits")
                clip_mask_path = fs.join(self.old_masks_path, name + "_" + clip_suffix + ".fits")
                softening_mask_path = fs.join(self.old_masks_path, name + "_" + softening_suffix + ".fits")
                deprojection_path = fs.join(self.old_deprojection_path, name + ".fits")
                deprojection_skirt_path = fs.join(self.old_deprojection_skirt_path, name + ".fits")
                edgeon_path = fs.join(self.old_edgeon_path, name + ".fits")

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(clip_mask_path): fs.remove_file(clip_mask_path)
                if fs.is_file(softening_mask_path): fs.remove_file(softening_mask_path)
                if fs.is_file(deprojection_path): fs.remove_file(deprojection_path)
                if fs.is_file(deprojection_skirt_path): fs.remove_file(deprojection_skirt_path)
                if fs.is_file(edgeon_path): fs.remove_file(edgeon_path)

            # Rerun all steps output
            if self.config.rerun_all_old:

                # Determine directory path
                path = self.old_steps_path_for_map(name)

                # Remove region files
                fs.remove_files_in_path(path, extension="reg")

                # Remove FITS files
                fs.remove_files_in_path(path, extension="fits")

                # Remove directories
                fs.remove_directories_in_path(path)

    # -----------------------------------------------------------------

    def set_rerun_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Setting rerun for the young stellar maps ...")

        # Debugging
        if self.has_rerun_steps_young: log.debug("Rerun steps: " + tostr(self.rerun_steps_young))

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Loop over the rerun steps
            for step in self.rerun_steps_young:

                # Debugging
                log.debug("Removing the intermediate output from the '" + step + "' step for the '" + name + "' young stellar map ...")

                # Determine the path for the intermediate map and potential mask
                map_path = self.young_step_path_for_map(name, step)
                mask_path = self.young_step_path_for_mask(name, step)

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(mask_path): fs.remove_file(mask_path)

            # Remove the end result
            if self.has_rerun_steps_young:

                # Debugging
                log.debug("Removing the end results for the '" + name + "' young stellar map ...")

                # Set the paths
                map_path = fs.join(self.young_component_maps_path, name + ".fits")
                clip_mask_path = fs.join(self.young_masks_path, name + "_" + clip_suffix + ".fits")
                softening_mask_path = fs.join(self.young_masks_path, name + "_" + softening_suffix + ".fits")
                deprojection_path = fs.join(self.young_deprojection_path, name + ".fits")
                deprojection_skirt_path = fs.join(self.young_deprojection_skirt_path, name + ".fits")
                edgeon_path = fs.join(self.young_edgeon_path, name + ".fits")

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(clip_mask_path): fs.remove_file(clip_mask_path)
                if fs.is_file(softening_mask_path): fs.remove_file(softening_mask_path)
                if fs.is_file(deprojection_path): fs.remove_file(deprojection_path)
                if fs.is_file(deprojection_skirt_path): fs.remove_file(deprojection_skirt_path)
                if fs.is_file(edgeon_path): fs.remove_file(edgeon_path)

            # Rerun all steps output
            if self.config.rerun_all_young:

                # Determine the directory path
                path = self.young_steps_path_for_map(name)

                # Remove region files
                fs.remove_files_in_path(path, extension="reg")

                # Remove FITS files
                fs.remove_files_in_path(path, extension="fits")

                # Remove directories
                fs.remove_directories_in_path(path)

    # -----------------------------------------------------------------

    def set_rerun_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting rerun for the ionizing stellar maps ...")

        # Debugging
        if self.has_rerun_steps_ionizing: log.debug("Rerun steps: " + tostr(self.rerun_steps_ionizing))

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Loop over the rerun steps
            for step in self.rerun_steps_ionizing:

                # Debugging
                log.debug("Removing the intermediate output from the '" + step + "' step for the '" + name + "' ionizing stellar map ...")

                # Determine the path for the intermediate map and potential mask
                map_path = self.ionizing_step_path_for_map(name, step)
                mask_path = self.ionizing_step_path_for_mask(name, step)

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(mask_path): fs.remove_file(mask_path)

            # Remove the end result
            if self.has_rerun_steps_ionizing:

                # Debugging
                log.debug("Removing the end results for the '" + name + "' ionizing stellar map ...")

                # Set the paths
                map_path = fs.join(self.ionizing_component_maps_path, name + ".fits")
                clip_mask_path = fs.join(self.ionizing_masks_path, name + "_" + clip_suffix + ".fits")
                softening_mask_path = fs.join(self.ionizing_masks_path, name + "_" + softening_suffix + ".fits")
                deprojection_path = fs.join(self.ionizing_deprojection_path, name + ".fits")
                deprojection_skirt_path = fs.join(self.ionizing_deprojection_skirt_path, name + ".fits")
                edgeon_path = fs.join(self.ionizing_edgeon_path, name + ".fits")

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(clip_mask_path): fs.remove_file(clip_mask_path)
                if fs.is_file(softening_mask_path): fs.remove_file(softening_mask_path)
                if fs.is_file(deprojection_path): fs.remove_file(deprojection_path)
                if fs.is_file(deprojection_skirt_path): fs.remove_file(deprojection_path)
                if fs.is_file(edgeon_path): fs.remove_file(edgeon_path)

            # Rerun all steps output
            if self.config.rerun_all_ionizing:

                # Determine the directory path
                path = self.ionizing_steps_path_for_map(name)

                # Remove region files
                fs.remove_files_in_path(path, extension="reg")

                # Remove FITS files
                fs.remove_files_in_path(path, extension="fits")

                # Remove directories
                fs.remove_directories_in_path(path)

    # -----------------------------------------------------------------

    def set_rerun_dust(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Setting rerun for the dust maps ...")

        # Debugging
        if self.has_rerun_steps_dust: log.debug("Rerun steps: " + tostr(self.rerun_steps_dust))

        # Loop over the dust maps
        for name in self.dust_selection:

            # Loop over the rerun steps
            for step in self.rerun_steps_dust:

                # Debugging
                log.debug("Removing the intermediate output from the '" + step + "' step for the '" + name + "' dust map ...")

                # Determine the path for the intermediate map and potential mask
                map_path = self.dust_step_path_for_map(name, step)
                mask_path = self.dust_step_path_for_mask(name, step)

                # Remove if present
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(mask_path): fs.remove_file(mask_path)

            # Remove the end result
            if self.has_rerun_steps_dust:

                # Debugging
                log.debug("Removing the end results for the '" + name + "' dust map ...")

                # Set the paths
                map_path = fs.join(self.dust_component_maps_path, name + ".fits")
                clip_mask_path = fs.join(self.dust_masks_path, name + "_" + clip_suffix + ".fits")
                softening_mask_path = fs.join(self.dust_masks_path, name + "_" + softening_suffix + ".fits")
                deprojection_path = fs.join(self.dust_deprojection_path, name + ".fits")
                deprojection_skirt_path = fs.join(self.dust_deprojection_skirt_path, name + ".fits")
                edgeon_path = fs.join(self.dust_edgeon_path, name + ".fits")

                # If present, remove
                if fs.is_file(map_path): fs.remove_file(map_path)
                if fs.is_file(clip_mask_path): fs.remove_file(clip_mask_path)
                if fs.is_file(softening_mask_path): fs.remove_file(softening_mask_path)
                if fs.is_file(deprojection_path): fs.remove_file(deprojection_path)
                if fs.is_file(deprojection_skirt_path): fs.remove_file(deprojection_skirt_path)
                if fs.is_file(edgeon_path): fs.remove_file(edgeon_path)

            # Rerun all steps output
            if self.config.rerun_all_dust:

                # Determine the directory path
                path = self.dust_steps_path_for_map(name)

                # Remove region files
                fs.remove_files_in_path(path, extension="reg")

                # Remove FITS files
                fs.remove_files_in_path(path, extension="fits")

                # Remove directories
                fs.remove_directories_in_path(path)

    # -----------------------------------------------------------------

    @property
    def redeproject(self):

        """
        Thsif unction ...
        :return:
        """

        return self.redeproject_old or self.redeproject_young or self.redeproject_ionizing or self.redeproject_dust

    # -----------------------------------------------------------------

    @property
    def redeproject_old(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_old

    # -----------------------------------------------------------------

    @property
    def redeproject_young(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_young

    # -----------------------------------------------------------------

    @property
    def redeproject_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_ionizing

    # -----------------------------------------------------------------

    @property
    def redeproject_dust(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.redeproject_dust

    # -----------------------------------------------------------------

    def set_redeproject(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting maps ...")

        # Old
        if self.redeproject_old: self.set_redeproject_old()

        # Young
        if self.redeproject_young: self.set_redeproject_young()

        # Ionizing
        if self.redeproject_ionizing: self.set_redeproject_ionizing()

        # Dust
        if self.redeproject_dust: self.set_redeproject_dust()

    # -----------------------------------------------------------------

    def set_redeproject_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting old stellar maps ...")

        # Loop over the old stellar maps
        for name in self.old_selection:

            # Determine file path
            map_path = fs.join(self.old_deprojection_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.old_deprojection_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting young stellar maps ...")

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Determine file path
            map_path = fs.join(self.young_deprojection_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.young_deprojection_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting ionizing stellar maps ...")

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Determine file path
            map_path = fs.join(self.ionizing_deprojection_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.ionizing_deprojection_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_dust(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting dust stellar maps ...")

        # Loop over the dust maps
        for name in self.dust_selection:

            # Determine file path
            map_path = fs.join(self.dust_deprojection_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.dust_deprojection_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    @property
    def redeproject_skirt(self):

        """
        This function ...
        :return:
        """

        return self.redeproject_skirt_old or self.redeproject_skirt_young or self.redeproject_skirt_ionizing or self.redeproject_skirt_dust

    # -----------------------------------------------------------------

    @property
    def redeproject_skirt_old(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_skirt_old

    # -----------------------------------------------------------------

    @property
    def redeproject_skirt_young(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_skirt_young

    # -----------------------------------------------------------------

    @property
    def redeproject_skirt_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_skirt_ionizing

    # -----------------------------------------------------------------

    @property
    def redeproject_skirt_dust(self):

        """
        This function ...
        :return:
        """

        return self.config.redeproject_skirt_dust

    # -----------------------------------------------------------------

    def set_redeproject_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting maps with SKIRT ...")

        # Old
        if self.redeproject_skirt_old: self.set_redeproject_skirt_old()

        # Young
        if self.redeproject_skirt_young: self.set_redeproject_skirt_young()

        # Ionizing
        if self.redeproject_skirt_ionizing: self.set_redeproject_skirt_ionizing()

        # Dust
        if self.redeproject_skirt_dust: self.set_redeproject_skirt_dust()

    # -----------------------------------------------------------------

    def set_redeproject_skirt_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting old stellar maps with SKIRT ....")

        # Loop over the old stellar maps
        for name in self.old_selection:

            # Determine file path
            map_path = fs.join(self.old_deprojection_skirt_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.old_deprojection_skirt_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_skirt_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting young stellar maps with SKIRT ...")

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Determine file path
            map_path = fs.join(self.young_deprojection_skirt_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.young_deprojection_skirt_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_skirt_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting ionizing stellar maps with SKIRT ...")

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Determine file path
            map_path = fs.join(self.ionizing_deprojection_skirt_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.ionizing_deprojection_skirt_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_redeproject_skirt_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for redeprojecting dust maps with SKIRT ...")

        # Loop over the dust maps
        for name in self.dust_selection:

            # Determine file path
            map_path = fs.join(self.dust_deprojection_skirt_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.dust_deprojection_skirt_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    @property
    def reproject(self):

        """
        Thsif unction ...
        :return:
        """

        return self.reproject_old or self.reproject_young or self.reproject_ionizing or self.reproject_dust

    # -----------------------------------------------------------------

    @property
    def reproject_old(self):

        """
        This function ...
        :return:
        """

        return self.config.reproject_old

    # -----------------------------------------------------------------

    @property
    def reproject_young(self):

        """
        This function ...
        :return:
        """

        return self.config.reproject_young

    # -----------------------------------------------------------------

    @property
    def reproject_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.reproject_ionizing

    # -----------------------------------------------------------------

    @property
    def reproject_dust(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.reproject_dust

    # -----------------------------------------------------------------

    def set_reproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for reprojecting maps ...")

        # Old stars
        if self.reproject_old: self.set_reproject_old()

        # Young stars
        if self.reproject_young: self.set_reproject_young()

        # Ionizing stars
        if self.reproject_ionizing: self.set_reproject_ionizing()

        # Dust
        if self.reproject_dust: self.set_reproject_dust()

    # -----------------------------------------------------------------

    def set_reproject_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for reprojecting old stellar maps ...")

        # Loop over the old stellar maps
        for name in self.old_selection:

            # Determine file path
            map_path = fs.join(self.old_edgeon_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.old_edgeon_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_reproject_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Settings things in order for reprojecting young stellar maps ...")

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Determine file path
            map_path = fs.join(self.young_edgeon_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.young_edgeon_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_reproject_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for reprojecting ionizing stellar maps ...")

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Determine file path
            map_path = fs.join(self.ionizing_edgeon_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.ionizing_edgeon_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.clear_directory(directory_path)

    # -----------------------------------------------------------------

    def set_reproject_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting things in order for reprojecting dust maps ...")

        # Loop overt he dust maps
        for name in self.dust_selection:

            # Determine file path
            map_path = fs.join(self.dust_edgeon_path, name + ".fits")

            # Determine directory path
            directory_path = fs.join(self.dust_edgeon_path, name)

            # Remove if present
            if fs.is_file(map_path): fs.remove_file(map_path)
            if fs.is_directory(directory_path): fs.remove_directory(directory_path)

    # -----------------------------------------------------------------

    @property
    def remove(self):

        """
        Thisf unction ...
        :return:
        """

        return self.remove_old or self.remove_young or self.remove_ionizing or self.remove_dust

    # -----------------------------------------------------------------

    @property
    def remove_old(self):

        """
        Thisfunction ...
        :return:
        """

        return self.config.remove_other_old

    # -----------------------------------------------------------------

    @property
    def remove_young(self):

        """
        This function ...
        :return:
        """

        return self.config.remove_other_young

    # -----------------------------------------------------------------

    @property
    def remove_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.remove_other_ionizing

    # -----------------------------------------------------------------

    @property
    def remove_dust(self):

        """
        Thisn function ...
        :return:
        """

        return self.config.remove_other_dust

    # -----------------------------------------------------------------

    def remove_other(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing maps other than those that are selected ...")

        # Old
        if self.remove_old: self.remove_other_old()

        # Young
        if self.remove_young: self.remove_other_young()

        # Ionizing
        if self.remove_ionizing: self.remove_other_ionizing()

        # Dust
        if self.remove_dust: self.remove_other_dust()

    # -----------------------------------------------------------------

    def remove_other_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing other old stellar maps ...")

        # Loop over the FITS files in the directory
        for filepath, filename in fs.files_in_path(self.old_component_maps_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.old_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the masks directory
        for filepath, filename in fs.files_in_path(self.old_masks_path, extension="fits", returns=["path", "name"]):

            # Determine map name
            if filename.endswith(clip_suffix): map_name = filename.split(clip_suffix)[0]
            elif filename.endswith(softening_suffix): map_name = filename.split(softening_suffix)[0]
            else: raise ValueError("Unrecognized file: '" + filename)

            # Check whether map name in selection
            if map_name in self.old_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the directories in the steps directory
        for dirpath, dirname in fs.directories_in_path(self.old_steps_path, returns=["path", "name"]):

            # Check whether filename is in the selection
            if dirname in self.old_selection: continue

            # otherwise, remove it
            fs.remove_directory(dirpath)

        # Loop over the FITS files in the deprojected directory
        for filepath, filename in fs.files_in_path(self.old_deprojection_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.old_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the deprojected_skirt directory
        for filepath, filename in fs.files_in_path(self.old_deprojection_skirt_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.old_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the edgeon directory
        for filepath, filename in fs.files_in_path(self.old_edgeon_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.old_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

    # -----------------------------------------------------------------

    def remove_other_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing other young stellar maps ...")

        # Loop over the FITS files in the directory
        for filepath, filename in fs.files_in_path(self.young_component_maps_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.young_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the masks directory
        for filepath, filename in fs.files_in_path(self.young_masks_path, extension="fits", returns=["path", "name"]):

            # Determine map name
            if filename.endswith(clip_suffix): map_name = filename.split(clip_suffix)[0]
            elif filename.endswith(softening_suffix): map_name = filename.split(softening_suffix)[0]
            else: raise ValueError("Unrecognized file: '" + filename)

            # Check whether map name in selection
            if map_name in self.young_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the directories in the steps directory
        for dirpath, dirname in fs.directories_in_path(self.young_steps_path, returns=["path", "name"]):

            # Check whether filename is in the selection
            if dirname in self.young_selection: continue

            # otherwise, remove it
            fs.remove_directory(dirpath)

        # Loop over the FITS files in the deprojected directory
        for filepath, filename in fs.files_in_path(self.young_deprojection_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.young_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the deprojected_skirt directory
        for filepath, filename in fs.files_in_path(self.young_deprojection_skirt_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.young_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the edgeon directory
        for filepath, filename in fs.files_in_path(self.young_edgeon_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.young_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

    #  -----------------------------------------------------------------

    def remove_other_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing other ionizing stellar maps ...")

        # Loop over the FITS files in the directory
        for filepath, filename in fs.files_in_path(self.ionizing_component_maps_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.ionizing_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the masks directory
        for filepath, filename in fs.files_in_path(self.ionizing_masks_path, extension="fits", returns=["path", "name"]):

            # Determine map name
            if filename.endswith(clip_suffix): map_name = filename.split(clip_suffix)[0]
            elif filename.endswith(softening_suffix): map_name = filename.split(softening_suffix)[0]
            else: raise ValueError("Unrecognized file: '" + filename)

            # Check whether map name in selection
            if map_name in self.ionizing_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the directories in the steps directory
        for dirpath, dirname in fs.directories_in_path(self.ionizing_steps_path, returns=["path", "name"]):

            # Check whether filename is in the selection
            if dirname in self.ionizing_selection: continue

            # otherwise, remove it
            fs.remove_directory(dirpath)

        # Loop over the FITS files in the deprojected directory
        for filepath, filename in fs.files_in_path(self.ionizing_deprojection_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.ionizing_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the deprojected_skirt directory
        for filepath, filename in fs.files_in_path(self.ionizing_deprojection_skirt_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.ionizing_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the edgeon directory
        for filepath, filename in fs.files_in_path(self.ionizing_edgeon_path, extension="fits", returns=["path", "name"]):

            # Check wether the filename is in the selection
            if filename in self.ionizing_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

    # -----------------------------------------------------------------

    def remove_other_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Removing other young stellar maps ...")

        # Loop over the FITS files in the directory
        for filepath, filename in fs.files_in_path(self.dust_component_maps_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.dust_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the masks directory
        for filepath, filename in fs.files_in_path(self.dust_masks_path, extension="fits", returns=["path", "name"]):

            # Determine map name
            if filename.endswith(clip_suffix): map_name = filename.split(clip_suffix)[0]
            elif filename.endswith(softening_suffix): map_name = filename.split(softening_suffix)[0]
            else: raise ValueError("Unrecognized file: '" + filename)

            # Check whether map name in selection
            if map_name in self.dust_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the directories in the steps directory
        for dirpath, dirname in fs.directories_in_path(self.dust_steps_path, returns=["path", "name"]):

            # Check whether filename is in the selection
            if dirname in self.dust_selection: continue

            # otherwise, remove it
            fs.remove_directory(dirpath)

        # Loop over the FITS files in the deprojected directory
        for filepath, filename in fs.files_in_path(self.dust_deprojection_path, extension="fits", returns=["path", "name"]):

            # Check whether the filename is in the selection
            if filename in self.dust_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the deprojected_skirt directory
        for filepath, filename in fs.files_in_path(self.dust_deprojection_skirt_path, extension="fits", returns=["path", "name"]):

            # Check whether filename is in the selection
            if filename in self.dust_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

        # Loop over the FITS files in the edgeon directory
        for filepath, filename in fs.files_in_path(self.dust_edgeon_path, extension="fits", returns=["path", "name"]):

            # Check whether filename is in the selection
            if filename in self.dust_selection: continue

            # Otherwise, remove it
            fs.remove_file(filepath)

    # -----------------------------------------------------------------

    def load_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected maps ...")

        # Load the old stellar maps
        self.load_old()

        # Load the young stellar maps
        self.load_young()

        # Load the ionizing stellar maps
        self.load_ionizing()

        # Load the dust maps
        self.load_dust()

    # -----------------------------------------------------------------

    def load_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected old stellar maps ...")

        # Load
        for name in self.old_selection:

            # Debugging
            log.debug("Loading the '" + name + "' old stellar map ...")

            # Check whether intermediate result is present
            if self.has_old_steps(name):

                # Determine the latest step
                step, path = self.get_last_old_step(name)

                # Step is found
                if step is not None:

                    # Inform
                    log.success("Found a " + step + " '" + name + "' map")

                    # Load
                    self.old_maps[name] = Frame.from_file(path)

                    # Set metadata
                    for step in steps_before_and_including(step): self.old_maps[name].metadata[step] = True
                    for step in steps_after(step): self.old_maps[name].metadata[step] = False

                # Load raw map
                else: self.load_raw_old_map(name)

            # Load raw map
            else: self.load_raw_old_map(name)

            # Determine path to cropped negatives file
            cropped_negatives_path = self.old_extra_path_for_map(name, cropped_negatives_filename)

            # Cropped negatives mask has been created
            if self.is_cropped_old(name) and fs.is_file(cropped_negatives_path): self.old_negative_masks[name] = Mask.from_file(cropped_negatives_path)

            # Not yet cropped
            else:

                # Get negatives plane
                if negatives_plane_name in get_plane_names(self.old_map_paths[name]):
                    negatives = Mask.from_file(self.old_map_paths[name], plane=negatives_plane_name)
                    self.old_negative_masks[name] = negatives

    # -----------------------------------------------------------------

    def load_raw_old_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Load
        self.old_maps[name] = Frame.from_file(self.old_map_paths[name])

        # Set metadata
        for step in steps: self.old_maps[name].metadata[step] = False

    # -----------------------------------------------------------------

    def load_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected young stellar maps ...")

        # Load
        for name in self.young_selection:

            # Debugging
            log.debug("Loading the '" + name + "' young stellar map ...")

            # Check whether intermediate result is present
            if self.has_young_steps(name):

                # Determine the last step
                step, path = self.get_last_young_step(name)

                # Step is found
                if step is not None:

                    # Inform
                    log.success("Found a " + step + " '" + name + "' map")

                    # Load
                    self.young_maps[name] = Frame.from_file(path)

                    # Set metadata
                    for step in steps_before_and_including(step): self.young_maps[name].metadata[step] = True
                    for step in steps_after(step): self.young_maps[name].metadata[step] = False

                # Load raw map
                else: self.load_raw_young_map(name)

            # Load raw map
            else: self.load_raw_young_map(name)

            # Determine the path to the cropped negatives file
            cropped_negatives_path = self.young_extra_path_for_map(name, cropped_negatives_filename)

            # Cropped mask has been created
            if self.is_cropped_young(name) and fs.is_file(cropped_negatives_path): self.young_negative_masks[name] = Mask.from_file(cropped_negatives_path)

            # Not cropped yet
            else:

                # Get negatives plane
                if negatives_plane_name in get_plane_names(self.young_map_paths[name]):
                    negatives = Mask.from_file(self.young_map_paths[name], plane=negatives_plane_name)
                    self.young_negative_masks[name] = negatives

    # -----------------------------------------------------------------

    def load_raw_young_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Load
        self.young_maps[name] = Frame.from_file(self.young_map_paths[name])

        # Set metadata
        for step in steps: self.young_maps[name].metadata[step] = False

    # -----------------------------------------------------------------

    def load_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected ionizing stellar maps ...")

        # Load
        for name in self.ionizing_selection:

            # Debugging
            log.debug("Loading the '" + name + "' ionizing stellar map ...")

            # Check whether intermediate result is present
            if self.has_ionizing_steps(name):

                # Determine the last step
                step, path = self.get_last_ionizing_step(name)

                # Step is found
                if step is not None:

                    # Inform
                    log.success("Found a " + step + " '" + name + "' map")

                    # Load
                    self.ionizing_maps[name] = Frame.from_file(path)

                    # Set metadata
                    for step in steps_before_and_including(step): self.ionizing_maps[name].metadata[step] = True
                    for step in steps_after(step): self.ionizing_maps[name].metadata[step] = False

                # Load raw map
                else: self.load_raw_ionizing_map(name)

            # Load raw map
            else: self.load_raw_ionizing_map(name)

            # Get negatives plane
            #if negatives_plane_name in get_plane_names(self.ionizing_map_paths[name]):
            #    negatives = Mask.from_file(self.ionizing_map_paths[name], plane=negatives_plane_name)
            #    self.ionizing_negative_masks[name] = negatives

            # Determine the path to the cropped negatives file
            cropped_negatives_path = self.ionizing_extra_path_for_map(name, cropped_negatives_filename)

            # Cropped negatives mask has been created
            if self.is_cropped_ionizing(name) and fs.is_file(cropped_negatives_path): self.ionizing_negative_masks[name] = Mask.from_file(cropped_negatives_path)

            # Not yet cropped
            else:

                # Get negatives plane of the HOT DUST (not the negatives of the H-alpha map)
                hot_negatives_plane_name = "hot_negatives"
                if hot_negatives_plane_name in get_plane_names(self.ionizing_map_paths[name]):
                    negatives = Mask.from_file(self.ionizing_map_paths[name], plane=hot_negatives_plane_name)
                    self.ionizing_negative_masks[name] = negatives

    # -----------------------------------------------------------------

    def load_raw_ionizing_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Load
        self.ionizing_maps[name] = Frame.from_file(self.ionizing_map_paths[name])

        # Set metadata
        for step in steps: self.ionizing_maps[name].metadata[step] = False

    # -----------------------------------------------------------------

    def load_dust(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the selected dust maps ...")

        # Load
        for name in self.dust_selection:

            # Debugging
            log.debug("Loading the '" + name + "' dust map ...")

            # Check whether intermediate result is present
            if self.has_dust_steps(name):

                # Determine the last step
                step, path = self.get_last_dust_step(name)

                # Step is found
                if step is not None:

                    # Inform
                    log.success("Found a " + step + " '" + name + "' map")

                    # Load
                    self.dust_maps[name] = Frame.from_file(path)

                    # Set metadata
                    for step in steps_before_and_including(step): self.dust_maps[name].metadata[step] = True
                    for step in steps_after(step): self.dust_maps[name].metadata[step] = False

                # Load raw map
                else: self.load_raw_dust_map(name)

            # Load raw map
            else: self.load_raw_dust_map(name)

            # Determine the path to the cropped negatives file
            cropped_negatives_path = self.dust_extra_path_for_map(name, cropped_negatives_filename)

            # Cropped negatives mask has already been created
            if self.is_cropped_dust(name) and fs.is_file(cropped_negatives_path): self.dust_negative_masks[name] = Mask.from_file(cropped_negatives_path)

            # Not cropped yet
            else:

                # Get negatives plane
                if negatives_plane_name in get_plane_names(self.dust_map_paths[name]):
                    negatives = Mask.from_file(self.dust_map_paths[name], plane=negatives_plane_name)
                    self.dust_negative_masks[name] = negatives

    # -----------------------------------------------------------------

    def load_raw_dust_map(self, name):

        """
        Thisf ucntion ...
        :param name:
        :return:
        """

        # Load
        self.dust_maps[name] = Frame.from_file(self.dust_map_paths[name])

        # Set metadata
        for step in steps: self.dust_maps[name].metadata[step] = False

    # -----------------------------------------------------------------

    def old_clipping_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Determine the clipping directory path
        clipping_path = fs.join(map_path, "clipping")
        if create and not fs.is_directory(clipping_path): fs.create_directory(clipping_path)

        # Return the clipping path
        return clipping_path

    # -----------------------------------------------------------------

    def old_extra_directory_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def old_extra_path_for_map(self, name, filename):

        """
        This function ...
        :param name:
        :param filename:
        :return:
        """

        # Get directory path
        map_path = self.old_extra_directory_path_for_map(name)

        # Return the filepath
        return fs.join(map_path, filename)

    # -----------------------------------------------------------------

    def old_steps_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def old_step_path_for_map(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        # Determine the directory path
        map_path = self.old_steps_path_for_map(name)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the map path
        return fs.join(map_path, step + ".fits")

    # -----------------------------------------------------------------

    def old_step_path_for_mask(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.old_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the MASK path
        return fs.join(map_path, step + "_mask.fits")

    # -----------------------------------------------------------------

    def young_clipping_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Determine the clipping directory path
        clipping_path = fs.join(map_path, "clipping")
        if create and not fs.is_directory(clipping_path): fs.create_directory(clipping_path)

        # Return the clipping path
        return clipping_path

    # -----------------------------------------------------------------

    def young_extra_directory_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def young_extra_path_for_map(self, name, filename):

        """
        This function ...
        :param name:
        :param filename:
        :return:
        """

        # Get the directory path
        map_path = self.young_extra_directory_path_for_map(name)

        # Return the filepath
        return fs.join(map_path, filename)

    # -----------------------------------------------------------------

    def young_steps_path_for_map(self, name, create=True):

        """
        Thisf unction ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def young_step_path_for_map(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        # Determine the directory path
        map_path = self.young_steps_path_for_map(name)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the map path
        return fs.join(map_path, step + ".fits")

    # -----------------------------------------------------------------

    def young_step_path_for_mask(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the MASK path
        return fs.join(map_path, step + "_mask.fits")

    # -----------------------------------------------------------------

    def ionizing_clipping_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Determine the clipping directory path
        clipping_path = fs.join(map_path, "clipping")
        if create and not fs.is_directory(clipping_path): fs.create_directory(clipping_path)

        # Return the clipping path
        return clipping_path

    # -----------------------------------------------------------------

    def ionizing_extra_directory_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def ionizing_extra_path_for_map(self, name, filename):

        """
        This function ...
        :param name:
        :param filename:
        :return:
        """

        # Get the directory path
        map_path = self.ionizing_extra_directory_path_for_map(name)

        # Return the filepath
        return fs.join(map_path, filename)

    # -----------------------------------------------------------------

    def ionizing_steps_path_for_map(self, name, create=True):

        """
        Thisn function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def ionizing_step_path_for_map(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        # Determine the directory path
        map_path = self.ionizing_steps_path_for_map(name)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the map path
        return fs.join(map_path, step + ".fits")

    # -----------------------------------------------------------------

    def ionizing_step_path_for_mask(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the MASK path
        return fs.join(map_path, step + "_mask.fits")

    # -----------------------------------------------------------------

    def dust_clipping_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Determine the clipping directory path
        clipping_path = fs.join(map_path, "clipping")
        if create and not fs.is_directory(clipping_path): fs.create_directory(clipping_path)

        # Return the clipping path
        return clipping_path

    # -----------------------------------------------------------------

    def dust_extra_directory_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def dust_extra_path_for_map(self, name, filename):

        """
        This function ...
        :param name:
        :param filename:
        :return:
        """

        # Get the directory path
        map_path = self.dust_extra_directory_path_for_map(name)

        # Return the filepath
        return fs.join(map_path, filename)

    # -----------------------------------------------------------------

    def dust_steps_path_for_map(self, name, create=True):

        """
        This function ...
        :param name:
        :param create:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        if create and not fs.is_directory(map_path): fs.create_directory(map_path)
        return map_path

    # -----------------------------------------------------------------

    def dust_step_path_for_map(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        # Determine the directory path
        map_path = self.dust_steps_path_for_map(name)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the map path
        return fs.join(map_path, step + ".fits")

    # -----------------------------------------------------------------

    def dust_step_path_for_mask(self, name, step):

        """
        This function ...
        :param name:
        :param step:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        if not fs.is_directory(map_path): fs.create_directory(map_path)

        # Check the step
        if step not in steps: raise ValueError("Invalid step: '" + step + "'")

        # Return the MASK path
        return fs.join(map_path, step + "_mask.fits")

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Old
        self.load_old_masks()

        # Young
        self.load_young_masks()

        # Ionizing
        self.load_ionizing_masks()

        # Dust
        self.load_dust_masks()

    # -----------------------------------------------------------------

    def load_old_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading clipping and softening masks for the old stellar maps ...")

        # Loop over the maps
        for name in self.old_selection:

            # Debugging
            log.debug("Loading the '" + name + "' old stellar map masks ...")

            # Check whether the clipping step output is present
            if self.has_old_step(name, clip_step):

                # Inform
                log.success("Found a " + clip_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.old_step_path_for_mask(name, clip_step)

                # Load the mask
                self.old_clip_masks[name] = load_mask_or_alpha_mask(path)

            # Check whether the softening step output is present
            if self.has_old_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.old_step_path_for_mask(name, softened_step)

                # Load the mask
                self.old_softening_masks[name] = load_mask_or_alpha_mask(path)

    # -----------------------------------------------------------------

    def load_young_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading clipping and softening masks for the young stellar maps ...")

        # Loop over the maps
        for name in self.young_selection:

            # Debugging
            log.debug("Loading the '" + name + "' young stellar map masks ...")

            # Check whether the clipping step output is present
            if self.has_young_step(name, clip_step):

                # Inform
                log.success("Found a " + clip_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.young_step_path_for_mask(name, clip_step)

                # Load the mask
                self.young_clip_masks[name] = load_mask_or_alpha_mask(path)

            # Check whether the softening step output is present
            if self.has_young_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.young_step_path_for_mask(name, softened_step)

                # Load the mask
                self.young_softening_masks[name] = load_mask_or_alpha_mask(path)

    # -----------------------------------------------------------------

    def load_ionizing_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading clipping and softening masks for the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_selection:

            # Debugging
            log.debug("Loading the '" + name + "' ionizing stellar map masks ...")

            # Check whether the clipping step output is present
            if self.has_ionizing_step(name, clip_step):

                # Inform
                log.success("Found a " + clip_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.ionizing_step_path_for_mask(name, clip_step)

                # Load the mask
                self.ionizing_clip_masks[name] = load_mask_or_alpha_mask(path)

            # Check whether the softening step output is present
            if self.has_ionizing_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.ionizing_step_path_for_mask(name, softened_step)

                # Load the mask
                self.ionizing_softening_masks[name] = load_mask_or_alpha_mask(path)

    # -----------------------------------------------------------------

    def load_dust_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading clipping and softening masks for the dust maps ...")

        # Loop over the maps
        for name in self.dust_selection:

            # Debugging
            log.debug("Loading the '" + name + "' dust map masks ...")

            # Check whether the clipping step output is present
            if self.has_dust_step(name, clip_step):

                # Inform
                log.success("Found a " + clip_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.dust_step_path_for_mask(name, clip_step)

                # Load the mask
                self.dust_clip_masks[name] = load_mask_or_alpha_mask(path)

            # Check whether the softening step output is present
            if self.has_dust_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.dust_step_path_for_mask(name, softened_step)

                # Load the mask
                self.dust_softening_masks[name] = load_mask_or_alpha_mask(path)

    # -----------------------------------------------------------------

    @property
    def stop_after_correction(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.stop_after == correct_step

    # -----------------------------------------------------------------

    @property
    def stop_after_cropping(self):

        """
        This function ...
        :return:
        """

        return self.config.stop_after == crop_step

    # -----------------------------------------------------------------

    @property
    def stop_after_interpolation_negatives(self):

        """
        This function ...
        :return:
        """

        return self.config.stop_after == interpolate_negatives_step

    # -----------------------------------------------------------------

    @property
    def stop_after_interpolation(self):

        """
        This function ...
        :return:
        """

        return self.config.stop_after == interpolate_step

    # -----------------------------------------------------------------

    @property
    def stop_after_truncation(self):

        """
        This function ...
        :return:
        """

        return self.config.stop_after == truncate_step

    # -----------------------------------------------------------------

    @property
    def stop_after_clipping(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.stop_after == clip_step

    # -----------------------------------------------------------------

    @property
    def stop_after_softening(self):

        """
        Thisf unction ...
        :return:
        """

        return self.config.stop_after == softened_step

    # -----------------------------------------------------------------

    def write_old_info(self, name, step, **info):

        """
        This function ...
        :param name:
        :param step:
        :param info:
        :return:
        """

        # Determine the filename
        filename = step + ".dat"

        # Determine the path
        path = self.old_extra_path_for_map(name, filename)

        # Write the info
        write_dict(info, path)

    # -----------------------------------------------------------------

    def write_young_info(self, name, step, **info):

        """
        Thisn function ...
        :param step:
        :param info:
        :return:
        """

        # Determine the filename
        filename = step + ".dat"

        # Determine the path
        path = self.young_extra_path_for_map(name, filename)

        # Write the info
        write_dict(info, path)

    # -----------------------------------------------------------------

    def write_ionizing_info(self, name, step, **info):

        """
        This function ...
        :param name:
        :param step:
        :param info:
        :return:
        """

        # Determine the filename
        filename = step + ".dat"

        # Determine the path
        path = self.ionizing_extra_path_for_map(name, filename)

        # Write the info
        write_dict(info, path)

    # -----------------------------------------------------------------

    def write_dust_info(self, name, step, **info):

        """
        This function ...
        :param name:
        :param step:
        :param info:
        :return:
        """

        # Determine the filename
        filename = step + ".dat"

        # Determine the path
        path = self.dust_extra_path_for_map(name, filename)

        # Write the info
        write_dict(info, path)

    # -----------------------------------------------------------------

    def process_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the maps ...")

        # 1. Apply corrections
        self.correction()

        # 2. Apply cropping
        self.cropping()

        # 3. Apply interpolation of negatives
        self.interpolation_negatives()

        # 4. Apply interpolation
        self.interpolation()

        # 5. Truncation
        self.truncation()

        # 6. Clipping
        self.clipping()

        # 7. Softening
        self.softening()

        # Stop?
        if self.config.stop_after_all: exit()

    # -----------------------------------------------------------------

    def correction(self):

        """
        Thisf unction ...
        :return:
        """

        # Correct
        if self.config.correct: self.correct_maps()

        # Flag
        else: self.flag_corrected()

        # Stop?
        if self.stop_after_correction: exit()

    # -----------------------------------------------------------------

    def cropping(self):

        """
        This function ...
        :return:
        """

        # Crop
        if self.config.crop: self.crop_maps()

        # Flag
        else: self.flag_cropped()

        # Stop?
        if self.stop_after_cropping: exit()

    # -----------------------------------------------------------------

    def interpolation_negatives(self):

        """
        This function ..
        :return:
        """

        # Interpolate negatives
        if self.config.interpolate_negatives: self.interpolate_negatives_maps()

        # Flag
        else: self.flag_interpolated_negatives()

        # Stop?
        if self.stop_after_interpolation_negatives: exit()

    # -----------------------------------------------------------------

    def interpolation(self):

        """
        This function ...
        :return:
        """

        # Interpolate
        if self.config.interpolate: self.interpolate_maps()

        # Flag
        else: self.flag_interpolated()

        # Stop?
        if self.stop_after_interpolation: exit()

    # -----------------------------------------------------------------

    def truncation(self):

        """
        This function ...
        :return:
        """

        # Truncate
        if self.config.truncate: self.truncate_maps()

        # Flag
        else: self.flag_truncated()

        # Stop?
        if self.stop_after_truncation: exit()

    # -----------------------------------------------------------------

    def clipping(self):

        """
        This function ...
        :return:
        """

        # Clip
        if self.config.clip: self.clip_maps()

        # Flag
        else: self.flag_clipped()

        # Stop?
        if self.stop_after_clipping: exit()

    # -----------------------------------------------------------------

    def softening(self):

        """
        Thisf unction ...
        :return:
        """

        # Soften the edges
        if self.config.soften: self.soften_edges()

        # Flag
        else: self.flag_softened()

        # Stop?
        if self.stop_after_softening: exit()

    # -----------------------------------------------------------------

    def correct_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the maps ...")

        # Old
        self.correct_old_maps()

        # Young
        self.correct_young_maps()

        # Ionizing
        self.correct_ionizing_maps()

        # Dust
        self.correct_dust_maps()

    # -----------------------------------------------------------------

    def flag_corrected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging the maps as being corrected ...")

        # Old
        self.flag_corrected_old()

        # Young
        self.flag_corrected_young()

        # Ionizing
        self.flag_corrected_ionizing()

        # Dust
        self.flag_corrected_dust()

    # -----------------------------------------------------------------

    def flag_corrected_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[correct_step] = True

    # -----------------------------------------------------------------

    def flag_corrected_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[correct_step] = True

    # -----------------------------------------------------------------

    def flag_corrected_ionizing(self):

        """
        This fuction ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[correct_step] = True

    # -----------------------------------------------------------------

    def flag_corrected_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[correct_step] = True

    # -----------------------------------------------------------------

    def is_corrected_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[correct_step]

    # -----------------------------------------------------------------

    def is_corrected_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[correct_step]

    # -----------------------------------------------------------------

    def is_corrected_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[correct_step]

    # -----------------------------------------------------------------

    def is_corrected_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[correct_step]

    # -----------------------------------------------------------------

    def correct_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_corrected_old(name):
                log.success("The '" + name + "' old stellar map is already corrected")
                continue

            # Debugging
            log.debug("Correcting the '" + name + "' old stellar map ...")

            # Correct the map
            ninfs, nnans, nnegatives = self.correct_map(self.old_maps[name])

            # Set flag
            self.old_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, correct_step))

            # Save info
            if self.config.steps: self.write_old_info(name, correct_step, ninfs=ninfs, nnans=nnans, nnegatives=nnegatives)

    # -----------------------------------------------------------------

    def correct_young_maps(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_corrected_young(name):
                log.success("The '" + name + "' young stellar map is already corrected")
                continue

            # Debugging
            log.debug("Correcting the '" + name + "' young stellar map ...")

            # Correct the map
            ninfs, nnans, nnegatives = self.correct_map(self.young_maps[name])

            # Set flag
            self.young_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, correct_step))

            # Save info
            if self.config.steps: self.write_young_info(name, correct_step, ninfs=ninfs, nnans=nnans, nnegatives=nnegatives)

    # -----------------------------------------------------------------

    def correct_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_corrected_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already corrected")
                continue

            # Debugging
            log.debug("Correcting the '" + name + "' ionizing stellar map ...")

            # Correct the map
            ninfs, nnans, nnegatives = self.correct_map(self.ionizing_maps[name])

            # Set flag
            self.ionizing_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, correct_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, correct_step, ninfs=ninfs, nnans=nnans, nnegatives=nnegatives)

    # -----------------------------------------------------------------

    def correct_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_corrected_dust(name):
                log.success("The '" + name + "' dust map is already corrected")
                continue

            # Debugging
            log.debug("Correcting the '" + name + "' dust map ...")

            # Correct the map
            ninfs, nnans, nnegatives = self.correct_map(self.dust_maps[name])

            # Set flag
            self.dust_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.step: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, correct_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, correct_step, ninfs=ninfs, nnans=nnans, nnegatives=nnegatives)

    # -----------------------------------------------------------------

    @lazyproperty
    def negatives_central_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.negatives_central_ellipse_factor

    # -----------------------------------------------------------------

    def interpolate_negatives_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the negatives in the maps ...")

        # Old
        if self.config.interpolate_old_negatives: self.interpolate_negatives_old_maps()
        else: self.flag_interpolated_negatives_old()

        # Young
        if self.config.interpolate_young_negatives: self.interpolate_negatives_young_maps()
        else: self.flag_interpolated_negatives_young()

        # Ionizing
        if self.config.interpolate_ionizing_negatives: self.interpolate_negatives_ionizing_maps()
        else: self.flag_interpolated_negatives_ionizing()

        # Dust
        if self.config.interpolate_dust_negatives: self.interpolate_negatives_dust_maps()
        else: self.flag_interpolated_negatives_dust()

    # -----------------------------------------------------------------

    def flag_interpolated_negatives(self):

        """
        Thisfunction ...
        :return:
        """

        # Inform the user
        log.info("Flagging the maps as having negatives interpolated ...")

        # Old
        self.flag_interpolated_negatives_old()

        # Young
        self.flag_interpolated_negatives_young()

        # Ionizing
        self.flag_interpolated_negatives_ionizing()

        # Dust
        self.flag_interpolated_negatives_dust()

    # -----------------------------------------------------------------

    def flag_interpolated_negatives_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[interpolate_negatives_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_negatives_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[interpolate_negatives_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_negatives_ionizing(self):

        """
        Thisn function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[interpolate_negatives_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_negatives_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[interpolate_negatives_step] = True

    # -----------------------------------------------------------------

    def is_interpolated_negatives_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[interpolate_negatives_step]

    # -----------------------------------------------------------------

    def is_interpolated_negatives_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[interpolate_negatives_step]

    # -----------------------------------------------------------------

    def is_interpolated_negatives_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[interpolate_negatives_step]

    # -----------------------------------------------------------------

    def is_interpolated_negatives_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[interpolate_negatives_step]

    # -----------------------------------------------------------------

    def has_negatives_mask_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.old_negative_masks and self.old_negative_masks[name] is not None and self.old_negative_masks[name].has_masked

    # -----------------------------------------------------------------

    @property
    def interpolation_settings_old(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["in_cutout"] = self.config.interpolate_in_cutout
        settings["source_outer_factor"] = self.config.source_outer_factor
        settings["smoothing_factor"] = self.config.old_interpolation_smoothing_factor
        settings["sigma_clip"] = self.config.sigma_clip
        settings["interpolation_softening_range"] = self.interpolation_softening_range
        settings["plot"] = self.config.plot_interpolation_old
        settings["relative_softening_radius"] = self.config.relative_softening_radius_interpolation_old
        settings["return_background"] = True

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def get_negatives_dilation_radius_old(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        return the_map.xsize * self.config.old_negatives_relative_dilation_radius

    # -----------------------------------------------------------------

    def interpolate_negatives_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the negatives of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_interpolated_negatives_old(name):
                log.success("The negatives of the '" + name + "' old stellar map is already interpolated")
                continue

            # Check whether there is a negatives mask
            if not self.has_negatives_mask_old(name):
                log.debug("No negatives for the '" + name + "' old stellar map")
                self.old_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Debugging
            log.debug("Interpolating the negatives of '" + name + "' old stellar map ...")

            # Get the map
            the_map = self.old_maps[name]

            # Get the truncation mask
            central_ellipse = self.negatives_central_ellipse.to_pixel(the_map.wcs)
            central_mask = central_ellipse.to_mask(the_map.xsize, the_map.ysize)

            # Save the central ellipse
            if self.config.steps: central_ellipse.saveto(self.old_extra_path_for_map(name, central_ellipse_filename))

            # Get the mask
            mask = self.old_negative_masks[name]
            mask = intersection(mask, central_mask)

            # Check if any
            if mask.all_unmasked:
                log.debug("No negatives for the '" + name + "' old stellar map within the central ellipse")
                self.old_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Save the mask
            if self.config.steps: mask.saveto(self.old_extra_path_for_map(name, negatives_mask_filename))

            # Dilate?
            if self.config.dilate_negatives_old:

                # Determine the dilation radius
                dilation_radius = self.get_negatives_dilation_radius_old(the_map)

                # Debugging
                log.debug("The dilation radius for the negatives mask is " + str(dilation_radius) + " pixels")

                # Debugging
                log.debug("Dilating the negatives mask ...")

                # Dilate the mask
                mask.disk_dilate(radius=dilation_radius)

                # Save the mask again
                if self.config.steps: mask.saveto(self.old_extra_path_for_map(name, negatives_dilated_filename))

            # No dilation
            else: dilation_radius = None

            # Fill holes
            mask.fill_holes()

            # Save the mask again
            if self.config.steps: mask.saveto(self.old_extra_path_for_map(name, negatives_filled_filename))

            # Plot?
            if self.config.plot_interpolation_old: plotting.plot_mask(mask, title="negatives mask")

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.old_maps[name], self.config.interpolation_method, mask, **self.interpolation_settings_old)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.old_extra_path_for_map(name, interpolation_negatives_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.old_extra_path_for_map(name, interpolation_negatives_background_filename))

            # Plot?
            if self.config.plot_interpolation_old: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.old_maps[name].metadata[interpolate_negatives_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, interpolate_negatives_step))

            # Save info
            if self.config.steps: self.write_old_info(name, interpolate_negatives_step, dilation_radius=dilation_radius,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_old)

    # -----------------------------------------------------------------

    def has_negatives_mask_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.young_negative_masks and self.young_negative_masks[name] is not None and self.young_negative_masks[name].has_masked

    # -----------------------------------------------------------------

    @property
    def interpolation_settings_young(self):

        """
        Thisf unction ...
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["in_cutout"] = self.config.interpolate_in_cutout
        settings["source_outer_factor"] = self.config.source_outer_factor
        settings["smoothing_factor"] = self.config.young_interpolation_smoothing_factor
        settings["sigma_clip"] = self.config.sigma_clip
        settings["interpolation_softening_range"] = self.interpolation_softening_range
        settings["plot"] = self.config.plot_interpolation_young
        settings["relative_softening_radius"] = self.config.relative_softening_radius_interpolation_young
        settings["return_background"] = True

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def get_negatives_dilation_radius_young(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        return the_map.xsize * self.config.young_negatives_relative_dilation_radius

    # -----------------------------------------------------------------

    def interpolate_negatives_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the negatives of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_interpolated_negatives_young(name):
                log.success("The negatives of the '" + name + "' young stellar map is already interpolated")
                continue

            # Check whether there is a negatives mask
            if not self.has_negatives_mask_young(name):
                log.debug("No negatives for the '" + name + "' young stellar map")
                self.young_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Debugging
            log.debug("Interpolating the negatives of '" + name + "' young stellar map ...")

            # Get the map
            the_map = self.young_maps[name]

            # Get the truncation mask
            central_ellipse = self.negatives_central_ellipse.to_pixel(the_map.wcs)
            truncation_mask = central_ellipse.to_mask(the_map.xsize, the_map.ysize)

            # Save the central ellipse
            if self.config.steps: central_ellipse.saveto(self.young_extra_path_for_map(name, central_ellipse_filename))

            # Get the mask
            mask = self.young_negative_masks[name]
            mask = intersection(mask, truncation_mask)

            # Check if any
            if mask.all_unmasked:
                log.debug("No negatives for the '" + name + "' young stellar map within the central ellipse")
                self.young_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Save the mask
            if self.config.steps: mask.saveto(self.young_extra_path_for_map(name, negatives_mask_filename))

            # Dilate?
            if self.config.dilate_negatives_young:

                # Determine the dilation radius
                dilation_radius = self.get_negatives_dilation_radius_young(the_map)

                # Debugging
                log.debug("The dilation radius for the negatives mask is " + str(dilation_radius) + " pixels")

                # Debugging
                log.debug("Dilating the negatives mask ...")

                # Dilate the mask
                mask.disk_dilate(radius=dilation_radius)

                # Save the mask again
                if self.config.steps: mask.saveto(self.young_extra_path_for_map(name, negatives_dilated_filename))

            # No dilation
            else: dilation_radius = None

            # Fill holes
            mask.fill_holes()

            # Save the mask again
            if self.config.steps: mask.saveto(self.young_extra_path_for_map(name, negatives_filled_filename))

            # Plot?
            if self.config.plot_interpolation_young: plotting.plot_mask(mask, title="negatives mask")

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.young_maps[name], self.config.interpolation_method, mask, **self.interpolation_settings_young)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.young_extra_path_for_map(name, interpolation_negatives_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.young_extra_path_for_map(name, interpolation_negatives_background_filename))

            # Plot?
            if self.config.plot_interpolation_young: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.young_maps[name].metadata[interpolate_negatives_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, interpolate_negatives_step))

            # Save info
            if self.config.steps: self.write_young_info(name, interpolate_negatives_step, dilation_radius=dilation_radius,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_young)

    # -----------------------------------------------------------------

    def has_negatives_mask_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.ionizing_negative_masks and self.ionizing_negative_masks[name] is not None and self.ionizing_negative_masks[name].has_masked

    # -----------------------------------------------------------------

    @property
    def interpolation_settings_ionizing(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["in_cutout"] = self.config.interpolate_in_cutout
        settings["source_outer_factor"] = self.config.source_outer_factor
        settings["smoothing_factor"] = self.config.ionizing_interpolation_smoothing_factor
        settings["sigma_clip"] = self.config.sigma_clip
        settings["interpolation_softening_range"] = self.interpolation_softening_range
        settings["plot"] = self.config.plot_interpolation_ionizing
        settings["relative_softening_radius"] = self.config.relative_softening_radius_interpolation_ionizing
        settings["return_background"] = True

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def get_negatives_dilation_radius_ionizing(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        return the_map.xsize * self.config.ionizing_negatives_relative_dilation_radius

    # -----------------------------------------------------------------

    def interpolate_negatives_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the negatives of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_interpolated_negatives_ionizing(name):
                log.success("The negatives of the '" + name + "' ionizing stellar map is already interpolated")
                continue

            # Check whether there is a negatives mask
            if not self.has_negatives_mask_ionizing(name):
                log.debug("No negatives for the '" + name + "' ionizing stellar map")
                self.ionizing_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Debugging
            log.debug("Interpolating the negatives of '" + name + "' ionizing stellar map ...")

            # Get the map
            the_map = self.ionizing_maps[name]

            # Get the truncation mask
            central_ellipse = self.negatives_central_ellipse.to_pixel(the_map.wcs)
            truncation_mask = central_ellipse.to_mask(the_map.xsize, the_map.ysize)

            # Save the central ellipse
            if self.config.steps: central_ellipse.saveto(self.ionizing_extra_path_for_map(name, central_ellipse_filename))

            # Get the mask
            mask = self.ionizing_negative_masks[name]
            mask = intersection(mask, truncation_mask)

            # Check if any
            if mask.all_unmasked:
                log.debug("No negatives for the '" + name + "' ionizing stellar map within the central ellipse")
                self.ionizing_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Save the mask
            if self.config.steps: mask.saveto(self.ionizing_extra_path_for_map(name, negatives_mask_filename))

            # Dilate?
            if self.config.dilate_negatives_ionizing:

                # Determine the dilation radius
                dilation_radius = self.get_negatives_dilation_radius_ionizing(the_map)

                # Debugging
                log.debug("Dilating the negatives mask ...")

                # Debugging
                log.debug("The dilation radius for the negatives mask is " + str(dilation_radius) + " pixels")

                # Dilate the mask
                mask.disk_dilate(radius=dilation_radius)

                # Save the mask again
                if self.config.steps: mask.saveto(self.ionizing_extra_path_for_map(name, negatives_dilated_filename))

            # No dilation
            else: dilation_radius = None

            # Fill holes
            mask.fill_holes()

            # Save the mask again
            if self.config.steps: mask.saveto(self.ionizing_extra_path_for_map(name, negatives_filled_filename))

            # Plot?
            if self.config.plot_interpolation_ionizing: plotting.plot_mask(mask, title="negatives mask")

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.ionizing_maps[name], self.config.interpolation_method, mask, **self.interpolation_settings_ionizing)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.ionizing_extra_path_for_map(name, interpolation_negatives_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.ionizing_extra_path_for_map(name, interpolation_negatives_background_filename))

            # Plot?
            if self.config.plot_interpolation_ionizing: plotting.plot_box(interpolation_mask, title="interpolation mask")

            # Set flag
            self.ionizing_maps[name].metadata[interpolate_negatives_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, interpolate_negatives_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, interpolate_negatives_step, dilation_radius=dilation_radius,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_ionizing)

    # -----------------------------------------------------------------

    def has_negatives_mask_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.dust_negative_masks and self.dust_negative_masks[name] is not None and self.dust_negative_masks[name].has_masked

    # -----------------------------------------------------------------

    @property
    def interpolation_settings_dust(self):

        """
        This function ...
        :return:
        """

        # Initialize a dictionary
        settings = dict()

        # Add settings
        settings["in_cutout"] = self.config.interpolate_in_cutout
        settings["source_outer_factor"] = self.config.source_outer_factor
        settings["smoothing_factor"] = self.config.dust_interpolation_smoothing_factor
        settings["sigma_clip"] = self.config.sigma_clip
        settings["interpolation_softening_range"] = self.interpolation_softening_range
        settings["plot"] = self.config.plot_interpolation_dust
        settings["relative_softening_radius"] = self.config.relative_softening_radius_interpolation_dust
        settings["return_background"] = True

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def get_negatives_dilation_radius_dust(self, the_map):

        """
        This function ...
        :param the_map:
        :return:
        """

        return the_map.xsize * self.config.dust_negatives_relative_dilation_radius

    # -----------------------------------------------------------------

    def interpolate_negatives_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the negatives of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_interpolated_negatives_dust(name):
                log.success("The negatives of the '" + name + "' dust map is already interpolated")
                continue

            # Check whether there is a negatives mask
            if not self.has_negatives_mask_dust(name):
                log.debug("No negatives for the '" + name + "' dust map")
                self.dust_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Debugging
            log.debug("Interpolating the negatives of '" + name + "' dust stellar map ...")

            # Get the map
            the_map = self.dust_maps[name]

            # Get the truncation mask
            central_ellipse = self.negatives_central_ellipse.to_pixel(the_map.wcs)
            truncation_mask = central_ellipse.to_mask(the_map.xsize, the_map.ysize)

            # Save the central ellipse
            if self.config.steps: central_ellipse.saveto(self.dust_extra_path_for_map(name, central_ellipse_filename))

            # Get the mask
            mask = self.dust_negative_masks[name]
            mask = intersection(mask, truncation_mask)

            # Check if any
            if mask.all_unmasked:
                log.debug("No negatives for the '" + name + "' dust map within the central ellipse")
                self.dust_maps[name].metadata[interpolate_negatives_step] = True
                continue

            # Save the mask
            if self.config.steps: mask.saveto(self.dust_extra_path_for_map(name, negatives_mask_filename))

            # Dilate?
            if self.config.dilate_negatives_dust:

                # Determine the dilation radius
                dilation_radius = self.get_negatives_dilation_radius_dust(the_map)

                # Debugging
                log.debug("The dilation radius for the negatives mask is " + str(dilation_radius) + " pixels")

                # Debugging
                log.debug("Dilating the negatives mask ...")

                # Dilate the mask
                mask.disk_dilate(radius=dilation_radius)

                # Save the mask again
                if self.config.steps: mask.saveto(self.dust_extra_path_for_map(name, negatives_dilated_filename))

            # No dilation
            else: dilation_radius = None

            # Fill holes
            mask.fill_holes()

            # Save the mask again
            if self.config.steps: mask.saveto(self.dust_extra_path_for_map(name, negatives_filled_filename))

            # Plot?
            if self.config.plot_interpolation_dust: plotting.plot_mask(mask, title="negatives mask")

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.dust_maps[name], self.config.interpolation_method, mask, **self.interpolation_settings_dust)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.dust_extra_path_for_map(name, interpolation_negatives_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.dust_extra_path_for_map(name, interpolation_negatives_background_filename))

            # Plot?
            if self.config.plot_interpolation_dust: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.dust_maps[name].metadata[interpolate_negatives_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, interpolate_negatives_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, interpolate_negatives_step, dilation_radius=dilation_radius,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_dust)

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_softening_radius(self):

        """
        This function ...
        :return:
        """

        return numbers.geometric_mean(self.config.interpolation_softening_start, self.config.interpolation_softening_end)

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_softening_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(self.config.interpolation_softening_start / self.interpolation_softening_radius, self.config.interpolation_softening_end / self.interpolation_softening_radius)

    # -----------------------------------------------------------------

    def interpolate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the maps ...")

        # Old
        if self.config.interpolate_old: self.interpolate_old_maps()
        else: self.flag_interpolated_old()

        # Young
        if self.config.interpolate_young: self.interpolate_young_maps()
        else: self.flag_interpolated_young()

        # Ionizing
        if self.config.interpolate_ionizing: self.interpolate_ionizing_maps()
        else: self.flag_interpolated_ionizing()

        # Dust
        if self.config.interpolate_dust: self.interpolate_dust_maps()
        else: self.flag_interpolated_dust()

    # -----------------------------------------------------------------

    def flag_interpolated(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging the maps as being interpolated ...")

        # Old
        self.flag_interpolated_old()

        # Young
        self.flag_interpolated_young()

        # Ionizing
        self.flag_interpolated_ionizing()

        # Dust
        self.flag_interpolated_dust()

    # -----------------------------------------------------------------

    def flag_interpolated_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[interpolate_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[interpolate_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_ionizing(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[interpolate_step] = True

    # -----------------------------------------------------------------

    def flag_interpolated_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[interpolate_step] = True

    # -----------------------------------------------------------------

    def is_interpolated_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[interpolate_step]

    # -----------------------------------------------------------------

    def is_interpolated_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[interpolate_step]

    # -----------------------------------------------------------------

    def is_interpolated_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[interpolate_step]

    # -----------------------------------------------------------------

    def is_interpolated_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[interpolate_step]

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_old(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.old_core_region_factor
        ellipse.angle += self.config.interpolation_angle_offset_old
        return ellipse

    # -----------------------------------------------------------------

    def interpolate_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_interpolated_old(name):
                log.success("The '" + name + "' old stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' old stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_old.to_pixel(self.old_maps[name].wcs)

            # Save the interpolation ellipse
            if self.config.steps: ellipse.saveto(self.old_extra_path_for_map(name, interpolation_ellipse_filename))

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.old_maps[name], self.config.interpolation_method, ellipse, **self.interpolation_settings_old)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.old_extra_path_for_map(name, interpolation_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.old_extra_path_for_map(name, interpolation_background_filename))

            # Plot?
            if self.config.plot_interpolation_old: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.old_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, interpolate_step))

            # Save info
            if self.config.steps: self.write_old_info(name, interpolate_step,
                                                      core_region_factor=self.config.old_core_region_factor,
                                                      region_angle_offset = self.config.interpolation_angle_offset_old,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_old)

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_young(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.young_core_region_factor
        ellipse.angle += self.config.interpolation_angle_offset_young
        return ellipse

    # -----------------------------------------------------------------

    def interpolate_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_interpolated_young(name):
                log.success("The '" + name + "' young stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' young stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_young.to_pixel(self.young_maps[name].wcs)

            # Save the interpolation ellipse
            if self.config.steps: ellipse.saveto(self.young_extra_path_for_map(name, interpolation_ellipse_filename))

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.young_maps[name], self.config.interpolation_method, ellipse, **self.interpolation_settings_young)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.young_extra_path_for_map(name, interpolation_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.young_extra_path_for_map(name, interpolation_background_filename))

            # Plot?
            if self.config.plot_interpolation_young: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.young_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, interpolate_step))

            # Save info
            if self.config.steps: self.write_young_info(name, interpolate_step,
                                                      core_region_factor=self.config.young_core_region_factor,
                                                      region_angle_offset=self.config.interpolation_angle_offset_young,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_young)

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_ionizing(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.ionizing_core_region_factor
        ellipse.angle += self.config.interpolation_angle_offset_ionizing
        return ellipse

    # -----------------------------------------------------------------

    def interpolate_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_interpolated_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' ionizing stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_ionizing.to_pixel(self.ionizing_maps[name].wcs)

            # Save the interpolation ellipse
            if self.config.steps: ellipse.saveto(self.ionizing_extra_path_for_map(name, interpolation_ellipse_filename))

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.ionizing_maps[name], self.config.interpolation_method, ellipse, **self.interpolation_settings_ionizing)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.ionizing_extra_path_for_map(name, interpolation_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.ionizing_extra_path_for_map(name, interpolation_background_filename))

            # Plot?
            if self.config.plot_interpolation_ionizing: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.ionizing_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, interpolate_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, interpolate_step,
                                                      core_region_factor=self.config.ionizing_core_region_factor,
                                                      region_angle_offset=self.config.interpolation_angle_offset_ionizing,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_ionizing)

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_dust(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.dust_core_region_factor
        ellipse.angle += self.config.interpolation_angle_offset_dust
        return ellipse

    # -----------------------------------------------------------------

    def interpolate_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Interpolating the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_interpolated_dust(name):
                log.success("The '" + name + "' dust map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' dust map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_dust.to_pixel(self.dust_maps[name].wcs)

            # Save the interpolation ellipse
            if self.config.steps: ellipse.saveto(self.dust_extra_path_for_map(name, interpolation_ellipse_filename))

            # Interpolate
            interpolation_mask, interpolation_background = interpolate_map(self.dust_maps[name], self.config.interpolation_method, ellipse, **self.interpolation_settings_dust)

            # Save the interpolation mask
            if self.config.steps: interpolation_mask.saveto(self.dust_extra_path_for_map(name, interpolation_mask_filename))

            # Save the interpolation background
            if self.config.steps and interpolation_background is not None: interpolation_background.saveto(self.dust_extra_path_for_map(name, interpolation_background_filename))

            # Plot?
            if self.config.plot_interpolation_dust: plotting.plot_mask(interpolation_mask, title="interpolation mask")

            # Set flag
            self.dust_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, interpolate_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, interpolate_step,
                                                      core_region_factor=self.config.dust_core_region_factor,
                                                      region_angle_offset=self.config.interpolation_angle_offset_dust,
                                                      method=self.config.interpolation_method,
                                                      **self.interpolation_settings_dust)

    # -----------------------------------------------------------------

    def truncate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the maps ...")

        # Old
        if self.config.truncate_old: self.truncate_old_maps()
        else: self.flag_truncated_old()

        # Young
        if self.config.truncate_young: self.truncate_young_maps()
        else: self.flag_truncated_young()

        # Ionizing
        if self.config.truncate_ionizing: self.truncate_ionizing_maps()
        else: self.flag_truncated_ionizing()

        # Dust
        if self.config.truncate_dust: self.truncate_dust_maps()
        else: self.flag_truncated_dust()

    # -----------------------------------------------------------------

    def flag_truncated(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging the maps as being truncated ...")

        # Old
        self.flag_truncated_old()

        # Young
        self.flag_truncated_young()

        # Ionizing
        self.flag_truncated_ionizing()

        # Dust
        self.flag_truncated_dust()

    # -----------------------------------------------------------------

    def flag_truncated_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[truncate_step] = True

    # -----------------------------------------------------------------

    def flag_truncated_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[truncate_step] = True

    # -----------------------------------------------------------------

    def flag_truncated_ionizing(self):

        """
        Thisf unction ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[truncate_step] = True

    # -----------------------------------------------------------------

    def flag_truncated_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[truncate_step] = True

    # -----------------------------------------------------------------

    def is_truncated_old(self, name):

        """
        This function ..
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[truncate_step]

    # -----------------------------------------------------------------

    def is_truncated_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[truncate_step]

    # -----------------------------------------------------------------

    def is_truncated_ionizing(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[truncate_step]

    # -----------------------------------------------------------------

    def is_truncated_dust(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[truncate_step]

    # -----------------------------------------------------------------

    def truncate_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_truncated_old(name):
                log.success("The '" + name + "' old stellar map is already truncated")
                continue

            # Debugging
            log.debug("Truncating the '" + name + "' old stellar map ...")

            # Truncate the map
            npixels = self.truncate_map(self.old_maps[name])

            # Set flag
            self.old_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, truncate_step))

            # Save info
            if self.config.steps: self.write_old_info(name, truncate_step, npixels=npixels)

    # -----------------------------------------------------------------

    def truncate_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_truncated_young(name):
                log.success("The '" + name + "' young stellar map is already truncated")
                continue

            # Debugging
            log.debug("Truncating the '" + name + "' young stellar map ...")

            # Truncate the map
            npixels = self.truncate_map(self.young_maps[name])

            # Set flag
            self.young_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, truncate_step))

            # Save info
            if self.config.steps: self.write_young_info(name, truncate_step, npixels=npixels)

    # -----------------------------------------------------------------

    def truncate_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_truncated_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already truncated")
                continue

            # Debugging
            log.debug("Truncating the '" + name + "' ionizing stellar map ...")

            # Truncate the map
            npixels = self.truncate_map(self.ionizing_maps[name])

            # Set flag
            self.ionizing_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, truncate_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, truncate_step, npixels=npixels)

    # -----------------------------------------------------------------

    def truncate_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_truncated_dust(name):
                log.success("The '" + name + "' dust map is already truncated")
                continue

            # Debugging
            log.debug("Truncating the '" + name + "' dust map ...")

            # Truncate the map
            npixels = self.truncate_map(self.dust_maps[name])

            # Set flag
            self.dust_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, truncate_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, truncate_step, npixels=npixels)

    # -----------------------------------------------------------------

    def is_cropped_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[crop_step]

    # -----------------------------------------------------------------

    def is_cropped_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[crop_step]

    # -----------------------------------------------------------------

    def is_cropped_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[crop_step]

    # -----------------------------------------------------------------

    def is_cropped_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[crop_step]

    # -----------------------------------------------------------------

    def crop_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the maps ...")

        # Old
        if self.config.crop_old: self.crop_old_maps()
        else: self.flag_cropped_old()

        # Young
        if self.config.crop_young: self.crop_young_maps()
        else: self.flag_cropped_young()

        # Ionizing
        if self.config.crop_ionizing: self.crop_ionizing_maps()
        else: self.flag_cropped_ionizing()

        # Dust
        if self.config.crop_dust: self.crop_dust_maps()
        else: self.flag_cropped_dust()

    # -----------------------------------------------------------------

    def flag_cropped(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging the maps as being cropped ...")

        # Old
        self.flag_cropped_old()

        # Young
        self.flag_cropped_young()

        # Ionizing
        self.flag_cropped_ionizing()

        # Dust
        self.flag_cropped_dust()

    # -----------------------------------------------------------------

    def flag_cropped_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[crop_step] = True

    # -----------------------------------------------------------------

    def flag_cropped_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[crop_step] = True

    # -----------------------------------------------------------------

    def flag_cropped_ionizing(self):

        """
        Thisf ucntion ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[crop_step] = True

    # -----------------------------------------------------------------

    def flag_cropped_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[crop_step] = True

    # -----------------------------------------------------------------

    def crop_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_cropped_old(name):
                log.success("The '" + name + "' old stellar map is already cropped")
                continue

            # Debugging
            log.debug("Cropping the '" + name + "' old stellar map ...")

            # Crop the map
            mask = self.old_negative_masks[name] if self.has_negatives_mask_old(name) else None
            self.crop_map(self.old_maps[name], self.config.cropping_factor, mask=mask)

            # Save the cropped mask
            if mask is not None and self.config.steps: self.old_negative_masks[name].saveto(self.old_extra_path_for_map(name, cropped_negatives_filename))

            # Set flag
            self.old_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, crop_step))

            # Save info
            if self.config.steps: self.write_old_info(name, crop_step, cropping_factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def crop_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_cropped_young(name):
                log.success("The '" + name + "' young stellar map is already cropped")
                continue

            # Debugging
            log.debug("Cropping the '" + name + "' young stellar map ...")

            # Crop the map
            mask = self.young_negative_masks[name] if self.has_negatives_mask_young(name) else None
            self.crop_map(self.young_maps[name], self.config.cropping_factor, mask=mask)

            # Save the cropped mask
            if mask is not None and self.config.steps: self.young_negative_masks[name].saveto(self.young_extra_path_for_map(name, cropped_negatives_filename))

            # Set flag
            self.young_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, crop_step))

            # Save info
            if self.config.steps: self.write_young_info(name, crop_step, cropping_factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def crop_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_cropped_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already cropped")
                continue

            # Debugging
            log.debug("Cropping the '" + name + "' ionizing stellar map ...")

            # Crop the map
            mask = self.ionizing_negative_masks[name] if self.has_negatives_mask_ionizing(name) else None
            self.crop_map(self.ionizing_maps[name], self.config.cropping_factor, mask=mask)

            # Save the cropped mask
            if mask is not None and self.config.steps: self.ionizing_negative_masks[name].saveto(self.ionizing_extra_path_for_map(name, cropped_negatives_filename))

            # Set flag
            self.ionizing_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, crop_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, crop_step, cropping_factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def crop_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Cropping the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_cropped_dust(name):
                log.success("The '" + name + "' dust map is already cropped")
                continue

            # Debugging
            log.debug("Cropping the '" + name + "' dust map ...")

            # Crop the map
            mask = self.dust_negative_masks[name] if self.has_negatives_mask_dust(name) else None
            self.crop_map(self.dust_maps[name], self.config.cropping_factor, mask=mask)

            # Save the cropped mask
            if mask is not None and self.config.steps: self.dust_negative_masks[name].saveto(self.dust_extra_path_for_map(name, cropped_negatives_filename))

            # Set flag
            self.dust_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, crop_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, crop_step, cropping_factor=self.config.cropping_factor)

    # -----------------------------------------------------------------

    def is_clipped_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[clip_step]

    # -----------------------------------------------------------------

    def is_clipped_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[clip_step]

    # -----------------------------------------------------------------

    def is_clipped_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[clip_step]

    # -----------------------------------------------------------------

    def is_clipped_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[clip_step]

    # -----------------------------------------------------------------

    def clip_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the maps ...")

        # Old
        if self.config.clip_old: self.clip_old_maps()
        else: self.flag_clipped_old()

        # Young
        if self.config.clip_young: self.clip_young_maps()
        else: self.flag_clipped_young()

        # Ionizing
        if self.config.clip_ionizing: self.clip_ionizing_maps()
        else: self.flag_clipped_ionizing()

        # Dust
        if self.config.clip_dust: self.clip_dust_maps()
        else: self.flag_clipped_dust()

    # -----------------------------------------------------------------

    def flag_clipped(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging maps as being clipped ...")

        # Old
        self.flag_clipped_old()

        # Young
        self.flag_clipped_young()

        # Ionizing
        self.flag_clipped_ionizing()

        # Dust
        self.flag_clipped_dust()

    # -----------------------------------------------------------------

    def flag_clipped_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[clip_step] = True

    # -----------------------------------------------------------------

    def flag_clipped_young(self):

        """
        Thisn function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[clip_step] = True

    # -----------------------------------------------------------------

    def flag_clipped_ionizing(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[clip_step] = True

    # -----------------------------------------------------------------

    def flag_clipped_dust(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[clip_step] = True

    # -----------------------------------------------------------------

    def has_old_clip_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.old_masks_path, name + "_clip.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.old_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_old_settings(self, name):

        """
        Thisf unction ...
        :return:
        """

        # Initialize dictionary
        settings = dict()

        # Basic settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["output_path"] = self.old_clipping_path_for_map(name)

        # NEW
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_old
        settings["boundary"] = self.old_maps_boundary
        settings["plot"] = self.config.plot_clipping_old
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_old

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_clipped_old(name):
                log.success("The '" + name + "' old stellar map is already clipped")
                if name not in self.old_clip_masks and not self.has_old_clip_mask(name): log.warning("Yet the clip mask is not present: performing clipping step again ...")
                else: continue

            # Debugging
            log.debug("Clipping the '" + name + "' old stellar map ...")

            # Get the origins
            origins = self.old_map_origins[name]

            # Remove filters to be ignored
            if self.config.ignore_filters_clipping is not None: origins = sequences.removed(origins, self.config.ignore_filters_clipping)

            # Clip the map (returns the mask)
            settings = self.get_clip_old_settings(name)
            self.old_clip_masks[name] = self.clip_map(self.old_maps[name], origins, **settings)

            # Set flag
            self.old_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.old_clip_masks[name].saveto(self.old_step_path_for_mask(name, clip_step))

            # Save info
            settings.pop("remote")
            settings.pop("boundary")
            settings["ignore_filters_clipping"] = self.config.ignore_filters_clipping
            if self.config.steps: self.write_old_info(name, clip_step, origins=origins, **settings)

    # -----------------------------------------------------------------

    def has_young_clip_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.young_masks_path, name + "_clip.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def young_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.young_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_young_settings(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Initialize dictionary
        settings = dict()

        # Basic settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["output_path"] = self.young_clipping_path_for_map(name)

        # NEW
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_young
        settings["boundary"] = self.young_maps_boundary
        settings["plot"] = self.config.plot_clipping_young
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_young

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_clipped_young(name):
                log.success("The '" + name + "' young stellar map is already clipped")
                if name not in self.young_clip_masks and not self.has_young_clip_mask(name): log.warning("Yet the clip mask is not present: performing clipping step again ...")
                else: continue

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            origins = self.young_map_origins[name]

            # Remove filters to be ignored
            if self.config.ignore_filters_clipping is not None: origins = sequences.removed(origins, self.config.ignore_filters_clipping)

            # Clip the map (returns the mask)
            settings = self.get_clip_young_settings(name)
            self.young_clip_masks[name] = self.clip_map(self.young_maps[name], origins, **settings)

            # Set flag
            self.young_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.young_clip_masks[name].saveto(self.young_step_path_for_mask(name, clip_step))

            # Save info
            settings.pop("remote")
            settings.pop("boundary")
            settings["ignore_filters_clipping"] = self.config.ignore_filters_clipping
            if self.config.steps: self.write_young_info(name, clip_step, origins=origins, **settings)

    # -----------------------------------------------------------------

    def has_ionizing_clip_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.ionizing_masks_path, name + "_clip.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.ionizing_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_ionizing_settings(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Initialize dictionary
        settings = dict()

        # Basic settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["output_path"] = self.ionizing_clipping_path_for_map(name)

        # NEW
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_ionizing
        settings["boundary"] = self.ionizing_maps_boundary
        settings["plot"] = self.config.plot_clipping_ionizing
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_ionizing

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_clipped_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already clipped")
                if name not in self.ionizing_clip_masks and not self.has_ionizing_clip_mask(name): log.warning("Yet the clip mask is not present: performing clipping step again ...")
                else: continue

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Remove filters to be ignored
            if self.config.ignore_filters_clipping is not None: origins = sequences.removed(origins, self.config.ignore_filters_clipping)

            # Clip the map (returns the mask)
            settings = self.get_clip_ionizing_settings(name)
            self.ionizing_clip_masks[name] = self.clip_map(self.ionizing_maps[name], origins, **settings)

            # Set flag
            self.ionizing_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.ionizing_clip_masks[name].saveto(self.ionizing_step_path_for_mask(name, clip_step))

            # Save info
            settings.pop("remote")
            settings.pop("boundary")
            settings["ignore_filters_clipping"] = self.config.ignore_filters_clipping
            if self.config.steps: self.write_ionizing_info(name, clip_step, origins=origins, **settings)

    # -----------------------------------------------------------------

    def has_dust_clip_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.dust_masks_path, name + "_clip.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_maps_boundary(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.config.dust_compactness_factor

    # -----------------------------------------------------------------

    def get_clip_dust_settings(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Initialize dictionary
        settings = dict()

        # Basic settings
        settings["convolve"] = self.config.convolve
        settings["remote"] = self.remote
        settings["npixels"] = self.config.min_npixels
        settings["connectivity"] = self.config.connectivity
        settings["rebin_remote_threshold"] = self.config.rebin_remote_threshold
        settings["fuzzy"] = self.config.fuzzy_mask
        settings["fuzziness"] = self.config.fuzziness
        settings["fuzziness_offset"] = self.config.fuzzy_min_significance_offset
        settings["output_path"] = self.dust_clipping_path_for_map(name)

        # NEW
        settings["dilate"] = self.config.dilate_masks
        settings["dilate_fuzzy"] = self.config.dilate_fuzzy_masks
        settings["soften"] = self.config.soften_masks
        settings["relative_softening_radius"] = self.config.relative_softening_radius_dust
        settings["boundary"] = self.dust_maps_boundary
        settings["plot"] = self.config.plot_clipping_dust
        settings["relative_dilation_radius"] = self.config.relative_dilation_radius_dust

        # Return the settings
        return settings

    # -----------------------------------------------------------------

    def clip_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clipping the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_clipped_dust(name):
                log.success("The '" + name + "' dust map is already clipped")
                if name not in self.dust_clip_masks and not self.has_dust_clip_mask(name): log.warning("Yet the clip mask is not present: performing clipping step again ...")
                else: continue

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            origins = self.dust_map_origins[name]

            # Remove filters to be ignored
            if self.config.ignore_filters_clipping is not None: origins = sequences.removed(origins, self.config.ignore_filters_clipping)

            # Clip the map (returns the mask)
            settings = self.get_clip_dust_settings(name)
            self.dust_clip_masks[name] = self.clip_map(self.dust_maps[name], origins, **settings)

            # Set flag
            self.dust_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.dust_clip_masks[name].saveto(self.dust_step_path_for_mask(name, clip_step))

            # Save info
            settings.pop("remote")
            settings.pop("boundary")
            settings["ignore_filters_clipping"] = self.config.ignore_filters_clipping
            if self.config.steps: self.write_dust_info(name, clip_step, origins=origins, **settings)

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_ellipse(self):

        """
        This function ...
        :return:
        """

        return self.truncation_ellipse * self.softening_radius

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_radius(self):

        """
        This function ...
        :return:
        """

        return numbers.geometric_mean(self.config.softening_start, 1.)

    # -----------------------------------------------------------------

    @lazyproperty
    def softening_range(self):

        """
        This function ...
        :return:
        """

        return RealRange(self.config.softening_start / self.softening_radius, 1. / self.softening_radius)

    # -----------------------------------------------------------------

    def is_softened_old(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.old_maps[name].metadata[softened_step]

    # -----------------------------------------------------------------

    def is_softened_young(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.young_maps[name].metadata[softened_step]

    # -----------------------------------------------------------------

    def is_softened_ionizing(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.ionizing_maps[name].metadata[softened_step]

    # -----------------------------------------------------------------

    def is_softened_dust(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.dust_maps[name].metadata[softened_step]

    # -----------------------------------------------------------------

    def soften_edges(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the map edges ...")

        # Old
        if self.config.soften_old: self.soften_edges_old()
        else: self.flag_softened_old()

        # Young
        if self.config.soften_young: self.soften_edges_young()
        else: self.flag_softened_young()

        # Ionizing
        if self.config.soften_ionizing: self.soften_edges_ionizing()
        else: self.flag_softened_ionizing()

        # Dust
        if self.config.soften_dust: self.soften_edges_dust()
        else: self.flag_softened_dust()

    # -----------------------------------------------------------------

    def flag_softened(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Flagging maps as being softened ...")

        # Old
        self.flag_softened_old()

        # Young
        self.flag_softened_young()

        # Ionizing
        self.flag_softened_ionizing()

        # Dust
        self.flag_softened_dust()

    # -----------------------------------------------------------------

    def flag_softened_old(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.old_maps:

            # Set flag
            self.old_maps[name].metadata[softened_step] = True

    # -----------------------------------------------------------------

    def flag_softened_young(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.young_maps:

            # Set flag
            self.young_maps[name].metadata[softened_step] = True

    # -----------------------------------------------------------------

    def flag_softened_ionizing(self):

        """
        This function ...
        :return:
        """

        # Loop over the maps
        for name in self.ionizing_maps:

            # Set flag
            self.ionizing_maps[name].metadata[softened_step] = True

    # -----------------------------------------------------------------

    def flag_softened_dust(self):

        """
        Thisf nuction ...
        :return:
        """

        # Loop over the maps
        for name in self.dust_maps:

            # Set flag
            self.dust_maps[name].metadata[softened_step] = True

    # -----------------------------------------------------------------

    def has_old_softening_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.old_masks_path, name + "_softening.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def has_young_softening_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.young_masks_path, name + "_softening.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def has_ionizing_softening_mask(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.ionizing_masks_path, name + "_softening.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def has_dust_softening_mask(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        # Determine the path
        path = fs.join(self.dust_masks_path, name + "_softening.fits")
        return fs.is_file(path)

    # -----------------------------------------------------------------

    def soften_edges_old(self):

        """
        Thisn function ...
        :return:
        """

        # Infomr the user
        log.info("Softening the edges of the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Check
            if self.is_softened_old(name):
                log.success("The '" + name + "' old stellar map is already softened")
                if name not in self.old_softening_masks and not self.has_old_softening_mask(name):
                    log.warning("Yet the softening mask is not present: performing softening step again ...")
                else: continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' old stellar map ...")

            # Soften the map, get the mask
            alpha_mask = self.soften_map(self.old_maps[name], self.softening_ellipse, self.softening_range)

            # Set the mask
            self.old_softening_masks[name] = alpha_mask

            # Set flag
            self.old_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, softened_step))

            # Save the mask
            if self.config.steps: self.old_softening_masks[name].saveto(self.old_step_path_for_mask(name, softened_step))

            # Save info
            if self.config.steps: self.write_old_info(name, softened_step, softening_start=self.config.softening_start,
                                                      softening_range=self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Check
            if self.is_softened_young(name):
                log.success("The '" + name + "' young stellar map is already softened")
                if name not in self.young_softening_masks and not self.has_young_softening_mask(name):
                    log.warning("Yet the softening mask is not present: performing softening step again ...")
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' young stellar map ...")

            # Soften
            alpha_mask = self.soften_map(self.young_maps[name], self.softening_ellipse, self.softening_range)

            # Set the mask
            self.young_softening_masks[name] = alpha_mask

            # Set flag
            self.young_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, softened_step))

            # Save the mask
            if self.config.steps: self.young_softening_masks[name].saveto(self.young_step_path_for_mask(name, softened_step))

            # Save info
            if self.config.steps: self.write_young_info(name, softened_step, softening_start=self.config.softening_start,
                                                      softening_range=self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Check
            if self.is_softened_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already softened")
                if name not in self.ionizing_softening_masks and not self.has_ionizing_softening_mask(name):
                    log.warning("Yet the softening mask is not present: performing softening step again ...")
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' ionizing stellar map ...")

            # Soften the map, get the mask
            alpha_mask = self.soften_map(self.ionizing_maps[name], self.softening_ellipse, self.softening_range)

            # Set the mask
            self.ionizing_softening_masks[name] = alpha_mask

            # Set flag
            self.ionizing_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, softened_step))

            # Save the mask
            if self.config.steps: self.ionizing_softening_masks[name].saveto(self.ionizing_step_path_for_mask(name, softened_step))

            # Save info
            if self.config.steps: self.write_ionizing_info(name, softened_step, softening_start=self.config.softening_start,
                                                      softening_range=self.softening_range)

    # -----------------------------------------------------------------

    def soften_edges_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Softening the edges of the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Check
            if self.is_softened_dust(name):
                log.success("The '" + name + "' dust map is already softened")
                if name not in self.dust_softening_masks and not self.has_dust_softening_mask(name):
                    log.warning("Yet the softening mask is not present: performing softening step again ...")
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' dust map ...")

            # Soften the map, get the alpha mask
            alpha_mask = self.soften_map(self.dust_maps[name], self.softening_ellipse, self.softening_range)

            # Set the mask
            self.dust_softening_masks[name] = alpha_mask

            # Set flag
            self.dust_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, softened_step))

            # Save the mask
            if self.config.steps: self.dust_softening_masks[name].saveto(self.dust_step_path_for_mask(name, softened_step))

            # Save info
            if self.config.steps: self.write_dust_info(name, softened_step, softening_start=self.config.softening_start,
                                                      softening_range=self.softening_range)

    # -----------------------------------------------------------------

    @lazyproperty
    def old_scaleheight(self):

        """
        This function ...
        :return:
        """

        scale_height = self.disk2d_model.scalelength / self.config.scalelength_to_scaleheight
        return scale_height

    # -----------------------------------------------------------------

    @lazyproperty
    def young_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.config.young_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_scaleheight(self):

        """
        This function ...
        :return:
        """

        return self.config.ionizing_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_scaleheight(self):

        """
        This fucntion ...
        :return:
        """

        return self.config.dust_scaleheight_ratio * self.old_scaleheight

    # -----------------------------------------------------------------

    def deproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the maps ...")

        # Old
        self.deproject_old()

        # Young
        self.deproject_young()

        # Ionizing
        self.deproject_ionizing()

        # Dust
        self.deproject_dust()

    # -----------------------------------------------------------------

    @property
    def has_old_deprojected(self):

        """
        Thisf unction ...
        :return:
        """

        # Loop over all old stellar maps
        for name in self.old_selection:

            # Determine path
            path = fs.join(self.old_deprojection_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_old_deprojected_skirt(self):

        """
        This function ...
        :return: 
        """

        # Loop over all old stellar maps
        for name in self.old_selection:

            # Determine path
            path = fs.join(self.old_deprojection_skirt_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_old_deprojections(self):

        """
        This function ...
        :return:
        """

        if self.old_deprojections is None: return False
        for name in self.old_selection:
            if name not in self.old_deprojections: return False
        return True

    # -----------------------------------------------------------------

    def set_old_deprojection_aliases(self):

        """
        This function ...
        :return:
        """

        # Loop over the old stellar deprojection models
        for name in self.old_deprojections:

            # Set the alias
            #self.old_deprojections[name].filename = "old___" + name + ".fits"
            self.old_deprojections[name].filename = name

            # Set the map
            self.old_deprojections[name].map = self.old_maps[name]

    # -----------------------------------------------------------------

    def deproject_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the old stellar maps ...")

        # Create deprojected old stellar maps
        if not self.has_old_deprojected:
            self.old_deprojections, self.old_deprojected = self.deproject_maps(self.old_maps, scale_height=self.old_scaleheight, root_path=self.old_deprojection_path, return_deprojections=True, downsample_factor=self.config.downsample_factor)
            #self.set_old_deprojection_aliases()

        # Create deprojected old stellar maps with SKIRT
        if not self.has_old_deprojected_skirt:
            if self.has_old_deprojections:
                self.old_deprojected_skirt = self.deproject_models(self.old_deprojections, root_path=self.old_deprojection_skirt_path, method="skirt", downsample_factor=self.config.downsample_factor)
            else:
                self.old_deprojections, self.old_deprojected_skirt = self.deproject_maps(self.old_maps, scale_height=self.old_scaleheight, root_path=self.old_deprojection_skirt_path, return_deprojections=True, downsample_factor=self.config.downsample_factor, method="skirt")
                self.set_old_deprojection_aliases()

    # -----------------------------------------------------------------

    @property
    def has_young_deprojected(self):

        """
        This function ...
        :return:
        """

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Determine the path
            path = fs.join(self.young_deprojection_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_young_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Loop over the young stellar maps
        for name in self.young_selection:

            # Determine the path
            path = fs.join(self.young_deprojection_skirt_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_young_deprojections(self):

        """
        This function ...
        :return:
        """

        if self.young_deprojections is None: return False
        for name in self.young_selection:
            if name not in self.young_deprojections: return False
        return True

    # -----------------------------------------------------------------

    def set_young_deprojection_aliases(self):

        """
        This function ...
        :return:
        """

        # Loop over the young stellar deprojection models
        for name in self.young_deprojections:

            # Set the alias
            #self.young_deprojections[name].filename = "young___" + name + ".fits"
            self.young_deprojections[name].filename = name

            # Set the map
            self.young_deprojections[name].map = self.young_maps[name]

    # -----------------------------------------------------------------

    def deproject_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the young stellar maps ...")

        # Create deprojected young stellar maps
        if not self.has_young_deprojected:
            self.young_deprojections, self.young_deprojected = self.deproject_maps(self.young_maps, scale_height=self.young_scaleheight, root_path=self.young_deprojection_path, return_deprojections=True, downsample_factor=self.config.downsample_factor)
            #self.set_young_deprojection_aliases()

        # Create deprojected young stellar maps with SKIRT
        if not self.has_young_deprojected_skirt:
            if self.has_young_deprojections:
                self.young_deprojected_skirt = self.deproject_models(self.young_deprojections, root_path=self.young_deprojection_skirt_path, method="skirt", downsample_factor=self.config.downsample_factor)
            else:
                self.young_deprojections, self.young_deprojected_skirt = self.deproject_maps(self.young_maps, scale_height=self.young_scaleheight, root_path=self.young_deprojection_skirt_path, method="skirt", return_deprojections=True, downsample_factor=self.config.downsample_factor)
                self.set_young_deprojection_aliases()

    # -----------------------------------------------------------------

    @property
    def has_ionizing_deprojected(self):

        """
        Thisf unction ...
        :return:
        """

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Detrmine path
            path = fs.join(self.ionizing_deprojection_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_ionizing_deprojected_skirt(self):

        """
        Thisfunction ...
        :return:
        """

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:

            # Determine the path
            path = fs.join(self.ionizing_deprojection_skirt_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_ionizing_deprojections(self):

        """
        This function ...
        :return:
        """

        if self.ionizing_deprojections is None: return False
        for name in self.ionizing_selection:
            if name not in self.ionizing_deprojections: return False
        return True

    # -----------------------------------------------------------------

    def set_ionizing_deprojection_aliases(self):

        """
        This function ...
        :return:
        """

        # Loop over the ionizing stellar deprojection models
        for name in self.ionizing_deprojections:

            # Set alias
            #self.ionizing_deprojections[name].filename = "ionizing___" + name + ".fits"
            self.ionizing_deprojections[name].filename = name

            # Set map
            self.ionizing_deprojections[name].map = self.ionizing_maps[name]

    # -----------------------------------------------------------------

    def deproject_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the ionizing stellar maps ...")

        # Create deprojected ionizing stellar maps
        if not self.has_ionizing_deprojected:
            self.ionizing_deprojections, self.ionizing_deprojected = self.deproject_maps(self.ionizing_maps, scale_height=self.ionizing_scaleheight, root_path=self.ionizing_deprojection_path, return_deprojections=True, downsample_factor=self.config.downsample_factor)
            #self.set_ionizing_deprojection_aliases()

        # Create deprojected ionizing stellar maps with SKIRT
        if not self.has_ionizing_deprojected_skirt:
            if self.has_ionizing_deprojections:
                self.ionizing_deprojected_skirt = self.deproject_models(self.ionizing_deprojections, root_path=self.ionizing_deprojection_skirt_path, method="skirt", downsample_factor=self.config.downsample_factor)
            else:
                self.ionizing_deprojections, self.ionizing_deprojected_skirt = self.deproject_maps(self.ionizing_maps, scale_height=self.ionizing_scaleheight, root_path=self.ionizing_deprojection_skirt_path, method="skirt", return_deprojections=True, downsample_factor=self.config.downsample_factor)
                self.set_ionizing_deprojection_aliases()

    # -----------------------------------------------------------------

    @property
    def has_dust_deprojected(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust maps
        for name in self.dust_selection:

            # Detemrine path
            path = fs.join(self.dust_deprojection_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All cheks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_dust_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust maps
        for name in self.dust_selection:

            # Determine path
            path = fs.join(self.dust_deprojection_skirt_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @property
    def has_dust_deprojections(self):

        """
        This function ...
        :return:
        """

        if self.dust_deprojections is None: return False
        for name in self.dust_selection:
            if name not in self.dust_deprojections: return False
        return True

    # -----------------------------------------------------------------

    def set_dust_deprojection_aliases(self):

        """
        This function ...
        :return:
        """

        # Loop over the dust deprojection models
        for name in self.dust_deprojections:

            # Set the alias
            #self.dust_deprojections[name].filename = "dust___" + name + ".fits"
            self.dust_deprojections[name].filename = name

            # Set the map
            self.dust_deprojections[name].map = self.dust_maps[name]

    # -----------------------------------------------------------------

    def deproject_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the dust maps ...")

        # Create deproejcted dust maps
        if not self.has_dust_deprojected:
            self.dust_deprojections, self.dust_deprojected = self.deproject_maps(self.dust_maps, scale_height=self.dust_scaleheight, root_path=self.dust_deprojection_path, return_deprojections=True, downsample_factor=self.config.downsample_factor)
            #self.set_dust_deprojection_aliases()

        # Create deprojected dust maps with SKIRT
        if not self.has_dust_deprojected_skirt:
            if self.has_dust_deprojections:
                self.dust_deprojected_skirt = self.deproject_models(self.dust_deprojections, root_path=self.dust_deprojection_skirt_path, method="skirt", downsample_factor=self.config.downsample_factor)
            else:
                self.dust_deprojections, self.dust_deprojected_skirt = self.deproject_maps(self.dust_maps, scale_height=self.dust_scaleheight, root_path=self.dust_deprojection_skirt_path, method="skirt", return_deprojections=True, downsample_factor=self.config.downsample_factor)
                self.set_dust_deprojection_aliases()

    # -----------------------------------------------------------------

    @property
    def max_old_x_range(self):

        """
        This function ...
        :return:
        """

        x_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.old_deprojections: x_range.adjust(self.old_deprojections[name].x_range)
        return x_range

    # -----------------------------------------------------------------

    @property
    def max_old_y_range(self):

        """
        This function ...
        :return:
        """

        y_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.old_deprojections: y_range.adjust(self.old_deprojections[name].y_range)
        return y_range

    # -----------------------------------------------------------------

    @property
    def max_young_x_range(self):

        """
        This function ...
        :return:
        """

        x_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.young_deprojections: x_range.adjust(self.young_deprojections[name].y_range)
        return x_range

    # -----------------------------------------------------------------

    @property
    def max_young_y_range(self):

        """
        This function ...
        :return:
        """

        y_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.young_deprojections: y_range.adjust(self.young_deprojections[name].y_range)
        return y_range

    # -----------------------------------------------------------------

    @property
    def max_ionizing_x_range(self):

        """
        Thisn function ...
        :return:
        """

        x_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.ionizing_deprojections: x_range.adjust(self.ionizing_deprojections[name].x_range)
        return x_range

    # -----------------------------------------------------------------

    @property
    def max_ionizing_y_range(self):

        """
        This function ...
        :return:
        """

        y_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.ionizing_deprojections: y_range.adjust(self.ionizing_deprojections[name].y_range)
        return y_range

    # -----------------------------------------------------------------

    @property
    def max_dust_x_range(self):

        """
        This function ...
        :return: 
        """

        x_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.dust_deprojections: x_range.adjust(self.dust_deprojections[name].x_range)
        return x_range

    # -----------------------------------------------------------------

    @property
    def max_dust_y_range(self):


        """
        This function ...
        :return:
        """

        y_range = QuantityRange.infinitesimal(0.0 * u("pc"))
        for name in self.dust_deprojections: y_range.adjust(self.dust_deprojections[name].y_range)
        return y_range

    # -----------------------------------------------------------------

    def project(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the maps ...")

        # Old
        self.project_old()

        # Young
        self.project_young()

        # Ionizing
        self.project_ionizing()

        # Dust
        self.project_dust()

    # -----------------------------------------------------------------

    @property
    def has_old_projected(self):

        """
        This function ...
        :return:
        """

        # Loop over all old stellar maps
        for name in self.old_selection:

            # Determine path
            path = fs.join(self.old_edgeon_path, name + ".fits")

            # Check
            if not fs.is_file(path): return False

        # All checks passed
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def old_projection_npixels(self):

        """
        Thisnfunction ...
        :return:
        """

        radial_extent = max(self.max_old_x_range.span, self.max_old_y_range.span)
        vertical_extent = 2. * self.old_scaleheight * self.config.scale_heights

        # Determine number of pixels
        nx = int(round(radial_extent / self.old_projection_physical_pixelscale))
        nz = int(round(vertical_extent / self.old_projection_physical_pixelscale))

        # Return the pixel shape
        return PixelShape.from_xy(nx, nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_old_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None
        for name in self.old_maps:
            if pixelscale is None or self.old_maps[name].average_pixelscale < pixelscale: pixelscale = self.old_maps[name].average_pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def old_projection_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.min_old_pixelscale * self.config.downsample_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def old_projection_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        return (abs(self.old_projection_pixelscale) * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

    # -----------------------------------------------------------------

    def project_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the old stellar maps ...")

        # Check
        if not self.has_old_deprojections: raise RuntimeError("Deprojections for the old stellar maps are not created. Run again with --redeproject_old or --redeproject_skirt_old")

        # Create projected old stellar maps
        if not self.has_old_projected: self.old_edgeon_skirt = self.project_models_edgeon(self.old_deprojections, self.old_projection_npixels, self.old_projection_pixelscale, root_path=self.old_edgeon_path, map_aliases=self.old_maps)

    # -----------------------------------------------------------------

    @property
    def has_young_projected(self):
        
        """
        This function ...
        :return: 
        """
        
        # Loop over the young stellar maps
        for name in self.young_selection:
            
            # Determine path
            path = fs.join(self.young_edgeon_path, name + ".fits")
            
            # Check
            if not fs.is_file(path): return False
        
        # All checks passed
        return True
        
    # -----------------------------------------------------------------

    @lazyproperty
    def young_projection_npixels(self):

        """
        Thisnfunction ...
        :return:
        """

        radial_extent = max(self.max_young_x_range.span, self.max_young_y_range.span)
        z_radius = 2. * self.young_scaleheight * self.config.scale_heights

        # Determine number of pixels
        nx = int(round(radial_extent / self.young_projection_physical_pixelscale))
        nz = int(round(z_radius / self.young_projection_physical_pixelscale))

        # Return the pixel shape
        return PixelShape.from_xy(nx, nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_young_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None
        for name in self.young_maps:
            if pixelscale is None or self.young_maps[name].average_pixelscale < pixelscale: pixelscale = self.young_maps[name].average_pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def young_projection_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.min_young_pixelscale * self.config.downsample_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def young_projection_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        return (abs(self.young_projection_pixelscale) * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

    # -----------------------------------------------------------------

    def project_young(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the young stellar maps ...")

        # Check
        if not self.has_young_deprojections: raise RuntimeError("Deprojections for the young stellar maps are not created. Run again with --redeproject_young or --redeproject_skirt_young")

        # Create projected young stellar maps
        if not self.has_young_projected: self.young_edgeon_skirt = self.project_models_edgeon(self.young_deprojections, self.young_projection_npixels, self.young_projection_pixelscale, root_path=self.young_edgeon_path, map_aliases=self.young_maps)

    # -----------------------------------------------------------------

    @property
    def has_ionizing_projected(self):
        
        """
        This function ...
        :return: 
        """

        # Loop over the ionizing stellar maps
        for name in self.ionizing_selection:
            
            # Determine the path
            path = fs.join(self.ionizing_edgeon_path, name + ".fits")
            
            # Check
            if not fs.is_file(path): return False
            
        # All checks passed
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_projection_npixels(self):

        """
        Thisnfunction ...
        :return:
        """

        radial_extent = max(self.max_ionizing_x_range.span, self.max_ionizing_y_range.span)
        z_radius = 2. * self.ionizing_scaleheight * self.config.scale_heights

        # Determine number of pixels
        nx = int(round(radial_extent / self.ionizing_projection_physical_pixelscale))
        nz = int(round(z_radius / self.ionizing_projection_physical_pixelscale))

        # Return the pixel shape
        return PixelShape.from_xy(nx, nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_ionizing_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None
        for name in self.ionizing_maps:
            if pixelscale is None or self.ionizing_maps[name].average_pixelscale < pixelscale: pixelscale = self.ionizing_maps[name].average_pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_projection_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.min_ionizing_pixelscale * self.config.downsample_factor

    # -----------------------------------------------------------------

    @lazyproperty
    def ionizing_projection_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        return (abs(self.ionizing_projection_pixelscale) * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

    # -----------------------------------------------------------------

    def project_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the ionizing stellar maps ...")

        # Check
        if not self.has_ionizing_deprojections: raise RuntimeError("Deprojections for the ionizing stellar maps are not created. Run again with --redeproject_ionizing or --redeproject_skirt_ionizing")

        # Create projected ionizing stellar maps
        if not self.has_ionizing_projected: self.ionizing_edgeon_skirt = self.project_models_edgeon(self.ionizing_deprojections, self.ionizing_projection_npixels, self.ionizing_projection_pixelscale, root_path=self.ionizing_edgeon_path, map_aliases=self.ionizing_maps)

    # -----------------------------------------------------------------

    @property
    def has_dust_projected(self):
        
        """
        This function ...
        :return: 
        """
        
        # Loop over the dust maps
        for name in self.dust_selection:
            
            # Determine the path
            path = fs.join(self.dust_edgeon_path, name + ".fits")
            
            # Check
            if not fs.is_file(path): return False
            
        # All checks passed
        return True

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_projection_npixels(self):

        """
        This function ...
        :return:
        """

        radial_extent = max(self.max_dust_x_range.span, self.max_dust_y_range.span)
        z_radius = 2. * self.dust_scaleheight * self.config.scale_heights

        # Determine number of pixels
        nx = int(round(radial_extent / self.dust_projection_physical_pixelscale))
        nz = int(round(z_radius / self.dust_projection_physical_pixelscale))

        # Return the pixel shape
        return PixelShape.from_xy(nx, nz)

    # -----------------------------------------------------------------

    @lazyproperty
    def min_dust_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None
        for name in self.dust_maps:
            if pixelscale is None or self.dust_maps[name].average_pixelscale < pixelscale: pixelscale = self.dust_maps[name].average_pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_projection_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.min_dust_pixelscale * self.config.downsample_factor

    # -----------------------------------------------------------------

    @property
    def dust_projection_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        return (abs(self.dust_projection_pixelscale) * self.galaxy_distance).to("pc", equivalencies=dimensionless_angles())

    # -----------------------------------------------------------------

    def project_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Projecting the dust maps ...")

        # Check
        if not self.has_dust_deprojections: raise RuntimeError("Deprojections for the dust maps are not created. Run again with --redeproject_dust or --redeproject_skirt_dust")

        # Create projected dust maps
        if not self.has_dust_projected: self.dust_edgeon_skirt = self.project_models_edgeon(self.dust_deprojections, self.dust_projection_npixels, self.dust_projection_pixelscale, root_path=self.dust_edgeon_path, map_aliases=self.dust_maps)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the maps
        self.write_maps()

        # Write the clip masks
        self.write_clip_masks()

        # Write the softening masks
        self.write_softening_masks()

        # Write the deprojected maps
        self.write_deprojected()

        # Write deprojected maps with SKIRT
        self.write_deprojected_skirt()

        # Write edgeon maps
        self.write_edgeon()

    # -----------------------------------------------------------------

    def write_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps ...")

        # Old
        self.write_old_maps()

        # Young
        self.write_young_maps()

        # Ionizing
        self.write_ionizing_maps()

        # Dust
        self.write_dust_maps()

    # -----------------------------------------------------------------

    def write_old_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the old stellar maps ...")

        # Loop over the maps
        for name in self.old_maps:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar map ...")

            # Determine the path
            path = fs.join(self.old_component_maps_path, name + ".fits")

            # Save the map
            self.old_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the young stellar maps ...")

        # Loop over the maps
        for name in self.young_maps:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar map ...")

            # Determine the path
            path = fs.join(self.young_component_maps_path, name + ".fits")

            # Save the map
            self.young_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_maps:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar map ...")

            # Determine the path
            path = fs.join(self.ionizing_component_maps_path, name + ".fits")

            # Save the map
            self.ionizing_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the dust maps ...")

        # Loop over the maps
        for name in self.dust_maps:

            # Debugging
            log.debug("Writing the '" + name + "' dust map ...")

            # Determine the path
            path = fs.join(self.dust_component_maps_path, name + ".fits")

            # Save the map
            self.dust_maps[name].saveto(path)

    # -----------------------------------------------------------------

    def write_clip_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks ...")

        # Old
        self.write_old_clip_masks()

        # Young
        self.write_young_clip_masks()

        # Ionizing
        self.write_ionizing_clip_masks()

        # Dust
        self.write_dust_clip_masks()

    # -----------------------------------------------------------------

    def write_old_clip_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks of the old stellar maps ...")

        # Loop over the masks
        for name in self.old_clip_masks:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar mask ...")

            # Determine the path
            path = fs.join(self.old_masks_path, name + "_clip.fits")

            # Write
            self.old_clip_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_clip_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks of the young stellar maps ...")

        # Loop over the masks
        for name in self.young_clip_masks:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar mask ...")

            # Determine the path
            path = fs.join(self.young_masks_path, name + "_clip.fits")

            # Write
            self.young_clip_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_clip_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks of the ionizing stellar maps ...")

        # Loop over the masks
        for name in self.ionizing_clip_masks:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar mask ...")

            # Determine the path
            path = fs.join(self.ionizing_masks_path, name + "_clip.fits")

            # Write
            self.ionizing_clip_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_clip_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks of the dust maps ...")

        # Loop over the masks
        for name in self.dust_clip_masks:

            # Debugging
            log.debug("Writing the '" + name + "' dust mask ...")

            # Determine the path
            path = fs.join(self.dust_masks_path, name + "_clip.fits")

            # Write
            self.dust_clip_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_softening_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the softening masks ...")

        # Old
        self.write_old_softening_masks()

        # Young
        self.write_young_softening_masks()

        # Ionizing
        self.write_ionizing_softening_masks()

        # Dust
        self.write_dust_softening_masks()

    # -----------------------------------------------------------------

    def write_old_softening_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the softening masks of the old stellar maps ...")

        # Loop over the masks
        for name in self.old_softening_masks:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar mask ...")

            # Determine the path
            path = fs.join(self.old_masks_path, name + "_softening.fits")

            # Write
            self.old_softening_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_softening_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the softening masks of the young stellar maps ...")

        # Loop over the masks
        for name in self.young_softening_masks:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar mask ...")

            # Determine the path
            path = fs.join(self.young_masks_path, name + "_softening.fits")

            # Write
            self.young_softening_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_softening_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the softening masks of the ionizing stellar maps ...")

        # Loop over the masks
        for name in self.ionizing_softening_masks:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar mask ...")

            # Determine the path
            path = fs.join(self.ionizing_masks_path, name + "_softening.fits")

            # Write
            self.ionizing_softening_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_softening_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the softening masks of the dust maps ...")

        # Loop over the masks
        for name in self.dust_softening_masks:

            # Debugging
            log.debug("Writing the '" + name + "' dust mask ...")

            # Determine the filepath
            path = fs.join(self.dust_masks_path, name + "_softening.fits")

            # Write
            self.dust_softening_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected maps ...")

        # Old
        self.write_old_deprojected()

        # Young
        self.write_young_deprojected()

        # Ionizing
        self.write_ionizing_deprojected()

        # Dust
        self.write_dust_deprojected()

    # -----------------------------------------------------------------

    def write_old_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected old stellar maps ...")

        # Loop over the maps
        for name in self.old_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected old stellar map ...")

            # Determine path
            path = fs.join(self.old_deprojection_path, name + ".fits")

            # Write
            self.old_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected young stellar maps ...")

        # Loop over the maps
        for name in self.young_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected young stellar map ...")

            # Determine path
            path = fs.join(self.young_deprojection_path, name + ".fits")

            # Write
            self.young_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected ionizing stellar maps ...")

        # Loop over the maps
        for name in self.ionizing_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected ionizing stellar map ...")

            # Determine the path
            path = fs.join(self.ionizing_deprojection_path, name + ".fits")

            # Write
            self.ionizing_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected dust maps ...")

        # Loop over the maps
        for name in self.dust_deprojected:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected dust map ...")

            # Determine the path
            path = fs.join(self.dust_deprojection_path, name + ".fits")

            # Write
            self.dust_deprojected[name].saveto(path)

    # -----------------------------------------------------------------

    def write_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the maps deprojected with SKIRT ...")

        # Old
        self.write_old_deprojected_skirt()

        # Young
        self.write_young_deprojected_skirt()

        # Ionizing
        self.write_ionizing_deprojected_skirt()

        # Dust
        self.write_dust_deprojected_skirt()

    # -----------------------------------------------------------------

    def write_old_deprojected_skirt(self):

        """
        Thisf unctino ...
        :return:
        """

        # Inform the user
        log.info("Writing old stellar maps deprojected with SKIRT ...")

        # Loop over the maps
        for name in self.old_deprojected_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected with SKIRT old stellar map ...")

            # Determine the path
            path = fs.join(self.old_deprojection_skirt_path, name + ".fits")

            # Write
            self.old_deprojected_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing young stellar maps deprojected with SKIRT ...")

        # Loop over the maps
        for name in self.young_deprojected_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected with SKIRT young stellar map ...")

            # Determine the path
            path = fs.join(self.young_deprojection_skirt_path, name + ".fits")

            # Write
            self.young_deprojected_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ionizing stellar maps deprojected with SKIRT ...")

        # Loop over the maps
        for name in self.ionizing_deprojected_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected with SKIRT ionizing stellar map ...")

            # Determine the path
            path = fs.join(self.ionizing_deprojection_skirt_path, name + ".fits")

            # Write
            self.ionizing_deprojected_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_deprojected_skirt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing dust maps deprojected with SKIRT ...")

        # Loop over the maps
        for name in self.dust_deprojected_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' deprojected with SKIRT dust map ...")

            # Determine the path
            path = fs.join(self.dust_deprojection_skirt_path, name + ".fits")

            # Write
            self.dust_deprojected_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the edgeon maps ...")

        # Old
        self.write_old_edgeon()

        # Young
        self.write_young_edgeon()

        # Ionizing
        self.write_ionizing_edgeon()

        # Dust
        self.write_dust_edgeon()

    # -----------------------------------------------------------------

    def write_old_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing old stellar edgeon maps ...")

        # Loop over the maps
        for name in self.old_edgeon_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' edgeon map ...")

            # Determine the path
            path = fs.join(self.old_edgeon_path, name + ".fits")

            # Write
            self.old_edgeon_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing young stellar edgeon maps ...")

        # Loop over the maps
        for name in self.young_edgeon_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' edgeon map ...")

            # Determine the path
            path = fs.join(self.young_edgeon_path, name + ".fits")

            # Write
            self.young_edgeon_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ionizing stellar edgeon maps ...")

        # Loop over the maps
        for name in self.ionizing_edgeon_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' edgeon map ...")

            # Determine the path
            path = fs.join(self.ionizing_edgeon_path, name + ".fits")

            # Write
            self.ionizing_edgeon_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_edgeon(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing dust edgeon maps ...")

        # Loop over the maps
        for name in self.dust_edgeon_skirt:

            # Debugging
            log.debug("Writing the '" + name + "' edgeon map ...")

            # Determine the path
            path = fs.join(self.dust_edgeon_path, name + ".fits")

            # Write
            self.dust_edgeon_skirt[name].saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Maps
        self.plot_components_maps(format=self.config.plotting_format)

        # Radial profiles
        self.plot_components_radial_profiles(format=self.config.plotting_format)

        # Masks
        self.plot_components_masks(format=self.config.plotting_format)

        # Deprojected
        self.plot_components_deprojected(format=self.config.plotting_format)

        # Deprojected with SKIRT
        self.plot_components_deprojected_skirt(format=self.config.plotting_format)

        # Edgeon
        self.plot_components_edgeon(format=self.config.plotting_format)

    # -----------------------------------------------------------------

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

# -----------------------------------------------------------------

def find_matching_colour_map_name(colour, map_names):

    """
    This function ...
    :param colour:
    :param map_names:
    :return:
    """

    # Loop over the map names
    for name in map_names:

        # Select
        if colours.same_colour(colour, name): return name

    # Nothing found
    return None

# -----------------------------------------------------------------

def find_factor_max_nnegatives(nnegatives, max_nnegatives):

    """
    This function ...
    :param nnegatives:
    :param max_nnegatives:
    :return:
    """

    # Loop over the factors in reverse order
    for factor in reversed(nnegatives.keys()):

        # Get the number of negatives
        negatives = nnegatives[factor]

        # Succes?
        if negatives < max_nnegatives: return factor

    # Error
    raise ValueError("None of the maps have a relative number of negatives lower than the limit of " + str(max_nnegatives*100) + "%")

# -----------------------------------------------------------------

def interpolate_map(the_map, method, region_or_mask, smoothing_factor=None, source_outer_factor=1.4, sigma_clip=True,
                    interpolation_softening_range=None, plot=False, in_cutout=True, relative_softening_radius=None,
                    return_background=False):

    """
    This function ...
    :param the_map:
    :param method:
    :param region_or_mask:
    :param smoothing_factor:
    :param source_outer_factor:
    :param sigma_clip:
    :param interpolation_softening_range:
    :param plot:
    :param in_cutout:
    :param relative_softening_radius:
    :param return_background:
    :return:
    """

    # Debuggging
    log.debug("Interpolating the map ...")

    # In cutout?
    if in_cutout: return interpolate_map_in_cutout(the_map, method, region_or_mask, source_outer_factor=source_outer_factor,
                                            smoothing_factor=smoothing_factor, sigma_clip=sigma_clip,
                                            interpolation_softening_range=interpolation_softening_range, plot=plot,
                                            relative_softening_radius=relative_softening_radius, return_background=return_background)

    # Not in cutout
    else:

        # Warnings
        log.warning("When 'in_cutout' is disabled, interpolation may be slower")
        log.warning("When 'in_cutout' is disabled, interpolation region edges cannot be softened")

        # INTERPOLATE IN FRAME
        return interpolate_map_in_frame(the_map, method, region_or_mask, smoothing_factor=smoothing_factor, plot=plot, return_background=return_background)

# -----------------------------------------------------------------

def interpolate_map_in_frame(the_map, method, region_or_mask, smoothing_factor=None, plot=False, return_background=False):

    """
    Thisf unction ...
    :param the_map:
    :param method:
    :param region_or_mask:
    :param smoothing_factor:
    :param plot:
    :param return_background:
    :return:
    """

    # Inform the user
    log.info("Interpolating the map as a frame ...")

    # Check interpolation method
    if method != "kernel": raise ValueError("The interpolation method can only be 'kernel'")

    # Plot?
    if plot: plotting.plot_frame(the_map, title="before interpolation")

    # Get the map and interpolate in the ellipse region
    mask = the_map.interpolate(region_or_mask, max_iterations=None, smoothing_factor=smoothing_factor)

    # Plot mask
    if plot: plotting.plot_mask(mask, title="interpolation mask")

    # Plot?
    if plot: plotting.plot_frame(the_map, title="after interpolation")

    # Return mask and background
    if return_background: return mask, None

    # Return the mask
    else: return mask

# -----------------------------------------------------------------

def interpolate_map_in_cutout(the_map, method, region_or_mask, source_outer_factor=1.4, smoothing_factor=None,
                              sigma_clip=True, interpolation_softening_range=None, plot=False,
                              relative_softening_radius=None, return_background=False):

    """
    This function ...
    :param the_map:
    :param method:
    :param region_or_mask:
    :param source_outer_factor:
    :param smoothing_factor:
    :param sigma_clip:
    :param interpolation_softening_range:
    :param plot:
    :param relative_softening_radius:
    :param return_background:
    :return:
    """

    # Debugging
    log.debug("Interpolating the map within a cutout ...")

    # Create source
    if isinstance(region_or_mask, MaskBase): source = Detection.from_mask(the_map, region_or_mask, source_outer_factor)
    elif isinstance(region_or_mask, PixelRegion): source = Detection.from_shape(the_map, region_or_mask, source_outer_factor)
    else: raise ValueError("Invalid value for 'region_or_mask'")

    # Plot the source
    if plot: source.plot(title="before background estimation")

    # Determine the FWHM for interpolation
    fwhm = the_map.fwhm_pix
    if smoothing_factor is not None:
        if method != "kernel": raise ValueError("Cannot specify smoothing factor for interpolation method other than 'kernel'")
        fwhm = fwhm * smoothing_factor

    # Debugging
    log.debug("The interpolation method is '" + method + "'")

    # Estimate the background, specifying the FWHM for 'kernel' method
    interpolation_mask = source.estimate_background(method, sigma_clip=sigma_clip, fwhm=fwhm, plot_interpolation=plot)

    # Plot the interpolation mask
    if plot: plotting.plot_mask(interpolation_mask, title="interpolation mask")

    # Plot the background
    if plot: plotting.plot_box(source.background, title="estimated background")

    # Plot
    if plot: source.plot(title="after background estimation")

    # Replace with background but with soft edge
    if interpolation_softening_range is not None:

        # Initialize an alpha mask for the pixels that will be replaced by the estimated background
        mask = AlphaMask.empty_like(the_map)

        # Debugging
        log.debug("Creating alpha mask ...")

        # A mask is passed: create alpha mask from this mask by softening instead of taking the source's shape
        if isinstance(region_or_mask, MaskBase):

            # Create alpha mask by softening the edges of the passed mask
            # min_radius = interpolation_softening_range.min * source.shape.semimajor
            # max_radius = interpolation_softening_range.max * source.shape.semimajor
            # radius = max_radius - min_radius
            # start = min_radius - radius
            # nbins = 10
            #
            # # Debugging
            # log.debug("Creating alpha mask ...")
            # log.debug("Min radius: " + str(min_radius))
            # log.debug("Max radius: " + str(max_radius))
            # log.debug("Radius: " + str(radius))
            # log.debug("Start: " + str(start))

            softening_radius = relative_softening_radius * source.xsize
            nbins = 10
            start = 0

            # Debugging

            log.debug("The softening radius is " + str(softening_radius) + " pixels")

            # Get the source mask
            source_mask = Mask(source.mask)

            # Create softened mask
            alpha_mask = source_mask.softened(radius=softening_radius, min_radius=start, nbins=nbins, per_level=False)

        # A pixel region is passed
        elif isinstance(region_or_mask, PixelRegion):

            # Create alpha mask from sources's shape, which would be the shifted version of the passed region
            region = source.shape
            shape = (source.ysize, source.xsize)
            alpha_mask = AlphaMask.from_ellipse(region, shape, interpolation_softening_range, dynamic_max=True)

        # Invalid value
        else: raise ValueError("Invalid value for 'region_or_mask'")

        # Plot the alpha mask
        if plot: plotting.plot_mask(alpha_mask, title="softening mask")

        # Replace the pixels by the background
        source.background.replace(the_map, where=alpha_mask)

        # Set alpha mask pixels
        mask.data[source.y_slice, source.x_slice] = alpha_mask.data

    # Replace in regular mask patch
    else:

        # Replace
        source.background.replace(the_map)

        # Initialize a mask for the pixels that will be replaced by the estimated background
        mask = Mask.empty_like(the_map)
        mask.data[source.y_slice, source.x_slice] = source.mask

    # Plot
    if plot: plotting.plot_box(the_map[source.y_slice, source.x_slice], title="after interpolation")

    # Return mask and background
    if return_background:

        # Create background frame
        background = Frame.zeros_like(the_map)
        background[source.y_slice, source.x_slice] = source.background

        # Return mask and backgroud
        return mask, background

    # Return the mask in which the interpolation happened
    else: return mask

# -----------------------------------------------------------------
