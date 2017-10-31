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
from ...magic.core.mask import Mask
from ...magic.core.alpha import AlphaMask
from ...magic.core.detection import Detection
from ...core.tools.stringify import tostr
from ...core.basics.range import QuantityRange
from ...core.units.parsing import parse_unit as u
from ...magic.basics.vector import PixelShape

# -----------------------------------------------------------------

steps_name = "steps"
masks_name = "masks"
deprojected_name = "deprojected"
deprojected_skirt_name = "deprojected_skirt"
edgeon_name = "edgeon"

# -----------------------------------------------------------------

correct_step = "corrected"
interpolate_step = "interpolated"
truncate_step = "truncated"
crop_step = "cropped"
clip_step = "clipped"
softened_step = "softened"
steps = [correct_step, interpolate_step, truncate_step, crop_step, clip_step, softened_step]

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

        # The deprojected maps
        self.old_deprojected = dict()
        self.young_deprojected = dict()
        self.ionizing_deprojected = dict()
        self.dust_deprojected = dict()

        # The deprojections
        self.old_deprojections = None
        self.young_deprojections = None
        self.ionizing_deprojections = None
        self.dust_deprojections = None

        # The deprojected maps with SKIRT
        self.old_deprojected_skirt = dict()
        self.young_deprojected_skirt = dict()
        self.ionizing_deprojected_skirt = dict()
        self.dust_deprojected_skirt = dict()

        # Edgeon maps with SKIRT
        self.old_edgeon_skirt = dict()
        self.young_edgeon_skirt = dict()
        self.ionizing_edgeon_skirt = dict()
        self.dust_edgeon_skirt = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # Auto-select
        if self.config.auto: self.auto_select()

        # 2. Prompt
        self.prompt()

        # 3. Set rerun
        if self.rerun: self.set_rerun()

        # 4. Set redeproject
        if self.redeproject: self.set_redeproject()

        # 5. Set redeproject SKIRT
        if self.redeproject_skirt: self.set_redeproject_skirt()

        # Set reproject
        if self.reproject: self.set_reproject()

        # 6. Remove other
        if self.remove: self.remove_other()

        # 7. Load the maps
        self.load_maps()

        # 8. Load the masks
        self.load_masks()

        # 9. Process the maps
        self.process_maps()

        # 10. Deproject the maps
        self.deproject()

        # 11. Project the maps to the edge-on view
        self.project()

        # 12. Writing
        self.write()

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
        if self.config.clear_all:
            fs.clear_directory(self.old_component_maps_path)
            fs.clear_directory(self.young_component_maps_path)
            fs.clear_directory(self.ionizing_component_maps_path)
            fs.clear_directory(self.dust_component_maps_path)

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

        # Set random
        if self.config.random: self.config.random_old = self.config.random_young = self.config.random_ionizing = self.config.random_dust = self.config.random

        # Set all
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # Make selections
        self.old_selection = sequences.make_selection(self.old_map_names, self.config.old, self.config.not_old, nrandom=self.config.random_old, all=self.config.all_old, indices=self.config.old_indices, not_indices=self.config.not_old_indices)
        self.young_selection = sequences.make_selection(self.young_map_names, self.config.young, self.config.not_young, nrandom=self.config.random_young, all=self.config.all_young, indices=self.config.young_indices, not_indices=self.config.not_young_indices)
        self.ionizing_selection = sequences.make_selection(self.ionizing_map_names, self.config.ionizing, self.config.not_ionizing, nrandom=self.config.random_ionizing, all=self.config.all_ionizing, indices=self.config.ionizing_indices, not_indices=self.config.not_ionizing_indices)
        self.dust_selection = sequences.make_selection(self.dust_map_names, self.config.dust, self.config.not_dust, nrandom=self.config.random_dust, all=self.config.all_dust, indices=self.config.dust_indices, not_indices=self.config.not_dust_indices)

        # Levels
        if self.config.levels is not None: self.levels = self.config.levels

        # Set rerun options
        if self.config.rerun is not None:
            self.config.rerun_old = self.config.rerun
            self.config.rerun_young = self.config.rerun
            self.config.rerun_ionizing = self.config.rerun
            self.config.rerun_dust = self.config.rerun

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
        log.info("Automatically making maps selections ...")

        # Select
        self.has_ssfr_maps

        print(self.old_selection)
        print(self.young_selection)
        print(self.ionizing_selection)
        print(self.dust_selection)

    # -----------------------------------------------------------------

    def prompt(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Prompting for user input ...")

        # Old
        if not self.has_old_selection: self.prompt_old()

        # Young
        if not self.has_young_selection: self.prompt_young()

        # Ionizing
        if not self.has_ionizing_selection: self.prompt_ionizing()

        # Dust
        if not self.has_dust_selection: self.prompt_dust()

        # Levels
        if not self.has_levels: self.prompt_levels()
        else: self.check_levels()

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

            # Prompt for the level
            if not self.config.all_levels: level = prompt_real(id + "_level", "sigma level for the " + filter_name + " image", default=self.config.default_level)
            else: level = self.config.default_level

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

        map_path = fs.join(self.old_steps_path, name)
        if not fs.is_directory(map_path) or fs.is_empty(map_path): return None

        for step in reversed(steps):
            filepath = fs.join(map_path, step + ".fits")
            if fs.is_file(filepath): return step, filepath

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    def get_last_young_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.young_steps_path, name)
        if not fs.is_directory(map_path) or fs.is_empty(map_path): return None

        for step in reversed(steps):
            filepath = fs.join(map_path, step + ".fits")
            if fs.is_file(filepath): return step, filepath

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    def get_last_ionizing_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.ionizing_steps_path, name)
        if not fs.is_directory(map_path) or fs.is_empty(map_path): return None

        for step in reversed(steps):
            filepath = fs.join(map_path, step + ".fits")
            if fs.is_file(filepath): return step, filepath

        raise RuntimeError("We shouldn't get here")

    # -----------------------------------------------------------------

    def get_last_dust_step(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        map_path = fs.join(self.dust_steps_path, name)
        if not fs.is_directory(map_path) or fs.is_empty(map_path): return None

        for step in reversed(steps):
            filepath = fs.join(map_path, step + ".fits")
            if fs.is_file(filepath): return step, filepath

        raise RuntimeError("We shouldn't get here")

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
    def rerun_old(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_old is not None

    # -----------------------------------------------------------------

    @property
    def rerun_young(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_young is not None

    # -----------------------------------------------------------------

    @property
    def rerun_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_ionizing is not None

    # -----------------------------------------------------------------

    @property
    def rerun_dust(self):

        """
        This function ...
        :return:
        """

        return self.config.rerun_dust is not None

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

                # Inform
                log.success("Found a " + step + " '" + name + "' map")

                # Load
                self.old_maps[name] = Frame.from_file(path)

                # Set metadata
                for step in steps_before_and_including(step): self.old_maps[name].metadata[step] = True
                for step in steps_after(step): self.old_maps[name].metadata[step] = False

            # Load raw map
            else:

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

                # Inform
                log.success("Found a " + step + " '" + name + "' map")

                # Load
                self.young_maps[name] = Frame.from_file(path)

                # Set metadata
                for step in steps_before_and_including(step): self.young_maps[name].metadata[step] = True
                for step in steps_after(step): self.young_maps[name].metadata[step] = False

            # Load raw map
            else:

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

                # Inform
                log.success("Found a " + step + " '" + name + "' map")

                # Load
                self.ionizing_maps[name] = Frame.from_file(path)

                # Set metadata
                for step in steps_before_and_including(step): self.ionizing_maps[name].metadata[step] = True
                for step in steps_after(step): self.ionizing_maps[name].metadata[step] = False

            # Load raw map
            else:

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

                # Inform
                log.success("Found a " + step + " '" + name + "' map")

                # Load
                self.dust_maps[name] = Frame.from_file(path)

                # Set metadata
                for step in steps_before_and_including(step): self.dust_maps[name].metadata[step] = True
                for step in steps_after(step): self.dust_maps[name].metadata[step] = False

            # Load raw map
            else:

                # Load
                self.dust_maps[name] = Frame.from_file(self.dust_map_paths[name])

                # Set metadata
                for step in steps: self.dust_maps[name].metadata[step] = False

    # -----------------------------------------------------------------

    def old_step_path_for_map(self, name, step):

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

    def young_step_path_for_map(self, name, step):

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

    def ionizing_step_path_for_map(self, name, step):

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

    def dust_step_path_for_map(self, name, step):

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
                self.old_clip_masks[name] = Mask.from_file(path)

            # Check whether the softening step output is present
            if self.has_old_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.old_step_path_for_mask(name, softened_step)

                # Load the mask
                self.old_softening_masks[name] = Mask.from_file(path)

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
                self.young_clip_masks[name] = Mask.from_file(path)

            # Check whether the softening step output is present
            if self.has_young_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.young_step_path_for_mask(name, softened_step)

                # Load the mask
                self.young_softening_masks[name] = Mask.from_file(path)

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
                self.ionizing_clip_masks[name] = Mask.from_file(path)

            # Check whether the softening step output is present
            if self.has_ionizing_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.ionizing_step_path_for_mask(name, softened_step)

                # Load the mask
                self.ionizing_softening_masks[name] = Mask.from_file(path)

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
                self.dust_clip_masks[name] = Mask.from_file(path)

            # Check whether the softening step output is present
            if self.has_dust_step(name, softened_step):

                # Inform
                log.success("Found a " + softened_step + " mask for the '" + name + "' map")

                # Determine mask path
                path = self.dust_step_path_for_mask(name, softened_step)

                # Load the mask
                self.dust_softening_masks[name] = Mask.from_file(path)

    # -----------------------------------------------------------------

    def process_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Processing the maps ...")

        # 1. Correct
        self.correct_maps()

        # 2. Interpolate the cores
        self.interpolate_maps()

        # 2. Truncate
        self.truncate_maps()

        # 3. Crop
        self.crop_maps()

        # 4. Clip
        self.clip_maps()

        # 5. Soften the edges
        self.soften_edges()

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
            self.correct_map(self.old_maps[name])

            # Set flag
            self.old_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, correct_step))

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
            self.correct_map(self.young_maps[name])

            # Set flag
            self.young_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, correct_step))

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
            self.correct_map(self.ionizing_maps[name])

            # Set flag
            self.ionizing_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, correct_step))

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
            self.correct_map(self.dust_maps[name])

            # Set flag
            self.dust_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.step: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, correct_step))

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
        self.interpolate_old_maps()

        # Young
        self.interpolate_young_maps()

        # Ionizing
        self.interpolate_ionizing_maps()

        # Dust
        self.interpolate_dust_maps()

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

        ellipse = self.truncation_ellipse * self.config.interpolate_old
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

            # No interpolation?
            if self.config.interpolate_old is None:
                self.old_maps[name].metadata[interpolate_step] = True
                continue

            # Check
            if self.is_interpolated_old(name):
                log.success("The '" + name + "' old stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' old stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_old.to_pixel(self.old_maps[name].wcs)

            # Create a source
            source = Detection.from_shape(self.old_maps[name], ellipse, self.config.source_outer_factor)

            # Estimate the background
            source.estimate_background(self.config.interpolation_method, sigma_clip=self.config.sigma_clip)

            # Create alpha mask
            alpha_mask = AlphaMask.from_ellipse(ellipse - source.shift, (source.ysize, source.xsize), self.interpolation_softening_range, wcs=self.old_maps[name])

            # Replace the pixels by the background
            source.background.replace(self.old_maps[name], where=alpha_mask)

            # Set flag
            self.old_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, interpolate_step))

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_young(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.interpolate_young
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

            # No interpolation?
            if self.config.interpolate_young is None:
                self.young_maps[name].metadata[interpolate_step] = True
                continue

            # Check
            if self.is_interpolated_young(name):
                log.success("The '" + name + "' young stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' young stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_young.to_pixel(self.young_maps[name].wcs)

            # Create a source
            source = Detection.from_shape(self.young_maps[name], ellipse, self.config.source_outer_factor)

            # Estimate the background
            source.estimate_background(self.config.interpolation_method, sigma_clip=self.config.sigma_clip)

            # Create alpha mask
            alpha_mask = AlphaMask.from_ellipse(ellipse - source.shift, (source.ysize, source.xsize), self.interpolation_softening_range, wcs=self.old_maps[name])

            # Replace the pixels by the background
            source.background.replace(self.young_maps[name], where=alpha_mask)

            # Set flag
            self.young_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, interpolate_step))

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_ionizing(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.interpolate_ionizing
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

            # No interpolation?
            if self.config.interpolate_ionizing is None:
                self.ionizing_maps[name].metadata[interpolate_step] = True
                continue

            # Check
            if self.is_interpolated_ionizing(name):
                log.success("The '" + name + "' ionizing stellar map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' ionizing stellar map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_ionizing.to_pixel(self.ionizing_maps[name].wcs)

            # Create a source
            source = Detection.from_shape(self.ionizing_maps[name], ellipse, self.config.source_outer_factor)

            # Estimate the background
            source.estimate_background(self.config.interpolation_method, sigma_clip=self.config.sigma_clip)

            # Create alpha mask
            alpha_mask = AlphaMask.from_ellipse(ellipse - source.shift, (source.ysize, source.xsize), self.interpolation_softening_range, wcs=self.old_maps[name])

            # Replace the pixels by the background
            source.background.replace(self.ionizing_maps[name], where=alpha_mask)

            # Set flag
            self.ionizing_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, interpolate_step))

    # -----------------------------------------------------------------

    @lazyproperty
    def interpolation_ellipse_dust(self):

        """
        This function ...
        :return:
        """

        ellipse = self.truncation_ellipse * self.config.interpolate_dust
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

            # No interpolation?
            if self.config.interpolate_dust is None:
                self.dust_maps[name].metadata[interpolate_step] = True
                continue

            # Check
            if self.is_interpolated_dust(name):
                log.success("The '" + name + "' dust map is already interpolated")
                continue

            # Debugging
            log.debug("Interpolating the '" + name + "' dust map ...")

            # Create interpolation ellipse
            ellipse = self.interpolation_ellipse_dust.to_pixel(self.dust_maps[name].wcs)

            # Create a source
            source = Detection.from_shape(self.dust_maps[name], ellipse, self.config.source_outer_factor)

            # Estimate the background
            source.estimate_background(self.config.interpolation_method, sigma_clip=self.config.sigma_clip)

            # Create alpha mask
            alpha_mask = AlphaMask.from_ellipse(ellipse - source.shift, (source.ysize, source.xsize), self.interpolation_softening_range, wcs=self.old_maps[name])

            # Replace the pixels by the background
            source.background.replace(self.dust_maps[name], where=alpha_mask)

            # Set flag
            self.dust_maps[name].metadata[interpolate_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, interpolate_step))

    # -----------------------------------------------------------------

    def truncate_maps(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Truncating the maps ...")

        # Old
        self.truncate_old_maps()

        # Young
        self.truncate_young_maps()

        # Ionizing
        self.truncate_ionizing_maps()

        # Dust
        self.truncate_dust_maps()

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
            self.truncate_map(self.old_maps[name])

            # Set flag
            self.old_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, truncate_step))

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
            self.truncate_map(self.young_maps[name])

            # Set flag
            self.young_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, truncate_step))

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
            self.truncate_map(self.ionizing_maps[name])

            # Set flag
            self.ionizing_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, truncate_step))

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
            self.truncate_map(self.dust_maps[name])

            # Set flag
            self.dust_maps[name].metadata[truncate_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, truncate_step))

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
        self.crop_old_maps()

        # Young
        self.crop_young_maps()

        # Ionizing
        self.crop_ionizing_maps()

        # Dust
        self.crop_dust_maps()

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
            self.crop_map(self.old_maps[name], self.config.cropping_factor)

            # Set flag
            self.old_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, crop_step))

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
            self.crop_map(self.young_maps[name], self.config.cropping_factor)

            # Set flag
            self.young_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, crop_step))

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
            self.crop_map(self.ionizing_maps[name], self.config.cropping_factor)

            # Set flag
            self.ionizing_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, crop_step))

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
            self.crop_map(self.dust_maps[name], self.config.cropping_factor)

            # Set flag
            self.dust_maps[name].metadata[crop_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, crop_step))

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
        self.clip_old_maps()

        # Young
        self.clip_young_maps()

        # Ionizing
        self.clip_ionizing_maps()

        # Dust
        self.clip_dust_maps()

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

            # Clip the map (returns the mask)
            self.old_clip_masks[name] = self.clip_map(self.old_maps[name], origins, convolve=self.config.convolve,
                                                 remote=self.remote, npixels=self.config.min_npixels,
                                                 connectivity=self.config.connectivity,
                                                 rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                 fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                 fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.old_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.old_clip_masks[name].saveto(self.old_step_path_for_mask(name, clip_step))

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

            # Clip the map (returns the mask)
            self.young_clip_masks[name] = self.clip_map(self.young_maps[name], origins, convolve=self.config.convolve,
                                                   remote=self.remote, npixels=self.config.min_npixels,
                                                   connectivity=self.config.connectivity,
                                                   rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                   fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                   fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.young_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.young_clip_masks[name].saveto(self.young_step_path_for_mask(name, clip_step))

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

            # Clip the map (returns the mask)
            self.ionizing_clip_masks[name] = self.clip_map(self.ionizing_maps[name], origins, convolve=self.config.convolve,
                                                      remote=self.remote, npixels=self.config.min_npixels,
                                                      connectivity=self.config.connectivity,
                                                      rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                      fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                      fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.ionizing_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.ionizing_clip_masks[name].saveto(self.ionizing_step_path_for_mask(name, clip_step))

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

            # Clip the map (returns the mask)
            self.dust_clip_masks[name] = self.clip_map(self.dust_maps[name], origins, convolve=self.config.convolve,
                                                  remote=self.remote, npixels=self.config.min_npixels,
                                                  connectivity=self.config.connectivity,
                                                  rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                  fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                  fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.dust_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, clip_step))

            # Save the mask
            if self.config.steps: self.dust_clip_masks[name].saveto(self.dust_step_path_for_mask(name, clip_step))

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
        self.soften_edges_old()

        # Young
        self.soften_edges_young()

        # Ionizing
        self.soften_edges_ionizing()

        # Dust
        self.soften_edges_dust()

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
            self.set_old_deprojection_aliases()

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
            self.set_young_deprojection_aliases()

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
            self.set_ionizing_deprojection_aliases()

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
            self.set_dust_deprojection_aliases()

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
        z_radius = 2. * self.old_scaleheight * self.config.scale_heights

        # Determine number of pixels
        nx = int(round(radial_extent / self.old_projection_physical_pixelscale))
        nz = int(round(z_radius / self.old_projection_physical_pixelscale))

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

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

# -----------------------------------------------------------------
