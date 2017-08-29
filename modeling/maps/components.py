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
from ..misc.deprojector import Deprojector
from ...core.remote.remote import Remote

# -----------------------------------------------------------------

steps_name = "steps"
masks_name = "masks"
deprojected_name = "deprojected"

# -----------------------------------------------------------------

correct_step = "corrected"
truncate_step = "truncated"
crop_step = "cropped"
clip_step = "clipped"
softened_step = "softened"
steps = [correct_step, truncate_step, crop_step, clip_step, softened_step]

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
    return reversed(after)

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

        # Paths
        self.old_masks_path = None
        self.young_masks_path = None
        self.ionizing_masks_path = None
        self.dust_masks_path = None
        self.old_deprojection_path = None
        self.young_deprojection_path = None
        self.ionizing_deprojection_path = None
        self.dust_deprojection_path = None

        # The clip masks
        self.old_masks = dict()
        self.young_masks = dict()
        self.ionizing_masks = dict()
        self.dust_masks = dict()

        # The deprojected maps
        self.old_deprojected = dict()
        self.young_deprojected = dict()
        self.ionizing_deprojected = dict()
        self.dust_deprojected = dict()

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Prompt
        self.prompt()

        # 3. Load the maps
        self.load_maps()

        # 4. Process the maps
        self.process_maps()

        # 5. Deproject the maps
        self.deproject()

        # Writing
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

        # Return the map path
        return fs.join(map_path, step + ".fits")

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

        # Return the map path
        return fs.join(map_path, step + ".fits")

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

        # Return the map path
        return fs.join(map_path, step + ".fits")

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

        # Return the map path
        return fs.join(map_path, step + ".fits")

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

            # Debugging
            log.debug("Correcting the '" + name + "' dust map ...")

            # Correct the map
            self.correct_map(self.dust_maps[name])

            # Set flag
            self.dust_maps[name].metadata[correct_step] = True

            # Save intermediate result
            if self.config.step: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, correct_step))

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
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' old stellar map ...")

            # Get the origins
            origins = self.old_map_origins[name]

            # Clip the map (returns the mask)
            self.old_masks[name] = self.clip_map(self.old_maps[name], origins, convolve=self.config.convolve,
                                                 remote=self.remote, npixels=self.config.min_npixels,
                                                 connectivity=self.config.connectivity,
                                                 rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                 fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                 fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.old_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, clip_step))

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
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' young stellar map ...")

            # Get the origins
            origins = self.young_map_origins[name]

            # Clip the map (returns the mask)
            self.young_masks[name] = self.clip_map(self.young_maps[name], origins, convolve=self.config.convolve,
                                                   remote=self.remote, npixels=self.config.min_npixels,
                                                   connectivity=self.config.connectivity,
                                                   rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                   fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                   fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.young_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, clip_step))

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
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' ionizing stellar map ...")

            # Get the origins
            origins = self.ionizing_map_origins[name]

            # Clip the map (returns the mask)
            self.ionizing_masks[name] = self.clip_map(self.ionizing_maps[name], origins, convolve=self.config.convolve,
                                                      remote=self.remote, npixels=self.config.min_npixels,
                                                      connectivity=self.config.connectivity,
                                                      rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                      fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                      fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.ionizing_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, clip_step))

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
                continue

            # Debugging
            log.debug("Clipping the '" + name + "' dust map ...")

            # Get the origins
            origins = self.dust_map_origins[name]

            # Clip the map (returns the mask)
            self.dust_masks[name] = self.clip_map(self.dust_maps[name], origins, convolve=self.config.convolve,
                                                  remote=self.remote, npixels=self.config.min_npixels,
                                                  connectivity=self.config.connectivity,
                                                  rebin_remote_threshold=self.config.rebin_remote_threshold,
                                                  fuzzy=self.config.fuzzy_mask, fuzziness=self.config.fuzziness,
                                                  fuzziness_offset=self.config.fuzzy_min_significance_offset)

            # Set flag
            self.dust_maps[name].metadata[clip_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, clip_step))

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
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' old stellar map ...")

            # Soften
            self.soften_map(self.old_maps[name], self.softening_ellipse, self.softening_range)

            # Set flag
            self.old_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.old_maps[name].saveto(self.old_step_path_for_map(name, softened_step))

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
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' young stellar map ...")

            # Soften
            self.soften_map(self.young_maps[name], self.softening_ellipse, self.softening_range)

            # Set flag
            self.young_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.young_maps[name].saveto(self.young_step_path_for_map(name, softened_step))

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
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' ionizing stellar map ...")

            # Soften
            self.soften_map(self.ionizing_maps[name], self.softening_ellipse, self.softening_range)

            # Set flag
            self.ionizing_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.ionizing_maps[name].saveto(self.ionizing_step_path_for_map(name, softened_step))

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
                continue

            # Debugging
            log.debug("Softening the edges of the '" + name + "' dust map ...")

            # Soften the map
            self.soften_map(self.dust_maps[name], self.softening_ellipse, self.softening_range)

            # Set flag
            self.dust_maps[name].metadata[softened_step] = True

            # Save intermediate result
            if self.config.steps: self.dust_maps[name].saveto(self.dust_step_path_for_map(name, softened_step))

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

    def deproject_old(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the old stellar maps ...")

        # Create the deprojector
        deprojector = Deprojector()

        # Set settings
        deprojector.config.method = "pts"
        deprojector.config.writing.deprojections = False

        # Run the deprojector
        deprojector.run(maps=self.old_maps, scale_height=self.old_scaleheight, root_path=self.old_deprojection_path)

        # Get the deprojected maps
        self.old_deprojected = deprojector.deprojected

    # -----------------------------------------------------------------

    def deproject_young(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the young stellar maps ...")

        # Create the deprojector
        deprojector = Deprojector()

        # Set settings
        deprojector.config.method = "pts"
        deprojector.config.writing.deprojections = False

        # Run the deprojector
        deprojector.run(maps=self.young_maps, scale_height=self.young_scaleheight, root_path=self.young_deprojection_path)

        # Get the deprojected maps
        self.young_deprojected = deprojector.deprojected

    # -----------------------------------------------------------------

    def deproject_ionizing(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the ionizing stellar maps ...")

        # Create the deprojector
        deprojector = Deprojector()

        # Set settings
        deprojector.config.method = "pts"
        deprojector.config.writing.deprojections = False

        # Run the deprojector
        deprojector.run(maps=self.ionizing_maps, scale_height=self.ionizing_scaleheight, root_path=self.ionizing_deprojection_path)

        # Get the deprojected maps
        self.ionizing_deprojected = deprojector.deprojected

    # -----------------------------------------------------------------

    def deproject_dust(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the dust maps ...")

        # Create the deprojector
        deprojector = Deprojector()

        # Set settings
        deprojector.config.method = "pts"
        deprojector.config.writing.deprojections = False

        # Run the deprojector
        deprojector.run(maps=self.dust_maps, scale_height=self.dust_scaleheight, root_path=self.dust_deprojection_path)

        # Get the deprojected maps
        self.dust_deprojected = deprojector.deprojected

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
        self.write_masks()

        # Write the deprojected maps
        self.write_deprojected()

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

    def write_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the clip masks ...")

        # Old
        self.write_old_masks()

        # Young
        self.write_young_masks()

        # Ionizing
        self.write_ionizing_masks()

        # Dust
        self.write_dust_masks()

    # -----------------------------------------------------------------

    def write_old_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the old stellar maps ...")

        # Loop over the masks
        for name in self.old_masks:

            # Debugging
            log.debug("Writing the '" + name + "' old stellar mask ...")

            # Determine the path
            path = fs.join(self.old_masks_path, name + ".fits")

            # Write
            self.old_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_young_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the young stellar maps ...")

        # Loop over the masks
        for name in self.young_masks:

            # Debugging
            log.debug("Writing the '" + name + "' young stellar mask ...")

            # Determine the path
            path = fs.join(self.young_masks_path, name + ".fits")

            # Write
            self.young_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_ionizing_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the ionizing stellar maps ...")

        # Loop over the masks
        for name in self.ionizing_masks:

            # Debugging
            log.debug("Writing the '" + name + "' ionizing stellar mask ...")

            # Determine the path
            path = fs.join(self.ionizing_masks_path, name + ".fits")

            # Write
            self.ionizing_masks[name].saveto(path)

    # -----------------------------------------------------------------

    def write_dust_masks(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the masks of the dust maps ...")

        # Loop over the masks
        for name in self.dust_masks:

            # Debugging
            log.debug("Writing the '" + name + "' dust mask ...")

            # Determine the path
            path = fs.join(self.dust_masks_path, name + ".fits")

            # Write
            self.dust_masks[name].saveto(path)

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

    @property
    def maps_sub_path(self):

        """
        This function ...
        :return: 
        """

        return None

# -----------------------------------------------------------------
