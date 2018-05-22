#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.maps.selector Contains the ComponentMapsSelector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from ...core.basics.log import log
from .selectioncomponent import MapsSelectionComponent
from ...core.tools.utils import lazyproperty
from ...core.tools import filesystem as fs
from ...core.tools import numbers
from ...core.tools import sequences
from ...magic.core.mask import Mask
from ...core.filter.filter import parse_filter
from ...core.tools.parsing import real
from ...magic.core.image import Image
from ...magic.tools import colours
from ...magic.basics.mask import Mask as oldMask
from ...core.basics.containers import ordered_by_key
from ...core.tools import formatting as fmt
from ...core.tools.serialization import write_dict

# -----------------------------------------------------------------

class ComponentMapsSelector(MapsSelectionComponent):

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
        super(ComponentMapsSelector, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Auto-select
        if self.config.auto: self.auto_select()

        # 3. Prompt
        self.prompt()

        # 4. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ComponentMapsSelector, self).setup(**kwargs)

        # Set random
        if self.config.random: self.config.random_old = self.config.random_young = self.config.random_ionizing = self.config.random_dust = self.config.random

        # Set all
        if self.config.all: self.config.all_old = self.config.all_young = self.config.all_ionizing = self.config.all_dust = True

        # Make selections
        self.old_selection = sequences.make_selection(self.old_map_names, self.config.old, self.config.not_old, nrandom=self.config.random_old, all=self.config.all_old, indices=self.config.old_indices, not_indices=self.config.not_old_indices)
        self.young_selection = sequences.make_selection(self.young_map_names, self.config.young, self.config.not_young, nrandom=self.config.random_young, all=self.config.all_young, indices=self.config.young_indices, not_indices=self.config.not_young_indices)
        self.ionizing_selection = sequences.make_selection(self.ionizing_map_names, self.config.ionizing, self.config.not_ionizing, nrandom=self.config.random_ionizing, all=self.config.all_ionizing, indices=self.config.ionizing_indices, not_indices=self.config.not_ionizing_indices)
        self.dust_selection = sequences.make_selection(self.dust_map_names, self.config.dust, self.config.not_dust, nrandom=self.config.random_dust, all=self.config.all_dust, indices=self.config.dust_indices, not_indices=self.config.not_dust_indices)

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
        preferred_colours = ["FUV-r", "FUV-i", "FUV-H", "FUV-g"]

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

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the selection
        self.write_selection()

    # -----------------------------------------------------------------

    def write_selection(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the selection ...")

        # Make single selection dictionary
        selection = dict()

        # Set the selected map names
        selection["old"] = self.old_selection
        selection["young"] = self.young_selection
        selection["ionizing"] = self.ionizing_selection
        selection["dust"] = self.dust_selection

        # Determine path for the selection file
        current_indices = fs.files_in_path(self.maps_components_path, extension="dat", returns="name", startswith="selection", convert=int, convert_split_pattern="_", convert_split_index=1)
        index = numbers.lowest_missing_integer(current_indices)
        selection_path = fs.join(self.maps_components_path, "selection_" + str(index) + ".dat")

        # Write the selection
        write_dict(selection, selection_path)

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
