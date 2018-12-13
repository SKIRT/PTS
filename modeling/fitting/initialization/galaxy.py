#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.initialization.galaxy Contains the GalaxyFittingInitializer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ....core.basics.log import log
from ...component.galaxy import GalaxyModelingComponent
from .base import FittingInitializerBase
from ....core.prep.wavelengthgrids import WavelengthGridGenerator, get_min_wavelength, get_max_wavelength
from ....core.basics.emissionlines import get_identifiers, important_lines, strong_lines
from ....core.tools.utils import lazyproperty
from ....core.tools import filesystem as fs
from ....core.simulation.wavelengthgrid import WavelengthGrid

# -----------------------------------------------------------------

class GalaxyFittingInitializer(FittingInitializerBase, GalaxyModelingComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        :return:
        """

        # Call the constructors of the base classes
        FittingInitializerBase.__init__(self, no_config=True)
        GalaxyModelingComponent.__init__(self, *args, **kwargs)

        # The wavelength grids
        self.basic_wavelength_grids = None
        self.refined_wavelength_grids = None
        self.highres_wavelength_grids = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load the ski file
        self.load_ski()

        # Load the model representation
        self.load_representation()

        # Set the stellar and dust components
        self.set_components()

        # Get the wavelength grids
        self.get_wavelength_grids()

        # Adjust the ski template
        self.adjust_ski()

        # Calculate the weight factor to give to each band
        self.calculate_weights()

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
        FittingInitializerBase.setup(self, **kwargs)
        GalaxyModelingComponent.setup(self, **kwargs)

    # -----------------------------------------------------------------

    def load_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the template ski file ...")

        # Load the ski template
        self.ski = self.fitting_run.ski_template

    # -----------------------------------------------------------------

    def set_components(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the stellar and dust components ...")

        # Add the components to the ski file and to the input map paths dictionary
        self.suite.add_model_components(self.model_name, self.ski, self.input_map_paths)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):
        return self.suite.get_model_definition(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def metallicity(self):
        return self.model_definition.metallicity

    # -----------------------------------------------------------------

    @property
    def old_age(self):
        return self.model_definition.old_stars_age

    # -----------------------------------------------------------------

    @property
    def young_age(self):
        return self.model_definition.young_stars_age

    # -----------------------------------------------------------------

    @property
    def fitting_filters(self):
        return self.fitting_run.fitting_filters

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_names(self):
        return [str(fltr) for fltr in self.fitting_filters]

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_filter_wavelengths(self):
        return [fltr.wavelength for fltr in self.fitting_filters]

    # -----------------------------------------------------------------

    def get_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the wavelength grids ...")

        # Basic wavelength grids
        if self.has_basic: self.load_basic_wavelength_grids()
        else: self.create_basic_wavelength_grids()

        # Refined wavelength grids
        if self.has_refined: self.load_refined_wavelength_grids()
        else: self.create_refined_wavelength_grids()

        # High-resolution wavelength grids
        if self.has_highres: self.load_highres_wavelength_grids()
        else: self.create_highres_wavelength_grids()

    # -----------------------------------------------------------------

    @property
    def has_basic(self):

        """
        This function ...
        :return:
        """

        has_any = False
        has_all = True

        # Loop over the target npoints of the basic grids
        for npoints in self.basic_npoints_list:

            # Determine path
            dirname = "basic_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Check
            if not fs.is_file(grid_path): has_all = False
            else: has_any = True

        # Clear
        if not has_all: fs.remove_directories_in_path(self.wavelength_grids_path, startswith="basic_")

        # Return
        return has_any and has_all

    # -----------------------------------------------------------------

    @property
    def has_refined(self):

        """
        This function ...
        :return:
        """

        has_any = False
        has_all = True

        # Loop over the target npoints of the refined grids
        for npoints in self.refined_npoints_list:

            # Determine path
            dirname = "refined_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Check
            if not fs.is_file(grid_path): has_all = False
            else: has_any = True

        # Clear
        if not has_all: fs.remove_directories_in_path(self.wavelength_grids_path, startswith="refined_")

        # Return
        return has_any and has_all

    # -----------------------------------------------------------------

    @property
    def has_highres(self):

        """
        This function ...
        :return:
        """

        has_any = False
        has_all = True

        # Loop over the target npoints of the high-res grids
        for npoints in self.highres_npoints_list:

            # Determine path
            dirname = "highres_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Check
            if not fs.is_file(grid_path): has_all = False
            else: has_any = True

        # Clear
        if not has_all: fs.remove_directories_in_path(self.wavelength_grids_path, startswith="highres_")

        # Return
        return has_any and has_all

    # -----------------------------------------------------------------

    @lazyproperty
    def template_sed_old(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.old_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def template_sed_young(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.young_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def template_sed_ionizing(self):

        """
        This function ...
        :return:
        """

        return self.model_definition.ionizing_sed

    # -----------------------------------------------------------------

    @lazyproperty
    def template_seds(self):

        """
        This function ...
        :return:
        """

        seds = OrderedDict()
        seds["old"] = self.template_sed_old
        seds["young"] = self.template_sed_young
        seds["ionizing"] = self.template_sed_ionizing
        return seds

    # -----------------------------------------------------------------

    @lazyproperty
    def min_wavelength(self):
        return get_min_wavelength(self.config.range.min, self.observed_filters_no_iras_planck)

    # -----------------------------------------------------------------

    @lazyproperty
    def max_wavelength(self):
        return get_max_wavelength(self.config.range.max, self.observed_filters_no_iras_planck)

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_npoints_list(self):
        return self.config.wg.npoints_range_basic.linear(self.config.wg.ngrids_basic)

    # -----------------------------------------------------------------

    @lazyproperty
    def refined_npoints_list(self):
        return self.config.wg.npoints_range_refined.linear(self.config.wg.ngrids_refined)

    # -----------------------------------------------------------------

    @lazyproperty
    def highres_npoints_list(self):
        return self.config.wg.npoints_range_highres.linear(self.config.wg.ngrids_highres)

    # -----------------------------------------------------------------

    @lazyproperty
    def basic_grid_paths(self):

        """
        This function ....
        :return:
        """

        paths = OrderedDict()
        for npoints in self.basic_npoints_list:

            # Determine path
            dirname = "basic_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)
            if fs.is_directory(path): fs.clear_directory(path)
            else: fs.create_directory(path)

            # Set path
            paths[npoints] = path

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def refined_grid_paths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        for npoints in self.refined_npoints_list:

            # Determine path
            dirname = "refined_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)
            if fs.is_directory(path): fs.clear_directory(path)
            else: fs.create_directory(path)

            # Set path
            paths[npoints] = path

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    @lazyproperty
    def highres_grid_paths(self):

        """
        This function ...
        :return:
        """

        paths = OrderedDict()
        for npoints in self.highres_npoints_list:

            # Determine path
            dirname = "highres_" + str(npoints)
            path = fs.join(self.wavelength_grids_path, dirname)
            if fs.is_directory(path): fs.clear_directory(path)
            else: fs.create_directory(path)

            # Set path
            paths[npoints] = path

        # Return the paths
        return paths

    # -----------------------------------------------------------------

    @property
    def normalization_filters(self):
        return self.model_definition.normalization_filters

    # -----------------------------------------------------------------

    @property
    def normalization_wavelengths(self):
        return self.model_definition.normalization_wavelengths

    # -----------------------------------------------------------------

    def create_basic_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating basic wavelength grids ...")

        # Generate the grids
        self.basic_wavelength_grids = create_basic_wavelength_grids(self.config.wg.ngrids_basic, self.config.wg.npoints_range_basic,
                                                                    self.config.wg.range, add_emission_lines=self.config.wg.add_emission_lines,
                                                                    filters=self.fitting_filters, fixed=self.normalization_wavelengths,
                                                                    check_filters=self.observed_filters_no_iras, plot_seds=self.template_seds,
                                                                    table=self.wg_table, out_paths=self.basic_grid_paths, plot_paths=self.basic_grid_paths)

    # -----------------------------------------------------------------

    def load_basic_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Initialize the dictionary
        self.basic_wavelength_grids = OrderedDict()

        # Loop over the target npoints of the basic grids
        for npoints in self.basic_npoints_list:

            # Determine name
            dirname = "basic_" + str(npoints)

            # Debugging
            log.debug("Loading the '" + dirname + "' wavelength grid ...")

            # Determine path
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Load grid
            grid = WavelengthGrid.from_file(grid_path)

            # Add the grid
            self.basic_wavelength_grids[npoints] = grid

    # -----------------------------------------------------------------

    def create_refined_wavelength_grids(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Creating refined wavelength grids ...")

        # Generate the wavelength grids
        self.refined_wavelength_grids = create_refined_wavelength_grids(self.config.wg.ngrids_refined, self.config.wg.npoints_range_refined,
                                                                        self.config.wg.range, add_emission_lines=self.config.wg.add_emission_lines,
                                                                        filters=self.fitting_filters, fixed=self.normalization_wavelengths,
                                                                        check_filters=self.observed_filters_no_iras, plot_seds=self.template_seds,
                                                                        table=self.wg_table, out_paths=self.refined_grid_paths, plot_paths=self.refined_grid_paths)

    # -----------------------------------------------------------------

    def load_refined_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Initialize the dictionary
        self.refined_wavelength_grids = OrderedDict()

        # Loop over the target npoints of the refined grids
        for npoints in self.refined_npoints_list:

            # Determine name
            dirname = "refined_" + str(npoints)

            # Debugging
            log.debug("Loading the '" + dirname + "' wavelength grid ...")

            # Determine path
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Load grid
            grid = WavelengthGrid.from_file(grid_path)

            # Add the grid
            self.refined_wavelength_grids[npoints] = grid

    # -----------------------------------------------------------------

    def create_highres_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating high-resolution wavelength grids ...")

        # Generate the wavelength grids
        self.highres_wavelength_grids = create_highres_wavelength_grids(self.config.wg.ngrids_highres, self.config.wg.npoints_range_highres,
                                                                        self.config.wg.range, add_emission_lines=self.config.wg.add_emission_lines,
                                                                        filters=self.fitting_filters, fixed=self.normalization_wavelengths,
                                                                        check_filters=self.observed_filters_no_iras, plot_seds=self.template_seds,
                                                                        table=self.wg_table, out_paths=self.highres_grid_paths, plot_paths=self.highres_grid_paths)

    # -----------------------------------------------------------------

    def load_highres_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Initialize the dictionary
        self.highres_wavelength_grids = OrderedDict()

        # Loop over the target npoints of the high-res grids
        for npoints in self.highres_npoints_list:

            # Determine name
            dirname = "highres_" + str(npoints)

            # Debugging
            log.debug("Loading the '" + dirname + "' wavelength grid ...")

            # Determine path
            path = fs.join(self.wavelength_grids_path, dirname)

            # Determine grid file path
            grid_path = fs.join(path, "grid.dat")

            # Load grid
            grid = WavelengthGrid.from_file(grid_path)

            # Add the grid
            self.highres_wavelength_grids[npoints] = grid

    # -----------------------------------------------------------------

    def adjust_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adjusting the ski file parameters ...")

        # Remove the existing instruments
        self.ski.remove_all_instruments()

        # Add the instrument
        self.ski.add_instrument("earth", self.representation.sed_instrument)

        # Set the number of photon packages
        self.ski.setpackages(self.config.npackages)

        # Set the name of the wavelength grid file
        self.ski.set_file_wavelength_grid("wavelengths.txt")

        # Set the dust emissivityex
        self.set_dust_emissivity()

        # Set the dust grid
        self.set_dust_grid()

        # Set all-cells dust library
        self.ski.set_allcells_dust_lib()

        # Set the dust selfabsorption
        self.set_selfabsorption()

        # Disable all writing options
        self.ski.disable_all_writing_options()

    # -----------------------------------------------------------------

    def set_dust_emissivity(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust emissivity ...")

        # Set dust emissivity if present
        if self.ski.has_dust_emissivity:

            # Enable or disable
            if self.config.transient_heating: self.ski.set_transient_dust_emissivity()
            else: self.ski.set_grey_body_dust_emissivity()

    # -----------------------------------------------------------------

    def set_selfabsorption(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the dust self-absorption ...")

        # Dust self-absorption
        if self.config.selfabsorption: self.ski.enable_selfabsorption()
        else: self.ski.disable_selfabsorption()

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # 1. Write the ski file
        self.write_ski()

        # 2. Write the weights table
        self.write_weights()

        # 3. Write the wavelength grids
        self.write_wavelength_grids()

        # 4. Write the wavelength grid table
        self.write_wavelength_grid_table()

        # 5. Write the paths to the input maps
        self.write_input_map_paths()

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file to " + self.fitting_run.template_ski_path + " ...")

        # Save the ski template file
        self.ski.save()

    # -----------------------------------------------------------------

    def write_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the wavelength grids ...")

        # Basic
        self.write_basic_wavelength_grids()

        # Refined
        self.write_refined_wavelength_grids()

        # Highres
        self.write_highres_wavelength_grids()

    # -----------------------------------------------------------------

    def write_basic_wavelength_grids(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the basic wavelength grids ...")

        # Loop over the grids
        for target_npoints in self.basic_wavelength_grids:

            # Debugging
            log.debug("Writing the " + str(target_npoints) + " points basic wavelength grid ...")

            # Determine the filepath
            filename = "basic_" + str(target_npoints) + ".dat"
            filepath = fs.join(self.wavelength_grids_path, filename)

            # Save the wavelength grid
            self.basic_wavelength_grids[target_npoints].to_skirt_input(filepath)

    # -----------------------------------------------------------------

    def write_refined_wavelength_grids(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Writing the refined wavelength grids ...")

        # Loop over the grids
        for target_npoints in self.refined_wavelength_grids:

            # Debugging
            log.debug("Writing the " + str(target_npoints) + " points refined wavelength grid ...")

            # Determine the filepath
            filename = "refined_" + str(target_npoints) + ".dat"
            filepath = fs.join(self.wavelength_grids_path, filename)

            # Save the wavelength grid
            self.refined_wavelength_grids[target_npoints].to_skirt_input(filepath)

    # -----------------------------------------------------------------

    def write_highres_wavelength_grids(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Writing the high-resolution wavelength grids ...")

        # Loop over the grids
        for target_npoints in self.highres_wavelength_grids:

            # Debugging
            log.debug("Writing the " + str(target_npoints) + " points high-resolution wavelength grid ...")

            # Determine the filepath
            filename = "highres_" + str(target_npoints) + ".dat"
            filepath = fs.join(self.wavelength_grids_path, filename)

            # Save the wavelength grid
            self.highres_wavelength_grids[target_npoints].to_skirt_input(filepath)

# -----------------------------------------------------------------

def create_basic_wavelength_grids(ngrids, npoints_range, wavelength_range, filters=None, fixed=None, plot_seds=None,
                                  table=None, add_emission_lines=True, check_filters=None, out_paths=None, plot_paths=None):

    """
    This function ...
    :param ngrids:
    :param npoints_range:
    :param wavelength_range:
    :param filters:
    :param fixed:
    :param plot_seds:
    :param table:
    :param add_emission_lines:
    :param check_filters:
    :param out_paths:
    :param plot_paths:
    :return:
    """

    # Get list of strong emission lines
    strong_emission_lines = strong_lines

    # Get wavelengths from filters
    if filters is not None: wavelengths = [fltr.wavelength for fltr in filters]
    else: wavelengths = None

    # Create generator
    basic_generator = WavelengthGridGenerator()

    # Set basic settings
    basic_generator.config.ngrids = ngrids
    basic_generator.config.npoints_range = npoints_range
    basic_generator.config.range = wavelength_range

    # Set other
    basic_generator.config.add_emission_lines = add_emission_lines
    # basic_generator.config.emission_lines = self.important_emission_lines
    basic_generator.config.emission_lines = strong_emission_lines
    basic_generator.config.check_filters = check_filters
    basic_generator.config.adjust_to = wavelengths
    basic_generator.config.filters = None
    basic_generator.config.fixed = fixed
    basic_generator.config.plotting_filters = filters

    # Set other flags
    basic_generator.config.show = False
    basic_generator.config.write = True
    basic_generator.config.plot = True

    # Writing
    basic_generator.config.write_grids = True
    basic_generator.config.write_elements = True
    basic_generator.config.write_table = False

    # Advanced
    basic_generator.config.plot_resampled = True
    basic_generator.config.plot_residuals = True

    # Other
    basic_generator.config.label = "basic"

    # Generate the wavelength grids
    basic_generator.run(table=table, seds=plot_seds, out_paths=out_paths, plot_paths=plot_paths)

    # Initialize dictionary for the grids
    grids = OrderedDict()

    # Loop over the grids
    for index in range(basic_generator.ngrids):

        # Get target npoints
        target_npoints = basic_generator.get_target_npoints(index)

        # Get wavelength grid
        wavelength_grid = basic_generator.grids[index]

        # Set the grid
        grids[target_npoints] = wavelength_grid

    # Return the grids
    return grids

# -----------------------------------------------------------------

def create_refined_wavelength_grids(ngrids, npoints_range, wavelength_range, filters=None, fixed=None, plot_seds=None,
                                    table=None, add_emission_lines=True, check_filters=None, out_paths=None, plot_paths=None):

    """
    This function ...
    :param ngrids:
    :param npoints_range:
    :param wavelength_range:
    :param filters:
    :param fixed:
    :param plot_seds:
    :param table:
    :param add_emission_lines:
    :param check_filters:
    :param out_paths:
    :param plot_paths:
    :return:
    """

    # Get list of strong emission lines
    strong_emission_lines = strong_lines

    # Get wavelengths from filters
    if filters is not None: wavelengths = [fltr.wavelength for fltr in filters]
    else: wavelengths = None

    # Create generator
    refined_generator = WavelengthGridGenerator()

    # Set basic settings
    refined_generator.config.ngrids = ngrids
    refined_generator.config.npoints_range = npoints_range
    refined_generator.config.range = wavelength_range

    # Set other
    refined_generator.config.add_emission_lines = add_emission_lines
    # refined_generator.config.emission_lines = self.important_emission_lines
    refined_generator.config.emission_lines = strong_emission_lines
    refined_generator.config.check_filters = check_filters
    refined_generator.config.adjust_to = wavelengths
    refined_generator.config.filters = filters
    refined_generator.config.fixed = fixed
    refined_generator.config.plotting_filters = filters

    # Set other flags
    refined_generator.config.show = False
    refined_generator.config.write = True
    refined_generator.config.plot = True

    # Writing
    refined_generator.config.write_grids = True
    refined_generator.config.write_elements = True
    refined_generator.config.write_table = False

    # Advanced
    refined_generator.config.plot_resampled = True
    refined_generator.config.plot_residuals = True

    # Other
    refined_generator.config.label = "refined"

    # Generate the refined wavelength grids
    refined_generator.run(table=table, seds=plot_seds, out_paths=out_paths, plot_paths=plot_paths)

    # Initialize dictionary for the grids
    grids = OrderedDict()

    # Loop over the grids
    for index in range(refined_generator.ngrids):

        # Get target npoints
        target_npoints = refined_generator.get_target_npoints(index)

        # Get wavelength grid
        wavelength_grid = refined_generator.grids[index]

        # Set the grid
        grids[target_npoints] = wavelength_grid

    # Return the grids
    return grids

# -----------------------------------------------------------------

def create_highres_wavelength_grids(ngrids, npoints_range, wavelength_range, filters=None, fixed=None, plot_seds=None,
                                    table=None, add_emission_lines=True, check_filters=None, out_paths=None, plot_paths=None):

    """
    This function ...
    :param ngrids:
    :param npoints_range:
    :param wavelength_range:
    :param filters:
    :param fixed:
    :param plot_seds:
    :param table:
    :param add_emission_lines:
    :param check_filters:
    :param out_paths:
    :param plot_paths:
    :return:
    """

    # Get list of all emission lines
    all_emission_lines = get_identifiers()

    # Get wavelengths from filters
    if filters is not None: wavelengths = [fltr.wavelength for fltr in filters]
    else: wavelengths = None

    # Create generator
    highres_generator = WavelengthGridGenerator()

    # Set basic settings
    highres_generator.config.ngrids = ngrids
    highres_generator.config.npoints_range = npoints_range
    highres_generator.config.range = wavelength_range

    # Set other
    highres_generator.config.add_emission_lines = add_emission_lines
    highres_generator.config.emission_lines = all_emission_lines
    highres_generator.config.check_filters = check_filters
    highres_generator.config.adjust_to = wavelengths
    highres_generator.config.filters = filters
    highres_generator.config.fixed = fixed
    highres_generator.config.plotting_filters = filters

    # Set other flags
    highres_generator.config.show = False
    highres_generator.config.write = True
    highres_generator.config.plot = True

    # Writing
    highres_generator.config.write_grids = True
    highres_generator.config.write_elements = True
    highres_generator.config.write_table = False

    # Advanced
    highres_generator.config.plot_resampled = True
    highres_generator.config.plot_residuals = True

    # Other
    highres_generator.config.label = "highres"

    # Generate the high-resolution wavelength grids
    highres_generator.run(table=table, seds=plot_seds, out_paths=out_paths, plot_paths=plot_paths)

    # Initialize dictionary for the grids
    grids = OrderedDict()

    # Loop over the grids
    for index in range(highres_generator.ngrids):

        # Get target npoints
        target_npoints = highres_generator.get_target_npoints(index)

        # Get wavelength grid
        wavelength_grid = highres_generator.grids[index]

        # Set the grid
        grids[target_npoints] = wavelength_grid

    # Return the grids
    return grids

# -----------------------------------------------------------------
