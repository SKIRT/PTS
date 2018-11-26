#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.analysis.run Contains the AnalysisRun class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta, abstractproperty
from collections import OrderedDict

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import SkiFile
from ..core.model import RTModel
from ...core.tools import sequences
from ...core.basics.composite import SimplePropertyComposite
from ..fitting.run import FittingRun
from ...core.tools.utils import lazyproperty
from ...core.tools.serialization import load_dict
from ...core.simulation.tree import DustGridTree
from ...core.simulation.grids import FileTreeDustGrid, load_grid
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..basics.projection import GalaxyProjection, EdgeOnProjection, FaceOnProjection
from ..basics.instruments import FullInstrument, SimpleInstrument, SEDInstrument
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.basics.log import log
from ...core.remote.remote import load_remote
from ...core.basics.configuration import Configuration
from ...core.simulation.logfile import LogFile
from ...core.tools import strings
from ...core.extract.progress import ProgressTable
from ...core.extract.timeline import TimeLineTable
from ...core.extract.memory import MemoryUsageTable
from ...magic.core.dataset import DataSet
from ..core.environment import colours_name as colour_maps_name
from ..core.environment import ssfr_name as ssfr_maps_name
from ..core.environment import tir_name as tir_maps_name
from ..core.environment import attenuation_name as attenuation_maps_name
from ..core.environment import old_name as old_maps_name
from ..core.environment import young_name as young_maps_name
from ..core.environment import ionizing_name as ionizing_maps_name
from ..core.environment import dust_name as dust_maps_name
from ..core.environment import rgb_name as rgb_maps_name
from ...core.data.sed import ObservedSED, SED
from ...magic.core.datacube import DataCube
from ...core.simulation.tree import get_nleaves
from ..build.definition import ModelDefinition
from ..basics.properties import GalaxyProperties
from ...core.tools import tables
from ...core.filter.filter import parse_filter
from ...magic.core.frame import Frame
from ...core.basics.distribution import Distribution
from ...magic.region.list import SkyRegionList

# -----------------------------------------------------------------

wavelengths_filename = "wavelengths.txt"
dustgridtree_filename = "tree.dat"

# -----------------------------------------------------------------

# Set contribution nmes
total = "total"
bulge = "bulge"
disk = "disk"
old = "old"
young = "young"
ionizing = "ionizing"
unevolved = "unevolved"
extra = "extra"

# All contributions
contributions = [total, bulge, disk, old, young, ionizing, unevolved]

# -----------------------------------------------------------------

def get_analysis_run_cwd(name):

    """
    This function ...
    :param name:
    :return:
    """

    return get_analysis_run(fs.cwd(), name)

# -----------------------------------------------------------------

def get_analysis_run(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :param name:
    :return:
    """

    # Get runs object
    runs = AnalysisRuns(modeling_path)

    # Return analysis run
    return runs.load(name)

# -----------------------------------------------------------------

def get_analysis_model_cwd(name):

    """
    This function ...
    :param name:
    :return:
    """

    return get_analysis_model(fs.cwd(), name)

# -----------------------------------------------------------------

def get_analysis_model(modeling_path, name):

    """
    This function ...
    :param modeling_path:
    :return:
    """

    # Get analysis run
    run = get_analysis_run(modeling_path, name)

    # Return model
    return run.model

# -----------------------------------------------------------------

class AnalysisRunInfo(SimplePropertyComposite):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(AnalysisRunInfo, self).__init__()

        # Define properties
        self.add_string_property("name", "name of the analysis run")
        self.add_string_property("path", "path of the analysis run")
        self.add_string_property("fitting_run", "fitting run name")
        self.add_string_property("generation_name", "generation name")
        self.add_string_property("simulation_name", "simulation name")
        self.add_string_property("model_name", "name of the model")
        self.add_real_property("chi_squared", "chi squared value of the fitted model")
        self.add_string_property("reference_deprojection", "name of the deprojection model that is used for the creating the instruments")

        # Parameter values dictionary
        self.add_section("parameter_values", "parameter values", dynamic=True)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

# Various filenames
dust_grid_filename = "dust_grid.dg"
wavelength_grid_filename = "wavelength_grid.dat"
dust_grid_build_name = "dust grid"
info_filename = "info.dat"
config_filename = "config.cfg"
launch_config_filename = "launch_config.cfg"
input_filename = "input.dat"
dust_grid_tree_filename = "tree.dat"

# Directories
model_name = "model"
instruments_name = "instruments"
projections_name = "projections"
extract_name = "extr"
plot_name = "plot"
misc_name = "misc"
evaluation_name = "evaluation"
contributions_name = "contributions"

# Analysis directories
properties_name = "properties"
attenuation_name = "attenuation"
colours_name = "colours"
fluxes_name = "fluxes"
images_name = "images"
residuals_name = "residuals"
maps_name = "maps"
absorption_name = "absorption"
heating_name = "heating"
energy_name = "energy"
sfr_name = "sfr"
correlations_name = "correlations"

# Projection filenames
earth_projection_filename = "earth.proj"
faceon_projection_filename = "faceon.proj"
edgeon_projection_filename = "edgeon.proj"

# Instrument filenames
sed_earth_instrument_filename = "earth_sed.instr"
full_sed_earth_instrument_filename = "earth_full_sed.instr"
simple_earth_instrument_filename = "earth_simple.instr"
full_earth_instrument_filename = "earth_full.instr"
simple_faceon_instrument_filename = "faceon_simple.instr"
full_faceon_instrument_filename = "faceon_full.instr"
simple_edgeon_instrument_filename = "edgeon_simple.instr"
full_edgeon_instrument_filename = "edgeon_full.instr"

# -----------------------------------------------------------------

class AnalysisRunBase(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @property
    def from_fitting(self):
        return self.fitting_run_name is not None

    # -----------------------------------------------------------------

    @property
    def from_model(self):
        return self.fitting_run_name is None

    # -----------------------------------------------------------------

    @property
    def from_generation(self):
        # Otherwise: from initial guess
        return self.from_fitting and self.generation_name is not None

    # -----------------------------------------------------------------

    @property
    def from_initial_guess(self):
        # Otherwise: from best simulation of a certain generation
        return self.from_fitting and self.generation_name is None

    # -----------------------------------------------------------------

    @property
    def name(self):
        return self.info.name

    # -----------------------------------------------------------------

    @property
    def generation_name(self):
        return self.info.generation_name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):
        return self.info.simulation_name

    # -----------------------------------------------------------------

    @property
    def model_name(self):
        return self.info.model_name

    # -----------------------------------------------------------------

    @property
    def input_file_path(self):
        return fs.join(self.path, input_filename)

    # -----------------------------------------------------------------

    @property
    def ski_file_path(self):
        return fs.join(self.path, self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_path(self):
        # Set the path to the wavelength grid file
        return fs.join(self.path, wavelength_grid_filename)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):
        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @property
    def dust_grid_path(self):
        # Set the path to the dust grid file
        return fs.join(self.path, dust_grid_filename)

    # -----------------------------------------------------------------

    @property
    def info_path(self):
        # Set the path to the analysis run info file
        return fs.join(self.path, info_filename)

    # -----------------------------------------------------------------

    @property
    def config_path(self):
        return fs.join(self.path, config_filename)

    # -----------------------------------------------------------------

    @property
    def heating_config_path(self):
        return fs.join(self.heating_path, config_filename)

    # -----------------------------------------------------------------

    @property
    def launch_config_path(self):
        return fs.join(self.path, launch_config_filename)

    # -----------------------------------------------------------------

    @property
    def total_simulation_path(self):
        return self.simulation_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_output_path(self):
        return self.output_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def bulge_output_path(self):
        return self.output_path_for_contribution(bulge)

    # -----------------------------------------------------------------

    @property
    def disk_output_path(self):
        return self.output_path_for_contribution(disk)

    # -----------------------------------------------------------------

    @property
    def old_output_path(self):
        return self.output_path_for_contribution(old)

    # -----------------------------------------------------------------

    @property
    def young_output_path(self):
        return self.output_path_for_contribution(young)

    # -----------------------------------------------------------------

    @property
    def ionizing_output_path(self):
        return self.output_path_for_contribution(ionizing)

    # -----------------------------------------------------------------

    @property
    def unevolved_output_path(self):
        return self.output_path_for_contribution(unevolved)

    # -----------------------------------------------------------------

    @property
    def extra_output_path(self):
        path = self.output_path_for_contribution(extra, create=False)
        if not fs.is_directory(path): return None
        else: return path

    # -----------------------------------------------------------------

    @property
    def total_logfile_path(self):
        return self.logfile_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_extr_path(self):
        return self.extraction_path_for_contribution(total)

    # -----------------------------------------------------------------

    @property
    def total_extract_path(self):
        return self.total_extr_path

    # -----------------------------------------------------------------

    @property
    def progress_path(self):
        return fs.join(self.total_extr_path, "progress.dat")

    # -----------------------------------------------------------------

    @property
    def timeline_path(self):
        return fs.join(self.total_extr_path, "timeline.dat")

    #  -----------------------------------------------------------------

    @property
    def memory_path(self):
        return fs.join(self.total_extr_path, "memory.dat")

    # -----------------------------------------------------------------

    @property
    def total_plot_path(self):
        return fs.create_directory_in(self.total_simulation_path, plot_name)

    # -----------------------------------------------------------------

    @property
    def total_misc_path(self):
        return fs.create_directory_in(self.total_simulation_path, misc_name)

    # -----------------------------------------------------------------

    @property
    def evaluation_path(self):
        return fs.join(self.path, evaluation_name)

    # -----------------------------------------------------------------

    @property
    def contributions_path(self):
        return fs.join(self.path, contributions_name)

    # -----------------------------------------------------------------

    @property
    def properties_path(self):
        return fs.join(self.path, properties_name)

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):
        return fs.join(self.path, attenuation_name)

    # -----------------------------------------------------------------

    @property
    def colours_path(self):
        return fs.join(self.path, colours_name)

    # -----------------------------------------------------------------

    @property
    def colours_observed_path(self):
        return fs.join(self.colours_path, "observed")

    # -----------------------------------------------------------------

    @property
    def colours_simulated_path(self):
        return fs.join(self.colours_path, "simulated")

    # -----------------------------------------------------------------

    @property
    def colours_residuals_path(self):
        return fs.join(self.colours_path, "residuals")

    # -----------------------------------------------------------------

    @abstractproperty
    def colour_names(self):
        pass

    # -----------------------------------------------------------------

    @property
    def fluxes_path(self):
        return fs.join(self.path, fluxes_name)

    # -----------------------------------------------------------------

    @property
    def images_path(self):
        return fs.join(self.path, images_name)

    # -----------------------------------------------------------------

    @property
    def residuals_path(self):
        return fs.join(self.path, residuals_name)

    # -----------------------------------------------------------------

    @property
    def maps_path(self):
        return fs.join(self.path, maps_name)

    # -----------------------------------------------------------------

    @property
    def colour_maps_path(self):
        return fs.join(self.maps_path, colour_maps_name)

    # -----------------------------------------------------------------

    @property
    def colour_maps_name(self):
        return fs.name(self.colour_maps_path)

    # -----------------------------------------------------------------

    @property
    def ssfr_maps_path(self):
        return fs.join(self.maps_path, ssfr_maps_name)

    # -----------------------------------------------------------------

    @property
    def ssfr_maps_name(self):
        return fs.name(self.ssfr_maps_path)

    # -----------------------------------------------------------------

    @property
    def tir_maps_path(self):
        return fs.join(self.maps_path, tir_maps_name)

    # -----------------------------------------------------------------

    @property
    def tir_maps_name(self):
        return fs.name(self.tir_maps_path)

    # -----------------------------------------------------------------

    @property
    def attenuation_maps_path(self):
        return fs.join(self.maps_path, attenuation_maps_name)

    # -----------------------------------------------------------------

    @property
    def attenuation_maps_name(self):
        return fs.name(self.attenuation_maps_path)

    # -----------------------------------------------------------------

    @property
    def old_maps_path(self):
        return fs.join(self.maps_path, old_maps_name)

    # -----------------------------------------------------------------

    @property
    def old_maps_name(self):
        return fs.name(self.old_maps_path)

    # -----------------------------------------------------------------

    @property
    def dust_maps_path(self):
        return fs.join(self.maps_path, dust_maps_name)

    # -----------------------------------------------------------------

    @property
    def dust_maps_name(self):
        return fs.name(self.dust_maps_path)

    # -----------------------------------------------------------------

    @property
    def young_maps_path(self):
        return fs.join(self.maps_path, young_maps_name)

    # -----------------------------------------------------------------

    @property
    def young_maps_name(self):
        return fs.name(self.young_maps_path)

    # -----------------------------------------------------------------

    @property
    def ionizing_maps_path(self):
        return fs.join(self.maps_path, ionizing_maps_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_maps_name(self):
        return fs.name(self.ionizing_maps_path)

    # -----------------------------------------------------------------

    @property
    def rgb_maps_path(self):
        return fs.join(self.maps_path, rgb_maps_name)

    # -----------------------------------------------------------------

    @property
    def rgb_maps_name(self):
        return fs.name(self.rgb_maps_path)

    # -----------------------------------------------------------------

    @property
    def absorption_path(self):
        return fs.join(self.path, absorption_name)

    # -----------------------------------------------------------------

    @property
    def heating_path(self):
        return fs.join(self.path, heating_name)

    # -----------------------------------------------------------------

    @property
    def energy_path(self):
        return fs.join(self.path, energy_name)

    # -----------------------------------------------------------------

    @property
    def sfr_path(self):
        return fs.join(self.path, sfr_name)

    # -----------------------------------------------------------------

    @property
    def correlations_path(self):
        return fs.join(self.path, correlations_name)

    # -----------------------------------------------------------------

    @property
    def dust_grid_build_path(self):
        return fs.join(self.path, dust_grid_build_name)

    # -----------------------------------------------------------------

    @property
    def dust_grid_simulation_out_path(self):
        return fs.join(self.dust_grid_build_path, "out")

    # -----------------------------------------------------------------

    def simulation_path_for_contribution(self, contribution, create=True):

        """
        This function ...
        :param contribution:
        :param create:
        :return:
        """

        if create: return fs.create_directory_in(self.contributions_path, contribution)
        else: return fs.join(self.contributions_path, contribution)

    # -----------------------------------------------------------------

    def ski_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.simulation_path_for_contribution(contribution), self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    def output_path_for_contribution(self, contribution, create=True):

        """
        This function ...
        :param contribution:
        :param create:
        :return:
        """

        if create: return fs.create_directory_in(self.simulation_path_for_contribution(contribution, create=create), "out")
        else: return fs.join(self.simulation_path_for_contribution(contribution, create=create), "out")

    # -----------------------------------------------------------------

    def logfile_path_for_contribution(self, contribution):

        """
        Thisf unction ...
        :return:
        """

        return fs.join(self.output_path_for_contribution(contribution), self.galaxy_name + "_log.txt")

    # -----------------------------------------------------------------

    def extraction_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.contributions_path, contribution, extract_name)

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):
        return self.info.name

    # -----------------------------------------------------------------

    @property
    def dust_grid_simulation_logfile_path(self):

        """
        This function ...
        :return:
        """

        # Determine the output path of the dust grid simulation
        out_path = self.dust_grid_simulation_out_path

        # Determine the log file path
        logfile_path = fs.join(out_path, "dustgrid_log.txt")

        # Return the log file path
        return logfile_path

    # -----------------------------------------------------------------

    @abstractproperty
    def has_dust_grid_simulation_logfile(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def dust_grid_simulation_logfile(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def dust_grid_tree(self):
        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def ncells(self):

        """
        This function ...
        :return:
        """

        # Read the log file
        if self.has_dust_grid_simulation_logfile:

            # Debugging
            log.debug("Determining the number of dust cells by reading the dust grid simulation's log file ...")

            # Load the log file and get the number of dust cells
            return self.dust_grid_simulation_logfile.dust_cells

        # Log file cannot be found
        else:

            # Debugging
            log.debug("Determining the number of dust cells by reading the dust cell tree data file (this can take a while) ...")

            # Get the number of leave nodes
            #return self.dust_grid_tree.nleaves  # requires loading the entire tree file!
            return get_nleaves(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    def get_remote_script_input_paths_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        paths = []

        # Loop over the commands
        for command in self.get_remote_script_commands_for_host(host_id):

            input_path = command.split("-i ")[1]
            if strings.is_quote_character(input_path[0]): input_path = input_path[1:].split(input_path[0])[0]
            else: input_path = input_path.split(" ")[0]

            paths.append(input_path)

        # Return the list of paths
        return paths

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_young(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def nyoung_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_tir(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def ntir_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_ssfr(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def nssfr_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_old(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def nold_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_ionizing(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def nionizing_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_dust(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def ndust_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_colours(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def ncolour_maps(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def has_maps_attenuation(self):
        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def nattenuation_maps(self):
        pass

    # -----------------------------------------------------------------

    @property
    def has_maps(self):
        return self.has_maps_attenuation or self.has_maps_colours or self.has_maps_dust or self.has_maps_ionizing or self.has_maps_old or self.has_maps_ssfr or self.has_maps_tir or self.has_maps_young

    # -----------------------------------------------------------------

    @abstractproperty
    def has_heating(self):
        pass

    # -----------------------------------------------------------------

    @property
    def simulated_sed_path(self):
        return fs.join(self.total_output_path, self.simulation_prefix + "_earth_sed.dat")

    # -----------------------------------------------------------------

    @abstractproperty
    def has_simulated_sed(self):
        pass

    # -----------------------------------------------------------------

    @property
    def simulated_fluxes_path(self):
        return fs.join(self.total_misc_path, self.simulation_prefix + "_earth_fluxes.dat")

    # -----------------------------------------------------------------

    @abstractproperty
    def has_simulated_fluxes(self):
        pass

    # -----------------------------------------------------------------

    @property
    def simulation_prefix(self):
        return self.galaxy_name

    # -----------------------------------------------------------------

    @property
    def simulated_datacube_path(self):
        return fs.join(self.total_output_path, self.galaxy_name + "_earth_total.fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_fluxes_path(self):
        return fs.create_directory_in(self.evaluation_path, "fluxes")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_images_path(self):
        return fs.create_directory_in(self.evaluation_path, "images")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_image_names(self):
        return fs.files_in_path(self.evaluation_images_path, returns="name", extension="fits", not_contains="error")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_image_filters(self):
        return [parse_filter(name) for name in self.evaluation_image_names]

    # -----------------------------------------------------------------

    def has_evaluation_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_image_filters

    # -----------------------------------------------------------------

    def get_evaluation_image_path_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_images_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def get_evaluation_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_evaluation_image_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_images_path(self):
        return fs.create_directory_in(self.evaluation_path, "proper_images")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_image_names(self):
        return fs.files_in_path(self.evaluation_proper_images_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_image_filters(self):
        return [parse_filter(name) for name in self.evaluation_proper_image_names]

    # -----------------------------------------------------------------

    def has_evaluation_proper_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_proper_image_filters

    # -----------------------------------------------------------------

    def get_evaluation_proper_image_path_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_proper_images_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def get_evaluation_proper_image_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_evaluation_proper_image_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_residuals_path(self):
        return fs.create_directory_in(self.evaluation_path, "residuals")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_residuals_names(self):
        return fs.files_in_path(self.evaluation_residuals_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_residuals_filters(self):
        return [parse_filter(name) for name in self.evaluation_residuals_names]

    # -----------------------------------------------------------------

    def has_evaluation_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_residuals_filters

    # -----------------------------------------------------------------

    def get_evaluation_residuals_path_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_residuals_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def get_evaluation_residuals_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_evaluation_residuals_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_residuals_path(self):
        return fs.create_directory_in(self.evaluation_path, "proper_residuals")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_residuals_names(self):
        return fs.files_in_path(self.evaluation_proper_residuals_path, returns="name", extension="fits")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_residuals_filters(self):
        return [parse_filter(name) for name in self.evaluation_proper_residuals_names]

    # -----------------------------------------------------------------

    def has_evaluation_proper_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_proper_residuals_filters

    # -----------------------------------------------------------------

    def get_evaluation_proper_residuals_path_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_proper_residuals_path, str(fltr) + ".fits")

    # -----------------------------------------------------------------

    def get_evaluation_proper_residuals_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Frame.from_file(self.get_evaluation_proper_residuals_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_residuals_distribution_names(self):
        return fs.files_in_path(self.evaluation_proper_residuals_path, returns="name", extension="dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_residuals_distributions_filters(self):
        return [parse_filter(name.split("_distribution")[0]) for name in self.evaluation_residuals_distribution_names]

    # -----------------------------------------------------------------

    def has_evaluation_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_residuals_distributions_filters

    # -----------------------------------------------------------------

    def get_evaluation_residuals_distribution_path_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_proper_residuals_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def get_evaluation_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_evaluation_residuals_distribution_path_for_filter(fltr))

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_residuals_distribution_names(self):
        return fs.files_in_path(self.evaluation_proper_residuals_path, returns="name", extension="dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def evaluation_proper_residuals_distribution_filters(self):
        return [parse_filter(name.split("_distribution")[0]) for name in self.evaluation_proper_residuals_distribution_names]

    # -----------------------------------------------------------------

    def has_evaluation_proper_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return fltr in self.evaluation_proper_residuals_distribution_filters

    # -----------------------------------------------------------------

    def get_evaluation_proper_residuals_distribution_path_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return fs.join(self.evaluation_proper_residuals_path, str(fltr) + "_distribution.dat")

    # -----------------------------------------------------------------

    def get_evaluation_proper_residuals_distribution_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        return Distribution.from_file(self.get_evaluation_proper_residuals_distribution_path_for_filter(fltr))

# -----------------------------------------------------------------

class AnalysisRun(AnalysisRunBase):

    """
    This class ...
    """

    def __init__(self, galaxy_name=None, info=None, hubble_stage=None):

        """
        The constructor ...
        :param galaxy_name:
        :param info:
        :param hubble_stage:
        """

        # Set the analysis run info
        self.info = info

        # Set galaxy properties
        self.galaxy_name = galaxy_name
        self.hubble_stage = hubble_stage

        ## Create directories

        # The directory for the model
        if not fs.is_directory(self.model_path): fs.create_directory(self.model_path)

        # The directory for the projections and the instruments
        if not fs.is_directory(self.projections_path): fs.create_directory(self.projections_path)
        if not fs.is_directory(self.instruments_path): fs.create_directory(self.instruments_path)

        # The directory for the dust grid output
        if not fs.is_directory(self.dust_grid_build_path): fs.create_directory(self.dust_grid_build_path)
        if not fs.is_directory(self.dust_grid_simulation_out_path): fs.create_directory(self.dust_grid_simulation_out_path)

        # Contributions directory
        if not fs.is_directory(self.contributions_path): fs.create_directory(self.contributions_path)

        # Evaluation
        if not fs.is_directory(self.evaluation_path): fs.create_directory(self.evaluation_path)

        # Analysis directories
        if not fs.is_directory(self.fluxes_path): fs.create_directory(self.fluxes_path)
        if not fs.is_directory(self.images_path): fs.create_directory(self.images_path)
        if not fs.is_directory(self.residuals_path): fs.create_directory(self.residuals_path)
        if not fs.is_directory(self.properties_path): fs.create_directory(self.properties_path)
        if not fs.is_directory(self.attenuation_path): fs.create_directory(self.attenuation_path)
        if not fs.is_directory(self.colours_path): fs.create_directory(self.colours_path)
        if not fs.is_directory(self.maps_path): fs.create_directory(self.maps_path)
        if not fs.is_directory(self.absorption_path): fs.create_directory(self.absorption_path)
        if not fs.is_directory(self.heating_path): fs.create_directory(self.heating_path)
        if not fs.is_directory(self.energy_path): fs.create_directory(self.energy_path)
        if not fs.is_directory(self.sfr_path): fs.create_directory(self.sfr_path)
        if not fs.is_directory(self.correlations_path): fs.create_directory(self.correlations_path)

        # Maps subdirectories
        if not fs.is_directory(self.colour_maps_path): fs.create_directory(self.colour_maps_path)
        if not fs.is_directory(self.ssfr_maps_path): fs.create_directory(self.ssfr_maps_path)
        if not fs.is_directory(self.tir_maps_path): fs.create_directory(self.tir_maps_path)
        if not fs.is_directory(self.attenuation_maps_path): fs.create_directory(self.attenuation_maps_path)
        if not fs.is_directory(self.old_maps_path): fs.create_directory(self.old_maps_path)
        if not fs.is_directory(self.dust_maps_path): fs.create_directory(self.dust_maps_path)
        if not fs.is_directory(self.young_maps_path): fs.create_directory(self.young_maps_path)
        if not fs.is_directory(self.ionizing_maps_path): fs.create_directory(self.ionizing_maps_path)
        if not fs.is_directory(self.rgb_maps_path): fs.create_directory(self.rgb_maps_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_name(cls, modeling_path, name, hubble_stage=None):

        """
        This function ...
        :param modeling_path:
        :param name:
        :param hubble_stage:
        :return:
        """

        analysis_path = fs.join(modeling_path, "analysis")
        run_path = fs.join(analysis_path, name)
        return cls.from_path(run_path, hubble_stage=hubble_stage)

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, path, hubble_stage=None):

        """
        This function ...
        :param path:
        :param hubble_stage:
        :return:
        """

        # Determine the info path
        info_path = fs.join(path, info_filename)
        if not fs.is_file(info_path): raise IOError("Could not find the info file")
        else: return cls.from_info(info_path, hubble_stage=hubble_stage)

    # -----------------------------------------------------------------

    @classmethod
    def from_info(cls, info_path, hubble_stage=None):

        """
        This function ...
        :param info_path:
        :param hubble_stage:
        :return:
        """

        # Load the analysis run info
        info = AnalysisRunInfo.from_file(info_path)

        # Create the instance
        run = cls(info=info, hubble_stage=hubble_stage)

        # Set galaxy name
        modeling_path = fs.directory_of(fs.directory_of(run.info.path))
        run.galaxy_name = fs.name(modeling_path)

        # Return the analysis run object
        return run

    # -----------------------------------------------------------------

    @property
    def has_dust_grid_simulation_logfile(self):
        return fs.is_file(self.dust_grid_simulation_logfile_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_simulation_logfile(self):
        return LogFile.from_file(self.dust_grid_simulation_logfile_path)

    # -----------------------------------------------------------------

    @property
    def has_output(self):
        return fs.has_files_in_path(self.output_path)

    # -----------------------------------------------------------------

    @property
    def has_logfile(self):
        return fs.is_file(self.logfile_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def logfile(self):
        return LogFile.from_file(self.total_logfile_path)

    # -----------------------------------------------------------------

    @property
    def has_misc(self):
        return fs.has_files_in_path(self.total_misc_path)

    # -----------------------------------------------------------------

    @property
    def has_extracted(self):
        return fs.has_files_in_path(self.total_extr_path)

    # -----------------------------------------------------------------

    @property
    def has_progress(self):
        return fs.is_file(self.progress_path)

    # -----------------------------------------------------------------

    @property
    def has_timeline(self):
        return fs.is_file(self.timeline_path)

    # -----------------------------------------------------------------

    @property
    def has_memory(self):
        return fs.is_file(self.memory_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def progress(self):
        return ProgressTable.from_file(self.progress_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def timeline(self):
        return TimeLineTable.from_file(self.timeline_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def memory(self):
        return MemoryUsageTable.from_file(self.memory_path)

    # -----------------------------------------------------------------

    @property
    def has_plots(self):
        return fs.has_files_in_path(self.plot_path)

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):
        return fs.is_directory(self.attenuation_path) and not fs.is_empty(self.attenuation_path)

    # -----------------------------------------------------------------

    @property
    def has_colours(self):
        return fs.is_directory(self.colours_path) and not fs.is_empty(self.colours_path)

    # -----------------------------------------------------------------

    @property
    def colour_names(self):
        return fs.files_in_path(self.colours_simulated_path, extension="fits", returns="name")

    # -----------------------------------------------------------------

    @property
    def has_residuals(self):
        return fs.is_directory(self.residuals_path) and not fs.is_empty(self.residuals_path)

    # -----------------------------------------------------------------

    @property
    def residual_image_names(self):
        return fs.files_in_path(self.residuals_path, extension="fits", not_contains=["significance"], returns="name")

    # -----------------------------------------------------------------

    @property
    def has_maps_attenuation(self):
        return fs.is_directory(self.attenuation_maps_path) and not fs.is_empty(self.attenuation_maps_path)

    # -----------------------------------------------------------------

    @property
    def nattenuation_maps(self):
        if fs.has_files_in_path(self.attenuation_maps_path, extension="fits"): return fs.nfiles_in_path(self.attenuation_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.attenuation_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_colours(self):
        return fs.is_directory(self.colour_maps_path) and not fs.is_empty(self.colour_maps_path)

    # -----------------------------------------------------------------

    @property
    def ncolour_maps(self):
        if fs.has_files_in_path(self.colour_maps_path, extension="fits"): return fs.nfiles_in_path(self.colour_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.colour_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_dust(self):
        return fs.is_directory(self.dust_maps_path) and not fs.is_empty(self.dust_maps_path)

    # -----------------------------------------------------------------

    @property
    def ndust_maps(self):
        if fs.has_files_in_path(self.dust_maps_path, extension="fits"): return fs.nfiles_in_path(self.dust_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.dust_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_ionizing(self):
        return fs.is_directory(self.ionizing_maps_path) and not fs.is_empty(self.ionizing_maps_path)

    # -----------------------------------------------------------------

    @property
    def nionizing_maps(self):
        if fs.has_files_in_path(self.ionizing_maps_path, extension="fits"): return fs.nfiles_in_path(self.ionizing_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.ionizing_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_rgb(self):
        return fs.is_directory(self.rgb_maps_path) and not fs.is_empty(self.rgb_maps_path)

    # -----------------------------------------------------------------

    @property
    def has_maps_old(self):
        return fs.is_directory(self.old_maps_path) and not fs.is_empty(self.old_maps_path)

    # -----------------------------------------------------------------

    @property
    def nold_maps(self):
        if fs.has_files_in_path(self.old_maps_path, extension="fits"): return fs.nfiles_in_path(self.old_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.old_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_ssfr(self):
        return fs.is_directory(self.ssfr_maps_path) and not fs.is_empty(self.ssfr_maps_path)

    # -----------------------------------------------------------------

    @property
    def nssfr_maps(self):
        if fs.has_files_in_path(self.ssfr_maps_path, extension="fits"): return fs.nfiles_in_path(self.ssfr_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.ssfr_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_tir(self):
        return fs.is_directory(self.tir_maps_path) and not fs.is_empty(self.tir_maps_path)

    # -----------------------------------------------------------------

    @property
    def ntir_maps(self):
        if fs.has_files_in_path(self.tir_maps_path, extension="fits"): return fs.nfiles_in_path(self.tir_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.tir_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_young(self):
        return fs.is_directory(self.young_maps_path) and not fs.is_empty(self.young_maps_path)

    # -----------------------------------------------------------------

    @property
    def nyoung_maps(self):
        if fs.has_files_in_path(self.young_maps_path, extension="fits"): return fs.nfiles_in_path(self.young_maps_path, extension="fits")
        else: return fs.nfiles_in_path(self.young_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_heating(self):
        return fs.is_file(self.heating_config_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def nfiles(self):
        return fs.nfiles_in_path(self.path, recursive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_space(self):
        return fs.directory_size(self.path)

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):
        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):
        return fs.directory_of(self.analysis_path)

    # -----------------------------------------------------------------

    @property
    def path(self):
        return self.info.path

    # -----------------------------------------------------------------

    @lazyproperty
    def config(self):
        return Configuration.from_file(self.config_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_config(self):
        return Configuration.from_file(self.heating_config_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):
        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):
        return load_grid(self.dust_grid_path)

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):
        return self.info.path

    # -----------------------------------------------------------------

    @property
    def ski_file(self):
        return SkiFile(self.ski_file_path)

    # -----------------------------------------------------------------

    @property
    def input_paths(self):
        return load_dict(self.input_file_path)

    # -----------------------------------------------------------------

    @property
    def dust_grid_tree_path(self):
        return fs.join(self.dust_grid_build_path, dust_grid_tree_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):

        """
        This function ...
        :return:
        """

        # Give debug message
        log.debug("Loading the dust grid tree, this may take a while (depending on the number of nodes) ...")

        # Return the tree
        return DustGridTree.from_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    def create_file_tree_dust_grid(self, search_method="Neighbor", write=False):

        """
        This function ...
        :param search_method:
        :param write:
        :return:
        """

        grid = FileTreeDustGrid(filename=self.dust_grid_tree_path, search_method=search_method, write=write)
        return grid

    # -----------------------------------------------------------------

    @lazyproperty
    def has_dust_grid_tree(self):
        return fs.is_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @property
    def model_path(self):
        return fs.join(self.path, model_name)

    # -----------------------------------------------------------------

    @property
    def instruments_path(self):
        return fs.join(self.path, instruments_name)

    # -----------------------------------------------------------------

    @property
    def projections_path(self):
        return fs.join(self.path, projections_name)

    # -----------------------------------------------------------------

    @property
    def earth_projection_path(self):
        return fs.join(self.projections_path, earth_projection_filename)

    # -----------------------------------------------------------------

    @property
    def faceon_projection_path(self):
        return fs.join(self.projections_path, faceon_projection_filename)

    # -----------------------------------------------------------------

    @property
    def edgeon_projection_path(self):
        return fs.join(self.projections_path, edgeon_projection_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_projection(self):
        return GalaxyProjection.from_file(self.earth_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):
        return EdgeOnProjection.from_file(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):
        return FaceOnProjection.from_file(self.faceon_projection_path)

    # -----------------------------------------------------------------

    @property
    def sed_earth_instrument_path(self):
        return fs.join(self.instruments_path, sed_earth_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_sed_earth_instrument(self):
        return fs.is_file(self.sed_earth_instrument_path)

    # -----------------------------------------------------------------

    @property
    def full_sed_earth_instrument_path(self):
        return fs.join(self.instruments_path, full_sed_earth_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_full_sed_earth_instrument(self):
        return fs.is_file(self.full_sed_earth_instrument_path)

    # -----------------------------------------------------------------

    @property
    def simple_earth_instrument_path(self):
        return fs.join(self.instruments_path, simple_earth_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_simple_earth_instrument(self):
        return fs.is_file(self.simple_earth_instrument_path)

    # -----------------------------------------------------------------

    @property
    def full_earth_instrument_path(self):
        return fs.join(self.instruments_path, full_earth_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_full_earth_instrument(self):
        return fs.is_file(self.full_earth_instrument_path)

    # -----------------------------------------------------------------

    @property
    def simple_faceon_instrument_path(self):
        return fs.join(self.instruments_path, simple_faceon_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_simple_faceon_instrument(self):
        return fs.is_file(self.simple_faceon_instrument_path)

    # -----------------------------------------------------------------

    @property
    def full_faceon_instrument_path(self):
        return fs.join(self.instruments_path, full_faceon_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_full_faceon_instrument(self):
        return fs.is_file(self.full_faceon_instrument_path)

    # -----------------------------------------------------------------

    @property
    def simple_edgeon_instrument_path(self):
        return fs.join(self.instruments_path, simple_edgeon_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_simple_edgeon_instrument(self):
        return fs.is_file(self.simple_edgeon_instrument_path)

    # -----------------------------------------------------------------

    @property
    def full_edgeon_instrument_path(self):
        return fs.join(self.instruments_path, full_edgeon_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def has_full_edgeon_instrument(self):
        return fs.is_file(self.full_edgeon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def sed_earth_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_sed_earth_instrument:
            instrument = SEDInstrument.from_projection(self.earth_projection)
            instrument.saveto(self.sed_earth_instrument_path)

        return SEDInstrument.from_file(self.sed_earth_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_earth_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_simple_earth_instrument:
            instrument = SimpleInstrument.from_projection(self.earth_projection)
            instrument.saveto(self.sed_earth_instrument_path)

        return SimpleInstrument.from_file(self.simple_earth_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_earth_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_full_earth_instrument:
            instrument = FullInstrument.from_projection(self.earth_projection)
            instrument.saveto(self.full_earth_instrument_path)

        return FullInstrument.from_file(self.full_earth_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_faceon_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_simple_faceon_instrument:
            instrument = SimpleInstrument.from_projection(self.faceon_projection)
            instrument.saveto(self.simple_faceon_instrument_path)

        return SimpleInstrument.from_file(self.simple_faceon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_faceon_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_full_faceon_instrument:
            instrument = FullInstrument.from_projection(self.faceon_projection)
            instrument.saveto(self.full_faceon_instrument_path)

        return FullInstrument.from_file(self.full_faceon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simple_edgeon_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_simple_edgeon_instrument:
            instrument = SimpleInstrument.from_projection(self.edgeon_projection)
            instrument.saveto(self.simple_edgeon_instrument_path)

        return SimpleInstrument.from_file(self.simple_edgeon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def full_edgeon_instrument(self):

        """
        This function ...
        :return:
        """

        if not self.has_full_edgeon_instrument:
            instrument = FullInstrument.from_projection(self.edgeon_projection)
            instrument.saveto(self.full_edgeon_instrument_path)

        return FullInstrument.from_file(self.full_edgeon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties_path(self):
        from ..core.environment import properties_name, data_name
        return fs.join(self.modeling_path, data_name, properties_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_properties(self):
        return GalaxyProperties.from_file(self.galaxy_properties_path)

    # -----------------------------------------------------------------

    @property
    def galaxy_distance(self):
        return self.galaxy_properties.distance

    # -----------------------------------------------------------------

    @property
    def galaxy_center(self):
        return self.galaxy_properties.center

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_path(self):
        return fs.join(self.modeling_path, "truncated")

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse_path(self):
        return fs.join(self.truncation_path, "ellipse.reg")

    # -----------------------------------------------------------------

    @lazyproperty
    def truncation_ellipse(self):
        return SkyRegionList.from_file(self.truncation_ellipse_path)[0]

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run_name(self):
        return self.info.fitting_run

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):
        return FittingRun.from_name(self.modeling_path, self.fitting_run_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_suite(self):
        from ..build.suite import ModelSuite
        return ModelSuite.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def model_definition_path(self):
        return self.model_suite.get_model_path(self.model_name)

    # -----------------------------------------------------------------

    @property
    def model_definition_stellar_path(self):
        return fs.join(self.model_definition_path, "stellar")

    # -----------------------------------------------------------------

    @property
    def model_definition_dust_path(self):
        return fs.join(self.model_definition_path, "dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def model_stellar_path(self):
        return fs.create_directory_in(self.model_path, "stellar")

    # -----------------------------------------------------------------

    @lazyproperty
    def model_dust_path(self):
        return fs.create_directory_in(self.model_path, "dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def stellar_component_paths(self):
        # Return the paths as a dictionary based on component name
        return fs.directories_in_path(self.model_stellar_path, returns="dict")

    # -----------------------------------------------------------------

    @property
    def stellar_component_names(self):
        return self.stellar_component_paths.keys()

    # -----------------------------------------------------------------

    @property
    def nstellar_components(self):
        return len(self.stellar_component_names)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_component_paths(self):
        # Return the paths as a dictionary based on component name
        return fs.directories_in_path(self.model_dust_path, returns="dict")

    # -----------------------------------------------------------------

    @property
    def dust_component_names(self):
        return self.dust_component_paths.keys()

    # -----------------------------------------------------------------

    @property
    def ndust_components(self):
        return len(self.dust_component_names)

    # -----------------------------------------------------------------

    @property
    def has_model(self):

        """
        This function ...
        :return:
        """

        if self.nstellar_components > 0 and self.ndust_components > 0: return True
        elif self.nstellar_components == 0 and self.ndust_components == 0: return False
        else: raise RuntimeError("Inconsistent state of the model: there are stellar components but no dust components or vice versa")

    # -----------------------------------------------------------------

    @lazyproperty
    def model_definition(self):

        """
        This function ...
        :return:
        """

        from .initialization import create_model_definition_in_path

        # Has the model definition been created yet?
        if not self.has_model: return create_model_definition_in_path(self.model_name, self.model_path, self.model_definition_stellar_path, self.model_definition_dust_path, parameter_values=self.parameter_values)

        # Create the model definition and return
        else: return ModelDefinition(model_name, self.model_path, stellar_paths=self.stellar_component_paths, dust_paths=self.dust_component_paths)

    # -----------------------------------------------------------------

    @property
    def model_old_map_name(self):
        return self.model_suite.get_old_map_name_for_model(self.model_name)

    # -----------------------------------------------------------------

    @property
    def model_young_map_name(self):
        return self.model_suite.get_young_map_name_for_model(self.model_name)

    # -----------------------------------------------------------------

    @property
    def model_ionizing_map_name(self):
        return self.model_suite.get_ionizing_map_name_for_model(self.model_name)

    # -----------------------------------------------------------------

    @property
    def model_dust_map_name(self):
        return self.model_suite.get_dust_map_name_for_model(self.model_name)

    # -----------------------------------------------------------------

    @property
    def model_old_map(self):
        return self.model_suite.load_stellar_component_map(self.model_name, "old")

    # -----------------------------------------------------------------

    @property
    def model_young_map(self):
        return self.model_suite.load_stellar_component_map(self.model_name, "young")

    # -----------------------------------------------------------------

    @property
    def model_ionizing_map(self):
        return self.model_suite.load_stellar_component_map(self.model_name, "ionizing")

    # -----------------------------------------------------------------

    @property
    def model_dust_map(self):
        return self.model_suite.load_dust_component_map(self.model_name, "disk")

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_name(self):
        return self.info.generation_name

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_name(self):
        return self.info.simulation_name

    # -----------------------------------------------------------------

    @lazyproperty
    def parameter_values(self):

        """
        This function ...
        :return:
        """

        # Get the ski file
        ski = self.ski_file

        # Get the values of all the labeled parameters
        values = ski.get_labeled_values()

        # Return the parameter values
        return values

    # -----------------------------------------------------------------

    @property
    def free_parameter_labels(self):
        return self.parameter_values.keys()

    # -----------------------------------------------------------------

    @lazyproperty
    def chi_squared(self):
        return self.info.chi_squared

    # -----------------------------------------------------------------

    @lazyproperty
    def model(self):

        """
        This function ...
        :return:
        """

        # Create the model and return it
        return RTModel(self.model_definition, simulation_name=self.simulation_name, chi_squared=self.chi_squared,
                      free_parameter_labels=self.free_parameter_labels, wavelength_grid=self.wavelength_grid,
                       observed_total_output_path=self.total_output_path, observed_bulge_output_path=self.bulge_output_path,
                       observed_disk_output_path=self.disk_output_path, observed_old_output_path=self.old_output_path,
                       observed_young_output_path=self.young_output_path, observed_sfr_output_path=self.ionizing_output_path,
                       observed_unevolved_output_path=self.unevolved_output_path,
                       observed_extra_output_path=self.extra_output_path, center=self.galaxy_center,
                       galaxy_name=self.galaxy_name, hubble_stage=self.hubble_stage, earth_wcs=self.reference_wcs,
                       truncation_ellipse=self.truncation_ellipse)

    # -----------------------------------------------------------------

    @lazyproperty
    def mappings(self):
        return self.model.mappings

    # -----------------------------------------------------------------

    @lazyproperty
    def normalized_mappings(self):
        return self.model.normalized_mappings

    # -----------------------------------------------------------------

    @property
    def uses_grid_resolution(self):
        return self.info.reference_deprojection == "grid"

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection_component_name(self):
        if self.uses_grid_resolution: return None
        return self.info.reference_deprojection

    # -----------------------------------------------------------------

    @lazyproperty
    def is_stellar_reference_deprojection(self):
        if self.uses_grid_resolution: raise ValueError("This function shouldn't be called")
        return self.reference_deprojection_component_name in self.model_suite.get_stellar_component_names(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def is_dust_reference_deprojection(self):
        if self.uses_grid_resolution: raise ValueError("This function shouldn't be called")
        return self.reference_deprojection_component_name in self.model_suite.get_dust_component_names(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection_component(self):
        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component(self.model_name, self.reference_deprojection_component_name, add_map=False)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component(self.model_name, self.reference_deprojection_component_name, add_map=False)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection(self):
        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component_deprojection(self.model_name, self.reference_deprojection_component_name, add_map=False)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component_deprojection(self.model_name, self.reference_deprojection_component_name, add_map=False)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map(self):
        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component_map(self.model_name, self.reference_deprojection_component_name)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component_map(self.model_name, self.reference_deprojection_component_name)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map_path(self):
        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.get_stellar_component_map_path(self.model_name, self.reference_deprojection_component_name)
            elif self.is_dust_reference_deprojection: return self.model_suite.get_dust_component_map_path(self.model_name, self.reference_deprojection_component_name)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):
        if self.reference_map_path is None: return None
        else: return CoordinateSystem.from_file(self.reference_map_path)

    # -----------------------------------------------------------------

    @property
    def remote_script_paths(self):
        return fs.files_in_path(self.path, extension="sh")

    # -----------------------------------------------------------------

    def get_remote_script_commands(self):

        """
        This fucntion ...
        :return:
        """

        commands = dict()

        # Loop over the script paths
        for path in self.remote_script_paths:

            # Get host ID
            host_id = fs.strip_extension(fs.name(path))

            lines = []

            # Open the file
            for line in fs.read_lines(path):

                if line.startswith("#"): continue
                if not line.strip(): continue

                lines.append(line)

            # Set the commands
            commands[host_id] = lines

        # Return the commands
        return commands

    # -----------------------------------------------------------------

    def get_remote_script_commands_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        commands = self.get_remote_script_commands()
        if host_id in commands: return commands[host_id]
        else: return []

    # -----------------------------------------------------------------

    def get_ski_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return SkiFile(self.ski_path_for_contribution(contribution))

    # -----------------------------------------------------------------

    @property
    def has_simulated_sed(self):
        return fs.is_file(self.simulated_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_sed(self):
        return SED.from_skirt(self.simulated_sed_path)

    # -----------------------------------------------------------------

    @property
    def has_simulated_fluxes(self):
        return fs.is_file(self.simulated_fluxes_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_fluxes(self):
        return ObservedSED.from_file(self.simulated_fluxes_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_datacube(self):

        """
        This function ...
        :return:
        """

        # Load the datacube
        datacube = DataCube.from_file(self.simulated_datacube_path, self.wavelength_grid)

        # Set the wcs
        datacube.wcs = self.reference_wcs

        # Return the datacube
        return datacube

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_dataset(self):

        """
        This function ...
        :return:
        """

        #get_name_function = lambda filename: filename.split("__")[1]
        #return DataSet.from_directory(self.total_misc_path, get_name=get_name_function)
        return DataSet.from_directory(self.images_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_frame_list(self):
        return self.simulated_dataset.get_framelist(named=False)  # on filter

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_named_frame_list(self):
        return self.simulated_dataset.get_framelist(named=True)  # on name

    # -----------------------------------------------------------------

    def get_simulated_frame_for_filter(self, fltr):

        """
        THis function ...
        :param fltr:
        :return:
        """

        # Return the simulated frame
        return self.simulated_frame_list[fltr]

    # -----------------------------------------------------------------

    @lazyproperty
    def maps_collection(self):
        from ..maps.collection import MapsCollection
        return MapsCollection.from_modeling_path(self.modeling_path, analysis_run_name=self.name)

    # -----------------------------------------------------------------

    @lazyproperty
    def observation_maps_collection(self):
        from ..maps.collection import MapsCollection
        return MapsCollection.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @property
    def colours_methods(self):
        return self.observation_maps_collection.get_colours_methods(flatten=True)

    # -----------------------------------------------------------------

    @property
    def colours_origins(self):
        return self.observation_maps_collection.get_colours_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_methods(self):
        return self.observation_maps_collection.get_ssfr_methods(flatten=True)

    # -----------------------------------------------------------------

    @property
    def ssfr_origins(self):
        return self.observation_maps_collection.get_ssfr_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def tir_methods(self):
        return self.observation_maps_collection.get_tir_methods(flatten=True)

    # -----------------------------------------------------------------

    @property
    def tir_origins(self):
        return self.observation_maps_collection.get_tir_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_methods(self):
        return self.observation_maps_collection.get_attenuation_methods(flatten=True)

    # -----------------------------------------------------------------

    @property
    def attenuation_origins(self):
        return self.observation_maps_collection.get_attenuation_origins(flatten=True)

    # -----------------------------------------------------------------

    @property
    def old_methods(self):
        return self.observation_maps_collection.get_old_methods(flatten=False)

    # -----------------------------------------------------------------

    @property
    def old_map_methods(self):
        #return self.old_methods[self.model_old_map_name] # for flattened
        #return find_value_for_unique_key_nested(self.old_methods, self.model_old_map_name)
        return self.old_methods["disk"][self.model_old_map_name]

    # -----------------------------------------------------------------

    @property
    def old_origins(self):
        return self.observation_maps_collection.get_old_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def old_map_origins(self):
        #return self.old_origins[self.model_old_map_name] # for flattened
        #return find_value_for_unique_key_nested(self.old_origins, self.model_old_map_name)
        return self.old_origins["disk"][self.model_old_map_name]

    # -----------------------------------------------------------------

    @property
    def old_map_method_and_name(self):
        return "disk", self.model_old_map_name

    # -----------------------------------------------------------------

    @property
    def young_methods(self):
        return self.observation_maps_collection.get_young_methods(flatten=False)

    # -----------------------------------------------------------------

    @property
    def young_map_methods(self):
        #return self.young_methods[self.model_young_map_name]
        return find_value_for_unique_key_nested(self.young_methods, self.model_young_map_name)

    # -----------------------------------------------------------------

    @property
    def young_origins(self):
        return self.observation_maps_collection.get_young_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def young_map_origins(self):
        #return self.young_origins[self.model_young_map_name]
        return find_value_for_unique_key_nested(self.young_origins, self.model_young_map_name)

    # -----------------------------------------------------------------

    @property
    def young_map_method_and_name(self):

        """
        This function ...
        :return:
        """

        keys = find_keys_for_unique_key_nested(self.young_methods, self.model_young_map_name)
        if len(keys) == 1:
            method = None
            map_name = keys[0]
        elif len(keys) == 2:
            method = keys[0]
            map_name = keys[1]
        else: raise ValueError("Something is wrong")
        return method, map_name

    # -----------------------------------------------------------------

    @property
    def ionizing_methods(self):
        return self.observation_maps_collection.get_ionizing_methods(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_methods(self):
        #return self.ionizing_methods[self.model_ionizing_map_name]
        return find_value_for_unique_key_nested(self.ionizing_methods, self.model_ionizing_map_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_origins(self):
        return self.observation_maps_collection.get_ionizing_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_origins(self):
        #return self.ionizing_origins[self.model_ionizing_map_name]
        return find_value_for_unique_key_nested(self.ionizing_origins, self.model_ionizing_map_name)

    # -----------------------------------------------------------------

    @property
    def ionizing_map_method_and_name(self):

        """
        This function ...
        :return:
        """

        keys = find_keys_for_unique_key_nested(self.ionizing_methods, self.model_ionizing_map_name)
        if len(keys) == 1:
            method = None
            map_name = keys[0]
        elif len(keys) == 2:
            method = keys[0]
            map_name = keys[1]
        else: raise ValueError("Something is wrong")
        return method, map_name

    # -----------------------------------------------------------------

    @property
    def dust_methods(self):
        return self.observation_maps_collection.get_dust_methods(flatten=False)

    # -----------------------------------------------------------------

    @property
    def dust_map_methods(self):
        #return self.dust_methods[self.model_dust_map_name]
        try: return find_value_for_unique_key_nested(self.dust_methods, self.model_dust_map_name)
        except ValueError: return find_value_for_unique_key_nested(self.dust_methods, self.model_dust_map_name.split("_", 1)[1])

    # -----------------------------------------------------------------

    @property
    def dust_origins(self):
        return self.observation_maps_collection.get_dust_origins(flatten=False)

    # -----------------------------------------------------------------

    @property
    def dust_map_origins(self):

        """
        This function ...
        :return:
        """

        #return self.dust_origins[self.model_dust_map_name]
        try: return find_value_for_unique_key_nested(self.dust_origins, self.model_dust_map_name)
        except ValueError: return find_value_for_unique_key_nested(self.dust_origins, self.model_dust_map_name.split("_", 1)[1])

    # -----------------------------------------------------------------

    @property
    def dust_map_method_and_name(self):

        """
        This function ...
        :return:
        """

        try: keys = find_keys_for_unique_key_nested(self.dust_methods, self.model_dust_map_name)
        except ValueError: keys = find_keys_for_unique_key_nested(self.dust_methods, self.model_dust_map_name.split("_", 1)[1])

        if len(keys) == 1:
            method = None
            map_name = keys[0]
        elif len(keys) == 2:
            method = keys[0]
            map_name = keys[1]
        else: raise ValueError("Something is wrong")
        return method, map_name

# -----------------------------------------------------------------

class AnalysisRuns(object):

    """
    This function ...
    """

    def __init__(self, modeling_path):

        """
        This function ...
        :param modeling_path:
        """

        # Set the modeling path
        self.modeling_path = modeling_path

    # -----------------------------------------------------------------

    @property
    def modeling_data_path(self):
        return fs.join(self.modeling_path, "data")

    # -----------------------------------------------------------------

    @property
    def galaxy_info_path(self):
        #  Set the path to the galaxy info file
        return fs.join(self.modeling_data_path, "info.dat")

    # -----------------------------------------------------------------

    @lazyproperty
    def galaxy_info(self):

        """
        This function ...
        :return:
        """

        # Load the info table
        table = tables.from_file(self.galaxy_info_path)

        # To ordered dict
        info = OrderedDict()
        for name in table.colnames: info[name] = table[name][0]

        # Return the info
        return info

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_type(self):
        return self.galaxy_info["Hubble Type"]

    # -----------------------------------------------------------------

    @lazyproperty
    def hubble_stage(self):
        return self.galaxy_info["Hubble Stage"]

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):
        return fs.join(self.modeling_path, "analysis")

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):
        return fs.directories_in_path(self.analysis_path, returns="name")

    # -----------------------------------------------------------------

    @lazyproperty
    def paths(self):
        return fs.directories_in_path(self.analysis_path, returns="path")

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def empty(self):
        return sequences.is_empty(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single(self):
        return sequences.is_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_name(self):
        return sequences.get_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_path(self):
        return self.get_path(self.single_name)

    # -----------------------------------------------------------------

    def get_path(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return fs.join(self.analysis_path, name)

    # -----------------------------------------------------------------

    def load(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        analysis_run_path = self.get_path(name)
        if not fs.is_directory(analysis_run_path): raise ValueError("Analysis run '" + name + "' does not exist")
        return AnalysisRun.from_path(analysis_run_path, hubble_stage=self.hubble_stage)

    # -----------------------------------------------------------------

    @lazyproperty
    def single(self):
        return AnalysisRun.from_path(self.single_path, hubble_stage=self.hubble_stage)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_name(self):

        """
        This function ...
        :return:
        """

        #if self.empty: return None
        #if self.has_single: return self.single_name
        #return sorted(self.names)[-1]

        return fs.name(self.last_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_path(self):
        return fs.last_created_path(*self.paths)

    # -----------------------------------------------------------------

    @lazyproperty
    def last(self):
        return self.load(self.last_name)

# -----------------------------------------------------------------

class CachedAnalysisRun(AnalysisRunBase):

    """
    This class ...
    """

    def __init__(self, run_path, remote):

        """
        The constructor ...
        :param run_path:
        :param remote:
        """

        # Set attributes
        self.path = run_path
        self.remote = remote

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):
        return fs.name(self.original_modeling_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, path, remote):

        """
        This function ...
        :param path:
        :param remote:
        :return:
        """

        return cls(path, remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def info(self):
        return AnalysisRunInfo.from_remote_file(self.info_path, self.remote)

    # -----------------------------------------------------------------

    @property
    def original_path(self):
        return self.info.path

    # -----------------------------------------------------------------

    @property
    def original_analysis_path(self):
        return fs.directory_of(self.original_path)

    # -----------------------------------------------------------------

    @property
    def original_modeling_path(self):
        return fs.directory_of(self.original_analysis_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def config(self):
        return Configuration.from_remote_file(self.config_path, self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def heating_config(self):
        return Configuration.from_remote_file(self.heating_config_path, self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):
        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path, remote=self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):
        return load_grid(self.dust_grid_path, remote=self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def nfiles(self):
        return self.remote.nfiles_in_path(self.path, recursive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def disk_space(self):
        return self.remote.directory_size(self.path)

    # -----------------------------------------------------------------

    @property
    def has_output(self):
        return self.remote.has_files_in_path(self.output_path)

    # -----------------------------------------------------------------

    @property
    def has_logfile(self):
        return self.remote.is_file(self.logfile_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def logfile(self):
        return LogFile.from_remote_file(self.total_logfile_path, self.remote)

    # -----------------------------------------------------------------

    @property
    def has_misc(self):
        return self.remote.has_files_in_path(self.total_misc_path)

    # -----------------------------------------------------------------

    @property
    def has_extracted(self):
        return self.remote.has_files_in_path(self.total_extr_path)

    # -----------------------------------------------------------------

    @property
    def has_progress(self):
        return self.remote.is_file(self.progress_path)

    # -----------------------------------------------------------------

    @property
    def has_timeline(self):
        return self.remote.is_file(self.timeline_path)

    # -----------------------------------------------------------------

    @property
    def has_memory(self):
        return self.remote.is_file(self.memory_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def progress(self):
        return ProgressTable.from_remote_file(self.progress_path, remote=self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def timeline(self):
        return TimeLineTable.from_remote_file(self.timeline_path, remote=self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def memory(self):
        return MemoryUsageTable.from_remote_file(self.memory_path, remote=self.remote)

    # -----------------------------------------------------------------

    @property
    def has_plots(self):
        return self.remote.has_files_in_path(self.plot_path)

    # -----------------------------------------------------------------

    @property
    def has_attenuation(self):
        return self.remote.is_directory(self.attenuation_path) and not self.remote.is_empty(self.attenuation_path)

    # -----------------------------------------------------------------

    @property
    def has_colours(self):
        return self.remote.is_directory(self.colours_path) and not self.remote.is_empty(self.colours_path)

    # -----------------------------------------------------------------

    @property
    def colour_names(self):
        return self.remote.files_in_path(self.colours_simulated_path, extension="fits", returns="name")

    # -----------------------------------------------------------------

    @property
    def has_residuals(self):
        return self.remote.is_directory(self.residuals_path) and not self.remote.is_empty(self.residuals_path)

    # -----------------------------------------------------------------

    @property
    def residual_image_names(self):
        return self.remote.files_in_path(self.residuals_path, extension="fits", not_contains=["significance"], returns="name")

    # -----------------------------------------------------------------

    @property
    def has_maps_attenuation(self):
        return self.remote.is_directory(self.attenuation_maps_path) and not self.remote.is_empty(self.attenuation_maps_path)

    # -----------------------------------------------------------------

    @property
    def nattenuation_maps(self):
        if self.remote.has_files_in_path(self.attenuation_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.attenuation_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.attenuation_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_colours(self):
        return self.remote.is_directory(self.colour_maps_path) and not self.remote.is_empty(self.colour_maps_path)

    # -----------------------------------------------------------------

    @property
    def ncolour_maps(self):
        if self.remote.has_files_in_path(self.colour_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.colour_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.colour_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_dust(self):
        return self.remote.is_directory(self.dust_maps_path) and not self.remote.is_empty(self.dust_maps_path)

    # -----------------------------------------------------------------

    @property
    def ndust_maps(self):
        if self.remote.has_files_in_path(self.dust_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.dust_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.dust_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_ionizing(self):
        return self.remote.is_directory(self.ionizing_maps_path) and not self.remote.is_empty(self.ionizing_maps_path)

    # -----------------------------------------------------------------

    @property
    def nionizing_maps(self):
        if self.remote.has_files_in_path(self.ionizing_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.ionizing_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.ionizing_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_rgb(self):
        return self.remote.is_directory(self.rgb_maps_path) and not self.remote.is_empty(self.rgb_maps_path)

    # -----------------------------------------------------------------

    @property
    def has_maps_old(self):
        return self.remote.is_directory(self.old_maps_path) and not self.remote.is_empty(self.old_maps_path)

    # -----------------------------------------------------------------

    @property
    def nold_maps(self):
        if self.remote.has_files_in_path(self.old_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.old_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.old_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_ssfr(self):
        return self.remote.is_directory(self.ssfr_maps_path) and not self.remote.is_empty(self.ssfr_maps_path)

    # -----------------------------------------------------------------

    @property
    def nssfr_maps(self):
        if self.remote.has_files_in_path(self.ssfr_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.ssfr_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.ssfr_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_tir(self):
        return self.remote.is_directory(self.tir_maps_path) and not self.remote.is_empty(self.tir_maps_path)

    # -----------------------------------------------------------------

    @property
    def ntir_maps(self):
        if self.remote.has_files_in_path(self.tir_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.tir_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.tir_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_maps_young(self):
        return self.remote.is_directory(self.young_maps_path) and not self.remote.is_empty(self.young_maps_path)

    # -----------------------------------------------------------------

    @property
    def nyoung_maps(self):
        if self.remote.has_files_in_path(self.young_maps_path, extension="fits"): return self.remote.nfiles_in_path(self.young_maps_path, extension="fits")
        else: return self.remote.nfiles_in_path(self.young_maps_path, extension="fits", recursive=True, recursion_level=1)

    # -----------------------------------------------------------------

    @property
    def has_heating(self):
        return self.remote.is_file(self.heating_config_path)

    # -----------------------------------------------------------------

    @property
    def ski_file(self):
        return SkiFile.from_remote_file(self.ski_file_path, self.remote)

    # -----------------------------------------------------------------

    @property
    def has_dust_grid_simulation_logfile(self):
        return self.remote.is_file(self.dust_grid_simulation_logfile_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_simulation_logfile(self):
        return LogFile.from_remote_file(self.dust_grid_simulation_logfile_path, self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid_tree(self):

        """
        This function ...
        :return:
        """

        # Give debug message
        log.debug("Loading the dust grid tree, this may take a while (depending on the number of nodes) ...")

        # Return the tree
        return DustGridTree.from_remote_file(self.dust_grid_tree_path, self.remote)

    # -----------------------------------------------------------------

    @property
    def remote_script_paths(self):
        return self.remote.files_in_path(self.path, extension="sh")

    # -----------------------------------------------------------------

    def get_remote_script_commands(self):

        """
        This fucntion ...
        :return:
        """

        commands = dict()

        # Loop over the script paths
        for path in self.remote_script_paths:

            # Get host ID
            host_id = fs.strip_extension(fs.name(path))

            lines = []

            # Open the file
            for line in self.remote.read_lines(path):

                if line.startswith("#"): continue
                if not line.strip(): continue

                lines.append(line)

            # Set the commands
            commands[host_id] = lines

        # Return the commands
        return commands

    # -----------------------------------------------------------------

    def get_remote_script_commands_for_host(self, host_id):

        """
        This function ...
        :param host_id:
        :return:
        """

        commands = self.get_remote_script_commands()
        if host_id in commands: return commands[host_id]
        else: return []

    # -----------------------------------------------------------------

    def get_heating_ski_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return SkiFile.from_remote_file(self.ski_path_for_contribution(contribution), self.remote)

    # -----------------------------------------------------------------

    @property
    def has_simulated_sed(self):
        return self.remote.is_file(self.simulated_sed_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_sed(self):
        return SED.from_skirt(self.simulated_sed_path, remote=self.remote)

    # -----------------------------------------------------------------

    @property
    def has_simulated_fluxes(self):
        return self.remote.is_file(self.simulated_fluxes_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def simulated_fluxes(self):
        return ObservedSED.from_remote_file(self.simulated_fluxes_path, self.remote)

# -----------------------------------------------------------------

class CachedAnalysisRuns(AnalysisRunBase):

    """
    This class ...
    """

    def __init__(self, modeling_path, remote):

        """
        This function ...
        :param modeling_path:
        :param remote:
        """

        # Set attributes
        self.modeling_path = modeling_path
        self.remote = load_remote(remote, silent=True)

    # -----------------------------------------------------------------

    @property
    def galaxy_name(self):
        return fs.name(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_name(self):
        return self.galaxy_name + "_analysis"

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.remote.home_directory, self.cache_directory_name)
        if not self.remote.is_directory(path): self.remote.create_directory(path)
        return path

    # -----------------------------------------------------------------

    def get_cache_directory_path_run(self, run_name):

        """
        This function ...
        :param run_name:
        :return:
        """

        path = fs.join(self.cache_directory_path, run_name)
        #if not self.remote.is_directory(path): self.remote.create_directory(path)
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):
        return self.remote.directories_in_path(self.cache_directory_path, returns="name")

    # -----------------------------------------------------------------

    def __len__(self):

        """
        This function ...
        :return:
        """

        return len(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def empty(self):
        return sequences.is_empty(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single(self):
        return sequences.is_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_name(self):
        return sequences.get_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_path(self):
        return self.get_path(self.single_name)

    # -----------------------------------------------------------------

    def get_path(self, name):

        """
        Thisn function ...
        :param name:
        :return:
        """

        return self.get_cache_directory_path_run(name)

    # -----------------------------------------------------------------

    def load(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        analysis_run_path = self.get_path(name)
        if not self.remote.is_directory(analysis_run_path): raise ValueError("Analysis run '" + name + "' does not exist")
        return CachedAnalysisRun.from_path(analysis_run_path, self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def single(self):
        return CachedAnalysisRun.from_path(self.single_path, self.remote)

    # -----------------------------------------------------------------

    # ABSTRACT PROPERTIES FROM BASE CLASS?

    @property
    def colour_names(self):
        return None

    # -----------------------------------------------------------------

    @property
    def dust_grid_simulation_logfile(self):
        return None

    # -----------------------------------------------------------------

    @property
    def dust_grid_tree(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_dust_grid_simulation_logfile(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_heating(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_attenuation(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_colours(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_dust(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_ionizing(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_old(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_ssfr(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_tir(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_maps_young(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_simulated_fluxes(self):
        return None

    # -----------------------------------------------------------------

    @property
    def has_simulated_sed(self):
        return None

    # -----------------------------------------------------------------

    @property
    def nattenuation_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def ncolour_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def ndust_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def nionizing_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def nold_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def nssfr_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def ntir_maps(self):
        return None

    # -----------------------------------------------------------------

    @property
    def nyoung_maps(self):
        return None

# -----------------------------------------------------------------

def find_value_for_unique_key_nested(dictionary, key, allow_none=False):

    """
    This function ...
    :param dictionary:
    :param key:
    :param allow_none:
    :return:
    """

    values = []

    for key_i in dictionary:

        # Sub-dict
        if isinstance(dictionary[key_i], dict):
            value = find_value_for_unique_key_nested(dictionary[key_i], key, allow_none=True)
            if value is not None: values.append(value)

        # Matches
        elif key_i == key:
            value = dictionary[key_i]
            values.append(value)

    if len(values) == 0 and not allow_none: raise ValueError("Key not found")
    if len(values) > 1: raise ValueError("Not unique")

    # Return the only value
    if len(values) == 0: return None
    else: return values[0]

# -----------------------------------------------------------------

def find_keys_for_unique_key_nested(dictionary, key, allow_none=False):

    """
    This function ...
    :param dictionary:
    :param key:
    :param allow_none:
    :return:
    """

    keys = []

    for key_i in dictionary:

        # Sub-dict
        if isinstance(dictionary[key_i], dict):

            keys_subdict = find_keys_for_unique_key_nested(dictionary[key_i], key, allow_none=True)
            if keys_subdict is not None:
                keys.append([key_i] + keys_subdict)

        # Matches
        elif key_i == key: keys.append([key])

    if len(keys) == 0 and not allow_none: raise ValueError("Key not found")
    if len(keys) > 1: raise ValueError("Not unique")

    if len(keys) == 0: return None
    else: return keys[0]

# -----------------------------------------------------------------
