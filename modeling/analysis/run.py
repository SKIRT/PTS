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
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import LabeledSkiFile
from ..core.model import Model
from ...core.tools import sequences
from ...core.basics.composite import SimplePropertyComposite
from ..fitting.run import FittingRun
from pts.core.tools.utils import lazyproperty
from ...core.tools.serialization import load_dict
from ...core.simulation.tree import DustGridTree
from ...core.simulation.grids import FileTreeDustGrid, load_grid
from ...core.simulation.wavelengthgrid import WavelengthGrid
from ..basics.projection import GalaxyProjection, EdgeOnProjection, FaceOnProjection
from ..basics.instruments import FullInstrument
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...core.basics.log import log
from ...core.remote.remote import load_remote
from ...core.basics.configuration import Configuration

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

dust_grid_filename = "dust_grid.dg"
wavelength_grid_filename = "wavelength_grid.dat"
dust_grid_build_name = "dust grid"
info_filename = "info.dat"
config_filename = "config.cfg"
launch_config_filename = "launch_config.cfg"
input_filename = "input.dat"
instruments_name = "instruments"
projections_name = "projections"
extract_name = "extr"
plot_name = "plot"
misc_name = "misc"
attenuation_name = "attenuation"
colours_name = "colours"
residuals_name = "residuals"
weighed_residuals_name = "weighed_residuals"
heating_name = "heating"
dust_grid_tree_filename = "tree.dat"

# Projections
earth_projection_filename = "earth.proj"
faceon_projection_filename = "faceon.proj"
edgeon_projection_filename = "edgeon.proj"

# Instruments
earth_instrument_filename = "earth.instr"
faceon_instrument_filename = "faceon.instr"
edgeon_instrument_filename = "edgeon.instr"

# -----------------------------------------------------------------

class AnalysisRunBase(object):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    @property
    def from_fitting(self):

        """
        This function ...
        :return:
        """

        #return self.fitting_run is not None
        return self.fitting_run_name is not None

    # -----------------------------------------------------------------

    @property
    def from_model(self):

        """
        This function ...
        :return:
        """

        #return self.fitting_run is None
        return self.fitting_run_name is None

    # -----------------------------------------------------------------

    @property
    def from_generation(self):

        """
        This function ...
        :return:
        """

        # Otherwise: from initial guess

        return self.from_fitting and self.generation_name is not None

    # -----------------------------------------------------------------

    @property
    def from_initial_guess(self):

        """
        This function ...
        :return:
        """

        # Otherwise: from best simulation of a certain generation

        return self.from_fitting and self.generation_name is None

    # -----------------------------------------------------------------

    @property
    def name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

    # -----------------------------------------------------------------

    @property
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.info.generation_name

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):

        """
        This function ...
        :return:
        """

        return self.info.simulation_name

    # -----------------------------------------------------------------

    @property
    def model_name(self):

        """
        This function ...
        :return:
        """

        return self.info.model_name

    # -----------------------------------------------------------------

    @property
    def input_file_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, input_filename)

    # -----------------------------------------------------------------

    @property
    def ski_file_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        # Set the path to the wavelength grid file
        return fs.join(self.path, wavelength_grid_filename)

    # -----------------------------------------------------------------

    @property
    def nwavelengths(self):

        """
        This function ...
        :return:
        """

        return len(self.wavelength_grid)

    # -----------------------------------------------------------------

    @property
    def dust_grid_path(self):

        """
        This function ...
        :return:
        """

        # Set the path to the dust grid file
        return fs.join(self.path, dust_grid_filename)

    # -----------------------------------------------------------------

    @property
    def info_path(self):

        """
        This function ...
        :return:
        """

        # Set the path to the analysis run info file
        return fs.join(self.path, info_filename)

    # -----------------------------------------------------------------

    @property
    def config_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, config_filename)

    # -----------------------------------------------------------------

    @property
    def launch_config_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, launch_config_filename)

    # -----------------------------------------------------------------

    @property
    def out_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "out")

    # -----------------------------------------------------------------

    @property
    def output_path(self):

        """
        This function ...
        :return:
        """

        return self.out_path

    # -----------------------------------------------------------------

    @property
    def extr_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, extract_name)

    # -----------------------------------------------------------------

    @property
    def extract_path(self):

        """
        This function ...
        :return:
        """

        return self.extr_path

    # -----------------------------------------------------------------

    @property
    def plot_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, plot_name)

    # -----------------------------------------------------------------

    @property
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, misc_name)

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, attenuation_name)

    # -----------------------------------------------------------------

    @property
    def colours_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, colours_name)

    # -----------------------------------------------------------------

    @property
    def residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, residuals_name)

    # -----------------------------------------------------------------

    @property
    def weighed_residuals_path(self):

        """
        This function ...
        :return:
        """

        return self.join(self.path, weighed_residuals_name)

    # -----------------------------------------------------------------

    @property
    def heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, heating_name)

    # -----------------------------------------------------------------

    @property
    def heating_wavelength_grid_path(self):

        """
        This fucntion ...
        :return:
        """

        return fs.join(self.heating_path, wavelength_grid_filename)

    # -----------------------------------------------------------------

    @property
    def heating_instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.heating_path, instruments_name)

    # -----------------------------------------------------------------

    def heating_simulation_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_path, contribution)

    # -----------------------------------------------------------------

    def heating_ski_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_simulation_path_for_contribution(contribution), self.galaxy_name + ".ski")

    # -----------------------------------------------------------------

    def heating_output_path_for_contribution(self, contribution):

        """
        This function ...
        :param contribution:
        :return:
        """

        return fs.join(self.heating_simulation_path_for_contribution(contribution), "out")

    # -----------------------------------------------------------------

    @property
    def analysis_run_name(self):

        """
        This function ...
        :return:
        """

        return self.info.name

# -----------------------------------------------------------------

class AnalysisRun(AnalysisRunBase):

    """
    This class ...
    """

    def __init__(self, galaxy_name=None, info=None):

        """
        The constructor ...
        :param galaxy_name:
        :param info:
        """

        # Set attributes
        self.galaxy_name = galaxy_name
        self.info = info

        ## Create directories

        # The directory for the projections and the instruments
        if not fs.is_directory(self.projections_path): fs.create_directory(self.projections_path)
        if not fs.is_directory(self.instruments_path): fs.create_directory(self.instruments_path)

        # The directory for the dust grid output
        if not fs.is_directory(self.dust_grid_build_path): fs.create_directory(self.dust_grid_build_path)
        if not fs.is_directory(self.dust_grid_simulation_out_path): fs.create_directory(self.dust_grid_simulation_out_path)

        # Simulation directories
        if not fs.is_directory(self.output_path): fs.create_directory(self.output_path)
        if not fs.is_directory(self.extract_path): fs.create_directory(self.extract_path)
        if not fs.is_directory(self.plot_path): fs.create_directory(self.plot_path)
        if not fs.is_directory(self.misc_path): fs.create_directory(self.misc_path)

        # Analysis directories
        if not fs.is_directory(self.attenuation_path): fs.create_directory(self.attenuation_path)
        if not fs.is_directory(self.colours_path): fs.create_directory(self.colours_path)
        if not fs.is_directory(self.residuals_path): fs.create_directory(self.residuals_path)
        if not fs.is_directory(self.weighed_residuals_path): fs.create_directory(self.weighed_residuals_path)
        if not fs.is_directory(self.heating_path): fs.create_directory(self.heating_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_name(cls, modeling_path, name):

        """
        This function ...
        :param modeling_path:
        :param name:
        :return:
        """

        analysis_path = fs.join(modeling_path, "analysis")
        run_path = fs.join(analysis_path, name)
        return cls.from_path(run_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_path(cls, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Determine the info path
        info_path = fs.join(path, info_filename)
        if not fs.is_file(info_path): raise IOError("Could not find the info file")
        else: return cls.from_info(info_path)

    # -----------------------------------------------------------------

    @classmethod
    def from_info(cls, info_path):

        """
        This function ...
        :param info_path:
        :return:
        """

        # Load the analysis run info
        info = AnalysisRunInfo.from_file(info_path)

        # Create the instance
        run = cls(info=info)

        # Set galaxy name
        modeling_path = fs.directory_of(fs.directory_of(run.info.path))
        run.galaxy_name = fs.name(modeling_path)

        # Return the analysis run object
        return run

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.path)

    # -----------------------------------------------------------------

    @property
    def modeling_path(self):

        """
        This function ...
        :return:
        """

        return fs.directory_of(self.analysis_path)

    # -----------------------------------------------------------------

    @property
    def path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @lazyproperty
    def config(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_file(self.config_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):

        """
        This function ...
        :return:
        """

        return load_grid(self.dust_grid_path)

    # -----------------------------------------------------------------

    @property
    def analysis_run_path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @property
    def ski_file(self):

        """
        This function ...
        :return:
        """

        return LabeledSkiFile(self.ski_file_path)

    # -----------------------------------------------------------------

    @property
    def input_paths(self):

        """
        This function ...
        :return:
        """

        return load_dict(self.input_file_path)

    # -----------------------------------------------------------------

    @property
    def dust_grid_build_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, dust_grid_build_name)

    # -----------------------------------------------------------------

    @property
    def dust_grid_simulation_out_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.dust_grid_build_path, "out")

    # -----------------------------------------------------------------

    @property
    def dust_grid_tree_path(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return fs.is_file(self.dust_grid_tree_path)

    # -----------------------------------------------------------------

    @property
    def instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, instruments_name)

    # -----------------------------------------------------------------

    @property
    def projections_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, projections_name)

    # -----------------------------------------------------------------

    @property
    def earth_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, earth_projection_filename)

    # -----------------------------------------------------------------

    @property
    def faceon_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, faceon_projection_filename)

    # -----------------------------------------------------------------

    @property
    def edgeon_projection_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.projections_path, edgeon_projection_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_projection(self):

        """
        This function ...
        :return:
        """

        return GalaxyProjection.from_file(self.earth_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_projection(self):

        """
        This function ...
        :return:
        """

        return EdgeOnProjection.from_file(self.edgeon_projection_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_projection(self):

        """
        This function ...
        :return:
        """

        return FaceOnProjection.from_file(self.faceon_projection_path)

    # -----------------------------------------------------------------

    @property
    def earth_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, earth_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def faceon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, faceon_instrument_filename)

    # -----------------------------------------------------------------

    @property
    def edgeon_instrument_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.instruments_path, edgeon_instrument_filename)

    # -----------------------------------------------------------------

    @lazyproperty
    def earth_instrument(self):

        """
        This function ...
        :return:
        """

        return FullInstrument.from_file(self.earth_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def faceon_instrument(self):

        """
        This function ...
        :return:
        """

        return FullInstrument.from_file(self.faceon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def edgeon_instrument(self):

        """
        This function ...
        :return:
        """

        return FullInstrument.from_file(self.edgeon_instrument_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run_name(self):

        """
        This function ...
        :return:
        """

        return self.info.fitting_run

    # -----------------------------------------------------------------

    @lazyproperty
    def fitting_run(self):

        """
        This function ...
        :return:
        """

        return FittingRun.from_name(self.modeling_path, self.fitting_run_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def model_suite(self):

        """
        This function ...
        :return:
        """

        from ..build.suite import ModelSuite
        return ModelSuite.from_modeling_path(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_name(self):

        """
        This function ...
        :return:
        """

        return self.info.generation_name

    # -----------------------------------------------------------------

    @lazyproperty
    def simulation_name(self):

        """
        This function ...
        :return:
        """

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

    @lazyproperty
    def chi_squared(self):

        """
        This function ...
        :return:
        """

        return self.info.chi_squared

    # -----------------------------------------------------------------

    @lazyproperty
    def model(self):

        """
        This function ...
        :return:
        """

        # Create the model
        model = Model()

        # Set attributes
        model.simulation_name = self.simulation_name
        model.chi_squared = self.chi_squared
        model.parameter_values = self.parameter_values

        # Return the model
        return model

    # -----------------------------------------------------------------

    @property
    def uses_grid_resolution(self):

        """
        Thisf unction ...
        :return:
        """

        return self.info.reference_deprojection == "grid"

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection_component_name(self):

        """
        Thisf unction ...
        :return:
        """

        if self.uses_grid_resolution: return None
        return self.info.reference_deprojection

    # -----------------------------------------------------------------

    @lazyproperty
    def is_stellar_reference_deprojection(self):

        """
        This function ...
        :return:
        """

        if self.uses_grid_resolution: raise ValueError("This function shouldn't be called")
        return self.reference_deprojection_component_name in self.model_suite.get_stellar_component_names(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def is_dust_reference_deprojection(self):

        """
        This function ...
        :return:
        """

        if self.uses_grid_resolution: raise ValueError("This function shouldn't be called")
        return self.reference_deprojection_component_name in self.model_suite.get_dust_component_names(self.model_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection_component(self):

        """
        This function ...
        :return:
        """

        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component(self.model_name, self.reference_deprojection_component_name, add_map=False)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component(self.model_name, self.reference_deprojection_component_name, add_map=False)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_deprojection(self):

        """
        Thisf unction ...
        :return:
        """

        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component_deprojection(self.model_name, self.reference_deprojection_component_name, add_map=False)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component_deprojection(self.model_name, self.reference_deprojection_component_name, add_map=False)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map(self):

        """
        This function ...
        :return:
        """

        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.load_stellar_component_map(self.model_name, self.reference_deprojection_component_name)
            elif self.is_dust_reference_deprojection: return self.model_suite.load_dust_component_map(self.model_name, self.reference_deprojection_component_name)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_map_path(self):

        """
        This function ...
        :return:
        """

        if self.reference_deprojection_component_name is None: return None
        else:
            if self.is_stellar_reference_deprojection: return self.model_suite.get_stellar_component_map_path(self.model_name, self.reference_deprojection_component_name)
            elif self.is_dust_reference_deprojection: return self.model_suite.get_dust_component_map_path(self.model_name, self.reference_deprojection_component_name)
            else: raise ValueError("Reference deprojection component name '" + self.reference_deprojection_component_name + "' not recognized as either stellar or dust")

    # -----------------------------------------------------------------

    @lazyproperty
    def reference_wcs(self):

        """
        This function ...
        :return:
        """

        if self.reference_map_path is None: return None
        else: return CoordinateSystem.from_file(self.reference_map_path)

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

        self.modeling_path = modeling_path

    # -----------------------------------------------------------------

    @property
    def analysis_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.modeling_path, "analysis")

    # -----------------------------------------------------------------

    @lazyproperty
    def names(self):

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.analysis_path, returns="name")

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

        """
        This function ...
        :return:
        """

        return sequences.is_empty(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single(self):

        """
        This function ...
        :return:
        """

        return sequences.is_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_name(self):

        """
        This function ...
        :return:
        """

        return sequences.get_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_path(self):

        """
        This function ...
        :return:
        """

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
        return AnalysisRun.from_path(analysis_run_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def single(self):

        """
        This function ...
        :return:
        """

        return AnalysisRun.from_path(self.single_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def last_name(self):

        """
        This function ...
        :return:
        """

        if self.empty: return None
        if self.has_single: return self.single_name

        return sorted(self.names)[-1]

    # -----------------------------------------------------------------

    @lazyproperty
    def last_path(self):

        """
        This function ...
        :return:
        """

        return self.get_path(self.last_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def last(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return AnalysisRunInfo.from_remote_file(self.info_path, self.remote)

    # -----------------------------------------------------------------

    @property
    def original_path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

    # -----------------------------------------------------------------

    @lazyproperty
    def config(self):

        """
        This function ...
        :return:
        """

        return Configuration.from_remote_file(self.config_path, self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_skirt_input(self.wavelength_grid_path, remote=self.remote)

    # -----------------------------------------------------------------

    @lazyproperty
    def dust_grid(self):

        """
        This function ...
        :return:
        """

        return load_grid(self.dust_grid_path, remote=self.remote)

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

        """
        This function ...
        :return:
        """

        return fs.name(self.modeling_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def cache_directory_name(self):

        """
        Thisf unction ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return fs.directories_in_path(self.cache_directory_path, returns="name")

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

        """
        This function ...
        :return:
        """

        return sequences.is_empty(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_single(self):

        """
        This function ...
        :return:
        """

        return sequences.is_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_name(self):

        """
        This function ...
        :return:
        """

        return sequences.get_singleton(self.names)

    # -----------------------------------------------------------------

    @lazyproperty
    def single_path(self):

        """
        This function ...
        :return:
        """

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

        """
        This function ...
        :return:
        """

        return CachedAnalysisRun.from_path(self.single_path, self.remote)

# -----------------------------------------------------------------
