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

# Import the relevant PTS classes and modules
from ...core.tools import filesystem as fs
from ...core.simulation.skifile import LabeledSkiFile
from ..core.model import Model
from ...core.tools import sequences
from ...core.basics.composite import SimplePropertyComposite
from ..fitting.run import FittingRun
from pts.core.tools.utils import lazyproperty
from ...core.tools.serialization import load_dict

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
        self.add_string_property("simulation name", "simulation name")
        self.add_string_property("model_name", "name of the model")
        self.add_property("parameter_values", "integer_real_and_quantity_list", "parameter values")
        self.add_property("chi_squared", "real", "chi squared value")

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class AnalysisRun(object):

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
        info_path = fs.join(path, "info.dat")
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

        # Create the instance
        run = cls()

        # Set the analysis run info
        run.info = AnalysisRunInfo.from_file(info_path)

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
    def from_fitting(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run is not None

    # -----------------------------------------------------------------

    @property
    def from_model(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run is None

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
    def path(self):

        """
        This function ...
        :return:
        """

        return self.info.path

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

        return fs.join(self.path, "input.dat")

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
        return fs.join(self.path, "wavelength_grid.dat")

    # -----------------------------------------------------------------

    @property
    def dust_grid_path(self):

        """
        This function ...
        :return:
        """

        # Set the path to the dust grid file
        return fs.join(self.path, "dust_grid.dg")

    # -----------------------------------------------------------------

    @property
    def info_path(self):

        """
        This function ...
        :return:
        """

        # Set the path to the analysis run info file
        return fs.join(self.path, "info.dat")

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

        return fs.join(self.path, "extr")

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

        return fs.join(self.path, "plot")

    # -----------------------------------------------------------------

    @property
    def misc_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "misc")

    # -----------------------------------------------------------------

    @property
    def attenuation_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "attenuation")

    # -----------------------------------------------------------------

    @property
    def colours_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "colours")

    # -----------------------------------------------------------------

    @property
    def residuals_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "residuals")

    # -----------------------------------------------------------------

    @property
    def heating_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "heating")

    # -----------------------------------------------------------------

    @property
    def heating_wavelength_grid_path(self):

        """
        This fucntion ...
        :return:
        """

        return fs.join(self.heating_path, "wavelength_grid.dat")

    # -----------------------------------------------------------------

    @property
    def heating_instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.heating_path, "instruments")

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
    def dust_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "dust_grid.dg")

    # -----------------------------------------------------------------

    @property
    def wavelength_grid_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "wavelength_grid.dat")

    # -----------------------------------------------------------------

    @property
    def instruments_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.path, "instruments")

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
        model.parameter_values = self.parameter_values # set the parameter values

        # Return the model
        return model

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
