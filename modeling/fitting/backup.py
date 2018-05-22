#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.backup Contains the FitBackupper class

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import filesystem as fs

# -----------------------------------------------------------------

class FitBackupper(FittingComponent):

    """
    This function ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(FitBackupper, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

    # -----------------------------------------------------------------

    @lazyproperty
    def has_prob(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.fitting_run.prob_path, recursive=True)

    # -----------------------------------------------------------------

    @lazyproperty
    def has_best(self):

        """
        This function ...
        :return:
        """

        return fs.has_files_in_path(self.fitting_run.best_path, recursive=True)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Backup the fitting configuration
        self.backup_config()

        # 3. Backup the weights
        self.backup_weights()

        # 4. Backup the best parameters table
        self.backup_best_parameters()

        # 5. Backup the prob directory
        if self.has_prob: self.backup_prob()

        # 6. Backup the best directory
        if self.has_best: self.backup_best()

        # 7. Backup the fluxes
        self.backup_fluxes()

        # 8. Backup the fluxes plots
        self.backup_fluxes_plots()

        # 9. Backup the differences
        self.backup_differences()

        # 10. Backup the chi squared tables
        self.backup_chi_squared()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(FitBackupper, self).setup(**kwargs)

        # Load the fitting run
        if kwargs.get("fitting_run", None) is not None: self.fitting_run = kwargs.pop("fitting_run")
        else: self.fitting_run = self.load_fitting_run(self.config.fitting_run)

    # -----------------------------------------------------------------

    @property
    def refitting_path(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.refitting_path

    # -----------------------------------------------------------------

    @property
    def backup_name(self):

        """
        This function ...
        :return:
        """

        return self.config.name

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.refitting_path, self.backup_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_names(self):

        """
        This function ...
        :return:
        """

        if self.config.generations is not None: return self.config.generations
        else: return self.fitting_run.generation_names

    # -----------------------------------------------------------------

    @lazyproperty
    def generations(self):

        """
        This function ...
        :return:
        """

        gens = OrderedDict()
        for name in self.generation_names: gens[name] = self.fitting_run.get_generation(name)
        return gens

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.refitting_path, self.backup_name)

    # -----------------------------------------------------------------

    @lazyproperty
    def backup_generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.backup_path, "generations")

    # -----------------------------------------------------------------

    @memoize_method
    def backup_path_for_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.create_directory_in(self.backup_generations_path, generation_name)

    # -----------------------------------------------------------------

    @memoize_method
    def backup_path_for_simulation(self, generation_name, simulation_name):

        """
        This function ...
        :param generation_name:
        :param simulation_name:
        :return:
        """

        return fs.create_directory_in(self.backup_path_for_generation(generation_name), simulation_name)

    # -----------------------------------------------------------------

    def backup_config(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Making a backup of the fitting configuration ...")

        # Copy the configuration
        fs.copy_file(self.fitting_run.fitting_configuration_path, self.backup_path)

    # -----------------------------------------------------------------

    def backup_weights(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Making a backup of the weights ...")

        # Copy the weights
        fs.copy_file(self.fitting_run.weights_table_path, self.backup_path)

    # -----------------------------------------------------------------

    def backup_best_parameters(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Making a backup of the best parameters table ...")

        # Copy the file
        fs.copy_file(self.fitting_run.best_parameters_table_path, self.backup_path)

    # -----------------------------------------------------------------

    def backup_prob(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the model and parameter probabilities ...")

        # Copy the directory
        fs.copy_directory(self.fitting_run.prob_path, self.backup_path)

    # -----------------------------------------------------------------

    def backup_best(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the best simulations ...")

        # Copy the directory
        fs.copy_directory(self.fitting_run.best_path, self.backup_path)

    # -----------------------------------------------------------------

    def backup_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the mock fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Check whether the simulation has fluxes
                if not generation.has_mock_sed(simulation_name):
                    log.warning("No mock SED for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Creating backup for simulation '" + simulation_name + "' ...")

                # Copy the fluxes file
                fs.copy_file(generation.get_mock_sed_path(simulation_name), self.backup_path_for_simulation(generation_name, simulation_name))

    # -----------------------------------------------------------------

    def backup_fluxes_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the mock fluxes plots ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Check whether the simulation has fluxes plot
                if not generation.has_mock_sed_plot(simulation_name):
                    log.warning("No mock SED plot for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Creating backup for simulation '" + simulation_name + "' ...")

                # Copy the plot file
                fs.copy_file(generation.get_mock_sed_plot_path(simulation_name), self.backup_path_for_simulation(generation_name, simulation_name))

    # -----------------------------------------------------------------

    def backup_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the flux differences ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Check whether the simulation has differences
                if not generation.has_sed_differences(simulation_name):
                    log.warning("No differences table for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Creating backup for simulation '" + simulation_name + "' ...")

                # Copy the file
                fs.copy_file(generation.get_simulation_sed_differences_path(simulation_name), self.backup_path_for_simulation(generation_name, simulation_name))

    # -----------------------------------------------------------------

    def backup_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the chi squared tables ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Debugging
            log.debug("Creating backup for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Copy the file
            fs.copy_file(generation.chi_squared_table_path, self.backup_path_for_generation(generation_name))

# -----------------------------------------------------------------
