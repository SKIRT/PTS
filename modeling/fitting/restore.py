#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.restore Restore the backup of a fit to the fitting run directory.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import OrderedDict

# Import the relevant PTS classes and modules
from pts.modeling.fitting.backup import FitBackupper
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from .component import FittingComponent
from ...core.tools.utils import lazyproperty

# -----------------------------------------------------------------

class FitRestorer(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FitRestorer, self).__init__(*args, **kwargs)

        # The fitting run
        self.fitting_run = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Make backup
        if self.config.backup: self.backup()

        # 3. Restore
        self.restore()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(FitRestorer, self).setup(**kwargs)

        # Load the fitting run
        if kwargs.get("fitting_run", None) is not None: self.fitting_run = kwargs.pop("fitting_run")
        else: self.fitting_run = self.load_fitting_run(self.config.run)

    # -----------------------------------------------------------------

    def backup(self):

        """
        This function ...
        :return:
        """

        # Create the backupper
        backupper = FitBackupper()

        # Set the name for the backup
        backupper.config.name = self.config.backup_name

        # Set modeling path
        backupper.config.path = self.config.path

        # Run the backupper
        backupper.run(fitting_run=self.fitting_run)

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_path(self):

        """
        This function ...
        :return:
        """

        path = fs.join(self.fitting_run.refitting_path, self.config.name)
        if not fs.is_directory(path): raise ValueError("'" + self.config.name + "' is not a backup of a fit")
        return path

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_generations_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "generations")

    # -----------------------------------------------------------------

    def get_generation_restore_path(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.join(self.restore_generations_path, generation_name)

    # -----------------------------------------------------------------

    def has_generation(self, generation_name):

        """
        This function ...
        :param generation_name:
        :return:
        """

        return fs.is_directory(self.get_generation_restore_path(generation_name))

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_prob_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "prob")

    # -----------------------------------------------------------------

    @property
    def has_prob(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.restore_prob_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_best_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "best")

    # -----------------------------------------------------------------

    @property
    def has_best(self):

        """
        This function ...
        :return:
        """

        return fs.is_directory(self.restore_best_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_best_parameters_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "best_parameters.dat")

    # -----------------------------------------------------------------

    @property
    def has_best_parameters(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.restore_best_parameters_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_configuration_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "configuration.cfg")

    # -----------------------------------------------------------------

    @property
    def has_configuration(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.restore_configuration_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def restore_weights_path(self):

        """
        This function ...
        :return:
        """

        return fs.join(self.restore_path, "weights.dat")

    # -----------------------------------------------------------------

    @property
    def has_weights(self):

        """
        This function ...
        :return:
        """

        return fs.is_file(self.restore_weights_path)

    # -----------------------------------------------------------------

    def restore(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Restoring the previous fit ...")

        # Configuration
        if self.has_configuration: self.restore_config()

        # Weights
        if self.has_weights: self.restore_weights()

        # Best parameters
        if self.has_best_parameters: self.restore_best_parameters()

        # Prob
        if self.has_prob: self.restore_prob()

        # Best
        if self.has_best: self.restore_best()

        # Fluxes
        self.restore_fluxes()

        # Fluxes plots
        self.restore_fluxes_plots()

        # Differences
        self.restore_differences()

        # Chi squared
        self.restore_chi_squared()

    # -----------------------------------------------------------------

    def restore_config(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Restoring the fitting configuration ...")

        # Copy the configuration
        fs.copy_file(self.restore_configuration_path, self.fitting_run.path)

    # -----------------------------------------------------------------

    def restore_weights(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Restoring the weights ...")

        # Copy the weights
        fs.copy_file(self.restore_weights_path, self.fitting_run.path)

    # -----------------------------------------------------------------

    def restore_best_parameters(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Restoring the best parameters table ...")

        # Copy the file
        fs.copy_file(self.restore_best_parameters_path, self.fitting_run.path)

    # -----------------------------------------------------------------

    def restore_prob(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Restoring the probabilities ...")

        # Copy the directory
        fs.clear_directory(self.fitting_run.prob_path, recursive=True)
        fs.copy_directory(self.restore_prob_path, self.fitting_run.path)

    # -----------------------------------------------------------------

    def restore_best(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Restoring the best simulations ...")

        # Copy the directory
        fs.clear_directory(self.fitting_run.best_path, recursive=True)
        fs.copy_directory(self.restore_best_path, self.fitting_run.path)

    # -----------------------------------------------------------------

    @property
    def generation_names(self):

        """
        This function ...
        :return:
        """

        return self.fitting_run.generation_names

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

    def restore_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the mock fluxes ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Has generation?
            if not self.has_generation(generation_name): continue

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get restore path
            generation_path = self.get_generation_restore_path(generation_name)

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Set simulation path
                simulation_path = fs.join(generation_path, simulation_name)

                # Check whether fluxes are present
                filepath = fs.join(simulation_path, "earth_fluxes.dat")
                if not fs.is_file(filepath):
                    log.warning("No mock SED for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Restoring fluxes for simulation '" + simulation_name + "' ...")

                # Copy the fluxes file
                fs.copy_file(filepath, generation.get_mock_sed_path(simulation_name))

    # -----------------------------------------------------------------

    def restore_fluxes_plots(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the mock fluxes plots ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Has generation?
            if not self.has_generation(generation_name): continue

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get restore path
            generation_path = self.get_generation_restore_path(generation_name)

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Set simulation path
                simulation_path = fs.join(generation_path, simulation_name)

                # Check whether fluxes are present
                filepath = fs.join(simulation_path, "earth_fluxes.pdf")
                if not fs.is_file(filepath):
                    log.warning("No mock SED plot for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Restoring fluxes plot for simulation '" + simulation_name + "' ...")

                # Copy the plot file
                fs.copy_file(filepath, generation.get_mock_sed_plot_path(simulation_name))

    # -----------------------------------------------------------------

    def restore_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making a backup of the flux differences ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Has generation?
            if not self.has_generation(generation_name): continue

            # Debugging
            log.debug("Creating backups for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get restore path
            generation_path = self.get_generation_restore_path(generation_name)

            # Loop over the simulations
            for simulation_name in generation.simulation_names:

                # Set simulation path
                simulation_path = fs.join(generation_path, simulation_name)

                # Check whether the simulation has differences
                filepath = fs.join(simulation_path, "differences.dat")
                if not fs.is_file(filepath):
                    log.warning("No differences table for simulation '" + simulation_name + "' of generation '" + generation_name + "'")
                    continue

                # Debugging
                log.debug("Restoring differences table for simulation '" + simulation_name + "' ...")

                # Copy the file
                fs.copy_file(filepath, generation.get_simulation_sed_differences_path(simulation_name))

    # -----------------------------------------------------------------

    def restore_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Restoring the chi squared tables ...")

        # Loop over the generations
        for generation_name in self.generation_names:

            # Has generation?
            if not self.has_generation(generation_name): continue

            # Debugging
            log.debug("Creating backup of chi squared table for generation '" + generation_name + "' ...")

            # Get the generation
            generation = self.generations[generation_name]

            # Get the generation path
            generation_path = self.get_generation_restore_path(generation_name)

            # Determine the filepath
            filepath = fs.join(generation_path, "chi_squared.dat")

            # Copy the file
            fs.copy_file(filepath, generation.chi_squared_table_path)

# -----------------------------------------------------------------
