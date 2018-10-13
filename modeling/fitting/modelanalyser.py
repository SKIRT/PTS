#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.fitting.modelanalyser Contains the FitModelAnalyser class, used for analysing the goodness
#  of the radiative transfer model.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import FittingComponent
from ...core.basics.log import log
from ...core.tools import filesystem as fs
from ...core.tools import tables, time
from ...core.basics.table import SmartTable
from ...core.filter.filter import parse_filter_from_instrument_and_band
from ...core.simulation.remote import get_simulation_id, get_simulation_for_host
from ...core.launch.analyser import SimulationAnalyser
from .generation import Generation, GenerationInfo
from ...core.tools.utils import lazyproperty
from ...core.tools.stringify import tostr

# -----------------------------------------------------------------

earth_instrument_name = "earth"

# -----------------------------------------------------------------

class FluxDifferencesTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(FluxDifferencesTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Instrument", str, None, "Instrument")
        self.add_column_info("Band", str, None, "Band")
        self.add_column_info("Flux difference", float, "Jy", "Flux difference")
        self.add_column_info("Relative difference", float, None, "Relative flux difference")
        self.add_column_info("Chi squared term", float, None, "Chi squared term")

    # -----------------------------------------------------------------

    def instruments(self):

        """
        Thisf ucntion ...
        :return:
        """

        return self.get_column_values("Instrument")

    # -----------------------------------------------------------------

    def bands(self):

        """
        This function ...
        :return:
        """

        return self.get_column_values("Band")

    # -----------------------------------------------------------------

    def differences(self, relative=True, add_unit=True):

        """
        This function ...
        :param relative:
        :param add_unit:
        :return:
        """

        if relative: return self.get_column_values("Relative difference")
        else: return self.get_column_values("Flux difference", add_unit=add_unit)

    # -----------------------------------------------------------------

    def chi_squared_terms(self):

        """
        This function ...
        :return:
        """

        return self.get_column_values("Chi squared term")

    # -----------------------------------------------------------------

    def filters(self, sort=False):

        """
        This function ...
        :param sort:
        :return:
        """

        # Get the filter objects
        filters = [parse_filter_from_instrument_and_band(instrument, band) for instrument, band in zip(self.instruments(), self.bands())]

        # Sort
        if sort: return list(sorted(filters, key=lambda fltr: fltr.wavelength.to("micron").value))
        else: return filters

    # -----------------------------------------------------------------

    def wavelengths(self, sort=False, add_unit=True, unit="micron"):

        """
        This function ...
        :param sort:
        :param add_unit:
        :param unit:
        :return:
        """

        # Return
        if add_unit: return [fltr.wavelength.to(unit) for fltr in self.filters(sort=sort)]
        else: return [fltr.wavelength.to(unit).value for fltr in self.filters(sort=sort)]

    # -----------------------------------------------------------------

    def index_for_filter(self, fltr):

        """
        Thisf unction ...
        :param fltr:
        :return:
        """

        return tables.find_index(self, [fltr.instrument, fltr.band], ["Instrument", "Band"])

    # -----------------------------------------------------------------

    def filter_for_index(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        instrument = self.get_value("Instrument", index)
        band = self.get_value("Band", index)
        return parse_filter_from_instrument_and_band(instrument, band)

    # -----------------------------------------------------------------

    def wavelength_for_index(self, index):

        """
        This function ....
        :param index:
        :return:
        """

        return self.filter_for_index(index).wavelength

    # -----------------------------------------------------------------

    def difference_for_filter(self, fltr, add_unit=True):

        """
        This function ...
        :param fltr:
        :param add_unit:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the difference
        #return self["Flux difference"][index]
        return self.get_value("Flux difference", index, add_unit=add_unit)

    # -----------------------------------------------------------------

    def relative_difference_for_filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the difference
        return self["Relative difference"][index]

    # -----------------------------------------------------------------

    def chi_squared_term_for_filter(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        # Get the index
        index = self.index_for_filter(fltr)

        # Return the chi squared term
        return self["Chi squared term"][index] if not self["Chi squared term"].mask[index] else None

    # -----------------------------------------------------------------

    def add_entry(self, instrument, band, difference, relative_difference, chi_squared_term=None, sort=False):

        """
        This function ...
        :param instrument:
        :param band:
        :param difference:
        :param relative_difference:
        :param chi_squared_term:
        :param sort:
        :return:
        """

        self.add_row([instrument, band, difference, relative_difference, chi_squared_term])

        # Sort the table
        if sort: self.sort()

    # -----------------------------------------------------------------

    def sort(self):

        """
        This function ...
        :return:
        """

        # Add a column with the wavelengths
        self["Wavelength"] = self.wavelengths(add_unit=False)
        self["Wavelength"].unit = "micron"

        # Sort the table on the wavelength column
        super(FluxDifferencesTable, self).sort("Wavelength")

        # Delete the wavelength column
        self.remove_column("Wavelength")

    # -----------------------------------------------------------------

    def add_from_filter_and_fluxes(self, fltr, flux, reference_flux, sort=False):

        """
        This function ...
        :param fltr:
        :param flux:
        :param reference_flux:
        :param sort:
        :return:
        """

        # Calculate the difference
        difference = flux - reference_flux
        relative_difference = difference / reference_flux

        # Add
        self.add_entry(fltr.instrument, fltr.band, difference, relative_difference, sort=sort)

# -----------------------------------------------------------------

class SEDFitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(SEDFitModelAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The fitting run
        self.fitting_run = None

        # The flux differences table
        self.differences = None

        # The calculated chi squared value
        self.chi_squared = None

        # The mock observed SEDs
        self.mock_seds = dict()

        # The weights
        self._weights = None

    # -----------------------------------------------------------------

    @classmethod
    def for_simulation(cls, simulation, config=None, **kwargs):

        """
        This function ...
        :param simulation:
        :param config:
        :param kwargs:
        :return:
        """

        # Create the instance
        analyser = cls(config, **kwargs)

        # Set the modeling path as the working path for this class
        analyser.config.path = simulation.analysis.modeling_path

        # Set the task
        analyser.simulation = simulation

        # Return the instance
        return analyser

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Calculate the differences
        self.calculate_differences()

        # 3. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 4. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def clear(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Clearing the fit model analyser ...")

        # Set the attributes to default values
        self.simulation = None
        self.differences = None
        self.chi_squared = None

    # -----------------------------------------------------------------

    @property
    def has_simulation(self):
        return self.simulation is not None

    # -----------------------------------------------------------------

    @property
    def simulation_base_path(self):
        return self.simulation.base_path

    # -----------------------------------------------------------------

    @property
    def generation_path(self):
        return fs.directory_of(self.simulation_base_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation_info(self):
        return GenerationInfo.from_generation_path(self.generation_path)

    # -----------------------------------------------------------------

    @lazyproperty
    def generation(self):
        return Generation.from_path(self.generation_path)

    # -----------------------------------------------------------------

    @property
    def generations_path(self):
        return fs.directory_of(self.generation_path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_path(self):
        return fs.directory_of(self.generations_path)

    # -----------------------------------------------------------------

    @property
    def generation_name(self):
        return fs.name(self.generation_path)

    # -----------------------------------------------------------------

    @property
    def fitting_run_name(self):
        return fs.name(self.fitting_run_path)

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(SEDFitModelAnalyser, self).setup(**kwargs)

        # Get the simulation
        if "simulation" in kwargs: self.simulation = kwargs.pop("simulation")
        elif not self.has_simulation: self.load_simulation()

        # Get the mock sed(s)
        if "mock_seds" in kwargs: self.mock_seds = kwargs.pop("mock_seds")
        elif "mock_sed" in kwargs: self.mock_seds[earth_instrument_name] = kwargs.pop("mock_sed")
        else: self.get_mock_seds(**kwargs)

        # Load the fitting run
        if "fitting_run" in kwargs: self.fitting_run = kwargs.pop("fitting_run")
        else: self.fitting_run = self.load_fitting_run(self.fitting_run_name)

        # Get weights?
        if "weights" in kwargs: self._weights = kwargs.pop("weights")

        # Initialize the differences table
        self.differences = FluxDifferencesTable()

    # -----------------------------------------------------------------

    @property
    def weights(self):

        """
        This function ...
        :return:
        """

        if self._weights is not None: return self._weights
        else: return self.fitting_run.weights

    # -----------------------------------------------------------------

    @property
    def simulation_id(self):

        """
        This function ...
        :return:
        """

        # Determine simulation ID
        if self.config.name is not None:
            if self.config.id is not None: raise ValueError("Cannot specifiy both name and simulation ID")
            return get_simulation_id(self.config.remote, self.config.name)
        else: return self.config.id

    # -----------------------------------------------------------------

    def load_simulation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading simulation with ID '" + str(self.simulation_id) + "' from remote host '" + self.config.remote + "' ...")

        # Load
        self.simulation = get_simulation_for_host(self.config.remote, self.simulation_id)

    # -----------------------------------------------------------------

    def get_mock_seds(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Get a reference to the flux calculator
        simulationanalyser = self.get_analyser(**kwargs)

        # Get mock observed fluxes
        if simulationanalyser.basic_analyser.mock_seds_from_images is not None: self.mock_seds = simulationanalyser.basic_analyser.mock_seds_from_images
        elif simulationanalyser.basic_analyser.mock_seds is not None: self.mock_seds = simulationanalyser.basic_analyser.mock_seds
        else: raise ValueError("No mock observed fluxes could be obtained from the simulation analyser")

    # -----------------------------------------------------------------

    def get_analyser(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Load from kwargs
        if "simulation_analyser" in kwargs: return kwargs.pop("simulation_analyser")
        else: return self.load_analyser()

    # -----------------------------------------------------------------

    def load_analyser(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Loading the simulation analyser ...")

        # Create simulation analyser
        analyser = SimulationAnalyser()
        analyser.config.extra = False

        # Run the analyser on the simulation
        analyser.run(simulation=self.simulation)

        # Return the analyser
        return analyser

    # -----------------------------------------------------------------

    @lazyproperty
    def fit_sed(self):

        """
        This function ...
        :return:
        """

        # Load the appropriate observed SED
        if self.generation_info.fit_not_clipped: sed = self.truncated_sed.copy()
        else: sed = self.observed_sed.copy()

        # Add additional relative error?
        if self.config.additional_error is not None: sed.add_relative_error(self.config.additional_error)

        # Return the SED
        return sed

    # -----------------------------------------------------------------

    def has_weight(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Find the index of the current band in the weights table
        index = tables.find_index(self.weights, key=[fltr.instrument, fltr.band], column_name=["Instrument", "Band"])
        return index is not None

    # -----------------------------------------------------------------

    def get_weight(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Get index
        index = tables.find_index(self.weights, key=[fltr.instrument, fltr.band], column_name=["Instrument", "Band"])
        weight = self.weights["Weight"][index] # apparently, this is a string, so parsing the table went wrong ...
        return float(weight)

    # -----------------------------------------------------------------

    def has_observed_flux(self, fltr):

        """
        Thisn function ...
        :param fltr:
        :return:
        """

        observed_fluxdensity = self.fit_sed.photometry_for_filter(fltr)
        return observed_fluxdensity is not None

    # -----------------------------------------------------------------

    def get_observed_flux(self, fltr, unit="Jy"):

        """
        This function ...
        :param fltr:
        :param unit:
        :return:
        """

        # Find the corresponding flux in the SED derived from observation
        observed_fluxdensity = self.fit_sed.photometry_for_filter(fltr, unit=unit).value

        # Find the corresponding flux error in the SED derived from observation
        observed_fluxdensity_error = self.fit_sed.error_for_filter(fltr, unit=unit).average.to(unit).value

        # Return the values
        return observed_fluxdensity, observed_fluxdensity_error

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences between the observed and simulated SED ...")

        # In the flux-density tables derived from the simulation (created by the ObservedFluxCalculator object),
        # search the one corresponding to the "earth" instrument
        if earth_instrument_name not in self.mock_seds: raise RuntimeError("Could not find a mock observed SED for the 'earth' instrument")

        # Get the mock SED
        mock_sed = self.mock_seds[earth_instrument_name]

        # Loop over the entries in the fluxdensity table (SED) derived from the simulation
        for i in range(len(mock_sed)):

            # Get instrument, band and flux density
            instrument = mock_sed["Instrument"][i]
            band = mock_sed["Band"][i]
            fluxdensity = mock_sed["Photometry"][i]
            fltr = parse_filter_from_instrument_and_band(instrument, band)

            # No weight?
            if not self.has_weight(fltr):
                log.warning("A weight is not found for the '" + str(fltr) + "' filter: skipping ...")
                continue

            # Get the weight
            weight = self.get_weight(fltr)

            # Show the weight
            log.debug("The weight for the '" + str(fltr) + "' filter is " + tostr(weight))

            # No observed flux?
            if not self.has_observed_flux(fltr):
                log.warning("The observed flux density could not be found for the " + instrument + " '" + str(fltr) + "' filter")
                continue

            # Get fluxdensity and error
            observed_fluxdensity, observed_fluxdensity_error = self.get_observed_flux(fltr)

            # Debugging
            log.debug("Calculating relative difference and chi squared term for the '" + str(fltr) + "' filter ...")
            log.debug("Observed flux: " + tostr(observed_fluxdensity) + " ± " + tostr(observed_fluxdensity_error))

            # Calculate difference and relative difference
            difference = fluxdensity - observed_fluxdensity
            relative_difference = difference / observed_fluxdensity

            # Calculate the chi squared term
            chi_squared_term = weight * difference ** 2 / observed_fluxdensity_error ** 2

            # Add an entry to the differences table
            self.differences.add_entry(instrument, band, difference, relative_difference, chi_squared_term)

    # -----------------------------------------------------------------

    @property
    def ndifferences(self):
        return len(self.differences)

    # -----------------------------------------------------------------

    @property
    def nfree_parameters(self):
        return self.fitting_run.nfree_parameters

    # -----------------------------------------------------------------

    @property
    def ndof(self):
        return self.ndifferences - self.nfree_parameters - 1 # number of data points - number of fitted parameters - 1

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the chi squared value for this model ...")

        # The (reduced) chi squared value is the sum of all the terms (for each band),
        # divided by the number of degrees of freedom
        self.chi_squared = np.sum(self.differences["Chi squared term"]) / self.ndof

        # Debugging
        log.debug("Found a (reduced) chi squared value of " + str(self.chi_squared))

    # -----------------------------------------------------------------

    @lazyproperty
    def chi_squared_table(self):
        return self.fitting_run.chi_squared_table_for_generation(self.generation_name)

    # -----------------------------------------------------------------

    @property
    def simulation_name(self):
        return self.simulation.name

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Update the status of the generation if necessary
        self.update_generation()

        # Write the flux differences
        self.write_differences()

        # Write the chi-squared value
        self.write_chi_squared()

    # -----------------------------------------------------------------

    def update_generation(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.info("Updating the generation status ...")

        # Find the index in the table for this generation
        index = tables.find_index(self.fitting_run.generations_table, self.generation_name, "Generation name")

        # Get the number of simulations for this generation
        nsimulations = self.fitting_run.generations_table["Number of simulations"][index]

        # Get the number of entries in the chi squared table
        nfinished_simulations = len(self.chi_squared_table)

        # If this is the last simulation
        if nsimulations == nfinished_simulations + 1:

            # Update the generations table
            self.fitting_run.generations_table.set_finishing_time(self.generation_name, time.timestamp())
            self.fitting_run.generations_table.save()

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table with the flux density differences for the current model ...")

        # Determine the path to the differences table
        path = fs.join(self.simulation.analysis.misc.path, "differences.dat")

        # Debugging
        log.debug("Writing the differences to '" + path + "' ...")

        # Save the differences table
        self.differences.saveto(path)

    # -----------------------------------------------------------------

    def write_chi_squared(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding the chi squared value for the current model to the chi squared data file ...")

        # Add entry
        #self.chi_squared_table.add_entry(self.simulation.name, self.chi_squared)
        self.chi_squared_table.add_or_set_chi_squared(self.simulation_name, self.chi_squared)

        # Save the table
        self.chi_squared_table.save()

# -----------------------------------------------------------------

class ImageResidualsTable(SmartTable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImageResidualsTable, self).__init__(*args, **kwargs)

        # Add column info
        self.add_column_info("Instrument", str, None, "Instrument")
        self.add_column_info("Band", str, None, "Band")
        self.add_column_info("Mean relative difference", float, None, "mean relative residual")
        self.add_column_info("Median relative difference", float, None, "median relative residual")
        self.add_column_info("Standard deviation of relative difference", float, None, "standard deviation of relative residual")
        self.add_column_info("Chi squared term", float, None, "Chi squared term")

    # -----------------------------------------------------------------

    def add_entry(self, instrument, band, mean, median, standard_deviation, chi_squared_term):

        """
        This function ...
        :param instrument:
        :param band:
        :param mean:
        :param median:
        :param standard_deviation:
        :param chi_squared_term
        :return:
        """

        self.add_row([instrument, band, mean, median, standard_deviation, chi_squared_term])

# -----------------------------------------------------------------

class ImagesFitModelAnalyser(FittingComponent):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImagesFitModelAnalyser, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation object
        self.simulation = None

        # The fitting run
        self.fitting_run = None

        # The observed image maker
        self.image_maker = None

    # -----------------------------------------------------------------

    @classmethod
    def for_simulation(cls, simulation, config=None, **kwargs):

        """
        This function ...
        :param simulation:
        :param config:
        :param kwargs:
        :return:
        """

        # Create the instance
        analyser = cls(config, **kwargs)

        # Set the modeling path as the working path for this class
        analyser.config.path = simulation.analysis.modeling_path

        # Set the task
        analyser.simulation = simulation

        # Return the instance
        return analyser

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Calculate the residual maps
        self.calculate_residuals()

        # 3. Calculate the chi squared for this model
        self.calculate_chi_squared()

        # 4. Load the chi squared table
        self.load_chi_squared_table()

        # 5. Update the status of the generation if necessary
        self.update_generation()

        # 6. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImagesFitModelAnalyser, self).setup(**kwargs)

    # -----------------------------------------------------------------

    def calculate_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the residual maps of the observed and simulated images ...")

    # -----------------------------------------------------------------

    def calculate_chi_squared(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def load_chi_squared_table(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def update_generation(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def calculate_chi_squared_from_differences(differences, nfree_parameters):

    """
    This function ...
    :param differences:
    :param nfree_parameters:
    :return:
    """

    # Calculate the degrees of freedom
    ndifferences = len(differences)
    ndof = ndifferences - nfree_parameters - 1  # number of data points - number of fitted parameters - 1

    # The (reduced) chi squared value is the sum of all the terms (for each band),
    # divided by the number of degrees of freedom
    chi_squared = np.sum(differences["Chi squared term"]) / ndof

    # Debugging
    log.debug("Found a (reduced) chi squared value of " + str(chi_squared))

    # Return
    return chi_squared

# -----------------------------------------------------------------
