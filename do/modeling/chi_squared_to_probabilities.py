#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.chi_squared_to_probabilities Convert chi squared table to proabilities.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.fitting.tables import ChiSquaredTable, ParametersTable
from pts.modeling.fitting.sedfitting import chi_squared_table_to_probabilities
from pts.core.basics.log import log
from pts.core.basics.distribution import Distribution
from pts.core.plot.distribution import plot_distribution

# -----------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition()

# The fitting run
definition.add_required("chisquared", "file_path", "name of the chi squared table file")
definition.add_required("parameters", "file_path", "name of parameters table file")

# For plotting
definition.add_flag("plot_ranks", "plot chi squared as a function of rank")
definition.add_flag("plot_distribution", "plot distribution of chi squared values")
definition.add_optional("max_chisquared", "positive_real", "maximum chi squared to show")
definition.add_optional("min_probability", "positive_real", "minimum probability to show")

# Create the configuration
config = parse_arguments("chi_squared_to_probabilities", definition, "Convert chi squared table to proabilities")

# -----------------------------------------------------------------

# Load the table
chi_squared = ChiSquaredTable.from_file(config.chisquared)

# Get the parameters table
parameters_table = ParametersTable.from_file(config.parameters)

# -----------------------------------------------------------------

probabilities_table = chi_squared_table_to_probabilities(chi_squared, parameters_table)
probabilities = probabilities_table.probabilities
nsimulations = len(probabilities)

# -----------------------------------------------------------------

# Get the number of zero probabilities
nzeros = probabilities_table.nzeros

#nzeros = sequences.nzeros(probabilities_table["Probability"])
#nzeros2 = sequences.nzeros(list(probabilities_table["Probability"]))
#print(nzeros, nzeros2)

# Show number of zeros
if nzeros > 5: log.warning(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probabilities of zero")
else: log.debug(str(nzeros) + " out of " + str(nsimulations) + " simulations have a probability of zero")

# -----------------------------------------------------------------

# Plot ranks?
if config.plot_ranks:

    # Plot histogram of the chi squared values w.r.t. simulation rank
    distribution = Distribution.by_rank("Simulation rank", chi_squared.chi_squared_values, y_name="Chi squared")
    plot_distribution(distribution, statistics=False, )

    # Plot histogram of the probability values w.r.t. simulation rank
    distribution = Distribution.by_rank("Simulation rank", probabilities, y_name="Probability")
    plot_distribution(distribution, statistics=False, color="red", logfrequency=True)

# -----------------------------------------------------------------

# Plot distribution?
if config.plot_distribution:

    # Make distribution of chi squared values and plot histogram
    distribution2 = Distribution.from_values("Chi squared", chi_squared.chi_squared_values, nbins=50, clip_above=config.max_chisquared, density=False)
    plot_distribution(distribution2)

    # Make distribution of probability values and plot histogram
    distribution2 = Distribution.from_values("Probability", probabilities, nbins=50, ignore_value=0, clip_below=config.min_probability, density=False)
    plot_distribution(distribution2, color="red")

# -----------------------------------------------------------------
