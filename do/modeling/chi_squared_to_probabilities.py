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

# Create the configuration
config = parse_arguments("chi_squared_to_probabilities", definition, "Convert chi squared table to proabilities")

# -----------------------------------------------------------------

# Load the table
chi_squared = ChiSquaredTable.from_file(config.chisquared)

# Get the parameters table
parameters_table = ParametersTable.from_file(config.parameters)

# -----------------------------------------------------------------

probabilities_table = chi_squared_table_to_probabilities(chi_squared, parameters_table)
probabilities = probabilities_table
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

# Plot histogram of the chi squared values w.r.t. simulation rank
distribution = Distribution.by_rank("Simulation rank", chi_squared.chi_squared_values, y_name="Chi squared")
plot_distribution(distribution, statistics=False)

# -----------------------------------------------------------------

# Make distribution of chi squared values and plot histogram
distribution2 = Distribution.from_values("Chi squared", chi_squared.chi_squared_values, nbins=50)
plot_distribution(distribution2)

# -----------------------------------------------------------------
