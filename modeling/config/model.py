#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition
from pts.core.remote.host import find_host_ids

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add flags
definition.add_flag("check_hosts", "check the availability of the remote hosts", True)
definition.add_flag("deploy", "deploy SKIRT and PTS where necessary", False)
definition.add_flag("check_versions", "check versions of SKIRT and PTS where necessary", True)

# Advanced settings
definition.add_flag("local", "keep computationaly heavy computations local")
definition.add_optional("remotes", "string_list", "remote hosts for computationally heavy computations (overrule the modeling configuration)", choices=find_host_ids())
definition.add_flag("fitting_local", "launch the simulations as part of the fitting locally (overrule the modeling configuration")
definition.add_optional("fitting_remotes", "string_list", "remote hosts for the fitting (overrule the modeling configuration)", choices=find_host_ids())

# Genetic algorithm settings
definition.add_optional("nsimulations", "positive_integer", "number of simulations per generation")
definition.add_optional("ngenerations", "positive_integer", "number of generations to run in one run (ngenerations > 1 is only allowed for local execution)", 1)
definition.add_flag("finish", "don't launch a new generation but finish evaluating all previous runs")

# Model settings
definition.add_optional("npackages_factor", "positive_real", "the factor with which to increase the number of photon packages for the new batch of simulations", 5.)
definition.add_flag("increase_npackages", "increase the number of photon packages with a certain factor", False)
definition.add_flag("refine_wavelengths", "increase the resolution of the wavelength grid for the new batch of simulations", False)
definition.add_flag("refine_dust", "increase the resolution of the dust cell grid for the new batch of simulations", False)
definition.add_flag("selfabsorption", "dust self-absorption (None means the initial values are respected)", None)
definition.add_flag("transient_heating", "transient (non-LTE) dust heating", None)

# -----------------------------------------------------------------
