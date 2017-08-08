#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.modeling.replot Replot the SEDs of individual simulations.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.tools import filesystem as fs
from pts.core.basics.log import log
from pts.modeling.component.component import load_modeling_configuration
from pts.modeling.component.galaxy import get_observed_sed as get_sed_galaxy
from pts.modeling.component.sed import get_observed_sed as get_sed_other
from pts.modeling.welcome import welcome
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.modeling.fitting.component import get_generation_names, get_simulations
from pts.core.plot.sed import SEDPlotter
from pts.core.data.sed import SED
from pts.modeling.core.environment import verify_modeling_cwd

# -----------------------------------------------------------------

# Set the modeling path
modeling_path = verify_modeling_cwd()

# -----------------------------------------------------------------

# Get names of all generations
all_generations = get_generation_names(modeling_path)

# Create configuration definition
definition = ConfigurationDefinition()
definition.add_optional("generations", "string_list", "generations for which to replot the SEDs", default=all_generations, choices=all_generations)

# -----------------------------------------------------------------

# Create the configuration
config = parse_arguments("replot", definition, description="replot the SEDs of individual simulations")

# Welcome message
welcome()

# -----------------------------------------------------------------

# Load the modeling configuration
modeling_config = load_modeling_configuration(modeling_path)

# -----------------------------------------------------------------

# Load the observed SED

# Galaxy modeling
if modeling_config.modeling_type == "galaxy": observed_sed = get_sed_galaxy(modeling_path)

# Other (SED modeling)
elif modeling_config.modeling_type == "other": observed_sed = get_sed_other(modeling_path)

# Invalid
else: raise RuntimeError("Invalid modeling type: " + modeling_config.modeling_type)

# -----------------------------------------------------------------

# Loop over the generations
for generation_name in config.generations:

    # Loop over the simulations
    for simulation in get_simulations(modeling_path, generation_name):

        # Inform the user
        log.info("Plotting SED for simulation '" + simulation.name + "' of generation '" + generation_name + "' ...")

        # Determine plot path
        plot_dir_path = fs.join(simulation.base_path, "plot")
        plot_path = fs.join(plot_dir_path, "sed.pdf")

        # Create SED plotter
        plotter = SEDPlotter()

        # Add observed sed
        plotter.add_sed(observed_sed, "observation")

        # Add simulated sed(s)
        for path in simulation.seddatpaths():
            sed = SED.from_skirt(path)
            plotter.add_sed(sed, fs.name(path))

        # Run the plotter
        plotter.run(output=plot_path)

# -----------------------------------------------------------------
