#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.do.developer.jobscript_template Generate a job script template.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import the relevant PTS classes and modules
from pts.core.remote.host import find_host_ids
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.simulation.jobscript import SKIRTJobScript
from pts.core.launch.options import LoggingOptions
from pts.core.simulation.arguments import SkirtArguments
from pts.core.simulation.definition import SingleSimulationDefinition
from pts.core.remote.host import load_host
from pts.core.simulation.parallelization import Parallelization

# -----------------------------------------------------------------

# Create the configuration definition
definition = ConfigurationDefinition()

# Add required arguments
definition.add_required("remote", "string", "remote host ID", choices=find_host_ids())
definition.add_required("ski", "string", "the name/path of the ski file")

# Input and output
definition.add_optional("input", "string", "input directory for the simulation(s)", letter="i")
definition.add_optional("output", "string", "output directory for the simulation(s)", "~/SKIRT/run", letter="o", convert_default=True)

definition.add_optional("cluster", "string", "the name of the cluster", letter="c")
definition.add_optional("parallel", "integer_pair", "the parallelization scheme (processes, threads)", letter="p", default=(8,4))
definition.add_optional("walltime", "duration", "an estimate for the walltime of the simulation for the specified parallelization scheme", "12:00:00", convert_default=True)
definition.add_flag("data_parallel", "enable data parallelization", None)

# Logging options
definition.import_section_from_composite_class("logging", "logging options", LoggingOptions)

# Create the configuration
config = parse_arguments("installation_commands", definition)

# -----------------------------------------------------------------

host = load_host(config.remote, config.cluster)

modules = []
modules.append("iimpi/2016b")

ski_path = config.ski
output_path = config.output
input_path = config.input
definition = SingleSimulationDefinition(ski_path, output_path, input_path)

logging = LoggingOptions()
logging.set_options(config.logging)

cores = config.parallel[0] * config.parallel[1]
threads_per_core = 2
parallelization = Parallelization(cores, threads_per_core, config.parallel[0])

arguments = SkirtArguments(definition, logging, parallelization)

name = definition.prefix
skirt_path = "skirt"
mpi_command = "mpirun"
walltime = config.walltime
jobscript = SKIRTJobScript(name, arguments, host.id, host.cluster, skirt_path, mpi_command, walltime, modules, mail=False, bind_to_cores=False)

lines = jobscript.to_lines()

for line in lines: print(line)

# -----------------------------------------------------------------
