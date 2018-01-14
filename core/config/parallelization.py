#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

# Import the relevant PTS classes and modules
from pts.core.basics.configuration import ConfigurationDefinition

# -----------------------------------------------------------------

# Create the configuration
definition = ConfigurationDefinition()

# Add prequired
definition.add_required("ski", "file_path", "name or path of the ski file")

# Add optional argument
definition.add_optional("input", "directory_path", "path to the input directory")

# Add required
definition.add_required("nnodes", "integer", "number of nodes")
definition.add_required("nsockets", "integer", "number of sockets per node")
definition.add_required("ncores", "integer", "number of cores per socket")
definition.add_required("memory", "real", "available virtual memory per node")

# Add flags
definition.add_flag("mpi", "mpi available", True)
definition.add_flag("hyperthreading", "use hyperthreading", False)
definition.add_optional("threads_per_core", "integer", "number of hyperthreads per core")

# Maximum number of threads per process
definition.add_optional("max_nthreads", "positive_integer", "maximum allowed number of threads per process", 8)

# Add optional
definition.add_optional("ncells", "integer", "number of dust cells (relevant if ski file uses a tree dust grid)")
definition.add_optional("nwavelengths", "integer", "number of wavelengths (relevant if ski file uses input file)")

# Flag
definition.add_flag("data_parallel", "use data-parallelization (None means automatic)", None)

# Flags
definition.add_flag("show", "show the parallelization", True)

# -----------------------------------------------------------------

# Advanced
definition.add_optional("min_nwavelengths_per_process", "positive_integer", "minimum number of wavelengths per process for data parallelization", 10)

# -----------------------------------------------------------------
