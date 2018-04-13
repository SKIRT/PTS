# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package

# Import the relevant PTS classes and modules
from pts.modeling.core.environment import load_modeling_environment_cwd
from pts.core.basics.configuration import ConfigurationDefinition, parse_arguments
from pts.core.filter.filter import parse_filter
from pts.core.tools import filesystem as fs
from pts.magic.tools.headers import get_header, get_filter
from pts.magic.core.frame import Frame
from pts.magic.config.plot_residuals import definition as plot_definition
from pts.magic.plot.imagegrid import ResidualImageGridPlotter

# ------------------------------------------------------------------------------

# Load the fitting runs
environment = load_modeling_environment_cwd()
runs = environment.fitting_runs

# Get the galaxy center
center = environment.galaxy_center
ra = center.ra.to("deg").value
dec = center.dec.to("deg").value

# ------------------------------------------------------------------------------

light_theme = "light"
dark_theme = "dark"
themes = [light_theme, dark_theme]

# ------------------------------------------------------------------------------

default_filter_names = ["FUV", "SDSS r", "I1", "MIPS 24mu", "Pacs red"]
default_filters = [parse_filter(name) for name in default_filter_names]

#other_filter_names = ["NUV", "SDSS i", "I1", "Pacs blue", "SPIRE 250"]
#other_filters = [parse_filter(name) for name in other_filter_names]

# ------------------------------------------------------------------------------

default_cmap = "inferno"
default_residual_cmap = 'RdBu'

# ------------------------------------------------------------------------------

# Create definition
definition = ConfigurationDefinition(write_config=False)

# The fitting run for which to explore the parameter space
if runs.empty: raise RuntimeError("No fitting runs are present")
elif runs.has_single: definition.add_fixed("run", "name of the fitting run", runs.single_name)
else: definition.add_required("run", "string", "name of the fitting run", choices=runs.names)

# Generation name
definition.add_required("generation", "string", "generation name")

# Filters
definition.add_positional_optional("filters", "filter_list", "filters for which to plot the images", default_filters)

# Import plotting settings
definition.import_settings(plot_definition)

# ------------------------------------------------------------------------------

# Get configuration
config = parse_arguments("plot_best_images", definition, "plot images of best simulation")

# ------------------------------------------------------------------------------

# Load the generation
run = runs.load(config.run)
generation = run.get_generation(config.generation)

# Get the chi squared table
chi_squared = generation.chi_squared_table

# Get the name of the best simulation
simulation_name = chi_squared.best_simulation_name

# ------------------------------------------------------------------------------

# Get the mock images
image_paths = generation.get_image_paths_for_simulation(simulation_name)

nfilters = len(config.filters)

# Create a list for the images
mock_images = [None] * nfilters

# Loop over the image paths
for filepath in image_paths:

    # Get image name
    name = fs.strip_extension(fs.name(filepath))

    # Get header
    header = get_header(filepath)

    # Get the filter
    fltr = get_filter(name, header=header)

    # Check whether the filter is in the list of filters to be plotted
    if fltr not in config.filters: continue

    # Get the index for this filter
    index = config.filters.index(fltr)

    # Load the image
    frame = Frame.from_file(filepath)

    # Replace zeroes and negatives
    frame.replace_zeroes_by_nans()
    frame.replace_negatives_by_nans()

    # Set the image
    mock_images[index] = frame

# ------------------------------------------------------------------------------

# Get the observed images
observed_images = []
for fltr in config.filters:

    # Check
    if not environment.has_photometry_image_for_filter(fltr):
        observed_images.append(None)
        continue

    # Get the photometry image
    filepath = environment.photometry_image_paths_for_filters[fltr]

    # Load the image
    frame = Frame.from_file(filepath)

    # Replace zeroes and negatives
    frame.replace_zeroes_by_nans()
    frame.replace_negatives_by_nans()

    # Set the image
    observed_images.append(frame)

# ------------------------------------------------------------------------------

# Create the plotter
plotter = ResidualImageGridPlotter(config=config)

# Add the images
for index in range(config.filters):

    # Get the filter
    fltr = config.filters[index]
    name = str(fltr)

    # Get frames
    observed = observed_images[index]
    model = mock_images[index]

    # Add to the plotter
    plotter.add_image(name, observed, model)

# Run the plotter
plotter.run()

# ------------------------------------------------------------------------------
