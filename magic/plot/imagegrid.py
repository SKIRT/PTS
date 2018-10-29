# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.plot.imagegrid Contains the ImageGridPlotter classes.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import aplpy
from abc import ABCMeta, abstractproperty
import matplotlib.pyplot as plt
from matplotlib import cm
from collections import OrderedDict, defaultdict

# Import the relevant PTS classes and modules
from ..tools.plotting import get_vmin_vmax
from ...core.tools import filesystem as fs
from ..core.frame import Frame
from ...core.basics.log import log
from ...core.basics.configurable import Configurable
from ...core.tools.utils import lazyproperty, memoize_method
from ...core.tools import sequences
from ..core.image import Image
from ...core.basics.distribution import Distribution
from ...core.basics.plot import MPLFigure
from ...core.basics.composite import SimplePropertyComposite
from ...core.basics.plot import normal_colormaps
from ..core.list import uniformize
from ...core.tools import numbers
from ...core.tools import types

# ------------------------------------------------------------------------------

light_theme = "light"
dark_theme = "dark"
themes = [light_theme, dark_theme]

# ------------------------------------------------------------------------------

default_cmap = "inferno"
default_residual_cmap = 'RdBu'
default_absolute_residual_cmap = "OrRd"

# ------------------------------------------------------------------------------

# Initialize dictionary for light theme settings
light_theme_settings = OrderedDict()

# Set parameters
light_theme_settings['axes.facecolor'] = 'white'
light_theme_settings['savefig.facecolor'] = 'white'
light_theme_settings['axes.edgecolor'] = 'black'
light_theme_settings['xtick.color'] = 'black'
light_theme_settings['ytick.color'] = 'black'
light_theme_settings["axes.labelcolor"] = 'black'
light_theme_settings["text.color"] = 'black'
# light_theme_settings["axes.titlecolor"]='black'

# ------------------------------------------------------------------------------

# Initialize dictionary for dark theme settings
dark_theme_settings = OrderedDict()

# Set parameters
dark_theme_settings['axes.facecolor'] = 'black'
dark_theme_settings['savefig.facecolor'] = 'black'
dark_theme_settings['axes.edgecolor'] = 'white'
dark_theme_settings['xtick.color'] = 'white'
dark_theme_settings['ytick.color'] = 'white'
dark_theme_settings["axes.labelcolor"] ='white'
dark_theme_settings["text.color"] = 'white'
#plt.rcParams["axes.titlecolor"] = 'white'

# ------------------------------------------------------------------------------

class ImagePlotSettings(SimplePropertyComposite):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # ------------------------------------------------------------------------------

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ImagePlotSettings, self).__init__()

        # Define properties
        self.add_property("label", "string", "label for the image", None)
        self.add_property("vmin", "real", "plotting minimum")
        self.add_property("vmax", "real", "plotting maximum")
        self.add_boolean_property("soft_vmin", "soft vmin", False) #, None) # use None as default to use plotter config if not defined
        self.add_boolean_property("soft_vmax", "soft vmax", False) #, None) # use None as default to use plotter config if not defined
        self.add_property("cmap", "string", "colormap", choices=normal_colormaps)

# ------------------------------------------------------------------------------

class ImageGridPlotter(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(ImageGridPlotter, self).__init__(*args, **kwargs)

        # The figure
        self.figure = None

        # The grid
        self.grid = None

        # The plots
        self.plots = None

        # The settings
        self.settings = defaultdict(self.image_settings_class)

    # -----------------------------------------------------------------

    @abstractproperty
    def image_settings_class(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @abstractproperty
    def names(self):

        """
        This function ...
        :return:
        """

        pass

    # ------------------------------------------------------------------------------

    @property
    def light(self):
        return self.config.theme == light_theme

    # -----------------------------------------------------------------

    @property
    def dark(self):
        return self.config.theme == dark_theme

    # -----------------------------------------------------------------

    @lazyproperty
    def text_color(self):

        """
        This function ...
        :return:
        """

        # Set light theme
        if self.light: return "black"

        # Dark theme
        elif self.dark: return "white"

        # Invalid
        else: raise ValueError("Invalid theme")

    # -----------------------------------------------------------------

    @lazyproperty
    def frame_color(self):

        """
        This function ...
        :return:
        """

        # Set light theme
        if self.light: return "black"

        # Dark theme
        elif self.dark: return "white"

        # Invalid
        else: raise ValueError("Invalid theme")

    # -----------------------------------------------------------------

    @lazyproperty
    def background_color(self):

        """
        This function ...
        :return:
        """

        # Set light theme
        if self.light: return "white"

        # Dark theme
        elif self.dark: return "black"

        # Invalid
        else: raise ValueError("Invalid theme")

    # -----------------------------------------------------------------

    @abstractproperty
    def first_frame(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    @lazyproperty
    def center(self):

        """
        This function ...
        :return:
        """

        # Center coordinate is defined
        if self.config.center is not None: return self.config.center

        # Not defined?
        return self.first_frame.center_sky

    # -----------------------------------------------------------------

    @property
    def ra_center(self):
        return self.center.ra

    # ------------------------------------------------------------------------------

    @property
    def dec_center(self):
        return self.center.dec

    # ------------------------------------------------------------------------------

    @lazyproperty
    def ra_center_deg(self):
        return self.ra_center.to("deg").value

    # ------------------------------------------------------------------------------

    @lazyproperty
    def dec_center_deg(self):
        return self.dec_center.to("deg").value

    # ------------------------------------------------------------------------------

    @lazyproperty
    def spacing_deg(self):
        return self.config.spacing.to("deg").value

    # ------------------------------------------------------------------------------

    @lazyproperty
    def radius_deg(self):
        return self.config.radius.to("deg").value

    # ------------------------------------------------------------------------------

    @lazyproperty
    def colormap(self):
        return cm.get_cmap(self.config.cmap)

    # -----------------------------------------------------------------

    @lazyproperty
    def nan_color(self):
        if self.config.nan_color is not None: return self.config.nan_color
        else: return self.colormap(0)

    # -----------------------------------------------------------------

    @lazyproperty
    def theme_settings(self):
        if self.light: return light_theme_settings
        elif self.dark: return dark_theme_settings
        else: raise ValueError("Invalid theme")

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ImageGridPlotter, self).setup(**kwargs)

        # plt.rcParams.update({'font.size':20})
        plt.rcParams["axes.labelsize"] = self.config.axes_label_size  # 16 #default 20
        plt.rcParams["xtick.labelsize"] = self.config.ticks_label_size  # 10 #default 16
        plt.rcParams["ytick.labelsize"] = self.config.ticks_label_size  # 10 #default 16
        plt.rcParams["legend.fontsize"] = self.config.legend_fontsize  # 10 #default 14
        plt.rcParams["legend.markerscale"] = self.config.legend_markers_cale
        plt.rcParams["lines.markersize"] = self.config.lines_marker_size  # 4 #default 4
        plt.rcParams["axes.linewidth"] = self.config.linewidth

        # Set theme-specific settings
        for label in self.theme_settings: plt.rcParams[label] = self.theme_settings[label]

        # plt.rcParams['xtick.major.size'] = 5
        # plt.rcParams['xtick.major.width'] = 2
        # plt.rcParams['ytick.major.size'] = 5
        # plt.rcParams['ytick.major.width'] = 2

# ------------------------------------------------------------------------------

def plot_images(images, **kwargs):

    """
    This function ...
    :param images:
    :param kwargs:
    :return:
    """

    # Create the plotter
    plotter = StandardImageGridPlotter(**kwargs)

    # Run the plotter
    plotter.run(images=images)

# -----------------------------------------------------------------

class StandardImagePlotSettings(ImagePlotSettings):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        This function ...
        :param kwargs:
        """

        # Call the constructor of the base class
        super(StandardImagePlotSettings, self).__init__(**kwargs)

        # Set properties
        self.set_properties(kwargs)

# -----------------------------------------------------------------

class StandardImageGridPlotter(ImageGridPlotter):

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
        super(StandardImageGridPlotter, self).__init__(*args, **kwargs)

        # The image frames
        self.frames = OrderedDict()

        # The error frames
        self.errors = OrderedDict()

        # The masks
        self.masks = OrderedDict()

        # The regions
        self.regions = OrderedDict()

    # ------------------------------------------------------------------------------

    @property
    def image_settings_class(self):

        """
        This function ...
        :return:
        """

        return StandardImagePlotSettings

    # ------------------------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Show stuff
        if self.config.show: self.show()

        # Write
        self.write()

        # Plot
        self.plot()

    # ------------------------------------------------------------------------------

    @property
    def names(self):

        """
        This function ...
        :return:
        """

        return self.frames.keys()

    # ------------------------------------------------------------------------------

    def add_image(self, name, image, errors=None, mask=None, regions=None, replace=False, settings=None):

        """
        This function ...
        :param name:
        :param image:
        :param errors:
        :param mask:
        :param regions:
        :param replace:
        :param settings:
        :return:
        """

        # Check if name already exists
        if not replace and name in self.names: raise ValueError("Already an image with name '" + name + "' added")

        # Image is passed
        if isinstance(image, Image):

            # Get the frame
            frame = image.primary

            # Get errors?
            # Get mask?
            # Get regions?

        # Frame is passed
        elif isinstance(image, Frame): frame = image

        # Invalid
        else: raise ValueError("Invalid value for 'image': must be Frame or Image")

        # Add frame
        self.frames[name] = frame

        # Add errors
        if errors is not None: self.errors[name] = errors

        # Add regions
        if regions is not None: self.regions[name] = regions

        # Add mask
        if mask is not None: self.masks[name] = mask

        # Set settings
        if settings is not None: self.settings[name].set_properties(settings)

    # ------------------------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # ------------------------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Images
        if self.config.write_images: self.write_images()

        # Frames
        if self.config.write_frames: self.write_frames()

        # Masks
        if self.config.write_masks: self.write_masks()

        # Regions
        if self.config.write_regions: self.write_regions()

    # ------------------------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

    # ------------------------------------------------------------------------------

    def write_frames(self):

        """
        This function ...
        :return:
        """

    # ------------------------------------------------------------------------------

    def write_masks(self):

        """
        This function ...
        :return:
        """

    # ------------------------------------------------------------------------------

    def write_regions(self):

        """
        This function ...
        :return:
        """

    # ------------------------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

# ------------------------------------------------------------------------------

images_name = "images"
observations_name = "observations"
models_name = "models"
errors_name = "errors"
model_errors_name = "model_errors"
residuals_name = "residuals"
distributions_name = "distributions"
settings_name = "settings"

# ------------------------------------------------------------------------------

observation_name = "observation"
model_name = "model"
observation_or_model = [observation_name, model_name]

# ------------------------------------------------------------------------------

horizontal_mode, vertical_mode = "horizontal", "vertical"
default_direction = vertical_mode
directions = [horizontal_mode, vertical_mode]

# ------------------------------------------------------------------------------

class ResidualImagePlotSettings(ImagePlotSettings):

    """
    This class ...
    """

    def __init__(self, **kwargs):

        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(ResidualImagePlotSettings, self).__init__()

        # Define properties
        self.add_property("residual_amplitude", "percentage", "amplitude of the residual plots")
        self.add_boolean_property("soft_residual_amplitude", "soft residual amplitude", False) #, None) # use None as default to use plotter config if not defined
        self.add_property("residual_cmap", "string", "colormap for the residual plots") # no choices because can be absolute or not

        # Set properties
        self.set_properties(kwargs)

# ------------------------------------------------------------------------------

def plot_residuals(observations, models, **kwargs):

    """
    This function ...
    :param observations:
    :param models:
    :param kwargs:
    :return:
    """

    # Create the plotter
    plotter = ResidualImageGridPlotter(**kwargs)

    # Run the plotter
    plotter.run(observations=observations, models=models)

# -----------------------------------------------------------------

class ResidualImageGridPlotter(ImageGridPlotter):

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
        super(ResidualImageGridPlotter, self).__init__(*args, **kwargs)

        # The image frames
        self.observations = OrderedDict()
        self.errors = OrderedDict()
        self.models = OrderedDict()
        self.model_errors = OrderedDict()
        self.residuals = OrderedDict()

        # The residual distributions
        self.distributions = OrderedDict()

    # ------------------------------------------------------------------------------

    @property
    def image_settings_class(self):
        return ResidualImagePlotSettings

    # ------------------------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Create the residual frames
        self.create_residuals()

        # Create the residual distributions
        self.create_distributions()

        # Show stuff
        if self.config.show: self.show()

        # Write
        self.write()

        # Plot
        self.plot()

    # ------------------------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ResidualImageGridPlotter, self).setup(**kwargs)

        # Load the images
        if kwargs.get(images_name, None) is not None: self.add_images(kwargs.pop(images_name))
        if kwargs.get(observations_name, None) is not None: self.add_observations(kwargs.pop(observations_name))
        if kwargs.get(models_name, None) is not None: self.add_models(kwargs.pop(models_name))
        if kwargs.get(errors_name, None) is not None: self.add_error_maps(kwargs.pop(errors_name))
        if kwargs.get(residuals_name, None) is not None: self.add_residual_maps(kwargs.pop(residuals_name))

        # Nothing added
        if self.config.from_directory is not None: self.load_from_directory(self.config.from_directory)
        elif not self.has_images: self.load_from_directory(self.config.path)

        # Initialize the figure
        self.initialize_figure()

    # ------------------------------------------------------------------------------

    @property
    def figsize(self):
        return (15,10)

    # ------------------------------------------------------------------------------

    @property
    def horizontal(self):
        return self.config.direction == horizontal_mode

    # ------------------------------------------------------------------------------

    @property
    def vertical(self):
        return self.config.direction == vertical_mode

    # ------------------------------------------------------------------------------

    @lazyproperty
    def npanels(self):
        if self.config.distributions: return 4  # observation, model, residual, distribution
        else: return 3  # observation, model, residual

    # ------------------------------------------------------------------------------

    @lazyproperty
    def nrows(self):
        if self.horizontal: return self.npanels
        elif self.vertical: return self.nimages
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    @lazyproperty
    def ncolumns(self):
        if self.horizontal: return self.nimages
        elif self.vertical: return self.npanels
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    @property
    def share_x(self):
        return True

    # ------------------------------------------------------------------------------

    @property
    def share_y(self):
        return True

    # ------------------------------------------------------------------------------

    def initialize_figure(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Initializing the figure with size " + str(self.figsize) + " ...")

        # Create the plot
        self.figure = MPLFigure(size=self.figsize)

        # Create plots
        #self.plots = self.figure.create_grid(self.nrows, self.ncolumns, sharex=self.share_x, sharey=self.share_y)

        # Create grid
        self.grid = self.figure.create_gridspec(self.nrows, self.ncolumns, hspace=0.0, wspace=0.0)

        # Initialize structure to contain the plots
        #print("NCOLUMNS", self.ncolumns)
        #print("NROWS", self.nrows)
        self.plots = [[None for i in range(self.ncolumns)] for j in range(self.nrows)]

    # ------------------------------------------------------------------------------

    @property
    def all_names(self):
        return sequences.combine_unique(self.observation_names, self.model_names, self.errors_names, self.residuals_names)

    # ------------------------------------------------------------------------------

    @property
    def observation_names(self):
        return self.observations.keys()

    # ------------------------------------------------------------------------------

    def has_observation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.observation_names

    # ------------------------------------------------------------------------------

    @property
    def model_names(self):
        return self.models.keys()

    # ------------------------------------------------------------------------------

    def has_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.model_names

    # ------------------------------------------------------------------------------

    @property
    def errors_names(self):
        return self.errors.keys()

    # ------------------------------------------------------------------------------

    def has_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.errors

    # ------------------------------------------------------------------------------

    @property
    def model_errors_names(self):
        return self.model_errors.keys()

    # ------------------------------------------------------------------------------

    def has_model_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.model_errors

    # ------------------------------------------------------------------------------

    @property
    def residuals_names(self):
        return self.residuals.keys()

    # ------------------------------------------------------------------------------

    def has_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.residuals

    # ------------------------------------------------------------------------------

    @property
    def distribution_names(self):
        return self.distributions.keys()

    # ------------------------------------------------------------------------------

    def has_distribution(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.distributions

    # ------------------------------------------------------------------------------

    @property
    def settings_names(self):
        return self.settings.keys()

    # ------------------------------------------------------------------------------

    def has_settings(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.settings_names

    # ------------------------------------------------------------------------------

    @property
    def names(self):
        return self.observation_names

    # ------------------------------------------------------------------------------

    @property
    def first_name(self):
        return self.names[0]

    # ------------------------------------------------------------------------------

    @property
    def first_observation(self):
        return self.get_observation(self.first_name)

    # ------------------------------------------------------------------------------

    @property
    def first_frame(self):
        return self.first_observation

    # ------------------------------------------------------------------------------

    @property
    def nimages(self):
        return len(self.names)

    # ------------------------------------------------------------------------------

    @property
    def has_images(self):
        return self.nimages > 0

    # ------------------------------------------------------------------------------

    def add_image(self, name, observation, model=None, errors=None, model_errors=None, residuals=None, replace=False,
                  settings=None):

        """
        This function ...
        :param name:
        :param observation:
        :param model:
        :param errors:
        :param model_errors:
        :param residuals:
        :param replace:
        :param settings:
        :return:
        """

        # Check if name already exists
        if not replace and name in self.names: raise ValueError("Already an image with name '" + name + "' added")

        # Check type of the image
        if isinstance(observation, Image):

            # Get observation frame
            if observation_name in observation.frame_names: observation = observation.frames[observation_name]
            else: observation = observation.primary

            # Get model frame
            if model_name in observation.frame_names:
                if model is not None: raise ValueError("Cannot pass model frame if image contains model frame")
                model = observation.frames[model_name]

            # Get errors frame
            if errors_name in observation.frame_names:
                if errors is not None: raise ValueError("Cannot pass error map if image contains error map")
                errors = observation.frames[errors_name]

            # Get model errors frame
            if model_errors_name in observation.frame_names:
                if model_errors is not None: raise ValueError("Cannot pass model error map if image contains model error map")
                model_errors = observation.frames[model_errors_name]

            # Get residuals frame
            if residuals_name in observation.frame_names:
                if residuals is not None: raise ValueError("Cannot pass residual map if image contains residual map")
                residuals = observation.frames[residuals_name]

        # Check the type of the model image
        if model is not None and isinstance(model, Image):

            # Get the model frame
            if model_name in model.frame_names: model = model.frames[model_name]
            else: model = model.primary

            # Get the model errors frame
            if model_errors_name in model.frame_names:
                if errors_name in model.frame_names: raise ValueError("Model image contains both 'errors' and 'model_errors' frame")
                if model_errors is not None: raise ValueError("Cannot pass model error map if model image contains model error map")
                model_errors = model.frames[model_errors_name]
            elif errors_name in model.frame_names:
                if model_errors is not None: raise ValueError("Cannot pass model error map if model image contains error map")
                model_errors = model.frames[errors_name]

        # Add observation
        self.observations[name] = observation

        # Add model
        if model is not None: self.models[name] = model

        # Add errors
        if errors is not None: self.errors[name] = errors

        # Add model errors
        if model_errors is not None: self.model_errors[name] = model_errors

        # Add residuals
        if residuals is not None: self.residuals[name] = residuals

        # Set settings
        if settings is not None: self.settings[name].set_properties(settings)

    # ------------------------------------------------------------------------------

    def add_observation(self, name, frame, errors=None):

        """
        This function ...
        :param name:
        :param frame:
        :param errors:
        :return:
        """

        # Check the type of the image
        if isinstance(frame, Image):

            # Get observation frame
            if observation_name in frame.frame_names: frame = frame.frames[observation_name]
            else: frame = frame.primary

            # Get error map
            if errors_name in frame.frame_names:
                if errors is not None: raise ValueError("Cannot pass error map if image contains error map")
                errors = frame.frames[errors_name]

            # Check whether there are no other frames
            if sequences.contains_more(frame.frame_names, ["primary", observation_name, errors_name]): raise ValueError("Observation image contains too many frames")

        # Add observation frame
        self.observations[name] = frame

        # Add error map
        if errors is not None: self.errors[name] = errors

    # ------------------------------------------------------------------------------

    def add_model(self, name, frame, errors=None):

        """
        This function ...
        :param name:
        :param frame:
        :param errors:
        :return:
        """

        # Check the type of the image
        if isinstance(frame, Image):

            # Get model frame
            if model_name in frame.frame_names: frame = frame.frames[model_name]
            else: frame = frame.primary

            # Get error map
            if errors_name in frame.frame_names:
                if model_errors_name in frame.frame_names: raise ValueError("Model image contains both 'errors' and 'model_errors' frame")
                if errors is not None: raise ValueError("Cannot pass error map if image contains error map")
                errors = frame.frames[errors_name]
            elif model_errors_name in frame.frame_names:
                if errors is not None: raise ValueError("Cannot pass error map if image contains error map")
                errors = frame.frames[model_errors_name]

            # Check whether there are no other frames
            if sequences.contains_more(frame.frame_names, ["primary", model_name, errors_name, model_errors_name]): raise ValueError("Model image contains too many frames")

        # Add model frame
        self.models[name] = frame

        # Add error map
        if errors is not None: self.model_errors[name] = errors

    # ------------------------------------------------------------------------------

    def add_errors(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Add
        self.errors[name] = frame

    # ------------------------------------------------------------------------------

    def add_model_errors(self, name, frame):

        """
        Thisn function ...
        :param name:
        :param frame:
        :return:
        """

        # Add
        self.model_errors[name] = frame

    # ------------------------------------------------------------------------------

    def add_residuals(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Add
        self.residuals[name] = frame

    # ------------------------------------------------------------------------------

    def add_distribution(self, name, distribution):

        """
        This function ...
        :param name:
        :param distribution:
        :return:
        """

        # Add
        self.distributions[name] = distribution

    # -----------------------------------------------------------------

    def add_settings(self, name, **settings):

        """
        This function ...
        :param name:
        :param settings:
        :return:
        """

        # Set settings
        self.settings[name].set_properties(settings)

    # ------------------------------------------------------------------------------

    def set_settings(self, name, settings):

        """
        This function ...
        :param name:
        :param settings:
        :return:
        """

        # Set settings
        self.settings[name] = settings

    # ------------------------------------------------------------------------------

    def set_setting(self, name, setting_name, value):
        
        """
        This function ...
        :param name: 
        :param setting_name: 
        :param value: 
        :return: 
        """

        # Set
        self.settings[name][setting_name] = value

    # ------------------------------------------------------------------------------

    def add_images(self, images):

        """
        This function ...
        :param images:
        :return:
        """

        # Debugging
        log.debug("Adding images ...")

        # Loop over the images
        for name in images:

            # Get the image
            image = images[name]

            # Add
            self.add_image(name, image)

    # ------------------------------------------------------------------------------

    def add_observations(self, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Debugging
        log.debug("Adding observations ...")

        # Loop over the frames
        for name in frames:

            # Get the frames
            frame = frames[name]

            # Add
            self.add_observation(name, frame)

    # ------------------------------------------------------------------------------

    def add_models(self, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Debugging
        log.debug("Adding models ...")

        # Loop over the frames
        for name in frames:

            # Get the frames
            frame = frames[name]

            # Add
            self.add_model(name, frame)

    # ------------------------------------------------------------------------------

    def add_error_maps(self, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Debugging
        log.debug("Adding error maps ...")

        # Loop over the frames
        for name in frames:

            # Get the frame
            frame = frames[name]

            # Add
            self.add_errors(name, frame)

    # ------------------------------------------------------------------------------

    def add_model_error_maps(self, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Debugging
        log.debug("Adding model error maps ...")

        # Loop over the frames
        for name in frames:

            # Get the frame
            frame = frames[name]

            # Add
            self.add_model_errors(name, frame)

    # ------------------------------------------------------------------------------

    def add_residual_maps(self, frames):

        """
        This function ...
        :param frames:
        :return:
        """

        # Debugging
        log.debug("Adding residual maps ...")

        # Loop over the frames
        for name in frames:

            # Get the frame
            frame = frames[name]

            # Add
            self.add_residuals(name, frame)

    # ------------------------------------------------------------------------------

    def load_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Are there FITS files in the directory?
        if fs.has_files_in_path(path, extension="fits"): self.load_images_from_directory(path)

        # Are there subdirectories?
        elif fs.has_directories_in_path(path):

            # Determine paths
            images_path = fs.join(path, images_name)
            observations_path = fs.join(path, observations_name)
            models_path = fs.join(path, models_name)
            residuals_path = fs.join(path, residuals_name)
            settings_path = fs.join(path, settings_name)

            # Load observations
            if fs.is_directory(images_path): self.load_images_from_directory(path)
            if fs.is_directory(observations_path): self.load_observations_from_directory(path)
            if fs.is_directory(models_path): self.load_models_from_directory(path)
            if fs.is_directory(residuals_path): self.load_residuals_from_directory(path)
            if fs.is_directory(settings_path): self.load_settings_from_directory(path)

        # No FITS files nor subdirectories
        else: raise IOError("No image files nor subdirectories found in '" + path + "'")

    # ------------------------------------------------------------------------------

    def load_images_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading image files from '" + path + "' ...")

        # Loop over the FITS files
        for name, filepath in fs.files_in_path(path, extension="fits", returns=["name", "path"]):

            # Debugging
            log.debug("Loading '" + name + "' image ...")

            # Load the image
            image = Image.from_file(filepath, always_call_first_primary=False)

            # Add the image
            self.add_image(name, image)

    # ------------------------------------------------------------------------------

    def load_observations_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading observed image frames from '" + path + "' ...")

        # Loop over the FITS files
        for name, filepath in fs.files_in_path(path, extension="fits", returns=["name", "path"]):

            # Debugging
            log.debug("Loading the '" + name + "' observed image ...")

            # Get header
            #header = get_header(filepath)

            # Get the filter
            #fltr = get_filter(name, header=header)

            # Check whether the filter is in the list of filters to be plotted
            #if fltr not in config.filters: continue

            # Get the index for this filter
            #index = config.filters.index(fltr)

            # Load the image
            #frame = Frame.from_file(filepath)
            image = Image.from_file(filepath, always_call_first_primary=False)

            # Replace zeroes and negatives
            image.primary.replace_zeroes_by_nans()
            image.primary.replace_negatives_by_nans()

            # Add the image
            self.add_observation(name, image)

    # ------------------------------------------------------------------------------

    def load_models_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading model image frames from '" + path + "' ...")

        # Loop over the FITS files
        for name, filepath in fs.files_in_path(path, extension="fits", returns=["name", "name"]):

            # Debugging
            log.debug("Loading the '" + name + "' model image ...")

            # Load the image
            image = Image.from_file(filepath, always_call_first_primary=False)

            # Replace zeroes and negatives
            image.primary.replace_zeroes_by_nans()
            image.primary.replace_negatives_by_nans()

            # Add the image
            self.add_model(name, image)

    # ------------------------------------------------------------------------------

    def load_residuals_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading residual image frames from '" + path + "' ...")

        # Loop over the FITS files
        for name, filepath in fs.files_in_path(path, extension="fits", returns=["name", "path"]):

            # Debugging
            log.debug("Loading the '" + name + "' residual map ...")

            # Load the frame
            frame = Frame.from_file(filepath)

            # Add the map
            self.add_residuals(name, frame)

    # ------------------------------------------------------------------------------

    def load_settings_from_directory(self, path):

        """
        This function ...
        :param path:
        :return:
        """

        # Debugging
        log.debug("Loading plotting settings from '" + path + "' ...")

        # Loop over the dat files
        for name, filepath in fs.files_in_path(path, extension="dat", returns=["name", "path"]):

            # Debugging
            log.debug("Loading the '" + name + "' settings ...")

            # Load the settings
            settings = ImagePlotSettings.from_file(filepath)

            # Set the settings
            self.set_settings(name, settings)

    # ------------------------------------------------------------------------------

    def get_observation_or_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        if self.has_observation(name): return self.get_observation(name)
        elif self.has_model(name): return self.get_model(name)
        else: raise ValueError("Doesn't have observation or model for name '" + name + "'")

    # ------------------------------------------------------------------------------

    def get_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.get_observation_or_model(name).filter

    # ------------------------------------------------------------------------------

    def get_wcs(self, name):

        """
        Thisf unction ...
        :param name:
        :return:
        """

        return self.get_observation_or_model(name).wcs

    # ------------------------------------------------------------------------------

    def calculate_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Get the frames
        #observation = self.observations[name]
        #model = self.models[name]

        # Uniformize
        observation, model = uniformize(self.observations[name], self.models[name])

        # Error-weighed residuals
        if self.config.weighed:

            if self.config.weighing_reference == observation_name:
                if not self.has_errors(name): raise ValueError("No errors for the '" + name + "' image")
                errors = self.get_errors(name)
            elif self.config.weighing_reference == model_name:
                if not self.has_model_errors(name): raise ValueError("No model errors for the '" + name + "' image")
                errors = self.get_model_errors(name)
            else: raise ValueError("Invalid value for 'weighing_reference'")

            # Calculate
            res = Frame((model - observation) / errors, wcs=observation.wcs)

        # Relative residuals
        elif self.config.relative: res = Frame((model - observation) / observation, wcs=observation.wcs)

        # Absolute residuals
        else: res = Frame(model - observation, wcs=observation.wcs)

        # Take absolute values?
        if self.config.absolute: res = res.absolute

        # Return the residual
        return res

    # ------------------------------------------------------------------------------

    def create_residuals(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Inform the user
        log.info("Creating the residual frames ...")

        # Loop over the observed images
        for name in self.names:

            # Checks
            if not self.has_model(name): continue
            if self.has_residuals(name): continue

            # Debugging
            log.debug("Creating residual frame for the '" + name + "' image ...")

            # Create
            res = self.calculate_residuals(name)

            # Add the residuals frame
            self.residuals[name] = res

    # ------------------------------------------------------------------------------

    def create_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the residual distributions ...")

        # Loop over the residual maps
        for name in self.residuals_names:

            # Checks
            if self.has_distribution(name): continue

            # Debugging
            log.debug("Creating distribution for the '" + name + "' residuals ...")

            # Get the residual map
            residuals = self.get_residuals(name)

            # Create the distribution
            distribution = Distribution.from_data("Residual", residuals, sigma_clip=self.config.sigma_clip_distributions, sigma_level=self.config.sigma_clip_level)

            # Add the distribution
            self.distributions[name] = distribution

    # ------------------------------------------------------------------------------

    def get_observation(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.observations[name]

    # ------------------------------------------------------------------------------

    @memoize_method
    def get_observation_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Create image
        image = Image(name=name)

        # Add observation frame
        image.add_frame(self.get_observation(name), observation_name)

        # Add error map
        if self.has_errors(name): image.add_frame(self.get_errors(name), errors_name)

        # Return the image
        return image

    # ------------------------------------------------------------------------------

    def get_model(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.models[name]

    # ------------------------------------------------------------------------------

    @memoize_method
    def get_model_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Create image
        image = Image(name=name)

        # Add model frame
        image.add_frame(self.get_model(name), model_name)

        # Add error map
        if self.has_model_errors(name): image.add_frame(self.get_model_errors(name), errors_name)

        # Return the image
        return image

    # ------------------------------------------------------------------------------

    def get_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.errors[name]

    # ------------------------------------------------------------------------------

    def get_model_errors(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.model_errors[name]

    # ------------------------------------------------------------------------------

    def get_residuals(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.residuals[name]

    # ------------------------------------------------------------------------------

    def get_distribution(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.distributions[name]

    # ------------------------------------------------------------------------------

    @memoize_method
    def get_image(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Create the image
        image = Image(name=name)

        # Add the observation
        if self.has_observation(name): image.add_frame(self.get_observation(name), observation_name)

        # Add the model
        if self.has_model(name): image.add_frame(self.get_model(name), model_name)

        # Add the errors
        if self.has_errors(name): image.add_frame(self.get_errors(name), errors_name)

        # Add the model errors
        if self.has_model_errors(name): image.add_frame(self.get_model_errors(name), model_errors_name)

        # Add the residuals
        if self.has_residuals(name): image.add_frame(self.get_residuals(name), residuals_name)

        # Return the image
        return image

    # ------------------------------------------------------------------------------

    def get_settings(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.settings[name]

    # ------------------------------------------------------------------------------

    def show(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Showing ...")

    # ------------------------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write observations
        if self.config.write_observations: self.write_observations()

        # Write models
        if self.config.write_models: self.write_models()

        # Write residual frames
        if self.config.write_residuals: self.write_residuals()

        # Write the images
        if self.config.write_images: self.write_images()

        # Write the distributions
        if self.config.write_distributions: self.write_distributions()

        # Write the settings
        if self.config.write_settings: self.write_settings()

    # ------------------------------------------------------------------------------

    @lazyproperty
    def images_path(self):
        return self.output_path_directory(images_name)

    # ------------------------------------------------------------------------------

    @lazyproperty
    def observations_path(self):
        return self.output_path_directory(observations_name)

    # ------------------------------------------------------------------------------

    @lazyproperty
    def models_path(self):
        return self.output_path_directory(models_name)

    # ------------------------------------------------------------------------------

    @lazyproperty
    def residuals_path(self):
        return self.output_path_directory(residuals_name)

    # ------------------------------------------------------------------------------

    @lazyproperty
    def distributions_path(self):
        return self.output_path_directory(distributions_name)

    # ------------------------------------------------------------------------------

    @lazyproperty
    def settings_path(self):
        return self.output_path_directory(settings_name)

    # ------------------------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Loop over all images
        for name in self.all_names:

            # Determine path
            path = fs.join(self.images_path, name + ".fits")

            # Debugging
            log.debug("Writing the '" + name + "' image ...")

            # Get image
            image = self.get_image(name)

            # Save the image
            image.saveto(path)

    # ------------------------------------------------------------------------------

    def write_observations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the observed frames ...")

        # Loop over the observed images
        for name in self.observation_names:

            # Determine the path
            path = fs.join(self.observations_path, name + ".fits")

            # Debugging
            log.debug("Writing the '" + name + "' observed image ...")

            # Get the frame
            frame = self.get_observation_image(name)

            # Save the frame
            frame.saveto(path)

    # ------------------------------------------------------------------------------

    def write_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the model frames ...")

        # Loop over the model images
        for name in self.model_names:

            # Determine the path
            path = fs.join(self.models_path, name + ".fits")

            # Debugging
            log.debug("Writing the '" + name + "' model image ...")

            # Get the frame
            frame = self.get_model_image(name)

            # Save the frame
            frame.saveto(path)

    # ------------------------------------------------------------------------------

    def write_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual frames ...")

        # Loop over the residual maps
        for name in self.residuals_names:

            # Determine the path
            path = fs.join(self.residuals_path, name + ".fits")

            # Debugging
            log.debug("Writing the '" + name + "' residual frame ...")

            # Get the residual map
            frame = self.get_residuals(name)

            # Save the frame
            frame.saveto(path)

    # ------------------------------------------------------------------------------

    def write_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the residual distributions ...")

        # Loop over the distributions
        for name in self.distribution_names:

            # Determine the path
            path = fs.join(self.distributions_path, name + ".fits")

            # Debugging
            log.debug("Writing the '" + name + "' residual distribution ...")

            # Get the distribution
            distribution = self.get_distribution(name)

            # Save
            distribution.saveto(path)

    # ------------------------------------------------------------------------------

    def write_settings(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the plotting settings ...")

        # Loop over the settings
        for name in self.settings_names:

            # Determine the path
            path = fs.join(self.settings_path, name + ".dat")

            # Debugging
            log.debug("Writing the '" + name + "' plotting settings ...")

            # Get the settings
            settings = self.get_settings(name)

            # Save
            settings.saveto(path)

    # ------------------------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot observations
        self.plot_observations()

        # Plot models
        self.plot_models()

        # Plot residuals
        self.plot_residuals()

        # Plot distributions
        if self.config.distributions: self.plot_distributions()

        # Finish the plot
        self.finish()

    # ------------------------------------------------------------------------------

    def get_label(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # No settings?
        if not self.has_settings(name): return name

        # Get the settings
        settings = self.get_settings(name)

        # Return
        if settings.label is not None: return settings.label
        else: return name

    # ------------------------------------------------------------------------------

    def get_colormap(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # No settings?
        if not self.has_settings(name): return self.config.cmap

        # Get the settings
        settings = self.get_settings(name)

        # Return
        if settings.cmap is not None: return settings.cmap
        else: return self.config.cmap

    # ------------------------------------------------------------------------------

    @property
    def config_residual_cmap(self):

        """
        This function ...
        :return:
        """

        if self.config.absolute: return self.config.absolute_residual_cmap
        else: return self.config.residual_cmap

    # ------------------------------------------------------------------------------

    def get_residual_colormap(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # No settings
        if not self.has_settings(name): return self.config_residual_cmap

        # Get the settings
        settings = self.get_settings(name)

        # Return
        if settings.residual_cmap is not None: return settings.residual_cmap
        else: return self.config_residual_cmap

    # ------------------------------------------------------------------------------

    def get_limits(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # No settings
        if not self.has_settings(name): return self.config.vmin, self.config.vmax, False, False

        # Get the settings
        settings = self.get_settings(name)

        # Get limits
        vmin = settings.vmin if settings.vmin is not None else self.config.vmin
        vmax = settings.vmax if settings.vmax is not None else self.config.vmax

        # Get flags
        soft_vmin = settings.soft_vmin if settings.vmin is not None else False # don't use True flag if vmin is not set in settings
        soft_vmax = settings.soft_vmax if settings.vmax is not None else False # don't use True flag if vmax is not set in settings

        # Return
        return vmin, vmax, soft_vmin, soft_vmax

    # ------------------------------------------------------------------------------

    def get_residual_amplitude(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # No settings
        if not self.has_settings(name): return self.config.residual_amplitude, False

        # Get the settings
        settings = self.get_settings(name)

        # Get amplitude
        amplitude = settings.residual_amplitude if settings.residual_amplitude is not None else self.config.residual_amplitude

        # Get flag
        soft_amplitude = settings.soft_residual_amplitude if settings.residual_amplitude is not None else False # don't use True flag if amplitude is not set in settings

        # Return
        return amplitude, soft_amplitude

    # ------------------------------------------------------------------------------

    def set_limits(self, name, vmin, vmax, soft_vmin=None, soft_vmax=None):

        """
        This function ...
        :param name:
        :param vmin:
        :param vmax:
        :param soft_vmin:
        :param soft_vmax:
        :return:
        """

        # Set vmin and vmax
        self.add_settings(name, vmin=vmin, vmax=vmax)

        # Set flags
        if soft_vmin is not None: self.set_setting(name, "soft_vmin", soft_vmin)
        if soft_vmax is not None: self.set_setting(name, "soft_vmax", soft_vmax)

    # ------------------------------------------------------------------------------

    def get_vmin_vmax(self, frame, vmin=None, vmax=None, soft_vmin=False, soft_vmax=False):

        """
        This function ...
        :param frame:
        :param vmin:
        :param vmax:
        :param soft_vmin:
        :param soft_vmax:
        :return:
        """

        # Defined?
        has_vmin = vmin is not None
        has_vmax = vmax is not None

        # Vmin and vmax don't have to be calculated
        if has_vmin and has_vmax and (not soft_vmin) and (not soft_vmax): return vmin, vmax

        # Calculate vmin and or vmax
        return get_vmin_vmax(frame.data, interval=self.config.interval, zmin=vmin, zmax=vmax, soft_zmin=soft_vmin, soft_zmax=soft_vmax)

    # ------------------------------------------------------------------------------

    def get_residual_vmin_vmax(self, frame, amplitude=None, soft_amplitude=False):

        """
        This function ...
        :param frame:
        :param amplitude:
        :param soft_amplitude:
        :return:
        """

        # Defined?
        if amplitude is not None and not soft_amplitude:

            if self.config.absolute: return 0., amplitude
            else: return -amplitude, amplitude

        # Calculate vmin and or vmax
        if self.config.absolute: return get_vmin_vmax(frame.data, interval=self.config.residual_interval, zmin=0, zmax=amplitude, soft_zmin=False, soft_zmax=soft_amplitude)
        else:
            zmin = -amplitude if amplitude is not None else None
            zmax = amplitude
            return get_vmin_vmax(frame.data, interval=self.config.residual_interval, zmin=zmin, zmax=zmax, soft_zmin=soft_amplitude, soft_zmax=soft_amplitude, around_zero=True, symmetric=True)

    # ------------------------------------------------------------------------------

    def get_observation_row_col(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Horizontal
        #if self.horizontal: return index, 0
        if self.horizontal: return 0, index

        # Vertical
        #elif self.vertical: return 0, index
        elif self.vertical: return index, 0

        # Invalid
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    def get_model_row_col(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Horizontal
        #if self.horizontal: return index, 1
        if self.horizontal: return 1, index

        # Vertical
        #elif self.vertical: return 1, index
        elif self.vertical: return index, 1

        # Invalid
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    def get_residuals_row_col(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Horizontal
        #if self.horizontal: return index, 2
        if self.horizontal: return 2, index

        # Vertical
        #elif self.vertical: return 2, index
        elif self.vertical: return index, 2

        # Invalid
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    def get_distribution_row_col(self, index):

        """
        This function ...
        :param index:
        :return:
        """

        # Horizontal
        #if self.horizontal: return index, 3
        if self.horizontal: return 3, index

        # Vertical
        #elif self.vertical: return 3, index
        elif self.vertical: return index, 3

        # Invalid
        else: raise ValueError("Invalid direction")

    # ------------------------------------------------------------------------------

    def get_observation_spec(self, index, return_row_col=False):

        """
        This function ...
        :param index:
        :param return_row_col:
        :return:
        """

        # Get row and col
        row, col = self.get_observation_row_col(index)

        #print(self.grid.get_geometry())
        #print(self.grid.get_height_ratios())

        # Return the grid spec
        #if return_row_col: return self.grid[row, col], row, col
        #else: return self.grid[row, col]
        #if return_row_col: return self.grid[index], row, col
        #else: return self.grid[index]
        # No, no, this was a mistake with 'get_observation_row_col'
        #if return_row_col: return self.grid[col, row], row, col # WHY?
        #else: return self.grid[col, row] # WHY?

        # This was right after all
        if return_row_col: return self.grid[row, col], row, col
        else: return self.grid[row, col]

    # ------------------------------------------------------------------------------

    def get_model_spec(self, index, return_row_col=False):

        """
        This function ...
        :param index:
        :param return_row_col:
        :return:
        """

        # Get row and col
        row, col = self.get_model_row_col(index)

        # Return the grid spec
        if return_row_col: return self.grid[row, col], row, col
        else: return self.grid[row, col]

    # ------------------------------------------------------------------------------

    def get_residuals_spec(self, index, return_row_col=False):

        """
        This function ...
        :param index:
        :param return_row_col:
        :return:
        """

        # Get row and col
        row, col = self.get_residuals_row_col(index)

        # Return the grid spec
        if return_row_col: return self.grid[row, col], row, col
        else: return self.grid[row, col]

    # ------------------------------------------------------------------------------

    def get_distribution_spec(self, index, return_row_col=False):

        """
        This function ...
        :param index:
        :param return_row_col:
        :return:
        """

        # Get row and col
        row, col = self.get_distribution_row_col(index)

        # Return the grid spec
        if return_row_col: return self.grid[row, col], row, col
        else: return self.grid[row, col]

    # ------------------------------------------------------------------------------

    def create_observation_plot(self, index, frame):

        """
        This function ...
        :param index:
        :param frame:
        :return:
        """

        # Get the subplot spec
        spec, row, col = self.get_observation_spec(index, return_row_col=True)
        #print(spec)
        #print("ROW", row, "COL", col)

        # Get coordinates of the subplot
        #points = spec.get_position(self.figure.figure).get_points()
        bbox = spec.get_position(self.figure.figure)
        coordinates = [bbox.x0, bbox.y0, bbox.width, bbox.height]

        # Create the plot
        # needs [xmin, ymin, dx, dy]
        plot = aplpy.FITSFigure(frame.to_hdu(), figure=self.figure.figure, subplot=coordinates)

        # Add the plot
        self.plots[row][col] = plot

        # Return the plot
        return plot

    # ------------------------------------------------------------------------------

    def create_model_plot(self, index, frame):

        """
        This function ...
        :param index:
        :param frame:
        :return:
        """

        # Get the subplot spec
        spec, row, col = self.get_model_spec(index, return_row_col=True)

        bbox = spec.get_position(self.figure.figure)
        coordinates = [bbox.x0, bbox.y0, bbox.width, bbox.height]

        # Create the plot
        plot = aplpy.FITSFigure(frame.to_hdu(), figure=self.figure.figure, subplot=coordinates)

        # Add the plot
        self.plots[row][col] = plot

        # Return the plot
        return plot

    # ------------------------------------------------------------------------------

    def create_residuals_plot(self, index, frame):

        """
        This function ...
        :param index:
        :param frame:
        :return:
        """

        # Get the subplot spec
        spec, row, col = self.get_residuals_spec(index, return_row_col=True)

        bbox = spec.get_position(self.figure.figure)
        coordinates = [bbox.x0, bbox.y0, bbox.width, bbox.height]

        # Create the plot
        plot = aplpy.FITSFigure(frame.to_hdu(), figure=self.figure.figure, subplot=coordinates)

        # Add the plot
        self.plots[row][col] = plot

        # Return the plot
        return plot

    # ------------------------------------------------------------------------------

    def _plot_observation(self, index, frame, cmap, label=None, vmin=None, vmax=None, soft_vmin=False, soft_vmax=False):

        """
        This function ...
        :param index:
        :param frame:
        :param cmap:
        :param label:
        :param vmin:
        :param vmax:
        :param soft_vmin:
        :param soft_vmax:
        :return:
        """

        # Create the plot
        plot = self.create_observation_plot(index, frame)

        # Get vmin and vmax
        vmin, vmax = self.get_vmin_vmax(frame, vmin=vmin, vmax=vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax)

        # Set colorscale
        plot.show_colorscale(vmin=vmin, vmax=vmax, cmap=cmap, stretch=self.config.scale)

        # Set tick label font
        plot.tick_labels.set_font(size='small')

        # Set center, radius and spacing
        plot.recenter(self.ra_center_deg, self.dec_center_deg, radius=self.radius_deg)
        plot.ticks.set_xspacing(self.spacing_deg)

        # Set color or frame
        plot.frame.set_color(self.frame_color)

        # FOR FIRST
        #f1._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
        #f1._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

        # Tick settings
        plot._ax2.tick_params(direction='in', which='major', length=self.config.major_tick_length, top=True, right=True, bottom=True, left=True)
        plot._ax2.tick_params(direction='in', which='minor', length=self.config.minor_tick_length, top=True, right=True, bottom=True, left=True)

        # Set image background color
        plot.set_nan_color(self.nan_color)

        # FOR FIRST
        #f1._ax1.scatter(ra, dec, marker='.', label='Observation')

        # FOR FIRST
        #legend1 = f1._ax1.legend(loc='upper right', fontsize=12, fancybox=True, framealpha=0, numpoints=None)
        #plt.setp(legend1.get_texts(), color=config.text_color_in)

        # Set title
        if label is not None: plot._ax1.set_title(label, fontsize=self.config.label_fontsize)

        # Return the vmin and vmax
        return vmin, vmax

    # ------------------------------------------------------------------------------

    def _plot_model(self, index, frame, cmap, vmin=None, vmax=None, soft_vmin=None, soft_vmax=None):

        """
        This function ...
        :param index:
        :param frame:
        :param vmin:
        :param vmax:
        :param soft_vmin:
        :param soft_vmax:
        :return:
        """

        # Create the plot
        plot = self.create_model_plot(index, frame)

        # Get vmin and vmax
        vmin, vmax = self.get_vmin_vmax(frame, vmin=vmin, vmax=vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax)

        # Set colorscale
        plot.show_colorscale(vmin=vmin, vmax=vmax, cmap=cmap, stretch=self.config.scale)

        # Set tick label font
        plot.tick_labels.set_font(size='small')

        # Set center, radius and spacing
        plot.recenter(self.ra_center_deg, self.dec_center_deg, radius=self.radius_deg)
        plot.ticks.set_xspacing(self.spacing_deg)

        # Set color for frame
        plot.frame.set_color(self.frame_color)

        # Set ticks
        plot._ax1.tick_params(direction='in', which='major', length=self.config.major_tick_length, top=True, right=True, bottom=True, left=True)
        plot._ax1.tick_params(direction='in', which='minor', length=self.config.minor_tick_length, top=True, right=True, bottom=True, left=True)

        # FOR FIRST
        #f6._ax1.scatter(ra, dec, marker='.', label='Model')
        #legend6 = f6._ax1.legend(loc='upper right', fontsize=12, fancybox=False, framealpha=0, numpoints=None)
        #plt.setp(legend6.get_texts(), color=config.text_color_in)

        # Set image background color
        plot.set_nan_color(self.nan_color)

    # ------------------------------------------------------------------------------

    def _plot_residuals(self, index, frame, cmap, amplitude=None, soft_amplitude=False):

        """
        This function ...
        :param index:
        :param frame:
        :param cmap:
        :param amplitude:
        :param soft_amplitude:
        :return:
        """

        # Create the plot
        plot = self.create_residuals_plot(index, frame)

        # Get vmin and vmax
        vmin, vmax = self.get_residual_vmin_vmax(frame, amplitude=amplitude, soft_amplitude=soft_amplitude)

        # Set colorscale
        plot.show_colorscale(vmin=vmin, vmax=vmax, cmap=cmap)

        # Set tick label font
        plot.tick_labels.set_font(size='small')

        # Set center, radius and spacing
        plot.recenter(self.ra_center_deg, self.dec_center_deg, radius=self.radius_deg)
        plot.ticks.set_xspacing(self.spacing_deg)

        # Set color for frame
        plot.frame.set_color(self.frame_color)

        # Set ticks
        plot._ax1.tick_params(direction='in', which='major', length=self.config.major_tick_length, top=True, right=True, bottom=True, left=True)
        plot._ax1.tick_params(direction='in', which='minor', length=self.config.minor_tick_length, top=True, right=True, bottom=True, left=True)

        # FOR FIRST
        # f11._ax1.scatter(ra, dec, marker='.', label='Relative \nResidual')

        # FOR FIRST
        # Set legend
        #legend11 = f11._ax1.legend(loc='lower right', fontsize=12, fancybox=False, framealpha=0, numpoints=None)
        #plt.setp(legend11.get_texts(), color=config.text_color_in)

        # Set background color
        plot.set_nan_color(self.background_color)

    # ------------------------------------------------------------------------------

    def _plot_distribution(self, index, distribution):

        """
        This function ...
        :param index:
        :param distribution:
        :return:
        """

        pass

    # ------------------------------------------------------------------------------

    def plot_observations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the observed image frames ...")

        # Loop over the names
        #print(self.names)
        #print(self.nimages)
        #print(len(self.names))
        for index, name in enumerate(self.names):

            # Debugging
            log.debug("Plotting the observed frame of the '" + name + "' image (panel " + str(index+1) + " of " + str(self.nimages) + ") ...")

            # Get the observation
            frame = self.get_observation(name)

            # Get the label for this image
            label = self.get_label(name)

            # Get the colormap for this image
            cmap = self.get_colormap(name)

            # Get the limits
            vmin, vmax, soft_vmin, soft_vmax = self.get_limits(name)

            # Plot
            vmin, vmax = self._plot_observation(index, frame, cmap, label=label, vmin=vmin, vmax=vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax)

            # Set new vmin and vmax (for corresponding model)
            self.set_limits(name, vmin, vmax, soft_vmin=False, soft_vmax=False)

    # ------------------------------------------------------------------------------

    def plot_models(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the model image frames ...")

        # Loop over the names
        for index, name in enumerate(self.names):

            # Check
            if not self.has_model(name): continue

            # Debugging
            log.debug("Plotting the model frame of the '" + name + "' image (panel " + str(index+1) + " of " + str(self.nimages) + ") ...")

            # Get the model
            frame = self.get_model(name)

            # Get the colormap for this image
            cmap = self.get_colormap(name)

            # Get the limits
            vmin, vmax, soft_vmin, soft_vmax = self.get_limits(name)

            # Plot
            self._plot_model(index, frame, cmap, vmin=vmin, vmax=vmax, soft_vmin=soft_vmin, soft_vmax=soft_vmax)

    # ------------------------------------------------------------------------------

    def plot_residuals(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residual image frames ...")

        # Loop over the names
        for index, name in enumerate(self.names):

            # Check
            if not self.has_residuals(name): continue

            # Debugging
            log.debug("Plotting the residuals frame of the '" + name + "' image (panel " + str(index+1) + " of " + str(self.nimages) + ") ...")

            # Get the residuals
            frame = self.get_residuals(name)

            # Get the colormap for this residual map
            cmap = self.get_residual_colormap(name)

            # Get the amplitude
            amplitude, soft_amplitude = self.get_residual_amplitude(name)

            # Plot
            # index, frame, cmap, amplitude=None, soft_amplitude=False
            self._plot_residuals(index, frame, cmap, amplitude=amplitude, soft_amplitude=soft_amplitude)

    # ------------------------------------------------------------------------------

    def plot_distributions(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the residual distributions ...")

        # Loop over the names
        for index, name in enumerate(self.names):

            # Check
            if not self.has_distribution(name): continue

            # Debugging
            log.debug("Plotting the residual distribution of the '" + name + "' image (panel " + str(index+1) + " of " + str(self.nimages) + " ) ...")

            # Get the distribution
            distribution = self.get_distribution(name)

    # ------------------------------------------------------------------------------

    def finish(self):

        """
        This function ...
        :param self:
        :return:
        """

        # Draw
        self.figure.draw()

        # Save to file
        if self.config.path is not None: self.figure.figure.savefig(self.config.path, dpi=self.config.dpi)

        # Show
        else: plt.show()

        # Close
        #plt.close(fig)
        plt.close()

# ------------------------------------------------------------------------------

def plot_images_aplpy(frames, filepath=None, center=None, radius=None, xy_ratio=None, dark=False, scale="log",
                      colormap="inferno", nrows=None, ncols=None, orientation="horizontal", plotsize=3., distance=None,
                      share_scale=None, descriptions=None, minmax_scaling=0.5):

    """
    This function ...
    :param frames:
    :param filepath:
    :param center:
    :param radius:
    :param xy_ratio:
    :param dark:
    :param scale:
    :param colormap:
    :param nrows:
    :param ncols:
    :param orientation:
    :param plotsize:
    :param distance:
    :param share_scale:
    :param descriptions:
    :param minmax_scaling: 0.5
    :return:
    """

    import matplotlib.gridspec as gridspec
    #from matplotlib.colorbar import ColorbarBase
    #from matplotlib.colors import LinearSegmentedColormap
    #from matplotlib.colors import Normalize

    from pts.magic.tools import plotting

    # Set
    set_theme(dark=dark)
    nimages = len(frames)

    xsize = plotsize
    #if xy_ratio is None: ysize = 3.5
    #else: ysize = xsize / xy_ratio
    if xy_ratio is None: xy_ratio = 0.85
    ysize = xsize / xy_ratio

    #print("plotsize", xsize, ysize)

    # Determine the number of columns and rows
    if nrows is None and ncols is None:

        if orientation == "horizontal": ncols, nrows = nimages, 1
        elif orientation == "vertical": ncols, nrows = 1, nimages
        else: raise ValueError("Invalid orientation: '" + orientation + "'")

    # Nrows is none but ncols is not
    elif nrows is None: ncols = numbers.round_up_to_int(nimages/nrows)

    # Ncols is none but nrows is not
    elif ncols is None: nrows = numbers.round_up_to_int(nimages/ncols)

    # Set figure size
    figxsize = xsize * ncols
    figysize = ysize * nrows

    #print("figsize", figxsize, figysize)

    # Create figure with appropriate size
    fig = plt.figure(figsize=(figxsize, figysize))

    # Create grid
    gs1 = gridspec.GridSpec(nrows, ncols)  # nimages ROWS, 4 COLUMNS
    # gs1.update(wspace=0.01, hspace=0.3)
    gs1.update(wspace=0., hspace=0.)
    plot_idx = 0

    # Get frame labels
    if types.is_dictionary(frames):
        labels = frames.keys()
        frames = frames.values()
    else: labels = [frame.filter_name for frame in frames]

    # Set scale for each image
    scales = dict()
    if types.is_string_type(scale):
        for label in labels: scales[label] = scale
    elif types.is_sequence(scale):
        for label, scalei in zip(labels, scale): scales[label] = scalei
    elif types.is_dictionary(scale): scales = scale
    else: raise ValueError("Invalid type for 'scale'")

    # Initialize dict for intervals
    intervals = dict()

    # Set descriptions
    if descriptions is None:
        descriptions = dict()
        for label in labels: descriptions[label] = None
    elif types.is_sequence(descriptions):
        descrpts = descriptions
        descriptions = dict()
        for label, descr in zip(labels, descrpts): descriptions[label] = descr
    elif types.is_dictionary(descriptions): pass # OK
    else: raise ValueError("Invalid type for 'descriptions'")

    # Set minmax scaling
    if types.is_real_type(minmax_scaling):
        factor = minmax_scaling
        minmax_scaling = dict()
        for label in labels: minmax_scaling[label] = factor
    elif types.is_dictionary(minmax_scaling):
        minmax_scaling_orig = minmax_scaling
        minmax_scaling = dict()
        for label in labels:
            if label in minmax_scaling_orig: minmax_scaling[label] = minmax_scaling_orig[label]
            else: minmax_scaling[label] = 0.5
    elif types.is_sequence(minmax_scaling):
        minmax_scaling_orig = minmax_scaling
        minmax_scaling = dict()
        for label, factor in zip(labels, minmax_scaling_orig): minmax_scaling[label] = factor
    else: raise ValueError("Invalid type for 'minmax_scaling'")

    # Loop over the frames
    for label, frame, index in zip(labels, frames, range(nimages)):

        rowi = index // ncols
        coli = index % ncols

        is_first_row = rowi == 0
        is_last_row = rowi == nrows - 1

        is_first_col = coli == 0
        is_last_col = coli == ncols - 1

        #print("row", rowi)
        #print("col", coli)

        # IS FIRST OR LAST IMAGE?
        is_first = index == 0
        is_last = index == nimages - 1

        # Debugging
        log.debug("Plotting the '" + label + "' image ...")

        # Get HDU
        hdu = frame.to_hdu()

        # Get interval
        if share_scale is not None and label in share_scale:

            share_with = share_scale[label]
            vmin, vmax = intervals[share_with]
            scalei = scales[share_with]

        else:

            # Get scale
            scalei = scales[label]
            is_logscale = scalei == "log"

            #print(label, minmax_scaling[label])
            vmin, vmax = plotting.get_vmin_vmax(frame.data, logscale=is_logscale, minmax_scaling=minmax_scaling[label])

        # Set interval
        intervals[label] = (vmin, vmax,)

        # Set title
        if descriptions[label] is not None: title = descriptions[label]
        else: title = label.replace("_", "\_").replace("um", "$\mu$m")

        # Has sky coordinate system?
        has_wcs = frame.has_wcs and frame.wcs.is_sky

        # OBSERVATION
        figi = aplpy.FITSFigure(hdu, figure=fig, subplot=list(gs1[plot_idx].get_position(fig).bounds))
        setup_map_plot(figi, colormap, vmin=vmin, vmax=vmax, label=r'' + str(title), center=center, radius=radius, scale=scalei, has_wcs=has_wcs)
        set_ticks(figi, is_first_row, is_last_row)

        # FIRST COLUMN
        if is_first_col:

            figi.tick_labels.show_y()
            figi.axis_labels.show_y()

        # LAST ROW
        if is_last_row:

            figi.tick_labels.show_x()
            figi.axis_labels.show_x()

        # Increment
        plot_idx += 1

    # Save the figure
    if filepath is not None: plt.savefig(filepath, bbox_inches='tight', dpi=300)
    else: plt.show()

    # Close
    plt.close()

    # Reset
    reset_theme()

# ------------------------------------------------------------------------------

def plot_one_residual_aplpy(observation, model, residual=None, path=None, scale="log", plotsize=3., dark=False,
                            center=None, radius=None, xy_ratio=None, first_label="Observation", second_label="Model",
                            residual_label="Residual", filter_label=True):

    """
    This function ...
    :param observation:
    :param model:
    :param residual:
    :param path:
    :param scale:
    :param plotsize:
    :param dark:
    :param center:
    :param radius:
    :param xy_ratio:
    :param first_label:
    :param second_label:
    :param residual_label:
    :param filter_label:
    :return:
    """

    # Make residual?
    if residual is None: residual = (model - observation) / observation

    # Colormaps
    colormap = "inferno"
    residual_colormap = "RdBu"

    import matplotlib.gridspec as gridspec
    from pts.magic.tools import plotting

    # Set theme
    set_theme(dark=dark)

    nrows = 1
    ncols = 3

    xsize = plotsize
    if xy_ratio is None: xy_ratio = 0.85
    ysize = xsize / xy_ratio

    # Set figure size
    figxsize = xsize * ncols
    figysize = ysize * nrows

    # Create figure with appropriate size
    #fig = plt.figure(figsize=(figxsize, figysize))
    figure = MPLFigure(size=(figxsize,figysize))

    # Create grid
    gs1 = gridspec.GridSpec(1, 4)  # nimages ROWS, 4 COLUMNS
    gs1.update(wspace=0., hspace=0.)
    plot_idx = 0

    # Percentual residuals
    residual = residual * 100.

    # Set title
    if filter_label and observation.has_filter: title = str(observation.filter).replace("um", " $\mu$m")
    else: title = first_label

    # Create HDU's for Aplpy
    observation_hdu = observation.to_hdu()
    model_hdu = model.to_hdu()
    residual_hdu = residual.to_hdu()

    # Get interval
    vmin, vmax = plotting.get_vmin_vmax(observation.data, logscale=scale=="log")

    # OBSERVATION
    fig1 = aplpy.FITSFigure(observation_hdu, figure=figure.figure, subplot=list(gs1[plot_idx].get_position(figure.figure).bounds))
    setup_map_plot(fig1, colormap, vmin=vmin, vmax=vmax, label=r'' + str(title), center=center, radius=radius, scale=scale, has_wcs=observation.has_celestial_wcs)
    set_ticks(fig1, True, True)

    # Enable y ticks and axis labels BECAUSE OBSERVATION IS THE FIRST COLUMN
    fig1.tick_labels.show_y()
    fig1.axis_labels.show_y()

    # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
    fig1.tick_labels.show_x()
    fig1.axis_labels.show_x()

    # Increment
    plot_idx += 1

    # MODEL
    fig2 = aplpy.FITSFigure(model_hdu, figure=figure.figure, subplot=list(gs1[plot_idx].get_position(figure.figure).bounds))
    setup_map_plot(fig2, colormap, vmin=vmin, vmax=vmax, label=second_label, center=center, radius=radius, scale=scale, has_wcs=model.has_celestial_wcs)
    set_ticks(fig2, True, True)

    # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
    fig2.tick_labels.show_x()
    fig2.axis_labels.show_x()

    # Increment
    plot_idx += 1

    # RESIDUAL
    fig3 = aplpy.FITSFigure(residual_hdu, figure=figure.figure, subplot=list(gs1[plot_idx].get_position(figure.figure).bounds))
    setup_map_plot(fig3, residual_colormap, vmin=-100, vmax=100, label=residual_label + ' (\%)', center=center, radius=radius, has_wcs=residual.has_celestial_wcs)
    set_ticks(fig3, True, True)

    # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
    fig3.tick_labels.show_x()
    fig3.axis_labels.show_x()

    # Show or save
    if path is None: figure.show()
    else: figure.saveto(path)

    # Reset theme
    reset_theme()

# ------------------------------------------------------------------------------

def plot_residuals_aplpy(observations, models, residuals, filepath=None, center=None, radius=None, xy_ratio=None,
                         dark=False, scale="log", plotsize=3., distance=None, mask_simulated=False, masks=None):

    """
    This function ...
    :param observations:
    :param models:
    :param residuals:
    :param filepath:
    :param center:
    :param radius:
    :param xy_ratio:
    :param dark:
    :param scale:
    :param plotsize:
    :param distance:
    :param mask_simulated:
    :param masks: if passed, both observations, models and residuals are masked
    :return:
    """

    import numpy as np
    import matplotlib.gridspec as gridspec
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import Normalize
    import seaborn as sns

    # Set theme
    set_theme(dark=dark)

    nimages = len(observations)
    ncols = 4
    nrows = nimages

    # Colormaps
    colormap = "inferno"
    residual_colormap = "RdBu"

    # Set individual map plot size
    xsize = plotsize
    #if xy_ratio is None: ysize = 3.5
    #else: ysize = xsize / xy_ratio
    #print("individual size", xsize, ysize)
    if xy_ratio is None: xy_ratio = 0.85
    ysize = xsize / xy_ratio

    # Set figure size
    figxsize = xsize * ncols
    figysize = ysize * nrows
    #print("figure size", figxsize, figysize)

    # Create figure with appropriate size
    fig = plt.figure(figsize=(figxsize, figysize))

    # Create grid
    gs1 = gridspec.GridSpec(nimages, 4) # nimages ROWS, 4 COLUMNS
    #gs1.update(wspace=0.01, hspace=0.3)
    gs1.update(wspace=0., hspace=0.)
    plot_idx = 0

    # Loop over the filters
    if masks is None: masks = [None] * nimages
    for observation, model, residual, mask, index in zip(observations, models, residuals, masks, range(nimages)):

        #print("units:")
        #print(observation.unit)
        #print(model.unit)
        observation.convert_to("mJy/sr", distance=distance)
        model.convert_to("mJy/sr", distance=distance)

        # MASK MODEL
        if mask_simulated:
            model.rebin(observation.wcs)
            model.apply_mask_nans(observation.nans)

        # MASK ALL?
        if mask is not None:
            observation.apply_mask_nans(mask)
            model.apply_mask_nans(mask)
            residual.apply_mask_nans(mask)

        # IS FIRST OR LAST IMAGE?
        is_first = index == 0
        is_last = index == nimages - 1

        # Debugging
        log.debug("Plotting the observation, model and residuals for the " + str(observation.filter) + " filter ...")

        # Percentual residuals
        residual = residual * 100.

        # Set title
        title = str(observation.filter).replace("um", " $\mu$m")

        # Create HDU's for Aplpy
        observation_hdu = observation.to_hdu()
        model_hdu = model.to_hdu()
        residual_hdu = residual.to_hdu()

        from pts.magic.tools import plotting

        vmin, vmax = plotting.get_vmin_vmax(observation.data, logscale=scale=="log")
        #vmax = 0.7 * vmax

        #print("VMIN", vmin)
        #print("VMAX", vmax)

        # ------------------------------------------------------------------------------
        # Plot obs, model and residual
        # ------------------------------------------------------------------------------

        # OBSERVATION
        fig1 = aplpy.FITSFigure(observation_hdu, figure=fig, subplot=list(gs1[plot_idx].get_position(fig).bounds))
        setup_map_plot(fig1, colormap, vmin=vmin, vmax=vmax, label=r'' + str(title), center=center, radius=radius, scale=scale)
        set_ticks(fig1, is_first, is_last)

        # Enable y ticks and axis labels BECAUSE OBSERVATION IS THE FIRST COLUMN
        fig1.tick_labels.show_y()
        fig1.axis_labels.show_y()

        # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
        if is_last: fig1.tick_labels.show_x()
        if is_last: fig1.axis_labels.show_x()

        # Increment
        plot_idx += 1

        # ------------------------------------------------------------------------------

        # MODEL
        fig2 = aplpy.FITSFigure(model_hdu, figure=fig, subplot=list(gs1[plot_idx].get_position(fig).bounds))
        setup_map_plot(fig2, colormap, vmin=vmin, vmax=vmax, label='Model', center=center, radius=radius, scale=scale)
        set_ticks(fig2, is_first, is_last)

        # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
        if is_last: fig2.tick_labels.show_x()
        if is_last: fig2.axis_labels.show_x()

        # Increment
        plot_idx += 1

        # ------------------------------------------------------------------------------

        # RESIDUAL
        fig3 = aplpy.FITSFigure(residual_hdu, figure=fig, subplot=list(gs1[plot_idx].get_position(fig).bounds))
        setup_map_plot(fig3, residual_colormap, vmin=-100, vmax=100, label='Residual (\%)', center=center, radius=radius)
        set_ticks(fig3, is_first, is_last)

        # SHOW THE X TICK LABELS AND AXIS LABELS ONLY IF LAST ROW
        if is_last: fig3.tick_labels.show_x()
        if is_last: fig3.axis_labels.show_x()

        # ------------------------------------------------------------------------------

        # COLORBAR
        colorbar_start_x = gs1[plot_idx].get_position(fig).bounds[0] + 0.025
        colorbar_start_y = gs1[plot_idx].get_position(fig).bounds[1] + 0.085 / (nimages)
        colorbar_x_width = gs1[plot_idx].get_position(fig).bounds[2] - 0.05
        colorbar_y_height = gs1[plot_idx].get_position(fig).bounds[3]
        cb_ax = fig.add_axes([colorbar_start_x, colorbar_start_y, colorbar_x_width, (0.02 + 0.002) / (nimages + 1)])

        # Colourbar
        cb = ColorbarBase(cb_ax, cmap=residual_colormap, norm=Normalize(vmin=-100, vmax=100), orientation='horizontal')
        cb.ax.xaxis.set_ticks_position('bottom')
        cb.ax.xaxis.set_label_position('bottom')
        cb.ax.zorder = 99
        cb.ax.xaxis.set_tick_params(color='white')
        cb.outline.set_edgecolor('white')
        plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='white')
        plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='white')
        cb.set_ticks([-100, -50, 0, 50, 100])

        # Increment
        plot_idx += 1

        # ------------------------------------------------------------------------------

        # KDE Plot of residuals

        residual = residual_hdu.data
        fig4 = plt.subplot(gs1[plot_idx])
        residuals_to_kde = np.where((residual <= 200) & (residual >= -200))

        if dark:
            sns.kdeplot(residual[residuals_to_kde], bw='silverman', c='white', shade=True)
            fig4.axes.set_facecolor("black")
        else:
            sns.kdeplot(residual[residuals_to_kde], bw='silverman', c='k', shade=True)
            fig4.axes.set_facecolor("white")

        fig4.tick_params(labelleft='off')
        plt.xlim([-150, 150])

        fig4.tick_params(direction='in', which='major', length=7, top=True, right=False, bottom=True, left=False)
        fig4.tick_params(direction='in', which='major', length=7, top=True, right=False, bottom=True, left=False)

        # Hide tick labels except for the last (bottom) plot
        if not is_last: fig4.tick_params(labelbottom=False)

        if dark: plt.axvline(0, c='white', ls='--', lw=2)
        else: plt.axvline(0, c='k', ls='--', lw=2)

        # Label for kde
        plt.xlabel('Residual (\%)')

        # Increment
        plot_idx += 1

    # Save the figure
    if filepath is not None: plt.savefig(filepath, bbox_inches='tight', dpi=300)
    else: plt.show()

    # Close
    plt.close()

    # Reset theme
    reset_theme()

# ------------------------------------------------------------------------------

def setup_map_plot(figure, colormap, vmin, vmax, label, smooth=None,text_x=0.05, text_y=0.95, center=None,
                   radius=None, scale="linear", has_wcs=True):

    """
    This function ...
    :param figure:
    :param colormap:
    :param vmin:
    :param vmax:
    :param label:
    :param smooth:
    :param text_x:
    :param text_y:
    :param center:
    :param radius:
    :param scale:
    :param has_wcs:
    :return:
    """

    figure.show_colorscale(cmap=colormap, vmin=vmin, vmax=vmax, smooth=smooth, stretch=scale)
    #figure.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')

    if has_wcs:
        figure.tick_labels.set_xformat('hh:mm:ss')
        figure.tick_labels.set_yformat('dd:mm:ss')

    figure._ax1.set_facecolor('black')
    figure.set_nan_color('black')

    # RECENTER
    if center is not None:
        if radius is None: raise ValueError("Cannot specify center without radius")
        if has_wcs: figure.recenter(center.ra.to("deg").value, center.dec.to("deg").value, radius=radius.to("deg").value)
        else: figure.recenter(center.x, center.y, radius=radius)

    # Hide axes labels and tick labels by default (enable for y for first column and for x for last row)
    figure.axis_labels.hide()
    figure.tick_labels.hide()

    # Axes spines
    figure._ax1.spines['bottom'].set_color('white')
    figure._ax1.spines['top'].set_color('white')
    figure._ax1.spines["left"].set_color("white")
    figure._ax1.spines["right"].set_color("white")

    # TICKS
    #figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
    #figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)
    #figure._ax2.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
    #figure._ax2.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

    # SET LABEL
    figure.add_label(text_x, text_y, r'' + str(label), relative=True, size=13, weight='bold', color='white',
                     horizontalalignment='left', verticalalignment='top',
                     bbox=dict(facecolor='black', edgecolor='none', alpha=0.5))

# ------------------------------------------------------------------------------

def set_ticks(figure, is_first_row, is_last_row):

    """
    This function ...
    :param figure:
    :param is_first_row:
    :param is_last_row:
    :return:
    """

    # ONLY ROW?
    is_only_row = is_first_row and is_last_row

    # ONLY
    if is_only_row:

        # IN EVERYWHERE
        figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
        figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

    # FIRST
    elif is_first_row:

        # LEFT, RIGHT AND TOP
        figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=False, left=True)
        figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=False, left=True)

    # LAST
    elif is_last_row:

        # TOP
        figure._ax1.tick_params(direction='inout', which='major', length=14, top=True, right=False, bottom=False, left=False)
        figure._ax1.tick_params(direction='inout', which='minor', length=8, top=True, right=False, bottom=False, left=False)

        #figure._ax1.tick_params(direction='out', which='major', length=7, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='out', which='minor', length=4, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=False, bottom=False, left=False)

        # BOTTOM, LEFT AND RIGHT
        figure._ax1.tick_params(direction='in', which='major', length=7, right=True, bottom=True, left=True)
        figure._ax1.tick_params(direction='in', which='minor', length=4, right=True, bottom=True, left=True)

        #figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=True, bottom=True, left=True)
        #figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=True, bottom=True, left=True)

    # In between
    else:

        # TOP
        figure._ax1.tick_params(direction='inout', which='major', length=14, top=True, right=False, bottom=False, left=False)
        figure._ax1.tick_params(direction='inout', which='minor', length=8, top=True, right=False, bottom=False, left=False)

        #figure._ax1.tick_params(direction='out', which='major', length=7, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='out', which='minor', length=4, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='in', which='major', length=7, top=True, right=False, bottom=False, left=False)
        #figure._ax1.tick_params(direction='in', which='minor', length=4, top=True, right=False, bottom=False, left=False)

        # LEFT AND RIGHT
        figure._ax1.tick_params(direction='in', which='major', length=7, right=True, bottom=False, left=True)
        figure._ax1.tick_params(direction='in', which='minor', length=4, right=True, bottom=False, left=True)

        #figure._ax1.tick_params(direction='in', which='major', length=7, top=False, right=True, bottom=False, left=True)
        #figure._ax1.tick_params(direction='in', which='minor', length=4, top=False, right=True, bottom=False, left=True)

# ------------------------------------------------------------------------------

def set_theme(dark=False):

    """
    This function ...
    :param dark:
    :return:
    """

    # General settings
    plt.rcParams["axes.labelsize"] = 14  # 16 #default 20
    plt.rcParams["xtick.labelsize"] = 8  # 10 #default 16
    plt.rcParams["ytick.labelsize"] = 8  # 10 #default 16
    plt.rcParams["legend.fontsize"] = 14  # 10 #default 14
    plt.rcParams["legend.markerscale"] = 0
    plt.rcParams["lines.markersize"] = 2.5  # 4 #default 4
    plt.rcParams["axes.linewidth"] = 1

    # Colors
    if dark:
        plt.rcParams['axes.facecolor'] = 'black'
        plt.rcParams['savefig.facecolor'] = 'black'
        plt.rcParams['axes.edgecolor'] = 'white'
        plt.rcParams['xtick.color'] = 'white'
        plt.rcParams['ytick.color'] = 'white'
        plt.rcParams["axes.labelcolor"] = 'white'
        plt.rcParams["text.color"] = 'white'
    else:
        plt.rcParams['axes.facecolor'] = "white"
        plt.rcParams['savefig.facecolor'] = 'white'
        plt.rcParams['axes.edgecolor'] = 'black'
        plt.rcParams['xtick.color'] = 'black'
        plt.rcParams['ytick.color'] = 'black'
        plt.rcParams["axes.labelcolor"] = 'black'
        plt.rcParams["text.color"] = 'black'

# ------------------------------------------------------------------------------

def reset_theme():

    """
    This function ...
    :return:
    """

    # Back to original settings
    plt.rcParams.update(plt.rcParamsDefault)

# ------------------------------------------------------------------------------
