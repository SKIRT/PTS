#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.animations Contains the DataCubeAnimationsMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..basics.log import log
from .datacubes import DatacubesMiscMaker, get_datacube_instrument_name
from ..basics.animation import Animation

# -----------------------------------------------------------------

class DataCubeAnimationsMaker(DatacubesMiscMaker):

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
        super(DataCubeAnimationsMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The frames path
        self.frames_path = None

        # The plot paths per instrument
        self.plot_paths = dict()

        # The frames
        self.frames = defaultdict(list)

        # The animations
        self.animations = dict()

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 2. Make frames
        self.make_frames()

        # 3. Make the animations
        self.make_animations()

        # 4. Write
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(DataCubeAnimationsMaker, self).setup(**kwargs)

        # Get output directory
        output_path = kwargs.pop("output_path", None)
        self.config.output = output_path

        # Create frames directory?
        if self.config.write_frames: self.frames_path = self.output_path_directory("frames", create=True)
        else: self.frames_path = self.output_path_directory("frames", create=False)

    # -----------------------------------------------------------------

    def make_frames(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Making the animation frames ...")

        from ...magic.core.datacube import DataCube
        from ...magic.core.rgba import RGBAImage
        from ...magic.core.rgb import DamagedImageFileError

        # Loop over the total datacube paths
        for path in self.total_datacube_paths:

            # Get the name of the instrument
            instr_name = get_instrument_name(path, self.simulation_prefix)

            # Make for this instrument?
            if not self.make_for_instrument(instr_name): continue

            # Debugging
            log.debug("Making frames for the '" + instr_name + "' instrument ...")

            # Create a plot directory for this instrument
            plot_path = fs.join(self.frames_path, instr_name)
            if self.config.write_frames: fs.create_directory(plot_path)

            # Load the datacube
            datacube = DataCube.from_file(path, self.wavelength_grid)

            # Plot paths for this instrument
            plot_paths = []

            # Loop over the frames, create PNG image for each frame
            for index in range(datacube.nframes):

                # Determine frame plot path
                frame_plot_path = fs.join(plot_path, str(index) + ".png")

                # Check if already exists
                if fs.is_file(frame_plot_path):

                    log.success("Frame " + str(index + 1) + " has already been created: loading it from file ...")

                    # Add the path
                    plot_paths.append(frame_plot_path)

                    # Load
                    try:
                        image = RGBAImage.from_file(frame_plot_path)
                        self.frames[instr_name].append(image)
                        # Skip creating the frame
                        continue
                    except DamagedImageFileError as e:
                        log.warning("The image frame is damaged: removing it and creating it again ...")
                        fs.remove_file(frame_plot_path)

                # Debugging
                log.debug("Creating frame " + str(index + 1) + " of " + str(self.nwavelengths) + " (" + str(float(index+1)/self.nwavelengths*100) + "%) ...")

                # Get the frame
                frame = datacube.frames[index]

                # Create RGBA image
                image, vmin, vmax = frame.to_rgba(interval=self.config.interval, scale=self.config.scale, alpha=self.config.alpha_method,
                                                  peak_alpha=self.config.peak_alpha, colours=self.config.colours, normalize_in=None,
                                                  return_minmax=True)

                # Add the frame
                self.frames[instr_name].append(image)

                # Save as PNG?
                if plot_path is not None:

                    # Add the path
                    plot_paths.append(frame_plot_path)

                    # Save
                    image.saveto(frame_plot_path)

            # Set plot paths for this instrument
            self.plot_paths[instr_name] = plot_paths

    # -----------------------------------------------------------------

    def make_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the animations ...")

        # Loop over the instruments
        for instrument_name in self.frames:

            # Initialize the animation
            animation = Animation()
            animation.fps = self.config.fps

            # Add the frames
            for frame in self.frames[instrument_name]:
                array = frame.asarray()
                animation.add_frame(array)

            # Add the animation
            self.animations[instrument_name] = animation

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the animation
        self.write_animations()

    # -----------------------------------------------------------------

    def write_animations(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the animations ...")

        # Loop over the instruments
        for instrument_name in self.animations:

            # Debugging
            log.debug("Writing the animation for the '" + instrument_name + "' instrument ...")

            # Determine the path
            path = self.output_path_file("animation_" + instrument_name + "." + self.config.format)

            # Write the animation
            self.animations[instrument_name].saveto(path)

# -----------------------------------------------------------------
