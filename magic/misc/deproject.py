#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.misc.project Contains the Deprojector class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from abc import ABCMeta

# Import the relevant PTS classes and modules
from ...core.basics.configurable import Configurable
from ...core.basics.log import log
from ...core.tools.stringify import tostr
from ..core.frame import Frame
from ...core.simulation.definition import SingleSimulationDefinition
from ...core.tools.utils import lazyproperty
from ...core.tools.numbers import round_to_int
from ...core.prep.dustgrids import create_one_dust_grid_for_galaxy_from_deprojection
from ...core.prep.smile import SKIRTSmileSchema
from ...core.tools import filesystem as fs
from ...core.launch.launcher import SKIRTLauncher
from ...core.units.parsing import parse_unit as u
from ...modeling.plotting.model import xy
from ...core.tools import introspection

# -----------------------------------------------------------------

map_filename = "map.fits"

# -----------------------------------------------------------------

def deproject(deprojection, method="skirt", output_path=None, simulation_path=None, npackages=1e7, parallelization=None):

    """
    This function ...
    :param deprojection:
    :param method:
    :param output_path:
    :param simulation_path:
    :param npackages:
    :param parallelization:
    :return:
    """

    # Create deprojector
    if method == "skirt": deprojector = SKIRTDeprojector()
    elif method == "pts": deprojector = PTSDeprojector()
    else: raise ValueError("Invalid method: '" + method + "'")

    # Set output path
    if output_path is not None:
        deprojector.config.path = output_path
        deprojector.config.write = True
    else: deprojector.config.write = False

    # Set simulation path
    if method == "skirt" and simulation_path is not None: deprojector.config.simulation_path = simulation_path

    # Set number of photon packages
    if method == "skirt": deprojector.config.npackages = npackages

    # Set parallelization
    if method == "skirt" and parallelization is not None: deprojector.config.parallelization = parallelization

    # Deproject
    deprojector.run(deprojection=deprojection)

    # Run the deprojected map
    return deprojector.deprojected

# -----------------------------------------------------------------

class Deprojector(Configurable):

    """
    This class ...
    """

    __metaclass__ = ABCMeta

    # -----------------------------------------------------------------

    def __init__(self, *args, **kwargs):

        """
        This class ...
        :param args:
        :param kwargs:
        """

        # Call the constructor of the base class
        super(Deprojector, self).__init__(*args, **kwargs)

        # The deprojection
        self.deprojection = None

        # The deprojected map
        self.deprojected = None

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        Thisn function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(Deprojector, self).setup(**kwargs)

        # Get the deprojection
        self.deprojection = kwargs.pop("deprojection")

    # -----------------------------------------------------------------

    def write_deprojected(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the deprojected map ...")

        # Determine path
        path = self.output_path_file("deprojected.fits")

        # Save
        self.deprojected.saveto(path)

    # -----------------------------------------------------------------

    @lazyproperty
    def map(self):

        """
        This function ...
        :return:
        """

        return self.deprojection.map

# -----------------------------------------------------------------

class PTSDeprojector(Deprojector):

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
        super(PTSDeprojector, self).__init__(*args, **kwargs)

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Deproject
        self.deproject()

        # 3. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def deproject(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deprojecting the map with PTS ...")

        unit = "pc"

        # Determine the number of pixels
        # npixels = int(round(float(max(self.maps[name].xsize, self.maps[name].ysize)) / self.config.downsample_factor))
        # shape = (npixels, npixels)

        # Debugging
        # log.debug("Using " + str(npixels) + " x " + str(npixels) + " pixels for the deprojected map")

        # Get the x and y range of the model
        x_range_scalar = self.deprojection.x_range.to(unit).value * 2.
        y_range_scalar = self.deprojection.y_range.to(unit).value * 2.

        # Determine the number of pixels
        x_span = x_range_scalar.span
        y_span = y_range_scalar.span
        x_to_y_ratio = x_span / y_span
        ratio = x_to_y_ratio if x_to_y_ratio > 1 else 1. / x_to_y_ratio
        max_npixels = max(self.map.xsize, self.map.ysize)
        largest_dimension = "x" if x_to_y_ratio > 1 else "y"

        if largest_dimension == "x":

            nxpixels = round_to_int(max_npixels / self.config.downsample_factor)
            nypixels = round_to_int(max_npixels / ratio / self.config.downsample_factor)
            exact_ratio = float(nxpixels) / float(nypixels)

        elif largest_dimension == "y":

            nxpixels = round_to_int(max_npixels / ratio / self.config.downsample_factor)
            nypixels = round_to_int(max_npixels / self.config.downsample_factor)
            exact_ratio = float(nypixels) / float(nxpixels)

        else: raise RuntimeError("An error occurred")

        # Debugging
        log.debug("Using " + str(nxpixels) + " x " + str(nypixels) + " pixels for the deprojected map")

        # Adjust the physical ranges to the exact pixel ratios
        if largest_dimension == "x": y_range_scalar = x_range_scalar.compressed(exact_ratio)
        elif largest_dimension == "y": x_range_scalar = y_range_scalar.compressed(exact_ratio)
        else: raise RuntimeError("An error occurred")

        # Debugging
        log.debug("Using an x range of " + tostr(x_range_scalar) + " " + unit)
        log.debug("Using an y range of " + tostr(y_range_scalar) + " " + unit)

        # Determine the output pixel shape
        shape = (nxpixels, nypixels)

        # Set the model limits
        limits = [x_range_scalar.as_tuple(), y_range_scalar.as_tuple()]

        # Create coordinate data
        x, y = xy(shape=shape, limits=limits)

        # Calculate the surface density
        density = self.deprojection.surface_density_function(normalize=True)(x, y)

        # Set deprojected map
        self.deprojected = Frame(density)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Infomr the user
        log.info("Writing ...")

        # Deprojected map
        self.write_deprojected()

# -----------------------------------------------------------------

class SKIRTDeprojector(Deprojector):

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
        super(SKIRTDeprojector, self).__init__(*args, **kwargs)

        # The dust grid
        self.dust_grid = None

        # The ski file
        self.ski = None

        # The SKIRT launcher
        self.launcher = SKIRTLauncher()

        # The simulation path
        self.simulation_path = None

        # The simulation input path
        self.in_path = None

    # -----------------------------------------------------------------

    def _run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 2. Create the dust grid
        self.create_dust_grid()

        # 3. Set the input path
        self.set_input_path()

        # 4. Create the ski file
        self.create_ski()

        # 5. Write the ski file
        self.write_ski()

        # 6. Launch SKIRT
        self.launch()

        # 7. Write
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(SKIRTDeprojector, self).setup(**kwargs)

        # Set the simulation path
        if self.config.simulation_path is not None: self.simulation_path = self.config.simulation_path
        else: self.simulation_path = introspection.create_unique_temp_dir("deprojection")

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        Thisn function ...
        :return:
        """

        return self.deprojection.distance

    # -----------------------------------------------------------------

    def create_dust_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the dust grid ...")

        # Set minimum level
        if self.config.dg.grid_type == "bintree": min_level = self.config.dg.bintree_min_level
        elif self.config.dg.grid_type == "octtree": min_level = self.config.dg.octtree_min_level
        else: min_level = None

        # Set max ndivisions per pixel
        max_ndivisions_per_pixel = 1. / self.config.dg.rel_scale  # default 1/0.5 = 2 divisions along each direction per pixel

        # Create the dust grid
        truncation_ellipse = None
        # grid_type, deprojection, distance, sky_ellipse, min_level, max_mass_fraction, max_ndivisions_per_pixel=2, nscaleheights=10.
        dust_grid = create_one_dust_grid_for_galaxy_from_deprojection(self.config.dg.grid_type, self.deprojection,
                                                                      self.distance,
                                                                      truncation_ellipse,
                                                                      min_level, self.config.dg.max_mass_fraction,
                                                                      max_ndivisions_per_pixel,
                                                                      self.config.dg.scale_heights)

        # Set the dust grid
        self.dust_grid = dust_grid

    # -----------------------------------------------------------------

    @lazyproperty
    def smile(self):

        """
        This function ...
        :return:
        """

        return SKIRTSmileSchema()

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_template(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Creating the ski file template ...")

        # Create
        ski = self.smile.create_oligochromatic_template()

        # Remove the existing instruments
        ski.remove_all_instruments()

        # Remove the stellar system
        ski.remove_stellar_system()

        # Set the number of photon packages
        ski.setpackages(0)

        # Return the ski template
        return ski

    # -----------------------------------------------------------------

    def set_input_path(self):

        """
        This function ...
        :return:
        """

        # Informt the user
        log.info("Setting the input path ...")

        # No map?
        if not self.deprojection.has_map: raise ValueError("The map cannot be loaded")

        # .. because file is present?
        map_path = self.deprojection.filepath
        if fs.is_file(map_path):

            self.in_path = fs.directory_of(map_path)  # directory of input map
            filename = fs.name(map_path)

        # .. or because map is loaded
        elif self.deprojection.map_is_loaded:

            # Create input directory
            self.in_path = fs.create_directory_in(self.simulation_path, "in")

            # Save the map
            map_path = fs.join(self.in_path, map_filename)
            filename = map_filename
            self.deprojection.map.saveto(map_path)  # save the map

        # We shouldn't get here
        else: raise RuntimeError("We shouldn't get here")

        # Set the model map filename
        self.deprojection.filename = filename

    # -----------------------------------------------------------------

    def create_ski(self):

        """
        Thisf unction ...
        :return:
        """

        # Inform the user
        log.info("Creating the ski file ...")

        # Make copy
        ski = self.ski_template.copy()

        # USE DUST GEOMETRY

        # Set filename for deprojection model
        deprojection = self.deprojection.copy()
        deprojection.filename = fs.name(deprojection.filename)

        # Get title
        #title = name

        # Create dust component
        dust_mass = 1e7 * u("Msun")  # dummy value
        mix = "themis"
        #ski.create_new_dust_component(title, geometry=deprojection, normalization_value=dust_mass, mix=mix)
        ski.create_new_dust_component(geometry=deprojection, normalization_value=dust_mass, mix=mix)

        # Set the dust grid
        ski.set_dust_grid(self.dust_grid)

        # Enable writing options
        ski.enable_all_writing_options()

        # Add the ski file
        self.ski = ski

    # -----------------------------------------------------------------

    @lazyproperty
    def ski_path(self):

        """
        Thisn function ...
        :return:
        """

        # Determine path
        return fs.join(self.simulation_path, "deprojection.ski")

    # -----------------------------------------------------------------

    @lazyproperty
    def out_path(self):

        """
        This function ...
        :return:
        """

        return fs.create_directory_in(self.simulation_path, "out")

    # -----------------------------------------------------------------

    @lazyproperty
    def definition(self):

        """
        This function ...
        :return:
        """

        # Create simulation definition
        return SingleSimulationDefinition(self.ski_path, self.out_path, self.in_path)

    # -----------------------------------------------------------------

    def write_ski(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the ski file ...")

        # Save the ski file
        self.ski.saveto(self.ski_path, fix=True)

    # -----------------------------------------------------------------

    def launch(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Launching ...")

        # Set settings
        self.launcher.config.show_progress = True
        self.launcher.config.finish_after = "Writing dust cell properties"  # finish after this line has been printed (when the next one comes)

        # Run
        self.launcher.run(definition=self.definition, parallelization=self.config.parallelization)
        simulation = self.launcher.simulation

        # Determine output filenames
        gridxy_filename = simulation.prefix() + "_ds_grhoxy.fits"
        geometryxy_filename = simulation.prefix() + "_ds_trhoxy.fits"

        # grid_xy_path = fs.join(out_path, gridxy_filename)
        geometry_xy_path = fs.join(self.out_path, geometryxy_filename)

        # Open the output frame
        frame = Frame.from_file(geometry_xy_path)

        # Set the deprojected map
        self.deprojected = frame

        ### EDGEON CUT (NOT AN EDGEON VIEW!!)

        # gridxz_filename = simulation.prefix() + "_ds_grhoxz.fits"
        # geometryxz_filename = simulation.prefix() + "_ds_trhoxz.fits"
        #
        # geometry_xz_path = fs.join(out_path, geometryxz_filename)
        #
        # # Open the frame
        # edgeon = Frame.from_file(geometry_xz_path)
        #
        # # Set the edgeon map
        # self.edgeon[name] = edgeon

    # -----------------------------------------------------------------

    def write(self):

        """
        Thins function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the deprojected map
        self.write_deprojected()

# -----------------------------------------------------------------
