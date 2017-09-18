#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedImageMaker class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from collections import defaultdict

# Import the relevant PTS classes and modules
from ..basics.log import log
from ..tools import filesystem as fs
from ..filter.filter import parse_filter
from ...magic.core.kernel import ConvolutionKernel
from ...magic.core.datacube import DataCube
from ...magic.basics.coordinatesystem import CoordinateSystem
from ...magic.core.remote import RemoteDataCube
from ...magic.core import fits
from ..simulation.wavelengthgrid import WavelengthGrid
from ..basics.configurable import Configurable
from ..simulation.simulation import createsimulations
from ..tools.utils import lazyproperty
from ..tools import types
from ..remote.remote import Remote

# -----------------------------------------------------------------

class ObservedImageMaker(Configurable):

    """
    This class ...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param interactive:
        :return:
        """

        # Call the constructor of the base class
        super(ObservedImageMaker, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The paths to the 'total' FITS files produced by SKIRT
        self.fits_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # The wavelength grid of the simulation
        self.wavelength_grid = None

        # Filter names
        self.filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                             "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

        # The instrument names
        self.instrument_names = None

        # The dictionary containing the different SKIRT output datacubes
        self.datacubes = dict()

        # The dictionary containing the created observation images
        self.images = dict()

        # The coordinate systems of each instrument
        self.coordinate_systems = None

        # The kernel paths
        self.kernel_paths = None

        # The target unit
        self.unit = None

        # The host id
        self.host_id = None

        # Threshold for using remote datacubes
        self.remote_threshold = None

        # Thresholds (frame size) for remote rebinning and convolution
        self.remote_rebin_threshold = None
        self.remote_convolve_threshold = None

        # The path to the output data cubes
        self.paths = defaultdict(dict)

        # The rebin coordinate systems
        self.rebin_coordinate_systems = None

    # -----------------------------------------------------------------

    @property
    def has_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        return self.coordinate_systems is not None

    # -----------------------------------------------------------------

    @property
    def convolution(self):

        """
        Thisf unction ...
        :return:
        """

        return self.kernel_paths is not None

    # -----------------------------------------------------------------

    @property
    def rebinning(self):

        """
        Thisf unction ...
        :return:
        """

        return self.rebin_coordinate_systems is not None

    # -----------------------------------------------------------------

    @property
    def has_remote(self):

        """
        This function ...
        :return:
        """

        return self.host_id is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def remote(self):

        """
        This function ...
        :return:
        """

        return Remote(host_id=self.host_id)

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the wavelength grid
        self.create_wavelength_grid()

        # 3. Load the datacubes
        self.load_datacubes()

        # 4. Set the coordinate systems of the datacubes
        if self.has_coordinate_systems: self.set_coordinate_systems()

        # 5. Make the observed images
        self.make_images()

        # 6. Do convolutions
        if self.convolution: self.convolve()

        # 7. Rebin
        if self.rebinning: self.rebin()

        # 8. Add sky
        self.add_sky()

        # 9. Add stars
        self.add_stars()

        # 10. Do unit conversions
        if self.unit is not None: self.convert_units()

        # 11. Write the results
        if self.output_path is not None: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :return:
        """

        # Call the setup function of the base class
        super(ObservedImageMaker, self).setup(**kwargs)

        # simulation, output_path=None, filter_names=None, instrument_names=None, wcs_path=None,
        # kernel_paths=None, unit=None, host_id=None
        if "simulation" in kwargs: simulation = kwargs.pop("simulation")
        elif "simulation_output_path" in kwargs: simulation = createsimulations(kwargs.pop("simulation_output_path"), single=True)
        else: raise ValueError("Simulation or simulation output path must be specified")

        # Obtain the paths to the 'total' FITS files created by the simulation
        self.fits_paths = simulation.totalfitspaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # Get output directory
        output_path = kwargs.pop("output_path", None)
        self.config.output = output_path

        # Get filter names for which to create observed images
        self.get_filter_names(**kwargs)

        # Get instrument (datacube names)
        self.get_instrument_names(**kwargs)

        # Get coordinate systems of the datacubes
        self.get_coordinate_systems(**kwargs)

        # Get unit for the images
        self.get_unit(**kwargs)

        # Get kernels
        self.get_kernels(**kwargs)

        # Get rebin coordinate systems
        self.get_rebin_coordinate_systems(**kwargs)

        # Get remote host ID
        self.get_host_id(**kwargs)

    # -----------------------------------------------------------------

    def get_filter_names(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Filter names
        if kwargs.get("filter_names", None) is not None:

            # Check
            if "filters" in kwargs: raise ValueError("Cannot specify 'filters' and 'filter_names' simultaneously")

            # Set filter names
            self.filter_names = kwargs.pop("filter_names")

        # Filters
        elif kwargs.get("filters", None) is not None: self.filter_names = [str(fltr) for fltr in kwargs.pop("filters")]

    # -----------------------------------------------------------------

    def get_instrument_names(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the instrument names ...")

        # Instrument names
        if kwargs.get("instrument_names", None) is not None:

            # Check
            if "instruments" in kwargs: raise ValueError("Cannot specify 'instruments' and 'instrument_names' simultaneously")

            # Get names of the instruments (datacubes) of which to create observed images
            self.instrument_names = kwargs.pop("instrument_names")

        # Instruments
        elif kwargs.get("instruments", None) is not None: self.instruments = kwargs.pop("instruments")

    # -----------------------------------------------------------------

    def get_coordinate_systems(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the coordinate systems ...")

        # WCS
        if kwargs.get("wcs", None) is not None:

            # Check that wcs_instrument is defined
            wcs_instrument = kwargs.pop("wcs_instrument")

            # Get the wcs
            wcs = kwargs.pop("wcs")

            # Set the coordinate system
            self.coordinate_systems = dict()
            self.coordinate_systems[wcs_instrument] = wcs

        # WCS paths
        elif kwargs.get("wcs_paths", None) is not None:

            # Get the paths
            wcs_paths = kwargs.pop("wcs_paths")

            # Defined for each instrument
            if types.is_dictionary(wcs_paths):

                # Initialize
                self.coordinate_systems = dict()

                # Loop over the instruments
                for instrument_name in wcs_paths:

                    # Load wcs
                    wcs = CoordinateSystem.from_file(wcs_paths[instrument_name])

                    # Set wcs
                    self.coordinate_systems[instrument_name] = wcs

            # Invalid
            else: raise ValueError("Invalid option for 'wcs_path'")

        # Single WCS path is defined
        elif kwargs.get("wcs_path", None) is not None:

            # Check that wcs_instrument is defined
            wcs_instrument = kwargs.pop("wcs_instrument")

            # Get the wcs
            wcs_path = kwargs.pop("wcs_path")
            wcs = CoordinateSystem.from_file(wcs_path)

            # Set the coordinate system
            self.coordinate_systems = dict()
            self.coordinate_systems[wcs_instrument] = wcs

    # -----------------------------------------------------------------

    def get_unit(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting the target unit ...")

        # Get the unit
        self.unit = kwargs.pop("unit", None)

    # -----------------------------------------------------------------

    def get_kernels(self, **kwargs):
        
        """
        This function ...
        :param kwargs: 
        :return: 
        """

        # Debugging
        log.debug("Getting the kernel paths ...")

        # Checks
        auto_psfs = kwargs.pop("auto_psfs", False)
        if kwargs.get("kernel_paths", None) is not None and kwargs.get("psf_paths", None) is not None: raise ValueError("Cannot specify 'kernel_paths' and 'psf_paths' simultaneously")
        if kwargs.get("psf_paths", None) is not None and auto_psfs: raise ValueError("Cannot specify 'psf_paths' when 'auto_psfs' is enabled")
        if auto_psfs and kwargs.get("kernel_paths", None) is not None: raise ValueError("Cannot specify 'kernel_paths' when 'auto_psfs' is enabled")

        # Kernel paths
        if kwargs.get("kernel_paths", None) is not None: self.kernel_paths = kwargs.pop("kernel_paths")

        # PSF paths
        elif kwargs.pop("psf_paths", None) is not None: self.kernel_paths = kwargs.pop("psf_paths")

        # Automatic PSF determination
        elif kwargs.pop("auto_psfs", None) is not None: self.set_psf_kernels()

    # -----------------------------------------------------------------

    def set_psf_kernels(self):

        """
        This function ...
        :return:
        """

        # Debugging
        log.debug("Determining the PSF kernel automatically for each image filter ...")

        # Get Aniano kernels object
        from pts.magic.convolution.aniano import AnianoKernels
        aniano = AnianoKernels()

        # Initialieze the kernel paths dictionary
        self.kernel_paths = dict()

        # Loop over the filter names
        for filter_name in self.filter_names:

            # Get the psf path
            psf_path = aniano.get_psf_path(filter_name)

            # Set the PSF kernel path
            self.kernel_paths[filter_name] = psf_path

    # -----------------------------------------------------------------

    def get_rebin_coordinate_systems(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting rebin coordinate systems ...")

        # Rebin WCS paths
        if kwargs.get("rebin_wcs_paths", None) is not None:

            # Initialize dictionary
            self.rebin_coordinate_systems = dict()

            # Get the argument
            rebin_wcs_paths = kwargs.pop("rebin_wcs_paths")

            # WCS paths are defined per instrument
            if types.is_dictionary_of_dictionaries(rebin_wcs_paths):

                # Loop over the different instruments
                for instrument_name in rebin_wcs_paths:
                    wcs_dict = dict()
                    # Loop over the filter names
                    for filter_name in rebin_wcs_paths[instrument_name]:

                        # Load the wcs
                        wcs = CoordinateSystem.from_file(rebin_wcs_paths[instrument_name][filter_name])
                        wcs_dict[filter_name] = wcs

                    # Set the coordinate systems for this instrument
                    self.rebin_coordinate_systems[instrument_name] = wcs_dict

            # WCS paths are only defined per filter name
            elif types.is_dictionary(rebin_wcs_paths):

                # Check that rebin_instrument is specified
                rebin_instrument = kwargs.pop("rebin_instrument")

                # Initialize
                self.rebin_coordinate_systems = dict()
                self.rebin_coordinate_systems[rebin_instrument] = dict()

                # Load the coordinate systems
                for filter_name in self.filter_names:
                    wcs = CoordinateSystem.from_file(rebin_wcs_paths[filter_name])
                    self.rebin_coordinate_systems[rebin_instrument][filter_name] = wcs

        # Rebin WCS
        elif kwargs.get("rebin_wcs", None) is not None:

            # Check that rebin_instrument is specified
            rebin_instrument = kwargs.pop("rebin_instrument")

            # Initialize
            self.rebin_coordinate_systems = dict()
            self.rebin_coordinate_systems[rebin_instrument] = dict()

            # Load the coordinate systems
            rebin_wcs = kwargs.pop("rebin_wcs")
            for filter_name in self.filter_names:
                self.rebin_coordinate_systems[rebin_instrument][filter_name] = rebin_wcs

        # Rebin wcs path
        elif kwargs.get("rebin_wcs_path", None) is not None:

            # Check that rebin_instrument is specified
            rebin_instrument = kwargs.pop("rebin_instrument")

            # INitialize
            self.rebin_coordinate_systems = dict()
            self.rebin_coordinate_systems[rebin_instrument] = dict()

            # Load the wcs
            rebin_wcs_path = kwargs.pop("rebin_wcs_path")
            rebin_wcs = CoordinateSystem.from_file(rebin_wcs_path)

            # Set the coordinate systems
            for filter_name in self.filter_names:
                self.rebin_coordinate_systems[rebin_instrument][filter_name] = rebin_wcs

        # Rebin dataset
        elif kwargs.get("rebin_dataset", None) is not None:

            from ...magic.core.dataset import DataSet

            # Get the dataset
            dataset = kwargs.pop("rebin_dataset")
            if types.is_string_type(dataset): dataset = DataSet.from_file(dataset)

            # Check that rebin_instrument is specified
            rebin_instrument = kwargs.pop("rebin_instrument")

            # Initialize
            self.rebin_coordinate_systems = dict()
            self.rebin_coordinate_systems[rebin_instrument] = dict()

            # Loop over the filter names
            for filter_name in self.filter_names:

                # Get the coordinate system
                wcs = dataset.get_coordinate_system_for_filter(filter_name, return_none=True)
                if wcs is None:
                    log.warning("The coordinate system for the '" + filter_name + "' filter is not found in the dataset: skipping ...")
                    continue

                # Set the coordinate system
                self.rebin_coordinate_systems[rebin_instrument][filter_name] = wcs

    # -----------------------------------------------------------------

    def get_host_id(self, **kwargs):

        """
        Thisf unction ...
        :param kwargs:
        :return:
        """

        # Debugging
        log.debug("Getting remote host ...")

        # Get the host ID
        self.host_id = kwargs.pop("host_id", None)

        # Remote threshold
        self.remote_threshold = kwargs.pop("remote_threshold", None)

        # Seperate rebin and convolve thresholds
        self.remote_rebin_threshold = kwargs.pop("remote_rebin_threshold", None)
        self.remote_convolve_threshold = kwargs.pop("remote_convolve_threshold", None)

    # -----------------------------------------------------------------

    def create_wavelength_grid(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating the wavelength grid ...")

        # Construct the wavelength grid from the array of wavelengths
        self.wavelength_grid = WavelengthGrid.from_wavelengths(self.wavelengths, "micron")

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        return {filter_name: parse_filter(filter_name) for filter_name in self.filter_names}

    # -----------------------------------------------------------------

    def load_datacubes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the SKIRT output datacubes ...")

        # Loop over the different simulated images
        for path in self.fits_paths:

            # Get the name of the instrument
            instr_name = instrument_name(path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Get the name of the datacube (as given by SKIRT)
            datacube_name = fs.strip_extension(fs.name(path))

            # Debugging
            log.debug("Loading datacube from '" + datacube_name + ".fits' ...")

            ## LOAD AND CONVERT UNITS TO SPECTRAL (WAVELENGTH-DENSITY)

            # Load the datacube (locally or remotely)
            if self.host_id is not None:

                # Remote threshold not specified or file is smaller than threshold
                if self.remote_threshold is None or fs.file_size(path) < self.remote_threshold:
                    try: datacube = DataCube.from_file(path, self.wavelength_grid)
                    except fits.DamagedFITSFileError as e:
                        log.error("The datacube '" + path + "' is damaged: images cannot be created. Skipping this datacube ...")
                        continue

                # Remote threshold
                else:

                    try: datacube = RemoteDataCube.from_file(path, self.wavelength_grid, self.host_id)
                    except fits.DamagedFITSFileError as e:
                        log.error("The datacube '" + path + "' is damaged: images cannot be created. Skipping this datacube ...")
                        continue

            # No host specified: local datacube
            else:

                try: datacube = DataCube.from_file(path, self.wavelength_grid)
                except fits.DamagedFITSFileError as e:
                    log.error("The datacube '" + path + "' is damaged: images cannot be created. Skipping this datacube ...")
                    continue

            # Convert the datacube from neutral flux density to wavelength flux density
            datacube.to_wavelength_density("W / (m2 * arcsec2 * micron)", "micron")

            # Add the datacube to the dictionary
            self.datacubes[datacube_name] = datacube

    # -----------------------------------------------------------------

    def set_coordinate_systems(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Setting the WCS of the simulated images ...")

        # Loop over the different datacubes and set the WCS
        for datacube_name in self.datacubes:

            # Check whether coordinate system is defined for this datacube
            if datacube_name not in self.coordinate_systems: continue

            # Debugging
            log.debug("Setting the coordinate system of the " + datacube_name + "' datacube ...")

            # Set the coordinate system for this datacube
            self.datacubes[datacube_name].wcs = self.coordinate_systems[datacube_name]

    # -----------------------------------------------------------------

    def make_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed images (this may take a while) ...")

        # Loop over the datacubes
        for datacube_name in self.datacubes:

            # Debugging
            log.debug("Making the observed images for " + datacube_name + ".fits ...")

            # Create a list of the filter names
            filter_names = self.filters.keys()

            # Create the corresponding list of filters
            filters = self.filters.values()

            # Initialize a dictionary, indexed by the filter names, to contain the images
            images = dict()

            # Determine the number of processes
            if isinstance(self.datacubes[datacube_name], RemoteDataCube): nprocesses = self.config.nprocesses_remote
            elif isinstance(self.datacubes[datacube_name], DataCube): nprocesses = self.config.nprocesses_local
            else: raise ValueError("Invalid datacube object for '" + datacube_name + "' instrument")

            # Create the observed images from the current datacube (the frames get the correct unit, wcs, filter)
            frames = self.datacubes[datacube_name].frames_for_filters(filters, convolve=self.config.spectral_convolution, nprocesses=nprocesses)

            # Add the observed images to the dictionary
            for filter_name, frame in zip(filter_names, frames): images[filter_name] = frame # these frames can be RemoteFrames if the datacube was a RemoteDataCube

            # Add the observed image dictionary for this datacube to the total dictionary (with the datacube name as a key)
            self.images[datacube_name] = images

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving the images ...")

        from ...magic.core.remote import RemoteFrame
        from ...magic.core.frame import Frame

        session = None

        # Loop over the images
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Check if the name of the image filter is a key in the 'kernel_paths' dictionary. If not, don't convolve.
                if filter_name not in self.kernel_paths or self.kernel_paths[filter_name] is None: continue

                # Check whether the pixelscale is defined
                if self.images[datacube_name][filter_name].pixelscale is None: raise ValueError("Pixelscale of the '" + filter_name + "' image of the '" + datacube_name + "' datacube is not defined, convolution not possible")

                # Load the kernel
                kernel = ConvolutionKernel.from_file(self.kernel_paths[filter_name])

                # Debugging
                log.debug("Convolving the '" + filter_name + "' image of the '" + datacube_name + "' instrument ...")

                # Get the frame
                frame = self.images[datacube_name][filter_name]

                # Convert into remote frame if necessary
                if self.remote_convolve_threshold is not None and isinstance(frame, Frame) and frame.data_size > self.remote_convolve_threshold:

                    # Create session if necessary
                    if session is None:
                        # START SESSION
                        new_connection = False
                        session = self.remote.start_python_session(attached=True, new_connection_for_attached=new_connection)

                    # Convert into remote
                    self.images[datacube_name][filter_name] = RemoteFrame.from_local(frame, session)

                # Convolve the frame
                self.images[datacube_name][filter_name].convolve(kernel)

        # End the session
        if session is not None: del session

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning the images to the requested coordinate systems ...")

        from ...magic.core.remote import RemoteFrame
        from ...magic.core.frame import Frame

        session = None

        # Loop over the datacubes
        for datacube_name in self.images:

            # Check if the name of the datacube appears in the rebin_wcs dictionary
            if datacube_name not in self.rebin_coordinate_systems: continue

            # Loop over the filters
            for filter_name in self.images[datacube_name]:

                # Check if the name of the image appears in the rebin_wcs[datacube_name] sub-dictionary
                if filter_name not in self.rebin_coordinate_systems[datacube_name]: continue

                # Debugging
                log.debug("Rebinning the '" + filter_name + "' image of the '" + datacube_name + "' instrument ...")

                # Get the original unit
                original_unit = self.images[datacube_name][filter_name].unit
                converted = False

                # Check if this is a surface brightness or intensity
                if not (original_unit.is_intensity or original_unit.is_surface_brightness):

                    # Determine the new (surface brightness or intensity) unit
                    #new_unit = original_unit / u("sr")
                    new_unit = original_unit.corresponding_angular_area_unit

                    # Debugging
                    log.debug("Converting the unit from '" + str(original_unit) + "' to '" + str(new_unit) + "' in order to be able to perform rebinning ...")

                    # Convert
                    self.images[datacube_name][filter_name].convert_to(new_unit)
                    converted = True

                # Get frame and target WCS
                frame = self.images[datacube_name][filter_name]
                wcs = self.rebin_coordinate_systems[datacube_name][filter_name]

                # Convert to remote frame if necessary
                if self.remote_rebin_threshold is not None and isinstance(frame, Frame) and frame.data_size > self.remote_rebin_threshold:

                    # Create session if necessary
                    if session is None:
                        # START SESSION
                        new_connection = False
                        session = self.remote.start_python_session(attached=True, new_connection_for_attached=new_connection)

                    # Convert
                    self.images[datacube_name][filter_name] = RemoteFrame.from_local(frame, session)

                # Rebin
                self.images[datacube_name][filter_name].rebin(wcs)

                # Convert the unit back
                if converted: self.images[datacube_name][filter_name].convert_to(original_unit)

        # End the session
        if session is not None: del session

    # -----------------------------------------------------------------

    def add_sky(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding artificial sky contribution to the images ...")

    # -----------------------------------------------------------------

    def add_stars(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Adding artificial stars to the images ...")

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the units of the images to " + str(self.unit) + " ...")

        # Loop over the images
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Debugging
                log.debug("Converting the unit of the " + filter_name + " image of the '" + datacube_name + "' instrument ...")

                # Convert
                factor = self.images[datacube_name][filter_name].convert_to(self.unit)

                # Debugging
                log.debug("The conversion factor is '" + str(factor) + "'")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the images ...")

        # Write (grouped or not)
        if self.config.group: self.write_images_grouped()
        else: self.write_images()

    # -----------------------------------------------------------------

    def write_images_grouped(self):

        """
        This function ...
        :return:
        """

        # Loop over the different instruments (datacubes)
        for datacube_name in self.images:

            # Make directory for this datacube
            datacube_path = self.output_path_directory(datacube_name)

            # Loop over the images
            for filter_name in self.images[datacube_name]:

                # Determine path to the output FITS file
                path = fs.join(datacube_path, filter_name + ".fits")

                # Save the image
                self.images[datacube_name][filter_name].saveto(path)

                # Set the path
                self.paths[datacube_name][filter_name] = path

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Loop over the different images (self.images is a nested dictionary of dictionaries)
        for datacube_name in self.images:
            for filter_name in self.images[datacube_name]:

                # Determine the path to the output FITS file
                path = self.output_path_file(datacube_name + "__" + filter_name + ".fits")

                # Save the image
                self.images[datacube_name][filter_name].saveto(path)

                # Set the path
                self.paths[datacube_name][filter_name] = path

# -----------------------------------------------------------------

def instrument_name(datacube_path, prefix):

    """
    This function ...
    :param datacube_path:
    :param prefix:
    :return:
    """

    return fs.name(datacube_path).split("_total.fits")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------
