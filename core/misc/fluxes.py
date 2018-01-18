#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.core.misc.fluxes Contains the ObservedFluxCalculator class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np
from collections import defaultdict, OrderedDict
from scipy.interpolate import interp1d

# Import astronomical modules
from astropy.units import spectral_density, spectral, Unit

# Import the relevant PTS classes and modules
from ..tools import filesystem as fs
from ..basics.log import log
from ..filter.broad import BroadBandFilter
from ..filter.filter import parse_filter
from ..data.sed import SED
from ...magic.services.spire import SPIRE
from ..units.parsing import parse_unit as u
from ..data.sed import ObservedSED
from ..tools import sequences
from ..basics.configurable import Configurable
from ..simulation.simulation import createsimulations
from ..tools import parsing
from ..tools import types
from ..tools import numbers
from ..tools.utils import lazyproperty
from .images import get_datacube_instrument_name
from ...magic.core.datacube import DataCube
from ...magic.core.remote import RemoteDataCube
from ...magic.core import fits
from ..simulation.wavelengthgrid import WavelengthGrid
from ..plot.sed import SEDPlotter
from ..tools import formatting as fmt
from ..remote.remote import Remote
from ..prep.deploy import Deployer

# -----------------------------------------------------------------

default_filter_names = ["FUV", "NUV", "u", "g", "r", "i", "z", "H", "J", "Ks", "I1", "I2", "I3", "I4", "W1", "W2",
                        "W3", "W4", "Pacs 70", "Pacs 100", "Pacs 160", "SPIRE 250", "SPIRE 350", "SPIRE 500"]

# -----------------------------------------------------------------

class ObservedFluxCalculator(Configurable):

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
        super(ObservedFluxCalculator, self).__init__(*args, **kwargs)

        # -- Attributes --

        # The simulation prefix
        self.simulation_prefix = None

        # The ski file
        self.ski = None

        # The paths to the SED files produced by SKIRT
        self.sed_paths = None

        # The paths to the total datacube files produced by SKIRT
        self.datacube_paths = None

        # The wavelengths of the simulation
        self.wavelengths = None

        # Filter names
        self.filter_names = default_filter_names

        # Instrument names
        self.instrument_names = None

        # The errors for the different filters
        self.errors = None

        # No spectral convolution for certain filters?
        self.no_spectral_convolution_filters = []

        # The output observed SEDs, per instrument
        self.mock_seds = dict()

        # The SPIRE instance
        self.spire = SPIRE()

        # The SED paths
        self.paths = dict()

        # The masks
        self.masks = None

        # The coordinate systems
        self.coordinate_systems = None

        # The images, per instrument
        self.images = dict()

        # For plotting, reference SED
        self.reference_sed = None

        # The host id
        self.host_id = None

        # Remote options
        self.remote_images_spectral_convolution = False
        self.remote_threshold = None
        self.remote_npixels_threshold = None
        self.remote_rebin_threshold = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Create the images
        if self.config.from_images: self.create_images()

        # 3. Calculate the observed fluxes
        self.calculate()

        # 4. Write the results
        if self.config.output is not None: self.write()

        # 5. Plot
        if self.config.plot: self.plot()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(ObservedFluxCalculator, self).setup(**kwargs)

        # Get properties
        if "simulation" in kwargs: simulation = kwargs.pop("simulation")
        elif "simulation_output_path" in kwargs: simulation = createsimulations(kwargs.pop("simulation_output_path"), single=True)
        else: raise ValueError("Either 'simulation' or 'simulation_output_path' must be specified")
        output_path = kwargs.pop("output_path", None)
        filter_names = kwargs.pop("filter_names", None)
        instrument_names = kwargs.pop("instrument_names", None)
        errors = kwargs.pop("errors", None)

        # Obtain the paths to the datacubes created by the simulation
        if self.config.from_images: self.datacube_paths = simulation.totalfitspaths()

        # Obtain the paths to the SED files created by the simulation
        else: self.sed_paths = simulation.seddatpaths()

        # Get the list of wavelengths for the simulation
        self.wavelengths = simulation.wavelengths()

        # Get the simulation prefix
        self.simulation_prefix = simulation.prefix()

        # The ski file
        self.ski = simulation.ski_file

        # No spectral convolution for certain filters?
        self.no_spectral_convolution_filters = kwargs.pop("no_spectral_convolution_filters", [])

        # Set the filter names
        if filter_names is not None: self.filter_names = filter_names

        # Set the instrument names
        self.instrument_names = instrument_names

        # Set the errors
        if errors is not None:
            self.errors = dict()
            # Set the error bars from strings
            for filter_name in errors:
                error = errors[filter_name]
                if types.is_string_type(error): error = parsing.errorbar(error)
                self.errors[filter_name] = error

        # Set output path
        self.config.output = output_path

        # Get the masks
        self.masks = kwargs.pop("masks", None)

        # Get the coordinate systems
        self.coordinate_systems = kwargs.pop("coordinate_systems", None)

        # Reference SED
        self.reference_sed = kwargs.pop("reference_sed", None)

        # Get remote host ID
        self.get_host_id(**kwargs)

        # Update the remote
        if self.has_remote and self.config.deploy_pts: self.deploy_pts()

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

        # Remote spectral convolution flag
        self.remote_images_spectral_convolution = kwargs.pop("remote_images_spectral_convolution", False)

        # Threshold for remote
        self.remote_threshold = kwargs.pop("remote_threshold", None)
        self.remote_npixels_threshold = kwargs.pop("remote_npixels_threshold", None)
        self.remote_rebin_threshold = kwargs.pop("remote_rebin_threshold", None)

    # -----------------------------------------------------------------

    def deploy_pts(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Deploying PTS remotely ...")

        # Create the deployer
        deployer = Deployer()

        # Don't do anything locally
        deployer.config.local = False

        # Only deploy PTS
        deployer.config.skirt = False
        deployer.config.pts = True

        # Set the host ids
        deployer.config.host_ids = [self.host_id]

        # Check versions between local and remote
        deployer.config.check = self.config.check_versions

        # Update PTS dependencies
        deployer.config.update_dependencies = self.config.update_dependencies

        # Run the deployer
        deployer.run()

    # -----------------------------------------------------------------

    @property
    def has_reference_sed(self):

        """
        This function ...
        :return:
        """

        return self.reference_sed is not None

    # -----------------------------------------------------------------

    @lazyproperty
    def wavelength_grid(self):

        """
        This function ...
        :return:
        """

        return WavelengthGrid.from_wavelengths(self.wavelengths, "micron")

    # -----------------------------------------------------------------

    def has_mask(self, instr_name, fltr):

        """
        This function ...
        :param instr_name:
        :param fltr:
        :return:
        """

        return instr_name in self.masks and fltr in self.masks[instr_name] and self.masks[instr_name][fltr] is not None

    # -----------------------------------------------------------------

    def get_mask(self, instr_name, fltr):

        """
        This function ...
        :param instr_name:
        :param fltr:
        :return:
        """

        return self.masks[instr_name][fltr]

    # -----------------------------------------------------------------

    def create_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Creating images from the simulated datacubes ...")

        # Loop over the different SEDs
        for datacube_path in self.datacube_paths:

            # Get the name of the instrument
            instr_name = get_datacube_instrument_name(datacube_path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Check if already encountered this instrument
            if instr_name in self.images: raise ValueError("Already encountered datacube for '" + instr_name + "' instrument")

            # Debugging
            log.debug("Creating images from the '" + instr_name + "' instrument ...")

            # Try loading the datacube
            datacube = self.load_datacube(datacube_path, instr_name)
            if datacube is None: continue

            # Get the instrument distance
            distance = self.ski.get_instrument_distance(instr_name)

            # If the distance is defined, set the distance
            if distance is not None: datacube.distance = distance

            # Debugging
            log.debug("Setting the coordinate system of the '" + instr_name + "' instrument ...")

            # Set the coordinate system for this datacube
            datacube.wcs = self.coordinate_systems[instr_name]

            # Debugging
            log.debug("Checking the units of the image ...")

            # Convert the datacube from neutral flux density to wavelength flux density
            datacube.convert_to_corresponding_wavelength_density_unit() # why?

            # Convert to non- angular or intrinsic area unit
            if datacube.is_per_angular_or_intrinsic_area: datacube.convert_to_corresponding_non_angular_or_intrinsic_area_unit()

            #print(self.config.spectral_convolution)
            #print("spectral convolution filters", self.spectral_convolution_filters)

            spec_filters = []
            #print("filters", self.filters)

            # Create the observed images from the current datacube (the frames get the correct unit, wcs, filter)
            nprocesses = 1
            frames = datacube.frames_for_filters(self.filters, convolve=spec_filters, nprocesses=nprocesses, check_previous_sessions=True, as_dict=True)

            # Mask the images
            for fltr in frames:

                # Check if should be masked
                if not self.has_mask(instr_name, fltr): continue

                # Get the frame
                frame = frames[fltr]

                # Get mask
                mask = self.get_mask(instr_name, fltr)

                # Rebin the frame
                frame.rebin(mask.wcs, convert=True)

                # Mask
                frame.apply_mask_nans(mask)

            # Add the frames for this instrument
            self.images[instr_name] = frames

    # -----------------------------------------------------------------

    def load_datacube(self, path, instr_name):

        """
        This function ...
        :param path:
        :param instr_name:
        :return:
        """

        # Debugging
        log.debug("Loading total datacube of '" + instr_name + "' instrument from '" + path + "' ...")

        # Load datacube remotely
        if self.needs_remote(path): datacube = self.load_datacube_remote(path)

        # Load datacube locally
        else: datacube = self.load_datacube_local(path)

        # Return
        return datacube

    # -----------------------------------------------------------------

    def load_datacube_local(self, path):

        """
        Thisj function ...
        :param self:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Trying to load the datacube '" + path + "' locally ...")

        # Load
        try: datacube = DataCube.from_file(path, self.wavelength_grid)
        except fits.DamagedFITSFileError as e:
            log.error("The datacube '" + path + "' is damaged: images cannot be created. Skipping this datacube ...")
            datacube = None

        # Return the datacube
        return datacube

    # -----------------------------------------------------------------

    def load_datacube_remote(self, path):

        """
        This function ...
        :param self:
        :param path:
        :return:
        """

        # Debugging
        log.debug("Trying to load the datacube '" + path + "' remotely ...")

        # Load
        try: datacube = RemoteDataCube.from_file(path, self.wavelength_grid, self.session)
        except fits.DamagedFITSFileError as e:
            log.error("The datacube '" + path + "' is damaged: images cannot be created. Skipping this datacube ...")
            datacube = None

        # Return
        return datacube

    # -----------------------------------------------------------------

    @lazyproperty
    def nspectral_convolution_filters(self):

        """
        This function ...
        :return:
        """

        return len(self.spectral_convolution_filters)

    # -----------------------------------------------------------------

    @property
    def has_spectral_convolution_filters(self):

        """
        This function ...
        :return:
        """

        return self.nspectral_convolution_filters > 0

    # -----------------------------------------------------------------

    def needs_remote(self, path):

        """
        This function ...
        :return:
        """

        from ...magic.core.fits import get_npixels

        # No remote is set
        if self.host_id is None: return False

        # File size is exceeded
        if self.remote_threshold is not None and fs.file_size(path) > self.remote_threshold: return True

        # Number of pixels is exceeded
        if self.remote_npixels_threshold is not None and get_npixels(path) > self.remote_npixels_threshold: return True

        # Remote spectral convolution
        if self.has_spectral_convolution_filters and self.remote_images_spectral_convolution: return True

        # Not remote
        return False

    # -----------------------------------------------------------------

    @lazyproperty
    def filters(self):

        """
        This function ...
        :return:
        """

        return [parse_filter(filter_name) for filter_name in self.filter_names]

    # -----------------------------------------------------------------

    @lazyproperty
    def spectral_convolution_filters(self):

        """
        This function ...
        :return:
        """

        # No spectral convolution for any filter
        if not self.config.spectral_convolution: return []

        # Initialize list
        filters = []

        # Loop over the filters
        for fltr in self.filters:

            if fltr in self.no_spectral_convolution_filters: pass
            else: filters.append(fltr)

        # Return the filters
        return filters

    # -----------------------------------------------------------------

    def calculate(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the observed fluxes ...")

        # From images
        if self.config.from_images: self.calculate_from_images()

        # From SEDs
        else: self.calculate_from_seds()

    # -----------------------------------------------------------------

    def calculate_from_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating fluxes from images ...")

        # Loop over the different instruments
        for instr_name in self.images:

            # Initialize an SED
            # Create an observed SED for the mock fluxes
            mock_sed = ObservedSED(photometry_unit="Jy")

            # Loop over the filters
            for fltr in self.images[instr_name]:

                # Set filter name
                filter_name = str(fltr)

                # Get the image frame
                frame = self.images[instr_name][fltr]

                # Calculate proper flux
                flux = frame.sum(add_unit=True)

                # Is there an error?
                # Add a point to the mock SED
                error = self.errors[filter_name] if self.errors is not None and filter_name in self.errors else None

                # Add this entry to the SED
                mock_sed.add_point(fltr, flux, error)

            # Add the mock SED for this instrument
            self.mock_seds[instr_name] = mock_sed

    # -----------------------------------------------------------------

    def calculate_from_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating fluxes from SEDs ...")

        # Loop over the different SEDs
        for sed_path in self.sed_paths:

            # Get the name of the instrument
            instr_name = get_sed_instrument_name(sed_path, self.simulation_prefix)

            # If a list of instruments is defined an this instrument is not in this list, skip it
            if self.instrument_names is not None and instr_name not in self.instrument_names: continue

            # Set the conversion info
            conversion_info = dict()
            conversion_info["distance"] = self.ski.get_instrument_distance(instr_name)

            # Debugging
            log.debug("Loading the modelled SED ...")

            # Load the modelled SED
            model_sed = SED.from_skirt(sed_path)

            # Debugging
            log.debug("Calculating the observed fluxes for the " + instr_name + " SED ...")

            # Convert model sed into observed SED
            mock_sed = create_mock_sed(model_sed, self.filters, self.spire, spectral_convolution=self.spectral_convolution_filters, errors=self.errors)

            # Add the complete SED to the dictionary (with the SKIRT SED name as key)
            self.mock_seds[instr_name] = mock_sed

    # -----------------------------------------------------------------

    @property
    def has_simulated_seds(self):

        """
        This function ...
        :return:
        """

        return self.sed_paths is not None

    # -----------------------------------------------------------------

    def get_simulated_sed(self, instr_name):

        """
        This function ...
        :param instr_name:
        :return:
        """

        # Loop over the different SEDs
        for sed_path in self.sed_paths:

            # Get the name of the instrument
            instrument_name = get_sed_instrument_name(sed_path, self.simulation_prefix)

            # Load the SED
            if instrument_name == instr_name: return SED.from_skirt(sed_path)

        # Error
        raise ValueError("Cannot find simulated SED for the '" + instr_name + "' instrument")

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the mock SEDs
        self.write_mock_seds()

        # Write the images
        if self.config.from_images and self.config.write_images: self.write_images()

    # -----------------------------------------------------------------

    def write_mock_seds(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mock observed SED ...")

        # Loop over the different flux tables
        for instr_name in self.mock_seds:

            # Determine the path to the output flux table
            path = self.output_path_file(instr_name + "_fluxes.dat")

            # Debugging
            log.debug("Writing the mock SED '" + instr_name + "' to '" + path + "' ...")

            # Write out the flux table
            self.mock_seds[instr_name].saveto(path)

            # Set the path
            self.paths[instr_name] = path

    # -----------------------------------------------------------------

    @lazyproperty
    def images_output_path(self):

        """
        This function ...
        :return:
        """

        return self.output_path_directory("images", create=True)

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the mock observed images ...")

        # Loop over the instruments
        for instr_name in self.images:

            # Debugging
            log.debug("Writing mock observed images for the '" + instr_name + "' instrument ...")

            # Create directory for this instrument
            instr_path = fs.create_directory_in(self.images_output_path, instr_name)

            # Loop over the filters
            for fltr in self.images[instr_name]:

                # Debugging
                log.debug("Writing mock observed '" + str(fltr) + "' image ...")

                # Get the image frame
                frame = self.images[instr_name][fltr]

                # Determine path
                path = fs.join(instr_path, str(fltr) + ".fits")

                # Save the image
                frame.saveto(path)

    # -----------------------------------------------------------------

    def plot(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting ...")

        # Plot the mock observed SEDs
        if self.config.plot_seds: self.plot_seds()

        # Plot the images
        if self.config.plot_images: self.plot_images()

    # -----------------------------------------------------------------

    def plot_seds(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SEDs ...")

        # Loop over the different flux tables
        for instr_name in self.mock_seds:

            # Determine the path to the plot file
            path = self.output_path_file(instr_name + "_fluxes.pdf")

            # Debugging
            log.debug("Plotting the mock SED '" + instr_name + "' to '" + path + "' ...")

            # Plot
            self.plot_sed(instr_name, path)

    # -----------------------------------------------------------------

    def plot_sed(self, instr_name, path):

        """
        This function ...
        :param instr_name:
        :param path:
        :return:
        """

        # Initialize the plotter
        plotter = SEDPlotter()

        # Get the mock observed sed
        sed = self.mock_seds[instr_name]

        # Add the SED
        plotter.add_sed(sed, "Mock observation")

        # Add simulated SED?
        if self.has_simulated_seds:

            # Get the SED
            simulated_sed = self.get_simulated_sed(instr_name)

            # Add
            plotter.add_sed(simulated_sed, "Simulation")

        # Add reference?
        if self.has_reference_sed: plotter.add_sed(self.reference_sed, "Observation")

        # Run the plotter
        plotter.run(output=path)

    # -----------------------------------------------------------------

    def plot_images(self):

        """
        This function ...
        :return:
        """

        from ...magic.plot.imagegrid import StandardImageGridPlotter

        # Inform the user
        log.info("Plotting the images ...")

        # Loop over the instruments
        for instr_name in self.images:

            # Debugging
            log.debug("Plotting the images for the '" + instr_name + "' instrument ...")

            # Create directory for this instrument
            instr_path = fs.create_directory_in(self.images_output_path, instr_name)

            # Create plotter
            plotter = StandardImageGridPlotter()

            # Set output directory
            plotter.config.output = instr_path

            # Extra
            plotter.config.normalize = True
            # plotter.config.colormap =

            # Write data
            plotter.config.write = False

            # Rebin and crop
            # plotter.rebin_to =
            # plotter.crop_to =

            # Loop over the filters
            for fltr in self.images[instr_name]:

                # Get the image frame
                frame = self.images[instr_name][fltr]

                # Add to plot
                plotter.add_frame(frame)

            # Run the plotter
            plotter.run()

# -----------------------------------------------------------------

def get_sed_instrument_name(sed_path, prefix):

    """
    This function ...
    :param sed_path:
    :param prefix:
    :return:
    """

    return fs.name(sed_path).split("_sed.dat")[0].split(prefix + "_")[1]

# -----------------------------------------------------------------

def get_wavelength_and_fluxdensity_arrays(model_sed, conversion_info=None, wavelength_unit="micron", photometry_unit="W/(m2*micron)", return_units=False):

    """
    This function ...
    :param model_sed: 
    :param conversion_info:
    :param wavelength_unit:
    :param photometry_unit:
    :param return_units:
    :return: 
    """

    # Get arrays
    wavelengths = model_sed.wavelengths(wavelength_unit, asarray=True)
    fluxdensities = model_sed.photometry(photometry_unit, asarray=True, conversion_info=conversion_info)

    # Return the wavelengths and fluxdensities
    if return_units: return wavelengths, fluxdensities, u(wavelength_unit), u(photometry_unit)
    else: return wavelengths, fluxdensities

# -----------------------------------------------------------------

def get_wavelength_and_fluxdensity_arrays_old(model_sed):

    """
    This function ...
    :param model_sed: 
    :return: 
    """

    # Get the wavelengths and flux densities
    wavelengths = model_sed.wavelengths("micron", asarray=True)
    fluxdensities = []
    for wavelength, fluxdensity_jy in zip(model_sed.wavelengths("micron"), model_sed.photometry("Jy", conversion_info=conversion_info)):

        # 2 different ways should be the same:
        # fluxdensity_ = fluxdensity_jy.to("W / (m2 * micron)", equivalencies=spectral_density(wavelength))
        # fluxdensity = fluxdensity_jy.to("W / (m2 * Hz)").value * spectral_factor_hz_to_micron(wavelength) * u("W / (m2 * micron)")
        # print(fluxdensity_, fluxdensity) # IS OK!

        fluxdensity = fluxdensity_jy.to("W / (m2 * micron)", wavelength=wavelength)
        fluxdensities.append(fluxdensity.value)

        # Add to data points for black body (to get spectral index)
        if wavelength > 50. * u("micron"):

            bb_frequencies.append(wavelength.to("Hz", equivalencies=spectral()).value)
            bb_fluxdensities.append(fluxdensity_jy.to("W / (m2 * Hz)").value)  # must be frequency-fluxdensity

    # Convert to array
    fluxdensities = np.array(fluxdensities)  # in W / (m2 * micron)

# -----------------------------------------------------------------

def calculate_spectral_indices(model_sed, return_frequencies=False):

    """
    This function ...
    :param model_sed:
    :param return_frequencies:
    :return: 
    """

    # Define minimum wavelength
    min_wavelength = 50. * u("micron")

    # Get frequencies and flux densities
    frequencies = model_sed.frequencies("Hz", min_wavelength=min_wavelength, asarray=True)
    fluxdensities = model_sed.photometry("W / (m2 * Hz)", min_wavelength=min_wavelength, asarray=True)

    # Take logarithm
    log_frequencies = np.log10(frequencies)
    log_fluxdensities = np.log10(fluxdensities)

    # Calculate the derivative of the log(Flux) to the frequency, make a 'spectral index' function
    #gradients = np.gradient(log_frequencies, log_fluxdensities) # DOESN'T WORK ANYMORE?
    #gradients = np.diff(log_fluxdensities) / np.diff(log_frequencies)
    new_frequencies, gradients = numbers.derivatives(log_frequencies, log_fluxdensities)

    # Create function
    spectral_indices = interp1d(new_frequencies, gradients)

    # Return the spectral indices
    if return_frequencies: return spectral_indices, 10.**np.array(new_frequencies)
    else: return spectral_indices

# -----------------------------------------------------------------

def get_kbeam(fltr, spectral_indices, spire):

    """
    This function ...
    :param fltr:
    :param spectral_indices:
    :param spire:
    :return:
    """

    # Calculate the spectral index for the simulated SED at this filter
    # Use the central wavelength (frequency)
    central_frequency = fltr.center.to("Hz", equivalencies=spectral()).value
    central_log_frequency = np.log10(central_frequency)
    spectral_index_filter = spectral_indices(central_log_frequency)

    # Get the Kbeam factor
    kbeam = spire.get_kbeam_spectral(fltr, spectral_index_filter)

    # Return
    return kbeam

# -----------------------------------------------------------------

def calculate_fluxdensity_convolution(fltr, wavelengths, fluxdensities, wavelength_unit, fluxdensity_unit, spectral_indices, spire, errors=None,
                                      return_wavelength_grid=False, result_unit="Jy"):

    """
    This function ...
    :param fltr:
    :param wavelengths:
    :param fluxdensities:
    :param wavelength_unit:
    :param fluxdensity_unit:
    :param spectral_indices:
    :param spire:
    :param errors:
    :param return_wavelength_grid:
    :param result_unit:
    :return: 
    """

    # Debugging
    log.debug("Calculating the observed flux for the " + str(fltr) + " filter by convolving spectrally ...")

    # Check wavelength unit: wavelengths must be in MICRONS
    wavelength_conversion_factor = wavelength_unit.to("micron")
    wavelengths = wavelengths * wavelength_conversion_factor

    # Check fluxdensity unit
    if not fluxdensity_unit.is_wavelength_density: raise ValueError("Flux densities must be in wavelength spectral density unit")

    # Get result unit
    result_unit = u(result_unit)

    # Calculate the flux: flux densities must be per wavelength instead of per frequency!
    value, wavelength_grid = fltr.convolve(wavelengths, fluxdensities, return_grid=True)
    fluxdensity = float(value) * fluxdensity_unit
    fluxdensity_value = fluxdensity.to(result_unit, equivalencies=spectral_density(fltr.pivot)).value  # convert back to Jy

    # For SPIRE, also multiply with Kbeam correction factor
    if fltr.is_spire: fluxdensity_value *= get_kbeam(fltr, spectral_indices, spire)

    # Determine filter name
    filter_name = str(fltr)

    # Add a point to the mock SED
    error = errors[filter_name] if errors is not None and filter_name in errors else None

    # Return the fluxdensity and error
    if return_wavelength_grid: return fluxdensity_value * result_unit, error, wavelength_grid
    else: return fluxdensity_value * result_unit, error

# -----------------------------------------------------------------

def calculate_fluxdensity_closest(fltr, wavelengths, fluxdensities, wavelength_unit, fluxdensity_unit, errors=None,
                                  return_wavelength=False, result_unit="Jy"):

    """
    This function ...
    :param fltr: 
    :param wavelengths: 
    :param fluxdensities:
    :param wavelength_unit:
    :param fluxdensity_unit:
    :param errors:
    :param return_wavelength:
    :param result_unit:
    :return: 
    """

    # Debugging
    log.debug("Getting the observed flux for the " + str(fltr) + " filter ...")

    # Determine the filter wavelength
    filter_wavelength = fltr.wavelength.to(wavelength_unit).value  # fltr.pivot.to(wavelength_unit).value

    # Get the index of the wavelength closest to that of the filter
    index = sequences.find_closest_index(wavelengths, filter_wavelength)  # wavelengths are in micron

    # Get result unit
    result_unit = u(result_unit)

    # Get the flux density
    fluxdensity = fluxdensities[index] * fluxdensity_unit  # flux densities are in W/(m2 * micron)
    fluxdensity = fluxdensity.to(result_unit, equivalencies=spectral_density(fltr.pivot))  # now in Jy

    # Determine filter name
    filter_name = str(fltr)

    # Add a point to the mock SED
    error = errors[filter_name] if errors is not None and filter_name in errors else None

    # Get the actual wavelength
    wavelength = wavelengths[index] * wavelength_unit

    # Return the fluxdensity and error
    if return_wavelength: return fluxdensity, error, wavelength
    else: return fluxdensity, error

# -----------------------------------------------------------------

def create_mock_sed(model_sed, filters, spire, spectral_convolution=True, errors=None):

    """
    This function ...
    :param model_sed: 
    :param filters:
    :param spire:
    :param spectral_convolution:
    :param errors:
    :return: 
    """

    # Get the wavelengths and flux densities as arrays
    wavelengths, fluxdensities, wavelength_unit, fluxdensity_unit = get_wavelength_and_fluxdensity_arrays(model_sed, return_units=True)

    # Calculate the spectral indices for wavelengths in the dust regime
    spectral_indices, frequencies = calculate_spectral_indices(model_sed, return_frequencies=True)

    # Create an observed SED for the mock fluxes
    mock_sed = ObservedSED(photometry_unit="Jy")

    # Keep track of the wavelengths that have already been used to
    used_wavelengths = defaultdict(list)

    # Wavelengths used for each filter
    wavelengths_for_filters = OrderedDict()

    # Loop over the different filters
    for fltr in filters:

        # Check whether the filter is covered by the SED
        if not model_sed.covers(fltr.wavelength):
            log.warning("The wavelength of the '" + str(fltr) + "' is not covered by the modeled SED: skipping this filter ...")
            continue

        # Needs spectral convolution?
        if needs_spectral_convolution(fltr, spectral_convolution):

            # Calculate
            fluxdensity, error, wavelength_grid = calculate_fluxdensity_convolution(fltr, wavelengths, fluxdensities, wavelength_unit, fluxdensity_unit, spectral_indices, spire, errors=errors, return_wavelength_grid=True)

            # Create list of wavelengths
            filter_wavelengths = [value * Unit("micron") for value in wavelength_grid]

            # Add
            wavelengths_for_filters[fltr] = filter_wavelengths

        # No spectral convolution
        else:

            # Get the flux density
            fluxdensity, error, wavelength = calculate_fluxdensity_closest(fltr, wavelengths, fluxdensities, wavelength_unit, fluxdensity_unit, errors=errors, return_wavelength=True)

            # Check the wavelength
            wavelength_micron = wavelength.to("micron").value
            if wavelength_micron in used_wavelengths:

                filters = used_wavelengths[wavelength_micron]
                filter_names = [str(f) for f in filters]
                log.warning("The simulated flux for the wavelength '" + str(wavelength) + "' has already been used to create SED point(s) for the " + ", ".join(filter_names) + " filter(s)")

            # Add the filter for the wavelength
            used_wavelengths[wavelength_micron].append(fltr)

        # Add a data point to the mock SED
        mock_sed.add_point(fltr, fluxdensity, error)

    # Show which wavelengths are used to create filter frames
    log.debug("Used the following wavelengths of the simulated SED to create mock observed flux points without spectral convolution:")
    log.debug("")
    for wavelength_micron in used_wavelengths:
        filters = used_wavelengths[wavelength_micron]
        filter_names = [str(f) for f in filters]
        nfilters = len(filter_names)
        if nfilters == 1: log.debug(" - " + str(wavelength_micron) + " micron: " + filter_names[0])
        else: log.debug(" - " + str(wavelength_micron) + " micron: " + fmt.bold + ", ".join(filter_names) + fmt.reset)
    log.debug("")

    # Show which wavelengths are used to create filter frames
    log.debug("Used the following wavelengths for the spectral convolution for the other filters:")
    log.debug("")
    for fltr in wavelengths_for_filters:
        filter_name = str(fltr)
        wavelength_strings = [str(wavelength) for wavelength in wavelengths_for_filters[fltr]]
        log.debug(" - " + filter_name + ": " + ", ".join(wavelength_strings))
    log.debug("")

    # Return the mock observed SED
    return mock_sed

# -----------------------------------------------------------------

def needs_spectral_convolution(fltr, spectral_convolution):

    """
    This function ...
    :param fltr:
    :param spectral_convolution: flag or list of Filters
    :return:
    """

    # Broad band filter
    if isinstance(fltr, BroadBandFilter):

        # Single boolean
        if types.is_boolean_type(spectral_convolution): return spectral_convolution

        # Sequence of filters: return whether the filter is in it
        elif types.is_sequence(spectral_convolution): return fltr in spectral_convolution

        # Invalid
        else: raise ValueError("Invalid option for 'spectral_convolution'")

    # Narrow band filters: no spectral convolution
    else: return False

# -----------------------------------------------------------------
