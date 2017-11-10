#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.modeling.photometry.photometry Contains the PhotoMeter class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
import numpy as np

# Import the relevant PTS classes and modules
from .component import PhotometryComponent
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from ...core.data.sed import ObservedSED
from ...core.basics.errorbar import ErrorBar, sum_errorbars_quadratically
from ...core.plot.sed import SEDPlotter
from ...magic.convolution.aniano import AnianoKernels
from ...magic.convolution.psfs import HerschelPSFs
from ...magic.core.kernel import ConvolutionKernel
from ...core.launch.pts import PTSRemoteLauncher, launch_local
from ...magic.misc.calibration import CalibrationError
from ...magic.photometry.aperturenoise import ApertureNoiseCalculator
from ..preparation import unitconversion
from .tables import FluxErrorTable, FluxDifferencesTable
from ...core.basics.configuration import DictConfigurationSetter, ConfigurationDefinition
from ..preparation.preparer import has_statistics, load_statistics
from ...dustpedia.core.properties import has_calibration_error, get_calibration_error
from ...core.tools.utils import lazyproperty
from ...magic.core.mask import intersection
from ...magic.core.image import Image

# -----------------------------------------------------------------

# CHRIS:

# For SPIRE: Section 5.2.12 of the SPIRE Handbook (http://herschel.esac.esa.int/Docs/SPIRE/spire_handbook.pdf)
# provides the "official" way of performing aperture corrections; we specifically care about the "Extended source photometry
# (starting from extended source maps)" part. If you're using arbitrarily-large, non-circular apertures, you basically
# have to use maps of the beam profile to work out what fraction of the flux is outside your aperture.
# Those can be retrieved from the SPIRE calibration wiki
# (http://herschel.esac.esa.int/twiki/bin/view/Public/SpirePhotometerBeamProfile2).

# In the meantime, I'm working on sensible automated aperture-corrections - current plan is to assume that the underlying
# flux distribution follows a 2D-sersic profile, and so fit to each source a 2D-sersic convolved with the beam, and
# hence estimate the amount of flux that gets spread beyond the aperture. Hopefully it will be easily applicable
# to your work too.

# For PACS: The most up-to-date document on the PACS calibration wiki
# (http://herschel.esac.esa.int/twiki/pub/Public/PacsCalibrationWeb/bolopsf_22.pdf)
# and its accompanying tar.gz (ftp://ftp.sciops.esa.int/pub/hsc-calibration/PACS/PSF/PACSPSF_PICC-ME-TN-033_v2.2.tar.gz)
# give the most recent PACS beam profiles, and encircled energy fractions (EEFs) for different aperture sizes in the
# various scanning modes.

# PIETER:

# Zie bijlage voor de tabellen met de correcties voor de enclosed energy fraction (EEF) voor elke aperture size
# for PACS en SPIRE.
# Voor de aperture size moet ge mogelijks voor elke band een andere waarde gebruiken
# (door evt de beamsize in rekening te brengen ofzo).
# De flux in elke band moet dan gedeeld worden door de correction factor voor die band, gebruik
# makend van de aperture size in die band.

# De Growth_Curve_Final_XXmicron.dat geven de aperture size in arcsec en de correction factor voor de 2 PACS bands.

# De PSF_correction_HATLAS_SPIRE.dat geeft de aperture size in arcsec en dan de correction factors for
# 250um, 350um en 500um.

# Deze corrections zijn voor een centrale pointsource en zijn dus een soort minimum correctie voor een extended source.
# Deze minimum correcties worden doorgaands toegepast op extended sources omdat ze een goed genoege 1ste orde
# benadering zijn.

# -----------------------------------------------------------------

class PhotoMeter(PhotometryComponent):
    
    """
    This class...
    """

    def __init__(self, *args, **kwargs):

        """
        The constructor ...
        :param kwargs:
        :return:
        """

        # Call the constructor of the base class
        super(PhotoMeter, self).__init__(*args, **kwargs)

        # The list of images
        self.images = dict()

        # Masks
        self.truncation_masks = dict()
        self.clip_masks = dict()

        # The fluxes
        self.fluxes = dict()

        # The flux errors
        self.errors = dict()

        # The SED
        self.sed = None

        # The differences
        self.differences = None

        # Create the PTS remote launcher
        self.launcher = None

        # The instances that keep track of PSFs and convolution kernels
        self.aniano = AnianoKernels()
        self.herschel = HerschelPSFs()

        # The flux error table
        self.error_table = None

        # The flux differences table
        self.differences_table = None

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Load the images
        self.load_images()

        # Load masks
        self.load_masks()

        self.apply_masks()

        # 4. Calculate the fluxes
        self.calculate_fluxes()

        # 5. Calculate the errors
        self.calculate_errors()

        # 6. Make the SED
        self.make_sed()

        # 7.  Make the error table
        self.make_error_table()

        # 8. Calculate the differences between the calculated photometry and the reference SEDs
        self.calculate_differences()

        # 9. Writing
        self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(PhotoMeter, self).setup(**kwargs)

        # Create an observed SED
        self.sed = ObservedSED(photometry_unit="Jy")

        # Setup the remote PTS launcher
        if self.config.remote is not None:
            self.launcher = PTSRemoteLauncher()
            self.launcher.setup(self.config.remote)

        # Initialize the flux error table
        self.error_table = FluxErrorTable()

        # Initialize the flux differences table
        self.differences_table = FluxDifferencesTable(labels=self.reference_sed_labels)

    # -----------------------------------------------------------------

    def load_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Loading the images and error maps ...")

        # Load all data
        for name in self.dataset.names:

            # Debugging
            log.debug("Loading the " + name + " image ...")

            # Load the frame
            frame = self.dataset.get_frame(name)

            # Only broad band filters
            if not frame.is_broad_band: continue

            # Load the truncated frame
            #truncated_frame = self.dataset.get_frame(name)
            #plotting.plot_frame(frame)
            #plotting.plot_frame(truncated_frame)

            # Check whether preparation statistics can be found
            if not has_statistics(self.config.path, name): raise ValueError("Something went wrong in preparation: statistics not found")

            # Load the statistics
            statistics = load_statistics(self.config.path, name)

            mean_frame = statistics.mean_frame
            median_frame = statistics.median_frame
            stddev_frame = statistics.stddev_frame
            mean_frame_not_clipped = statistics.mean_frame_not_clipped
            median_frame_not_clipped = statistics.median_frame_not_clipped
            stddev_frame_not_clipped = statistics.stddev_frame_not_clipped
            mean_sky = statistics.mean_sky
            median_sky = statistics.median_sky
            mean_noise = statistics.mean_noise
            mean_subtracted = statistics.mean_subtracted
            median_subtracted = statistics.median_subtracted
            stddev_subtracted = statistics.stddev_subtracted

            #print(frame.unit)

            # Debugging
            log.debug("Checking the units of the image ...")

            # Convert to non- angular or intrinsic area unit
            if frame.is_per_angular_or_intrinsic_area: frame.convert_to_corresponding_non_angular_or_intrinsic_area_unit()

            # Add to the appropriate dictionary
            self.images[name] = frame #+ sky # add the sky to the frame

    # -----------------------------------------------------------------

    @lazyproperty
    def image_names(self):

        """
        This funtion ...
        :return:
        """

        return list(sorted(self.images.keys(), key=lambda name: self.images[name].filter.wavelength.to("micron").value))

    # -----------------------------------------------------------------

    def get_filter(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.images[name].filter

    # -----------------------------------------------------------------

    def get_wcs(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return self.images[name].wcs

    # -----------------------------------------------------------------

    def load_masks(self):

        """
        Thisn function ...
        :return:
        """

        # Inform the user
        log.info("Loading the masks ...")

        # Loop over the images
        for name in self.image_names:

            # Get frame and error map
            frame = self.images[name]
            errors = self.dataset.get_errormap(name)
            wcs = frame.wcs

            # Get the masks
            truncation_mask = self.get_truncation_mask(wcs)
            clip_mask = self.get_significance_mask(frame, errors)

            # Add the masks
            self.truncation_masks[name] = truncation_mask
            self.clip_masks[name] = clip_mask

    # -----------------------------------------------------------------

    def get_masks(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return [self.truncation_masks[name], self.clip_masks[name]]

    # -----------------------------------------------------------------

    def apply_masks(self):

        """
        Thisf unction ..
        :return:
        """

        for name in self.image_names:

            # Create combined mask
            mask = intersection(*self.get_masks(name))

            frame = self.images[name]

            # Set masked pixels to NaN
            frame.apply_mask_nans(mask)

            # Create image
            image = Image()
            image.add_frame(frame, "primary")
            image.add_mask(self.truncation_masks[name], "truncation")
            image.add_mask(self.clip_masks[name], "clip")

            # Set
            self.images[name] = image

    # -----------------------------------------------------------------

    def calculate_fluxes(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the aperture fluxes ...")

        # Loop over all the images
        for name in self.images:

            # Debugging
            log.debug("Calculating the total flux in the " + name + " image ...")

            # Calculate the total flux in Jansky
            flux = self.images[name].primary.sum() # * Unit("Jansky")

            # Apply correction for EEF of aperture
            #if "Pacs" in name or "SPIRE" in name: flux *= self.calculate_aperture_correction_factor(name)

            # Add the flux to the dictionary
            self.fluxes[name] = flux

    # -----------------------------------------------------------------

    def calculate_errors(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the errors on the aperture fluxes ...")

        # Loop over all the images
        for name in self.image_names:

            # Debugging
            log.debug("Calculating the flux error for the " + name + " image ...")

            # Calculate the calibration error
            #calibration_error = self.calculate_calibration_error(name)

            # Calculate the aperture noise error
            #aperture_noise = self.calculate_aperture_noise(name)

            # Set the error contributions
            #self.errors[name] = (calibration_error, aperture_noise)
            #self.errors[name] = calibration_error

            fltr = self.get_filter(name)

            # Get calibataion error as error bar
            error = get_calibration_error(fltr, self.fluxes[name], errorbar=True)

            # Add
            self.errors[name] = error

    # -----------------------------------------------------------------

    def make_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the observed SED ...")

        # Loop over all images
        for name in self.image_names:

            # Debugging
            log.debug("Adding SED entry for the " + name + " image ...")

            # The flux
            flux = self.fluxes[name]

            # The error contributions
            #calibration_error, aperture_noise = self.errors[name]

            # Add the error contributions
            #errorbar = sum_errorbars_quadratically(calibration_error, aperture_noise)

            errorbar = self.errors[name]

            # Add this entry to the SED
            self.sed.add_point(self.images[name].filter, flux, errorbar)

    # -----------------------------------------------------------------

    def make_error_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Making the flux error table ...")

        # Loop over the errors
        for name in self.images:

            # Debugging
            log.debug("Adding an entry to the error table for the " + name + " image ...")

            # Get the error contributions
            calibration_error, aperture_noise = self.errors[name]

            # Add the error contributions
            total_error = sum_errorbars_quadratically(calibration_error, aperture_noise)

            # Calculate the relative error
            total_relative_error = total_error / self.fluxes[name]

            # Add an entry to the flux error table
            self.error_table.add_entry(self.images[name].filter, calibration_error, aperture_noise, total_error, total_relative_error)

    # -----------------------------------------------------------------

    def calculate_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Calculating the differences with the reference SEDs ...")

        # Get list of instruments, bands and fluxes of the calculated SED
        filters = self.sed.filters()
        fluxes = self.sed.photometry(unit="Jy", add_unit=False)

        # The number of data points
        number_of_points = len(filters)

        # Loop over the different points in the calculated SED
        for i in range(number_of_points):

            # The dictionary with the flux differences for the different reference SEDs
            differences = dict()

            # The instrument and band
            instrument = filters[i].instrument
            band = filters[i].band

            # Loop over the different reference SEDs
            for label in self.reference_sed_labels:

                relative_difference = None

                # Loop over the data points in the reference SED
                for j in range(len(self.reference_seds[label].table["Wavelength"])):

                    if self.reference_seds[label].table["Instrument"][j] == instrument and self.reference_seds[label].table["Band"][j] == band:

                        difference = fluxes[i] - self.reference_seds[label].table["Flux"][j]
                        relative_difference = difference / fluxes[i] * 100.

                        # Break because a match has been found within this reference SED
                        break

                # Add percentage to the dictionary (or None if no match was found in this reference SED)
                differences[label] = relative_difference

            # Add entry to the table
            self.differences_table.add_entry(filters[i], fluxes[i], differences)

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the images
        self.write_images()

        # Write SED table
        self.write_sed()

        # Write the contributions to the flux errors
        self.write_error_table()

        # Write the differences
        self.write_differences()

        # Plot the SED
        self.plot_sed()

        # Plot the SED with references
        self.plot_sed_with_references()

    # -----------------------------------------------------------------

    def write_images(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the iamges ...")

        for name in self.image_names:

            # Detemrine path
            path = fs.join(self.phot_images_path, name + ".fits")

            # SAve
            self.images[name].saveto(path)

    # -----------------------------------------------------------------

    def write_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing SED to a data file ...")

        # Save the SED
        self.sed.saveto(self.observed_sed_path)

    # -----------------------------------------------------------------

    def write_error_table(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the table of the flux error contributions ...")

        # Write the table
        self.error_table.saveto(self.phot_errors_path)

    # -----------------------------------------------------------------

    def write_differences(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing the percentual differences with reference fluxes to a data file ...")

        # Save the differences table
        self.differences_table.saveto(self.phot_differences_path)

    # -----------------------------------------------------------------

    def plot_sed(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Add the SED
        plotter.add_sed(self.sed, "PTS")

        # Determine the full path to the plot file
        path = fs.join(self.phot_path, "sed.pdf")
        plotter.run(output=path, title=self.galaxy_name)

    # -----------------------------------------------------------------

    def plot_sed_with_references(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Plotting the SED with reference fluxes ...")

        # Create a new SEDPlotter instance
        plotter = SEDPlotter()

        # Add the SED
        plotter.add_sed(self.sed, "PTS")

        # Add the reference SEDs
        for label in self.reference_seds: plotter.add_sed(self.reference_seds[label], label)

        # Determine the full path to the plot file
        path = fs.join(self.phot_path, "sed_with_references.pdf")
        plotter.run(ouput=path, title=self.galaxy_name)

    # -----------------------------------------------------------------

    def calculate_aperture_noise(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.info("Calculating the aperture noise for " + name + " ...")

        # Get the frame (this is with the sky NOT subtracted)
        frame = self.images[name]

        truncation_ellipse_sky = self.truncation_ellipse
        truncation_ellipse_image = truncation_ellipse_sky.to_pixel(frame.wcs)

        # CONFIGURATION DICTIONARY
        config_dict = dict()
        plot_path = self.noise_path_for_image(name)
        fs.create_directory(plot_path)
        config_dict["plot_path"] = plot_path
        config_dict["debug_plotting"] = dict()
        config_dict["debug_plotting"]["intersection"] = False
        config_dict["debug_plotting"]["oversampled"] = False
        config_dict["debug_plotting"]["nans"] = False
        config_dict["debug_plotting"]["annulus_nans"] = True

        config_dict["method"] = self.config.noise_method

        # Configuration setter
        command_name = "calculate_aperture_noise"
        description = "calculate aperture noise"
        setter = DictConfigurationSetter(config_dict, command_name, description)

        # Create the configuration
        definition = ConfigurationDefinition()
        config = setter.run(definition)

        # Create the aperture noise calculator, configure it with the configuration
        calculator = ApertureNoiseCalculator(config=config)

        # Set the input
        input_dict = dict()
        input_dict["cutout"] = frame.data
        input_dict["band_name"] = name
        input_dict["adj_semimaj_pix"] = truncation_ellipse_image.radius.x
        input_dict["adj_axial_ratio"] = truncation_ellipse_image.radius.x / truncation_ellipse_image.radius.y
        input_dict["adj_angle"] = truncation_ellipse_image.angle.to("deg").value
        input_dict["centre_i"] = truncation_ellipse_image.center.y
        input_dict["centre_j"] = truncation_ellipse_image.center.x
        input_dict["downsample_factor"] = 1.0

        # Get the sky annulus ellipses
        annulus_inner = self.sky_annulus_inner(name)
        annulus_outer = self.sky_annulus_outer(name)

        annulus_inner_radius = annulus_inner.radius.x
        annulus_outer_radius = annulus_outer.radius.x

        # THIS IS JUST A FIX BECAUSE SOMETHING WENT WRONG WITH CREATING THE ANNULUS REGION (FOR SDSS e.g.)
        if annulus_inner_radius == annulus_outer_radius:
            annulus_inner_radius = 0.4 * annulus_outer_radius

        #input_dict["semimaj_pix_annulus_outer"] = annulus_outer_radius
        #input_dict["semimaj_pix_annulus_inner"] = annulus_inner_radius

        axratio_annulus_outer = annulus_outer.radius.x / annulus_outer.radius.y
        axratio_annulus_inner = annulus_inner.radius.x / annulus_inner.radius.y

        # Check
        if not np.isclose(axratio_annulus_outer, axratio_annulus_inner): log.error("DIFFERENCE AX RATIO", axratio_annulus_outer, axratio_annulus_inner)
        #input_dict["axial_ratio_annulus"] = axratio_annulus_outer

        annulus_angle_outer = annulus_outer.angle.to("deg").value
        annulus_angle_inner = annulus_inner.angle.to("deg").value

        # Check
        if not np.isclose(annulus_angle_outer, annulus_angle_inner): log.error("DIFFERENCE ANNULUS ANGLE", annulus_angle_outer, annulus_angle_inner)
        #input_dict["annulus_angle"] = annulus_angle_inner
        #input_dict["annulus_centre_i"] = annulus_outer.center.y
        #input_dict["annulus_centre_j"] = annulus_outer.center.x

        input_dict["annulus_inner_factor"] = 1.5
        input_dict["annulus_outer_factor"] = 2.0

        input_dict["plot_path"] = self.phot_temp_path

        # Run the aperture noise calculator
        calculator.run(**input_dict)

        # Get the noise, create the error bar
        error_bar = ErrorBar(calculator.noise)

        # Return the error bar
        return error_bar

    # -----------------------------------------------------------------

    def calculate_aperture_correction_factor(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        frame = self.images[name]

        filter_name = str(frame.filter)

        # Inform the user
        log.info("Calculating the aperture correction factor for " + filter_name + " ...")

        #####

        # INPUT DICTIONARY
        input_dict = dict()

        # CONFIGURATION DICTIONARY
        config_dict = dict()

        #####


        # Set cutout
        input_dict["cutout"] = frame


        config_dict["pix_arcsec"] = frame.average_pixelscale.to("arcsec").value


        truncation_ellipse_sky = self.truncation_ellipse
        truncation_ellipse_image = truncation_ellipse_sky.to_pixel(frame.wcs)


        # PACS BLUE
        if filter_name == "Pacs blue":

            ## USE psf.data !!
            psf_path = self.aniano.get_psf_path(self.reference_filter) # reference filter = pacs red filter
            psf = ConvolutionKernel.from_file(psf_path)

            annulus_inner = self.sky_annulus_inner("Pacs blue")
            annulus_outer = self.sky_annulus_outer("Pacs blue")

        # PACS GREEN
        elif filter_name == "Pacs green":

            psf_path = self.aniano.get_psf_path(self.reference_filter) # reference filter = pacs red filter
            psf = ConvolutionKernel.from_file(psf_path)

            annulus_inner = self.sky_annulus_inner("Pacs green")
            annulus_outer = self.sky_annulus_outer("Pacs green")

        # PACS RED
        elif filter_name == "Pacs red":

            psf_path = self.aniano.get_psf_path(self.reference_filter)  # reference filter = pacs red filter
            psf = ConvolutionKernel.from_file(psf_path)

            annulus_inner = self.sky_annulus_inner("Pacs red")
            annulus_outer = self.sky_annulus_outer("Pacs red")

        # SPIRE PSW
        elif filter_name == "SPIRE PSW":

            psf = self.herschel.get_spire_psf("PSW")

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # SPIRE PMW
        elif filter_name == "SPIRE PMW":

            psf = self.herschel.get_spire_psf("PMW")

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # SPIRE PLW
        elif filter_name == "SPIRE PLW":

            psf = self.herschel.get_spire_psf("PLW")

            annulus_inner = self.sky_annulus_inner(filter_name)
            annulus_outer = self.sky_annulus_outer(filter_name)

        # INVALID FILTER
        else: raise ValueError("Invalid filter: '" + filter_name + "'")

        # Debugging
        log.debug("Preparing the PSF kernel ...")

        # PREPARE THE CONVOLUTION KERNEL
        psf.prepare_for(frame)

        psf.saveto(fs.join(self.phot_temp_path, filter_name + "_psf.fits"))

        # SET INPUT DICT

        input_dict["psf"] = psf

        #annulus_inner_factor_x = annulus_inner.radius.x / truncation_ellipse_image.radius.x
        #annulus_outer_factor_x = annulus_outer.radius.x / truncation_ellipse_image.radius.x

        #annulus_inner_factor_y = annulus_inner.radius.y / truncation_ellipse_image.radius.y
        #annulus_outer_factor_y = annulus_outer.radius.y / truncation_ellipse_image.radius.y

        #if annulus_inner_factor_x != annulus_inner_factor_y: print("DIFFERENCE INNER", annulus_inner_factor_x, annulus_inner_factor_y)
        #if annulus_outer_factor_x != annulus_outer_factor_y: print("DIFFERENCE OUTER", annulus_outer_factor_x, annulus_outer_factor_y)

        config_dict["semimaj_pix"] = truncation_ellipse_image.radius.x
        config_dict["axial_ratio"] = truncation_ellipse_image.radius.x / truncation_ellipse_image.radius.y
        config_dict["angle"] = truncation_ellipse_image.angle.to("deg").value
        config_dict["centre_i"] = truncation_ellipse_image.center.y
        config_dict["centre_j"] = truncation_ellipse_image.center.x

        # ANNULUS PROPERTIES

        config_dict["semimaj_pix_annulus_outer"] = annulus_outer.radius.x
        config_dict["semimaj_pix_annulus_inner"] = annulus_inner.radius.x

        axratio_annulus_outer = annulus_outer.radius.x / annulus_outer.radius.y
        axratio_annulus_inner = annulus_inner.radius.x / annulus_inner.radius.y

        # Check
        if not np.isclose(axratio_annulus_outer, axratio_annulus_inner): print("DIFFERENCE AX RATIO", axratio_annulus_outer, axratio_annulus_inner)

        config_dict["axial_ratio_annulus"] = axratio_annulus_outer

        annulus_angle_outer = annulus_outer.angle.to("deg").value
        annulus_angle_inner = annulus_inner.angle.to("deg").value

        # Check
        if not np.isclose(annulus_angle_outer, annulus_angle_inner): print("DIFFERENCE ANNULUS ANGLE", annulus_angle_outer, annulus_angle_inner)

        config_dict["annulus_angle"] = annulus_angle_inner

        config_dict["annulus_centre_i"] = annulus_outer.center.y
        config_dict["annulus_centre_j"] = annulus_outer.center.x

        #config_dict["annulus_outer_factor"] = 2.0
        #config_dict["annulus_inner_factor"] = 1.5

        # SUBPIXEL FACTOR (with consideration of sub-pixels when factor > 1.0)
        config_dict["subpixel_factor"] = 1.0

        # Debugging
        log.debug("Performing the aperture correction calculation remotely ...")

        # Calculate the aperture correction factor
        if self.launcher is not None: factor = self.launcher.run_attached("aperture_correction", config_dict, input_dict, return_output_names=["factor"], unpack=True)
        else: factor = launch_local() # TODO

        # Debugging
        log.debug("The aperture correction factor for " + filter_name + " is " + repr(factor))

        # Return the correction factor
        return factor

# -----------------------------------------------------------------
