#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.prepare.batch Contains the BatchImagePreparer class.

# -----------------------------------------------------------------

# Ensure Python 3 compatibility
from __future__ import absolute_import, division, print_function

# Import standard modules
from multiprocessing import Pool

# Import astronomical modules
from astropy.units import Unit

# Import the relevant PTS classes and modules
from ...core.tools.logging import log
from ...core.tools import filesystem as fs
from ...core.basics.filter import Filter
from ...core.basics.configurable import Configurable
from ..misc.calibration import CalibrationError
from ..misc.extinction import GalacticExtinction
from ..core.frame import Frame
from ..core.dataset import DataSet
from ..region.list import SkyRegionList
from ...magic.misc.kernels import aniano_names, AnianoKernels

# -----------------------------------------------------------------

class BatchImagePreparer(Configurable):
    
    """
    This class ...
    """
    
    def __init__(self, config=None):
        
        """
        The constructor ...
        """

        # Call the constructor of the base class
        super(BatchImagePreparer, self).__init__(config)

        # The frames
        self.frames = dict()

        # The errormaps
        self.errormaps = dict()

        # The attenuations
        self.attenuations = dict()

        # Initialize the process pool
        self.pool = Pool(processes=self.config.nprocesses)

        # Rebinning wcs and convolution filter
        self.rebinning_wcs = None
        self.convolution_filter = None
        self.convolution_fwhm = None

        # The Aniano kernel service
        self.aniano = AnianoKernels()

        # Don't convolve or rebin
        self.dont_convolve = []
        self.dont_rebin = []

        # Set the principal ellipse and saturation region in sky coordinates
        #self.image_preparer.principal_ellipse_sky = regions.largest_ellipse(galaxy_region).to_sky(self.image.wcs)
        #self.image_preparer.saturation_region_sky = saturation_region.to_sky(
        #    self.image.wcs) if saturation_region is not None else None

        self.principal_ellipses_sky = dict()
        self.saturation_regions_sky = dict()

    # -----------------------------------------------------------------

    def add_frame(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Check if name not already used
        if name in self.frames: raise ValueError("Already a frame with the name " + name)

        # Set the frame
        self.frames[name] = frame

    # -----------------------------------------------------------------

    @property
    def min_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.frames:
            if fwhm is None or self.frames[name].fwhm < fwhm: fwhm = self.frames[name].fwhm

        # Return the minimum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_fwhm(self):

        """
        This function ...
        :return:
        """

        fwhm = None

        # Loop over the images
        for name in self.frames:
            if fwhm is None or self.frames[name].fwhm > fwhm: fwhm = self.frames[name].fwhm

        # Return the maximum FWHM
        return fwhm

    # -----------------------------------------------------------------

    @property
    def max_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.frames:

            wcs = self.frames[name].wcs
            if pixelscale is None or wcs.average_pixelscale > pixelscale: pixelscale = wcs.average_pixelscale

        # Return the maximum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def min_pixelscale(self):

        """
        This function ...
        :return:
        """

        pixelscale = None

        # Loop over the images
        for name in self.frames:

            wcs = self.frames[name].wcs
            if pixelscale is None or wcs.average_pixelscale < pixelscale: pixelscale = wcs.average_pixelscale

        # Return the minimum pixelscale
        return pixelscale

    # -----------------------------------------------------------------

    @property
    def bounding_box(self):

        """
        This function ...
        :return:
        """

        # Region of all the bounding boxes
        boxes_region = SkyRegionList()

        # Add the bounding boxes as sky rectangles
        for name in self.frames: boxes_region.append(self.frames[name].wcs.bounding_box)

        # Return the bounding box of the region of rectangles
        return boxes_region.bounding_box

    # -----------------------------------------------------------------

    @property
    def center_coordinate(self):

        """
        This function ...
        :return:
        """

        return self.bounding_box.center

    # -----------------------------------------------------------------

    def run(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # 1. Call the setup function
        self.setup(**kwargs)

        # 2. Get the galactic attenuation
        self.get_extinction()

        # 3. Prepare the images
        self.prepare()

        # 4. Writing
        if self.config.write: self.write()

    # -----------------------------------------------------------------

    def setup(self, **kwargs):

        """
        This function ...
        :param kwargs:
        :return:
        """

        # Call the setup function of the base class
        super(BatchImagePreparer, self).setup()

        # Load the frames (from config or input kwargs)
        if "frames" in kwargs: self.frames = kwargs.pop("frames")
        elif "dataset" in kwargs: self.frames = kwargs.pop("dataset").get_frames()
        else: self.load_frames()

        # Set rebinning wcs
        self.set_rebinning_wcs()

        # Set convolution filter
        self.set_convolution_filter()

    # -----------------------------------------------------------------

    def set_rebinning_wcs(self):

        """
        This function ...
        :return:
        """

        max_pixelscale_wcs = None

        # Loop over the frames
        for name in self.frames:

            # Get the pixelscale
            pixelscale = self.frames[name].average_pixelscale

            # Ignore this frame if its pixelscale is larger than the minimum pixelscale that is specified by the user
            if self.config.max_pixelscale is not None and pixelscale > self.config.max_pixelscale:
                self.dont_rebin.append(name)
                continue

            # Check if the pixelscale is greater than that of the previous frame (but still below the maximum)
            if max_pixelscale_wcs is None or pixelscale > max_pixelscale_wcs.average_pixelscale: max_pixelscale_wcs = self.frames[name].wcs

        # Now set the rebinning wcs
        self.rebinning_wcs = max_pixelscale_wcs

    # -----------------------------------------------------------------

    def set_convolution_filter(self):

        """
        This function ...
        :return:
        """

        max_fwhm = None
        max_fwhm_filter = None

        # Loop over the frames
        for name in self.frames:

            # Get the FWHM
            fwhm = self.frames[name].fwhm

            # Ignore this frame if its FWHM is higher than the maximum FWHM specified by the user
            if self.config.max_fwhm is not None and fwhm > self.config.max_fwhm:
                self.dont_convolve.append(name)
                continue

            # Check if the FWHM is greater than that of the previous frame (but still below the maximum)
            if max_fwhm is None or fwhm > max_fwhm:
                max_fwhm = fwhm
                max_fwhm_filter = self.frames[name].filter

        # Set the convolution filter and FWHM
        self.convolution_filter = max_fwhm_filter
        self.convolution_fwhm = max_fwhm

    # -----------------------------------------------------------------

    def load_frames(self):

        """
        This function ...
        :return:
        """

        # Create new dataset
        if self.config.dataset.endswith(".fits"):

            # Load the frame
            frame = Frame.from_file(self.config.dataset)

            # Determine the name for this image
            name = str(frame.filter)

            # Add the frame
            self.add_frame(name, frame)

        # Load dataset from file
        elif self.config.dataset.endswith(".dat"):

            # Get the dataset
            dataset = DataSet.from_file(self.config.dataset)

            # Get the frames
            self.frames = dataset.get_frames()

        # Invalid value for 'dataset'
        else: raise ValueError("Parameter 'dataset' must be filename of a dataset file (.dat) or a FITS file (.fits)")

    # -----------------------------------------------------------------

    def get_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Getting the galactic extinction ...")

        # Create the galactic extinction calculator
        extinction = GalacticExtinction(self.center_coordinate)

        # Loop over the frames
        for label in self.frames:

            # Get the filter
            fltr = self.frames[label].filter

            # Get the extinction
            self.attenuations[label] = extinction.extinction_for_filter(fltr)

    # -----------------------------------------------------------------

    def prepare(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Preparing the images ...")

        # 1. Extract stars and galaxies from the image
        if self.config.extract_sources: self.extract_sources()

        # 2. If requested, correct for galactic extinction
        if self.config.correct_for_extinction: self.correct_for_extinction()

        # 3. If requested, convert the unit
        if self.config.convert_units: self.convert_units()

        # 4. If requested, convolve
        if self.config.convolve: self.convolve()

        # 5. Rebin
        if self.config.rebin: self.rebin()

        # 6. Subtract the sky
        if self.config.subtract_sky: self.subtract_sky()

        # 7. Calculate the calibration uncertainties
        if self.config.cself.calculate_calibration_uncertainties()

        # 8. Calculate the error maps
        self.calculate_errormaps()

    # -----------------------------------------------------------------

    def extract_sources(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Extracting sources ...")

        # Loop over the images
        for label in self.frames:

            # Check if source extraction is required
            if not self._needs_extinction_correction(label): continue

            # Execute
            result = self.pool.apply_async(_extract_sources, args=(,))  # All simple types (strings) ?

            # Add the result
            results.append(result)

            # Get and set the dust masses

        #self.dust_masses = [result.get() for result in results]

        # Close and join the process pool
        self.pool.close()
        self.pool.join()

    # -----------------------------------------------------------------

    def __needs_source_extraction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return not self.frames[label].source_extracted

    # -----------------------------------------------------------------

    def correct_for_extinction(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Correcting for galactic extinction ...")

        # Loop over the images
        for label in self.frames:

            # Check whether extinction correction is required
            if not self._needs_extinction_correction(label): continue

            # Do the correction

            corrected_path = fs.join(output_path, "corrected_for_extinction.fits")

    # -----------------------------------------------------------------

    def _needs_extinction_correction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return not self.frames[label].extinction_corrected

    # -----------------------------------------------------------------

    def convert_units(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Converting the units ...")

        for label in self.frames:

            # Check whether unit conversion if necessary
            if not self._needs_unit_conversion(label): continue

            converted_path = fs.join(output_path, "converted_unit.fits")

            assert self.image.unit == Unit("Jy/pix")

            # Get pixelscale
            pixelscale = self.image.average_pixelscale

            # Conversion from Jy / pix to MJy / pix
            conversion_factor = 1e-6

            # Conversion from MJy / pix to MJy / sr
            conversion_factor *= (1.0 / pixelscale ** 2).to("pix2/srpa").value

            # Multiply the image with the conversion factor
            self.image *= conversion_factor

            assert self.config.unit_conversion.to_unit == Unit("MJy/sr")

            # Set the new unit
            self.image.unit = self.config.unit_conversion.to_unit

            # Write intermediate result
            if self.config.write_steps: self.write_intermediate_result("converted_unit.fits")

            # Inform the user
            log.success("Units converted")

    # -----------------------------------------------------------------

    def _needs_unit_conversion(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        return self.frames[label].unit != self.config.unit_conversion.to_unit

    # -----------------------------------------------------------------

    def convolve(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Convolving to a common resolution ...")

        # Loop over the frames
        for label in self.frames:

            #convolved_path = fs.join(output_path, "convolved.fits")

            # Check if convolution is necessary
            if not self._needs_convolution(label): continue

            # Filters
            from_filter = self.frames[label].filter
            to_filter = self.convolution_filter

            # FWHMs
            from_fwhm = self.frames[label].fwhm
            to_fwhm = self.convolution_fwhm

            # Get the kernel path
            kernel_file_path = self.aniano.get_kernel_path(from_filter, to_filter, from_fwhm=from_fwhm, to_fwhm=to_fwhm)

            # Set the kernel path and FWHM
            #self.image_preparer.config.convolution.kernel_path = kernel_file_path  # set kernel path
            #self.image_preparer.config.convolution.kernel_fwhm = self.reference_fwhm  # set kernel FWHM (is a quantity here)

    # -----------------------------------------------------------------

    def _needs_convolution(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # If this frame does not need to be convolved
        if label in self.dont_convolve: return False

        # Check the resolution
        if self.frames[label].fwhm == self.convolution_fwhm: return False

        # Else, return True
        return True

    # -----------------------------------------------------------------

    def rebin(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Rebinning to a common pixel grid ...")

        # 5. If requested, rebin
        # rebin

        for label in self.frames:

            #rebinned_path = fs.join(output_path, "rebinned.fits")

            # Check if rebinning is necessary
            if not self._needs_rebinning(label): continue

            #if label == self.config.rebinning_reference: continue



            #if fs.is_file(rebinned_path): continue

    # -----------------------------------------------------------------

    def _needs_rebinning(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        # If this frame does not need to be rebinned
        if label in self.dont_rebin: return False

        # Check the coordinate system
        if self.frames[label].wcs == self.rebinning_wcs: return False

        # Else, return True
        return True

    # -----------------------------------------------------------------

    def subtract_sky(self):

        """
        This function ...
        :return:
        """

        # 6. If requested, subtract the sky
        # subtract_sky

        # Loop over the frames
        for label in self.frames:

            # Check whether the image has to be sky subtracted
            if not self._needs_sky_subtraction(label): continue



            #subtracted_path = fs.join(output_path, "sky_subtracted.fits")

    # -----------------------------------------------------------------

    def _needs_sky_subtraction(self, label):

        """
        This function ...
        :param label:
        :return:
        """

        #if image.frames.primary.sky_subtracted:
        #    log.debug("The " + image.name + " image has already been sky subtracted")
        #    self.image_preparer.config.subtract_sky = False
        #else: self.image_preparer.config.subtract_sky = True  # Has yet to be sky subtracted

        return not self.frames[label].sky_subtracted

    # -----------------------------------------------------------------

    def calculate_calibration_uncertainties(self):

        """
        This function ...
        :return:
        """


        for label in self.frames:

            calibration_error = CalibrationError.from_filter(self.frames[label].filter)

        # 7. Calculate the calibration uncertainties
        # calculate_calibration_uncertainties

    # -----------------------------------------------------------------

    def calculate_errormaps(self):

        """
        This function ...
        :return:
        """

        # 8. If requested, set the uncertainties
        # set_uncertainties

    # -----------------------------------------------------------------

    def write(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Writing ...")

        # Write the frames
        self.write_frames()

        # Write the new dataset
        self.write_dataset()

    # -----------------------------------------------------------------

    def write_frames(self):

        """
        This function ...
        :return:
        """

        pass

    # -----------------------------------------------------------------

    def write_dataset(self):

        """
        This function ...
        :return:
        """

        pass

# -----------------------------------------------------------------

def _extract_sources():

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Extracting the sources ...")

    # Create an animation to show the result of this step
    if self.visualisation_path is not None:

        # Create an animation
        animation = ImageBlinkAnimation()
        animation.add_image(self.image.frames.primary)

        # Create an animation to show the progress of the SourceExtractor
        source_extractor_animation = SourceExtractionAnimation(self.image.frames.primary)
    else:
        animation = None
        source_extractor_animation = None

    # Run the extractor
    self.extractor.run(self.image.frames.primary, self.galaxy_region, self.star_region, self.saturation_region,
                       self.other_region, self.galaxy_segments, self.star_segments, self.other_segments,
                       source_extractor_animation)

    # Write the animation
    if self.visualisation_path is not None:

        # Determine the path to the animation
        path = fs.join(self.visualisation_path, time.unique_name(self.image.name + "_sourceextraction") + ".gif")

        # Save the animation
        source_extractor_animation.save(path)

        # ...
        animation.add_image(self.image.frames.primary)

        # Determine the path to the animation
        path = fs.join(self.visualisation_path,
                       time.unique_name(self.image.name + "_imagepreparation_extractsources") + ".gif")

        # Save the animation
        animation.save(path)

    # Add the sources mask to the image
    self.image.add_mask(self.extractor.mask, "sources")

    # Get the principal shape in sky coordinates
    self.principal_shape_sky = self.extractor.principal_shape.to_sky(self.image.wcs)

    # Get the saturation region in sky coordinates
    self.saturation_region_sky = self.saturation_region.to_sky(
        self.image.wcs) if self.saturation_region is not None else None

    # Write intermediate result
    if self.config.write_steps: self.write_intermediate_result("extracted.fits")

    # Inform the user
    log.success("Sources extracted")

# -----------------------------------------------------------------

def _convolve(self):

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Convolving the image with kernel " + self.config.convolution.kernel_path + " ...")

    # Create an animation to show the result of this step
    if self.visualisation_path is not None:

        # Create an animation
        animation = ImageBlinkAnimation()
        animation.add_image(self.image.frames.primary)

    else:
        animation = None

    # Open the convolution kernel
    kernel = ConvolutionKernel.from_file(self.config.convolution.kernel_path, fwhm=self.config.convolution.kernel_fwhm)

    # Prepare the kernel
    kernel.prepare_for(self.image)

    # Save the kernel frame for manual inspection
    if self.config.write_steps:
        kernel_path = self.full_output_path("kernel.fits")
        kernel.save(kernel_path)

    # Check whether the convolution has to be performed remotely
    if self.config.convolution.remote is not None:

        # Inform the user
        log.info("Convolution will be performed remotely on host '" + self.config.convolution.remote + "' ...")

        # Create remote image, convolve and make local again
        remote_image = RemoteImage.from_local(self.image, self.config.convolution.remote)
        remote_image.convolve(kernel, allow_huge=True)
        self.image = remote_image.to_local()

    # The convolution is performed locally. # Convolve the image (the primary and errors frame)
    else:
        self.image.convolve(kernel, allow_huge=True)

    # Save convolved frame
    if self.config.write_steps: self.write_intermediate_result("convolved.fits")

    # Write the animation
    if self.visualisation_path is not None:

        animation.add_image(self.image.frames.primary)

        # Determine the path to the animation
        path = fs.join(self.visualisation_path,
                       time.unique_name(self.image.name + "_imagepreparation_convolve") + ".gif")

        # Save the animation
        animation.save(path)

    # Inform the user
    log.success("Convolution finished")

# -----------------------------------------------------------------

def _subtract_sky():

    """
    This function ...
    :return:
    """

    # Inform the user
    log.info("Subtracting the sky ...")

    # Create an animation to show the result of this step
    if self.visualisation_path is not None:

        # Create an animation to show the result of the sky subtraction step
        animation = ImageBlinkAnimation()
        animation.add_image(self.image.frames.primary)

        # Create an animation to show the progress of the SkySubtractor
        skysubtractor_animation = Animation()
    else:
        animation = None
        skysubtractor_animation = None

    # Convert the principal ellipse in sky coordinates into pixel coordinates
    principal_shape = self.principal_shape_sky.to_pixel(self.image.wcs)

    # Convert the saturation region in sky coordinates into pixel coordinates
    if self.saturation_region_sky is not None:
        saturation_region = self.saturation_region_sky.to_pixel(self.image.wcs)
    else: saturation_region = None

    # Create the 'extra' mask (bad and padded pixels)
    extra_mask = None
    if "bad" in self.image.masks:
        # Combine with padded mask or just use bad mask
        if "padded" in self.image.masks:
            extra_mask = self.image.masks.bad + self.image.masks.padded
        else:
            extra_mask = self.image.masks.bad
    elif "padded" in self.image.masks: extra_mask = self.image.masks.padded  # just use padded mask

    # Run the sky subtractor
    self.sky_subtractor.run(self.image.frames.primary, principal_shape, self.image.masks.sources, extra_mask, saturation_region, skysubtractor_animation)

    # Add the sky frame to the image
    self.image.add_frame(self.sky_subtractor.sky_frame, "sky")

    # Add the mask that is used for the sky estimation
    self.image.add_mask(self.sky_subtractor.mask, "sky")

    # Add the sky noise frame
    self.image.add_frame(self.sky_subtractor.noise_frame, "sky_errors")

    # Write intermediate result if requested
    if self.config.write_steps: self.write_intermediate_result("sky_subtracted.fits")

    # Write sky annuli maps if requested
    if self.config.write_sky_apertures:

        # Write the sky region
        region_path = fs.join(self.config.sky_apertures_path, "annulus.reg")
        self.sky_subtractor.region.save(region_path)

        # Write the apertures frame
        apertures_frame_path = fs.join(self.config.sky_apertures_path, "apertures.fits")
        self.sky_subtractor.apertures_frame.save(apertures_frame_path)

        # Write the apertures mean frame
        apertures_mean_path = fs.join(self.config.sky_apertures_path, "apertures_mean.fits")
        self.sky_subtractor.apertures_mean_frame.save(apertures_mean_path)

        # Write the apertures noise frame
        apertures_noise_path = fs.join(self.config.sky_apertures_path, "apertures_noise.fits")
        self.sky_subtractor.apertures_noise_frame.save(apertures_noise_path)

    # Write the animation
    if self.visualisation_path is not None:

        # Determine the path to the animation
        path = fs.join(self.visualisation_path, time.unique_name(self.image.name + "_skysubtraction") + ".gif")

        # Save the animation
        skysubtractor_animation.save(path)

        # ...
        animation.add_image(self.image.frames.primary)

        # Determine the path to the animation
        path = fs.join(self.visualisation_path,
                       time.unique_name(self.image.name + "_imagepreparation_subtractsky") + ".gif")

        # Save the animation
        animation.save(path)

    # Inform the user
    log.success("Sky subtracted")

# -----------------------------------------------------------------
