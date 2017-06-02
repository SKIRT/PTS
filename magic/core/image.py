#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       PTS -- Python Toolkit for working with SKIRT          **
# **       © Astronomical Observatory, Ghent University          **
# *****************************************************************

## \package pts.magic.core.image Contains the Image class.

# -----------------------------------------------------------------

# Ensure Python 3 functionality
from __future__ import absolute_import, division, print_function

# Import standard modules
import copy

# Import the relevant PTS classes and modules
from ..basics.layers import Layers
from ..region.list import PixelRegionList
from ..basics.mask import Mask
from .mask import Mask as newMask
from ...core.tools import filesystem as fs
from ...core.tools.logging import log
from .frame import Frame, sum_frames
from ...core.tools import types

# -----------------------------------------------------------------

class Image(object):

    """
    This class ...
    """

    # Set the default extension
    default_extension = "fits"

    # -----------------------------------------------------------------

    def __init__(self, name="untitled"):

        """
        The constructor ...
        :param name:
        :return:
        """

        # Initialize a set of layers to represent image frames, masks, segmentation maps and regions
        self.frames = Layers()
        self.masks = Layers()
        self.segments = Layers()
        self.regions = Layers()

        # The image name
        self.name = name

        # The original image header
        self.original_header = None

        # The dictionary containing meta information
        self.metadata = dict()

        # Temporary fix because fwhm is sometimes not transferred to a new primary Frame and therefore fwhm information is lost on the complete image
        self._fwhm = None

        # The image path
        self.path = None

    # -----------------------------------------------------------------

    @classmethod
    def from_file(cls, path, name=None, always_call_first_primary=True, hdulist_index=0, no_filter=False):

        """
        This function ...
        :param path:
        :param name:
        :param always_call_first_primary:
        :param hdulist_index:
        :param no_filter:
        :return:
        """

        # If no name is given, determine the name from the file path
        if name is None: name = fs.strip_extension(fs.name(path))

        # Create a new image
        image = cls(name)

        # Set the image path
        image.path = path

        # Load the image frames
        image.load_frames(path, always_call_first_primary=always_call_first_primary, hdulist_index=hdulist_index, no_filter=no_filter)

        # Return the image
        return image

    # -----------------------------------------------------------------

    def load_image(self, image, replace=True):

        """
        This function loads the frames, masks, regions and segmentation maps from another image into this image
        :param image:
        :param replace:
        :return:
        """

        # Remove the frames, masks, regions and segmentation maps from this image
        #image.remove_all_frames()
        #image.remove_all_masks()
        #image.remove_all_regions()
        #image.remove_all_segments()

        # Add the frames
        for label in image.frames: self.add_frame(image.frames[label], label, overwrite=replace)

        # Add the masks
        for label in image.masks: self.add_mask(image.masks[label], label, overwrite=replace)

        # Add the regions
        for label in image.regions: self.add_region(image.regions[label], label, overwrite=replace)

        # Add the segmentation maps
        for label in image.segments: self.add_segments(image.segments[label], label, overwrite=replace)

    # -----------------------------------------------------------------

    def copy(self):

        """
        This function ...
        :return:
        """

        return copy.deepcopy(self)

    # -----------------------------------------------------------------

    @property
    def primary(self):

        """
        This function ...
        :return:
        """

        return self.frames[0]

    # -----------------------------------------------------------------

    @property
    def has_errors(self):

        """
        This function ...
        :return:
        """

        return self.errors_name is not None

    # -----------------------------------------------------------------

    @property
    def errors_name(self):

        """
        This function returns the name of the error map frame
        :return:
        """

        frame_names = self.frames.keys()

        possible_error_map_names = ["errors", "error", "uncertainties"]

        if contains_two_or_more(frame_names, possible_error_map_names): raise RuntimeError("Image contains multiple error maps, don't know which to choose")

        # Loop over possible names, return frame if present
        for name in possible_error_map_names:
            if name in frame_names: return name

        # No error map found
        return None

    # -----------------------------------------------------------------

    @property
    def errors(self):

        """
        This function ...
        :return:
        """

        frame_name = self.errors_name
        if frame_name is None: return None

        return self.frames[frame_name]

    # -----------------------------------------------------------------

    @property
    def has_frames(self):

        """
        This function ...
        :return:
        """

        return self.nframes > 0

    # -----------------------------------------------------------------

    @property
    def has_regions(self):

        """
        This function ...
        :return:
        """

        return self.nregions > 0

    # -----------------------------------------------------------------

    @property
    def has_masks(self):

        """
        This function ...
        :return:
        """

        return self.nmasks > 0

    # -----------------------------------------------------------------

    @property
    def has_segments(self):

        """
        This function ...
        :return:
        """

        return self.nsegments > 0

    # -----------------------------------------------------------------

    @property
    def nframes(self):

        """
        This function ...
        :return:
        """

        return len(self.frames)

    # -----------------------------------------------------------------

    @property
    def nmasks(self):

        """
        This function ...
        :return:
        """

        return len(self.masks)

    # -----------------------------------------------------------------

    @property
    def nregions(self):

        """
        This function ...
        :return:
        """

        return len(self.regions)

    # -----------------------------------------------------------------

    @property
    def nsegments(self):

        """
        This function ...
        :return:
        """

        return len(self.segments)

    # -----------------------------------------------------------------

    @property
    def shape(self):

        """
        This function ...
        :return:
        """

        return self.primary.shape if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def xsize(self):

        """
        This function ...
        :return:
        """

        return self.primary.xsize if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def ysize(self):

        """
        This function ...
        :return:
        """

        return self.primary.ysize if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def filter(self):

        """
        This function ...
        :return:
        """

        # Return the filter of the primary frame
        return self.primary.filter if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def filter_name(self):

        """
        This function ...
        :return: 
        """

        return str(self.filter)

    # -----------------------------------------------------------------

    @property
    def wavelength(self):

        """
        This function ...
        :return:
        """

        # Return the wavelength of the primary frame
        return self.primary.wavelength if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def unit(self):

        """
        This function ...
        :return:
        """

        # Return the unit of the primary frame
        return self.primary.unit if self.has_frames else None

    # -----------------------------------------------------------------

    @unit.setter
    def unit(self, value):

        """
        This function ...
        :return:
        """

        # Set the unit of all frames
        for frame_name in self.frames: self.frames[frame_name].unit = value

    # -----------------------------------------------------------------

    @property
    def pixelscale(self):

        """
        This function ...
        :return:
        """

        # Return the pixelscale of the primary frame
        return self.primary.pixelscale if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def average_pixelscale(self):

        """
        This function ...
        :return:
        """

        # Return the averaged pixelscale of the primary frame
        return self.primary.average_pixelscale if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def pixelarea(self):

        """
        This function ...
        :return:
        """

        # Return the pixelarea of the primary frame
        return self.primary.pixelarea if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def fwhm(self):

        """
        This function ...
        :return:
        """

        if self._fwhm is not None: return self._fwhm

        # Return the FWHM of the primary frame
        return self.primary.fwhm if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def fwhm_pix(self):

        """
        This function ...
        :return:
        """

        if not self.has_frames: return None

        # Get the FWHM in sky coordinates
        fwhm = self.fwhm

        # Convert into pixels
        return (fwhm / self.average_pixelscale).to("").value

    # -----------------------------------------------------------------

    @property
    def wcs(self):

        """
        This function ...
        :return:
        """

        # Return the wcs of the primary frame
        return self.primary.wcs if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def coordinate_range(self):

        """
        This function ...
        :return:
        """

        # Return the coordinate range of the primary frame
        return self.primary.coordinate_range if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def description(self):

        """
        This function ...
        :return:
        """

        return self.primary.description

    # -----------------------------------------------------------------

    @description.setter
    def description(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.primary.description = value

    # -----------------------------------------------------------------

    @property
    def zero_point(self):

        """
        This function ...
        :return:
        """

        return self.primary.zero_point

    # -----------------------------------------------------------------

    @property
    def source_extracted(self):

        """
        This function ...
        :return:
        """

        return self.primary.source_extracted

    # -----------------------------------------------------------------

    @source_extracted.setter
    def source_extracted(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        self.primary.source_extracted = value

    # -----------------------------------------------------------------

    @property
    def extinction_corrected(self):

        """
        This function ...
        :return:
        """

        return self.primary.extinction_corrected

    # -----------------------------------------------------------------

    @extinction_corrected.setter
    def extinction_corrected(self, value):

        """
        This function ...
        :return:
        """

        self.primary.extinction_corrected = value

    # -----------------------------------------------------------------

    @property
    def sky_subtracted(self):

        """
        This function ...
        :return:
        """

        return self.primary.sky_subtracted

    # -----------------------------------------------------------------

    @sky_subtracted.setter
    def sky_subtracted(self, value):

        """
        This function ...
        :return:
        """

        self.primary.sky_subtracted = value

    # -----------------------------------------------------------------

    def __repr__(self):

        """
        This function ...
        :return:
        """

        return "<" + self.__class__.__name__ + " '" + self.name + "' with " + str(self.nframes) + " frames, " \
               + str(self.nregions) + " regions, " + str(self.nmasks) + " masks, and " + str(self.nsegments) \
               + " segmentation maps>"

    # -----------------------------------------------------------------

    def save(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Saving the image ...")

        # Check whether the path is defined
        if self.path is None: raise RuntimeError("Path is not defined for this image")

        # Save the image
        self.saveto(self.path)

    # -----------------------------------------------------------------

    def saveto(self, path, add_metadata=False, origin=None, add_masks=True, add_segments=True, add_regions=False):

        """
        This function exports the image (frames and masks) as a datacube into FITS file.
        :param path:
        :param add_metadata:
        :param origin:
        :param add_masks:
        :param add_segments:
        :param add_regions:
        :return:
        """

        # Create an array to contain the data cube
        datacube = []

        # Names of the planes
        plane_names = []

        plane_index = 0

        # Create a header from the wcs
        if self.wcs is not None: header = self.wcs.to_header() # Create a header from the coordinate system
        else: header = fits.Header() # Construct a new header

        # Get the names of the frames
        frame_names = self.frames.keys()

        last_addition = None

        # Export all frames to the specified file
        for frame_name in frame_names:

            # Inform the user that this frame will be saved to the image file
            log.debug("Exporting the " + frame_name + " frame to " + path)

            # Check if the coordinate system of this frame matches that of the other frames ?

            # Add this frame to the data cube, if its coordinates match those of the primary frame
            datacube.append(self.frames[frame_name]._data)
            plane_name = frame_name + " [frame]"
            plane_names.append(plane_name)

            # Add the name of the frame to the header
            header["PLANE" + str(plane_index)] = plane_name

            # Increment the plane index
            plane_index += 1

            last_addition = "Frame"

        # Add the masks
        if add_masks:

            # Export all masks to the specified file
            for mask_name in self.masks:

                # Inform the user that this mask will be saved to the image file
                log.debug("Exporting the " + mask_name + " mask to " + path)

                # Add this mask to the data cube
                datacube.append(self.masks[mask_name].astype(int))
                plane_name = mask_name + " [mask]"
                plane_names.append(plane_name)

                # Add the name of the mask to the header
                header["PLANE" + str(plane_index)] = plane_name

                # Increment the plane index
                plane_index += 1

                last_addition = "Mask"

        # Add the segmentation maps
        if add_segments:

            # Export all segmentation maps
            for segments_name in self.segments:

                # Inform the user that this segmentation map will be saved to the image file
                log.debug("Exporting the " + segments_name + " segmentation map to " + path)

                # Add this segmentation map to the data cube
                datacube.append(self.segments[segments_name].data)
                plane_names = segments_name + " [segments]"
                plane_names.append(plane_name)

                # Add the name of the segmentation map to the header
                header["PLANE" + str(plane_index)] = plane_name

                # Increment the plane index
                plane_index += 1

                last_addition = "SegmentationMap"

        #if add_regions:

            # http://docs.astropy.org/en/stable/io/fits/
            #tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='target', format='20A', array=a1), fits.Column(name='V_mag', format='E', array=a2)])

        # Add the meta information to the header
        if add_metadata:
            for key in self.metadata:
                try: header[key] = self.metadata[key]
                except ValueError: pass # Some values in the header gives errors in Astropy when adding them again to this new header ... (e.g. ValueError: Illegal value: = 'created by T.H. Jarrett'.)

        # Set plane information
        if plane_index > 1:
            header["NAXIS"] = 3
            header["NAXIS3"] = plane_index
        else:
            header.remove("PLANE0")
            header["PTSCLS"] = last_addition  # if only one Frame or Mask or SegmentationMap has been added

        # Set unit, FWHM and filter description
        if self.unit is not None: header.set("SIGUNIT", str(self.unit), "Unit of the map")
        if self.fwhm is not None: header.set("FWHM", self.fwhm.to("arcsec").value, "[arcsec] FWHM of the PSF")
        if self.filter is not None: header.set("FILTER", str(self.filter), "Filter used for this observation")
        if self.wcs is None and self.pixelscale is not None:
            header.set("XPIXSIZE", repr(self.pixelscale.x.to("arcsec").value), "[arcsec] Pixelscale for x axis")
            header.set("YPIXSIZE", repr(self.pixelscale.y.to("arcsec").value), "[arcsec] Pixelscale for y axis")

        # Add origin description
        if origin is not None: header["ORIGIN"] = origin
        else: header["ORIGIN"] = "Image class of PTS package"

        # FITS format
        if path.endswith(".fits"):

            # Import
            from . import fits as pts_fits

            if len(datacube) == 1: datacube = datacube[0]

            # Write
            pts_fits.write_datacube(datacube, header, path)

        # ASDF format
        elif path.endswith(".asdf"):

            # Import
            from asdf import AsdfFile

            # Create the tree
            tree = dict()

            if len(datacube) == 1:

                plane = datacube[0]
                name = plane_names[0].split(" [")[0]

                # Set the plane
                tree[name] = plane

            else:

                for i in range(len(datacube)):
                    tree[plane_names[i]] = datacube[i]
                tree["header"] = header

            # Create the asdf file
            ff = AsdfFile(tree)

            # Write
            ff.write_to(path)

        # Only FITS or ASDF format is allowed
        else: raise ValueError("Only the FITS or ASDF filetypes are supported")

        # Update the path
        self.path = path

    # -----------------------------------------------------------------

    def show_layers(self):

        """
        This function ...
        :return:
        """

        if self.has_frames: log.info("Frames:")

        frame_index = 0
        for name in self.frames:

            # Show
            log.info("[" + str(frame_index) + "] " + name)

            # Increment
            frame_index += 1

        if self.has_regions: log.info("Regions:")

        region_index = 0
        for name in self.regions:

            # Show
            log.info("[" + str(region_index) + "] " + name)

            # Increment
            region_index += 1

        if self.has_masks: log.info("Masks:")

        mask_index = 0
        for name in self.masks:

            # Show
            log.info("[" + str(mask_index) + "] " + name)

            # Increment
            mask_index += 1

        if self.has_segments: log.info("Segmentation maps:")

        segments_index = 0
        for name in self.segments:

            # Show
            log.info("[" + str(segments_index) + "] " + name)

            # Increment
            segments_index += 1

    # -----------------------------------------------------------------

    def apply_mask(self, mask, fill=0.0):

        """
        This function ..
        :param mask:
        :param fill:
        """

        # Replace the masked pixels in all frames by the fill value
        for frame_name in self.frames: self.frames[frame_name][mask] = fill

    # -----------------------------------------------------------------

    def import_region(self, path, name, overwrite=False):

        """
        This function imports a new region from a DS9 region file
        :param path:
        :param name:
        :param overwrite:
        :return:
        """

        # Check if the region file exists
        if not fs.is_file(path): raise IOError("The region file does not exist")

        # Create an Region object from the regions file
        region = PixelRegionList.from_file(path)

        # Add the region to the set of regions
        self.add_region(region, name, overwrite)

    # -----------------------------------------------------------------

    def get_state(self):

        """
        This function ...
        :return:
        """

        # Create an empty dictionary to contain the state of the current image
        state = dict()

        # Loop over all frames, regions and masks and record whether they are selected
        for frame_name in self.frames: state["frames/"+frame_name] = self.frames[frame_name].selected
        for region_name in self.regions: state["regions/"+region_name] = self.regions[region_name].selected
        for mask_name in self.masks: state["masks/"+mask_name] = self.masks[mask_name].selected

        # Return the state dictionary
        return state

    # -----------------------------------------------------------------

    def set_state(self, state):

        """
        This function ...
        :param state:
        :return:
        """

        # Loop over the entries in the state dictionary
        for identifier, selected in state.items():

            # Split the layer identifier into the layer type and the actual name of that layer
            layer_type, name = identifier.split("/")

            # Set the appropriate flag
            if layer_type == "frames": self.frames[name].selected = selected
            elif layer_type == "regions": self.regions[name].selected = selected
            elif layer_type == "masks": self.masks[name].selected = selected
            else: raise ValueError("Invalid state dictionary")

    # -----------------------------------------------------------------

    @fwhm.setter
    def fwhm(self, fwhm):

        """
        This function ...
        :param fwhm:
        :return:
        """

        self._fwhm = fwhm

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Setting the FWHM of the " + frame_name + " frame to " + str(fwhm) + " ...")

            # Set the unit for this frame
            self.frames[frame_name].fwhm = fwhm

    # -----------------------------------------------------------------

    @filter.setter
    def filter(self, fltr):

        """
        This function ...
        :param fltr:
        :return:
        """

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Setting the filter of the " + frame_name + " frame to " + str(fltr) + " ...")

            # Set the filter for this frame
            self.frames[frame_name].filter = fltr

    # -----------------------------------------------------------------

    @wcs.setter
    def wcs(self, wcs):

        """
        This function ...
        :param wcs:
        :return:
        """

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Setting the coordinate system of the " + frame_name + " frame ...")

            # Set the wcs for this frame
            self.frames[frame_name].wcs = wcs

    # -----------------------------------------------------------------

    @property
    def has_wcs(self):

        """
        This function ...
        :return:
        """

        return self.wcs is not None

    # -----------------------------------------------------------------

    def convert_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False, brightness_strict=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        """

        # Inform the user
        log.debug("Converting the unit of the image from " + str(self.unit) + " to " + str(to_unit) + " ...")

        # Calculate the conversion factor
        factor = self.unit.conversion_factor(to_unit, density=density, fltr=self.filter, pixelscale=self.pixelscale,
                                             distance=distance, brightness=brightness, density_strict=density_strict,
                                             brightness_strict=brightness_strict)

        # Debug message
        log.debug("Conversion factor = " + str(factor))

        # Multiply the image with the conversion factor
        self.__imul__(factor)

        # Set the new unit
        self.unit = to_unit

        # Return the conversion factor
        return factor

    # -----------------------------------------------------------------

    def convolve(self, kernel, allow_huge=True, fft=True):

        """
        This function ...
        :param kernel: of type ConvolutionKernel
        :param allow_huge:
        :param fft:
        """

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Convolving the " + frame_name + " frame ...")

            # Convolve this frame
            self.frames[frame_name].convolve(kernel, allow_huge=allow_huge, fft=fft)

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        """

        # Check whether the image has a WCS
        if self.has_wcs: raise RuntimeError("Cannot rebin a frame without coordinate system")

        # Create a copy of the current wcs
        original_wcs = self.wcs.deepcopy()

        footprint = None

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Rebinning the " + frame_name + " frame ...")

            # Rebin this frame (the reference wcs is automatically set in the new frame)
            footprint = self.frames[frame_name].rebin(reference_wcs, exact=exact, parallel=parallel)

        # Loop over the masks
        for mask_name in self.masks:

            # Inform the user
            log.debug("Rebinning the " + mask_name + " mask ...")

            # Create a frame for the mask
            mask_frame = Frame(self.masks[mask_name].astype(float), wcs=original_wcs)

            # Rebin the mask frame
            footprint = mask_frame.rebin(reference_wcs, exact=exact, parallel=parallel)

            # Return the rebinned mask
            # data, name, description
            self.masks[mask_name] = Mask(mask_frame > 0.5, name=self.masks[mask_name].name, description=self.masks[mask_name].description)

        # If there was any frame or mask, we have footprint
        if footprint is not None:

            # Add mask for padded pixels after rebinning
            # this mask now covers pixels added to the frame after rebinning plus more (radius 10 pixels)
            padded = Mask(footprint < 0.9).disk_dilation(radius=10)
            self.add_mask(padded, "padded")

    # -----------------------------------------------------------------

    def crop(self, x_min, x_max, y_min, y_max):

        """
        This function ...
        :param x_min:
        :param x_max:
        :param y_min:
        :param y_max:
        """

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Cropping the " + frame_name + " frame ...")

            # Crop this frame
            self.frames[frame_name].crop(x_min, x_max, y_min, y_max)

        # Loop over all masks
        for mask_name in self.masks:

            # Inform the user
            log.debug("Cropping the " + mask_name + " mask ...")

            # Crop this mask
            self.masks[mask_name] = self.masks[mask_name][y_min:y_max, x_min:x_max]

        # Loop over all regions
        for region_name in self.regions:

            # Inform the user
            log.debug("Cropping the " + region_name + " region ...")

            # Crop the region
            self.regions[region_name] = self.regions[region_name].cropped(x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    def intersect_masks(self):

        """
        This function ...
        :return:
        """

        # Calculate the total mask
        mask = newMask.intersection(*self.masks.values())

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def unite_masks(self):

        """
        This function ...
        :return:
        """

        # Calculate the total mask
        mask = newMask.union(*self.masks.values())

        # Return the mask
        return mask

    # -----------------------------------------------------------------

    def delete_frame(self, frame_name):

        """
        This function removes the specified frame
        :param frame_name:
        :return:
        """

        # Inform the user
        log.debug("Deleting the " + frame_name + " frame ...")

        # Remove this frame from the frames dictionary
        del self.frames[frame_name]

    # -----------------------------------------------------------------

    def delete_region(self, region_name):

        """
        This function removes the specified region
        :param region_name:
        :return:
        """

        # Inform the user
        log.debug("Deleting the " + region_name + " region ...")

        # Remove this region from the regions dictionary
        del self.regions[region_name]

    # -----------------------------------------------------------------

    def delete_mask(self, mask_name):

        """
        This function removes the specified mask
        :param mask_name:
        :return:
        """

        # Inform the user
        log.debug("Deleting the " + mask_name + " mask ...")

        # Remove this mask from the masks dictionary
        del self.masks[mask_name]

    # -----------------------------------------------------------------

    def __mul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__imul__(value)

    # -----------------------------------------------------------------

    def __imul__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Multiplying the " + frame_name + " frame by a factor of " + str(value) + " ...")

            # Multiply the frame by the given factor
            self.frames[frame_name] *= value

        # Return a reference to this instance
        return self

    # -----------------------------------------------------------------

    def __div__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.copy().__idiv__(value)

    # -----------------------------------------------------------------

    def __idiv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            log.debug("Dividing the " + frame_name + " frame by a factor of " + str(value) + " ...")

            # Divide the frame by the given factor
            self.frames[frame_name] /= value

        # Return a reference to this instance
        return self

    # -----------------------------------------------------------------

    def __truediv__(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        return self.__div__(value)

    # -----------------------------------------------------------------

    def __itruediv__(self, factor):

        """
        This function ...
        :param factor:
        :return:
        """

        return self.__idiv__(factor)

    # -----------------------------------------------------------------

    def __abs__(self):

        """
        This function ...
        :return:
        """

        new = self.copy()

        # Loop over all frames
        for frame_name in self.frames:

            # Get absolute values of frame
            new.frames[frame_name] = abs(self.frames[frame_name])

        # Return the new image
        return new

    # -----------------------------------------------------------------

    def __setitem__(self, item, value):

        """
        This function ...
        :param item:
        :param value:
        :return:
        """

        # Loop over the frames, and set the values for each frame
        for frame_name in self.frames: self.frames[frame_name][item] = value

    # -----------------------------------------------------------------

    def sum_frames(self):

        """
        This function ...
        :return:
        """

        # Inform the user
        log.info("Summing all frames ...")

        frames = []

        # Loop over the frames and add them to the list
        for frame_name in self.frames:

            # Debugging
            log.debug("Adding the '" + frame_name + "' frame ...")
            frames.append(self.frames[frame_name])

        # Sum the frames
        result = sum_frames(*frames)

        # Return the resulting frame
        return result

    # -----------------------------------------------------------------

    def rename_frame(self, old_name, new_name, keep_position=False):

        """
        This function renames a frame
        :param old_name:
        :param new_name:
        :param keep_position:
        :return:
        """

        # Get the current position
        if keep_position: current_position = self.frames.keys().index(old_name)
        else: current_position = None

        # Remove the frame from the dictionary of frames and re-add it under a different key
        self.frames[new_name] = self.frames.pop(old_name)

        # Move the frame to the original position
        if keep_position: self.move_frame(new_name, current_position)

    # -----------------------------------------------------------------

    def rename_region(self, old_name, new_name):

        """
        This function renames a region
        :param old_name:
        :param new_name:
        :return:
        """

        # Remove the region of the dictionary of regions and re-add it under a different key
        self.regions[new_name] = self.regions.pop(old_name)

    # -----------------------------------------------------------------

    def rename_mask(self, old_name, new_name):

        """
        This function ...
        :param old_name:
        :param new_name:
        :return:
        """

        # Remove the mask from the dictionary of masks and re-add it under a different key
        self.masks[new_name] = self.masks.pop(old_name)

    # -----------------------------------------------------------------

    def load_frames(self, path, index=None, name=None, description=None, always_call_first_primary=True, rebin_to_wcs=False, hdulist_index=0, no_filter=False, silent=False):

        """
        This function ...
        :param path:
        :param index:
        :param name:
        :param description:
        :param always_call_first_primary:
        :param rebin_to_wcs:
        :param hdulist_index:
        :param no_filter:
        :param silent:
        :return:
        """

        # Check if the file exists
        if not fs.is_file(path): raise IOError("File '" + path + "' does not exist")

        # Show which image we are importing
        if not silent: log.debug("Reading in file '" + path + "' ...")

        # Load frames
        from . import fits as pts_fits
        frames, masks, segments, meta = pts_fits.load_frames(path, index, name, description, always_call_first_primary,
                                                       rebin_to_wcs, hdulist_index, no_filter)

        # Set frames, masks and meta information
        for frame_name in frames: self.add_frame(frames[frame_name], frame_name)
        for mask_name in masks: self.add_mask(masks[mask_name], mask_name)
        for segments_name in segments: self.add_segments(segments[segments_name], segments_name)
        for keyword in meta: self.metadata[keyword] = meta[keyword]

    # -----------------------------------------------------------------

    def add_frame(self, frame, name, overwrite=False):

        """
        This function ...
        :param frame:
        :param name:
        :param overwrite:
        :return:
        """

        # Inform the user
        log.debug("Adding '" + name + "' to the set of frames ...")

        # Check whether a frame with this name already exists
        if name in self.frames and not overwrite: raise RuntimeError("A frame with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if frame.shape != self.shape: raise ValueError("Frame does not have the correct shape for this image")

        # Set the WCS
        if self.wcs is not None: frame.wcs = self.wcs

        # Add the layer to the layers dictionary
        self.frames[name] = frame

    # -----------------------------------------------------------------

    def remove_frames_except(self, *names):

        """
        This function ...
        :param names:
        :return:
        """

        # Loop over all frames
        for frame_name in list(self.frames.keys()):

            # Don't remove the frame with the specified name
            if frame_name in names: continue

            # Remove all other frames
            self.remove_frame(frame_name)

    # -----------------------------------------------------------------

    def remove_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.debug("Removing the '" + name + "' frame ...")

        # Check whether a frame with this name exists
        if name not in self.frames: raise RuntimeError("A frame with this name does not exist")

        # Delete the frame
        del self.frames[name]

    # -----------------------------------------------------------------

    def remove_all_frames(self):

        """
        This function ...
        :return:
        """

        # Loop over all frames
        for frame_name in list(self.frames.keys()): self.remove_frame(frame_name)

    # -----------------------------------------------------------------

    def replace_frame(self, name, frame):

        """
        This function ...
        :param name:
        :param frame:
        :return:
        """

        # Inform the user
        log.debug("Replacing the '" + name + "' frame ...")

        # Check whether a frame with this name exists
        if name not in self.frames: raise RuntimeError("A frame with this name does not exist")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if frame.shape != self.shape: raise ValueError("Frame does not have the correct shape for this image")

        # Set the WCS
        if self.wcs is not None: frame.wcs = self.wcs

        # Set the frame
        self.frames[name] = frame

    # -----------------------------------------------------------------

    def move_frame(self, name, position):

        """
        This function ...
        :param name:
        :param position:
        :return:
        """

        # Inform the user
        log.debug("Moving frame '" + name + "' to position " + str(position) + " in the frame layers ...")

        # Check whether a frame with this name exists
        if name not in self.frames: raise RuntimeError("A frame with this name does not exist")

        # Get the current position of the frame
        current_position = self.frames.keys().index(name)

        # Check if the current position is already the same as the requested position
        if current_position == position: return

        # Make an iterable from the indices of all the other frames
        other_indices = iter([index for index in range(len(self.frames)) if index != current_position])

        # Create a new layers instance
        reordered_frames = Layers()

        # Loop over the other indices
        while len(reordered_frames) < len(self.frames):

            # Check the number of frames present in the new dictionary
            current_nframes = len(reordered_frames)

            # If the number of current frames corresponds to the requested frame position, add the frame
            if current_nframes == position: reordered_frames[name] = self.frames[name]

            # Fill the rest of the reordered frames up with the other frames, in the order defined by other_indices
            else:

                index = next(other_indices)
                this_name = self.frames.keys()[index]
                reordered_frames[this_name] = self.frames[this_name]

        # Set the frames layers
        self.frames = reordered_frames

    # -----------------------------------------------------------------

    def add_region(self, region, name, overwrite=False):

        """
        This function ...
        :param region:
        :param name:
        :param overwrite:
        :return:
        """

        # Inform the user
        log.debug("Adding '" + name + "' to the set of regions ...")

        # Check whether a region with this name already exists
        if name in self.regions and not overwrite: raise RuntimeError("A region with this name already exists")

        # Add the region to the set of regions
        self.regions[name] = region

    # -----------------------------------------------------------------

    def remove_regions_except(self, names):

        """
        This function ...
        :param names:
        :return:
        """

        if isinstance(names, basestring): names = [names]

        # Loop over all regions
        for region_name in list(self.regions.keys()):

            # Don't remove the region with the specified name
            if region_name in names: continue

            # Remove all other regions
            self.remove_region(region_name)

    # -----------------------------------------------------------------

    def remove_region(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.debug("Removing the '" + name + "' region ...")

        # Check whether a region with this name exists
        if name not in self.regions: raise RuntimeError("A region with this name does not exist")

        # Delete the region
        del self.regions[name]

    # -----------------------------------------------------------------

    def remove_all_regions(self):

        """
        This function ...
        :return:
        """

        # Loop over all regions
        for region_name in list(self.regions.keys()): self.remove_region(region_name)

    # -----------------------------------------------------------------

    def add_mask(self, mask, name, overwrite=False):

        """
        This function ...
        :param mask:
        :param name:
        :param overwrite:
        :return:
        """

        # Inform the user
        log.debug("Adding '" + name + "' to the set of masks ...")

        # Check whether a mask with this name already exists
        if name in self.masks and not overwrite: raise RuntimeError("A mask with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if mask.shape != self.shape: raise ValueError("Mask does not have the correct shape for this image")

        # Add the mask to the set of masks
        self.masks[name] = mask

    # -----------------------------------------------------------------

    def remove_masks_except(self, names):

        """
        This function ...
        :param names:
        :return:
        """

        if isinstance(names, basestring): names = [names]

        # Loop over all masks
        for mask_name in list(self.masks.keys()):

            # Don't remove the mask with the specified name
            if mask_name in names: continue

            # Remove all other masks
            self.remove_mask(mask_name)

    # -----------------------------------------------------------------

    def remove_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.debug("Removing the '" + name + "' mask ...")

        # Check whether a mask with this name exists
        if name not in self.masks: raise RuntimeError("A mask with this name does not exist")

        # Delete the mask
        del self.masks[name]

    # -----------------------------------------------------------------

    def remove_all_masks(self):

        """
        This function ...
        :return:
        """

        # Loop over all masks
        for mask_name in list(self.masks.keys()): self.remove_mask(mask_name)

    # -----------------------------------------------------------------

    def add_segments(self, segments, name, overwrite=False):

        """
        This function ...
        :param segments:
        :param name:
        :param overwrite:
        :return:
        """

        # Inform the user
        log.debug("Adding '" + name + "' to the set of segmentation maps ...")

        # Check whether a segmentation map with this name already exists
        if name in self.segments and not overwrite: raise RuntimeError("A segmentation map with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if segments.shape != self.shape: raise ValueError("Segmentation map does not have the correct shape for this image")

        # Add the segmentation map to the set of segmentation maps
        self.segments[name] = segments

    # -----------------------------------------------------------------

    def remove_segments_except(self, names):

        """
        This function ...
        :param names:
        :return:
        """

        if isinstance(names, basestring): names = [names]

        # Loop over all segmentation maps
        for segments_name in list(self.segments.keys()):

            # Don't remove the segmentation map with the specified name
            if segments_name in names: continue

            # Remove all other segmentation maps
            self.remove_segments(segments_name)

    # -----------------------------------------------------------------

    def remove_segments(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        # Inform the user
        log.debug("Removing the '" + name + "' segmentation map ...")

        # Check whether a segmentation map with this name exists
        if name not in self.segments: raise RuntimeError("A segmentation map with this name does not exist")

        # Delete the segmentation map
        del self.segments[name]

    # -----------------------------------------------------------------

    def remove_all_segments(self):

        """
        This function ...
        :return:
        """

        # Loop over all segmentation maps
        for segments_name in list(self.segments.keys()): self.remove_segments(segments_name)

# -----------------------------------------------------------------

def ordered_dict_prepend(dct, key, value, dict_setitem=dict.__setitem__):

    """
    This function ...
    :param dct:
    :param key:
    :param value:
    :param dict_setitem:
    :return:
    """

    root = dct._OrderedDict__root
    first = root[1]

    if key in dct:
        link = dct._OrderedDict__map[key]
        link_prev, link_next, _ = link
        link_prev[1] = link_next
        link_next[0] = link_prev
        link[0] = root
        link[1] = first
        root[1] = first[0] = link
    else:
        root[1] = first[0] = dct._OrderedDict__map[key] = [root, first, key]
        dict_setitem(dct, key, value)

# -----------------------------------------------------------------

def contains_two_or_more(in_list, from_list):

    """
    This function ...
    :param in_list:
    :param from_list:
    :return:
    """

    counter = 0

    for item in in_list:

        if item in from_list: counter += 1

        if counter == 2: return True

    return False

# -----------------------------------------------------------------
