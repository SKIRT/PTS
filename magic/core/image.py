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
import numpy as np

# Import astronomical modules
from astropy.io.fits import Header

# Import the relevant PTS classes and modules
from ..basics.layers import Layers
from ..region.list import PixelRegionList
from .mask import Mask, union, intersection
from ...core.tools import filesystem as fs
from ...core.basics.log import log
from .frame import sum_frames, Frame
from ...core.tools.stringify import tostr
from ...core.units.unit import PhotometricUnit
from ...core.units.stringify import represent_unit
from ...core.units.unit import parse_unit as u
from ...core.tools import introspection

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
    def from_file(cls, path, name=None, always_call_first_primary=True, hdulist_index=0, no_filter=False, density=False,
                  brightness=False, density_strict=False, brightness_strict=False, indices=None, absolute_index_names=True):

        """
        This function ...
        :param path:
        :param name:
        :param always_call_first_primary:
        :param hdulist_index:
        :param no_filter:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param indices:
        :param absolute_index_names:
        :return:
        """

        # If no name is given, determine the name from the file path
        if name is None: name = fs.strip_extension(fs.name(path))

        # Create a new image
        image = cls(name)

        # Set the image path
        image.path = path

        # Load the image frames
        image.load_frames(path, always_call_first_primary=always_call_first_primary, hdulist_index=hdulist_index,
                          no_filter=no_filter, density=density, brightness=brightness, density_strict=density_strict,
                          brightness_strict=brightness_strict, indices=indices, absolute_index_names=absolute_index_names)

        # Return the image
        return image

    # -----------------------------------------------------------------

    @classmethod
    def from_frame(cls, frame, name=None):

        """
        This function ...
        :param frame:
        :param name:
        :return:
        """

        # Check type
        if not isinstance(frame, Frame): raise ValueError("Argument is not a frame")

        # Set name
        if name is None: name = frame.name

        # Create image
        if name is None: image = cls()
        else: image = cls(name=name)

        # Add frame
        image.add_frame(frame, name="primary")

        # Return the image
        return image

    # -----------------------------------------------------------------

    @classmethod
    def from_frames(cls, *args, **kwargs):

        """
        This function ...
        :param args:
        :param kwargs:
        :return:
        """

        # Give warning
        if introspection.is_python2(): log.warning("In python 2, the order of the frames in the image is not determined by the order of the keyword arguments, but effectively random")

        # Check primary frame
        if len(args) == 0: image = Image()
        elif len(args) == 1: image = Image.from_frame(args[0])
        else: raise ValueError("Cannot pass more than one primary frame, specifiy other frames as keyword arguments")

        # Add other frames
        for name in kwargs:
            frame = kwargs[name]
            if not isinstance(frame, Frame): raise ValueError("Not a frame: '" + name + "'")
            image.add_frame(frame, name=name)

        # Return the image
        return image

    # -----------------------------------------------------------------

    @classmethod
    def from_3d_array(cls, array, split_axis=-1, name=None, wcs=None, names=None):

        """
        This function ...
        :param array:
        :param split_axis:
        :param name:
        :param wcs:
        :param names:
        :return:
        """

        # Move axis to first position
        array = np.moveaxis(array, split_axis, 0)
        nframes = array.shape[0]

        # Create image
        if name is None: image = cls()
        else: image = cls(name=name)

        # Add the frames
        for index in range(nframes):

            # Are names defined?
            if names is not None: frame_name = names[index]
            else: frame_name = None

            # Create the frame
            print(array[index])
            frame = Frame(array[index], wcs=wcs, name=frame_name)

            # Add the frame
            if frame_name is None: frame_name = "frame" + str(index)
            image.add_frame(frame, frame_name)

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
    def primary_name(self):

        """
        This function ...
        :return:
        """

        return self.frame_names[0]

    # -----------------------------------------------------------------

    def is_primary(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name == self.primary_name

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
    def frame_names(self):

        """
        This function ...
        :return:
        """

        return self.frames.keys()

    # -----------------------------------------------------------------

    def has_frame(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.frame_names

    # -----------------------------------------------------------------

    @property
    def mask_names(self):

        """
        This function ...
        :return:
        """

        return self.masks.keys()

    # -----------------------------------------------------------------

    def has_mask(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.mask_names

    # -----------------------------------------------------------------

    @property
    def regions_names(self):

        """
        This function ...
        :return:
        """

        return self.regions.keys()

    # -----------------------------------------------------------------

    def has_region_list(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.regions_names

    # -----------------------------------------------------------------

    @property
    def segments_names(self):

        """
        This function ...
        :return:
        """

        return self.segments.keys()

    # -----------------------------------------------------------------

    def has_segmentation_map(self, name):

        """
        This function ...
        :param name:
        :return:
        """

        return name in self.segments_names

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
    def has_unit(self):

        """
        This function ...
        :return:
        """

        return self.unit is not None

    # -----------------------------------------------------------------

    @property
    def is_photometric(self):

        """
        This function ...
        :return:
        """

        if not self.has_frames: raise ValueError("No frames: no unit")
        return self.primary.is_photometric

    # -----------------------------------------------------------------

    @property
    def distance(self):

        """
        This function ...
        :return:
        """

        return self.primary.distance if self.has_frames else None

    # -----------------------------------------------------------------

    @distance.setter
    def distance(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        for frame_name in self.frames: self.frames[frame_name].distance = value

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

    @pixelscale.setter
    def pixelscale(self, value):

        """
        This function ...
        :param value:
        :return:
        """

        # Set the pixelscale of all frames
        for frame_name in self.frames: self.frames[frame_name].pixelscale = value

    # -----------------------------------------------------------------

    @property
    def has_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.pixelscale is not None

    # -----------------------------------------------------------------

    @property
    def has_angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_frames: raise ValueError("No frames: no pixelscale")
        return self.primary.has_angular_pixelscale

    # -----------------------------------------------------------------

    @property
    def has_physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        if not self.has_frames: raise ValueError("No frames: no pixelscale")
        return self.primary.has_physical_pixelscale

    # -----------------------------------------------------------------

    @property
    def angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.primary.angular_pixelscale if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def physical_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.primary.physical_pixelscale if self.has_frames else None

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
    def average_angular_pixelscale(self):

        """
        This function ...
        :return:
        """

        return self.primary.average_angular_pixelscale if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def average_physical_pixelscale(self):

        """
        Thisf unction ...
        :return:
        """

        return self.primary.average_physical_pixelscale if self.has_frames else None

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
    def angular_pixelarea(self):

        """
        This function ...
        :return:
        """

        return self.primary.angular_pixelarea if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def pixel_solid_angle(self):

        """
        This function ...
        :return:
        """

        return self.angular_pixelarea

    # -----------------------------------------------------------------

    @property
    def physical_pixelarea(self):

        """
        This function ...
        :return:
        """

        return self.primary.physical_pixelarea if self.has_frames else None

    # -----------------------------------------------------------------

    @property
    def pixel_area(self):

        """
        This function ...
        :return:
        """

        return self.physical_pixelarea

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

        if self.fwhm is None: return None
        #if not self.has_frames: return None

        # Get the FWHM in sky coordinates
        fwhm = self.fwhm

        # Convert into pixels
        return (fwhm / self.average_pixelscale).to("").value

    # -----------------------------------------------------------------

    @property
    def sigma(self):

        """
        This function ...
        :return:
        """

        if self.fwhm is None: return None
        #if not self.has_frames: return None

        from ..tools import statistics
        return self.fwhm * statistics.fwhm_to_sigma if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def sigma_pix(self):

        """
        This function ...
        :return:
        """

        if self.fwhm is None: return None
        #if not self.has_frames: return None

        from ..tools import statistics
        return self.fwhm_pix * statistics.fwhm_to_sigma if self.fwhm is not None else None

    # -----------------------------------------------------------------

    @property
    def psf_filter(self):

        """
        This function ...
        :return:
        """

        return self.primary.psf_filter if self.has_frames else None

    # -----------------------------------------------------------------

    @psf_filter.setter
    def psf_filter(self, value):

        """
        Thisfunction ...
        :param value:
        :return:
        """

        # Loop over all frames
        for frame_name in self.frames:

            # Inform the user
            # TOO MUCH OUTPUT
            #log.debug("Setting the PSF filter of the " + frame_name + " frame ...")

            # Set the PSF filter for this frame
            self.frames[frame_name].psf_filter = value

    # -----------------------------------------------------------------

    @property
    def psf_filter_name(self):

        """
        This function ...
        :return:
        """

        return str(self.psf_filter) if self.psf_filter is not None else None

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

    @property
    def file_size(self):

        """
        This function ...
        :return:
        """

        if self.path is None: raise ValueError("Path is not defined")
        return fs.file_size(self.path)

    # -----------------------------------------------------------------

    @property
    def data_size(self):

        """
        This function ...
        :return:
        """

        total_nbytes = 0.

        # Loop over all currently selected frames
        for frame_name in self.frames: total_nbytes += self.frames[frame_name].data.nbytes

        # Loop over the masks
        for mask_name in self.masks: total_nbytes += self.masks[mask_name].data.nbytes

        # Loop over the segmentation maps
        for segments_name in self.segments: total_nbytes += self.segments[segments_name].data.nbytes

        # Return the size in GB
        return (total_nbytes * u("byte")).to("GB")

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

    def saveto(self, path, add_metadata=False, origin=None, add_masks=True, add_segments=True, add_regions=False,
               extra_info=None, add_plane_names=True, add_filter=True, add_psf_filter=True):

        """
        This function exports the image (frames and masks) as a datacube into FITS file.
        :param path:
        :param add_metadata:
        :param origin:
        :param add_masks:
        :param add_segments:
        :param add_regions:
        :param extra_info:
        :param add_plane_names:
        :param add_filter:
        :return:
        """

        # Create an array to contain the data cube
        datacube = []

        # Names of the planes
        plane_names = []

        plane_index = 0

        # Create a header from the wcs
        if self.wcs is not None: header = self.wcs.to_header() # Create a header from the coordinate system
        else: header = Header() # Construct a new header

        # Get the names of the frames
        frame_names = self.frames.keys()

        last_addition = None

        # Export all frames to the specified file
        for frame_name in frame_names:

            # Inform the user that this frame will be saved to the image file
            #log.debug("Exporting the " + frame_name + " frame to '" + path + "' ...")

            # Check if the coordinate system of this frame matches that of the other frames ?

            # Add this frame to the data cube, if its coordinates match those of the primary frame
            datacube.append(self.frames[frame_name]._data)
            plane_name = frame_name + " [frame]"
            plane_names.append(plane_name)

            # Add the name of the frame to the header
            if add_plane_names: header["PLANE" + str(plane_index)] = plane_name

            # Increment the plane index
            plane_index += 1

            last_addition = "Frame"

        # Add the masks
        if add_masks:

            # Export all masks to the specified file
            for mask_name in self.masks:

                # Inform the user that this mask will be saved to the image file
                #log.debug("Exporting the " + mask_name + " mask to '" + path + "' ...")

                # Add this mask to the data cube
                datacube.append(self.masks[mask_name].astype(int))
                plane_name = mask_name + " [mask]"
                plane_names.append(plane_name)

                # Add the name of the mask to the header
                if add_plane_names: header["PLANE" + str(plane_index)] = plane_name

                # Increment the plane index
                plane_index += 1

                last_addition = "Mask"

        # Add the segmentation maps
        if add_segments:

            # Export all segmentation maps
            for segments_name in self.segments:

                # Inform the user that this segmentation map will be saved to the image file
                #log.debug("Exporting the " + segments_name + " segmentation map to '" + path + "' ...")

                # Add this segmentation map to the data cube
                datacube.append(self.segments[segments_name].data)
                plane_name = segments_name + " [segments]"
                plane_names.append(plane_name)

                # Add the name of the segmentation map to the header
                if add_plane_names: header["PLANE" + str(plane_index)] = plane_name

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

        # Add extra info
        if extra_info is not None:
            for key in extra_info:
                header[key] = extra_info[key]

        # Set plane information
        if plane_index > 1:
            header["NAXIS"] = 3
            header["NAXIS3"] = plane_index
        else:
            header.remove("PLANE0")
            header["PTSCLS"] = last_addition  # if only one Frame or Mask or SegmentationMap has been added

        # Set unit, FWHM and filter description
        if self.unit is not None:
            header.set("SIGUNIT", represent_unit(self.unit), "Unit of the map")
            header.set("PHYSTYPE", self.unit.physical_type, "Physical type of the unit")
        if self.fwhm is not None: header.set("FWHM", self.fwhm.to("arcsec").value, "[arcsec] FWHM of the PSF")

        # Set PSF FILTER
        if add_psf_filter and self.psf_filter is not None: header.set("PSFFLTR", str(self.psf_filter), "Filter to which the PSF of the frame corresponds")

        # Set filter
        if add_filter:
            if self.filter is not None: header.set("FILTER", str(self.filter), "Filter used for this observation")
            else: header.set("FILTER", "n/a", "This image does not correspond to a certain observational filter")

        # Pixelscale
        if self.wcs is None and self.pixelscale is not None:

            # Angular
            if self.has_angular_pixelscale:

                header.set("XPIXSIZE", repr(self.pixelscale.x.to("arcsec").value), "[arcsec] Pixelscale for x axis")
                header.set("YPIXSIZE", repr(self.pixelscale.y.to("arcsec").value), "[arcsec] Pixelscale for y axis")

            # Physical
            elif self.has_physical_pixelscale:

                header.set("XPIXSIZE", repr(self.pixelscale.x.to("pc").value), "[pc] Pixelscale for x axis")
                header.set("YPIXSIZE", repr(self.pixelscale.y.to("pc").value), "[pc] Pixelscale for y axis")

            # Invalid pixelscale
            else: raise RuntimeError("We shouldn't get here")

        # Set distance
        if self.distance is not None: header.set("DISTANCE", repr(self.distance.to("Mpc").value), "[Mpc] Distance to the object")

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
            # TOO MUCH OUTPUT
            #log.debug("Setting the FWHM of the " + frame_name + " frame to " + str(fwhm) + " ...")

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
            # TOO MUCH OUTPUT
            #log.debug("Setting the filter of the " + frame_name + " frame to " + str(fltr) + " ...")

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
            # TOO MUCH OUTPUT
            #log.debug("Setting the coordinate system of the " + frame_name + " frame ...")

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

    def converted_to(self, to_unit, distance=None, density=False, brightness=False, density_strict=False, brightness_strict=False):

        """
        This function ...
        :param to_unit:
        :param distance:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :return:
        """

        new = self.copy()
        new.convert_to(to_unit, distance=distance, density=density, brightness=brightness, density_strict=density_strict, brightness_strict=brightness_strict)
        return new

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

        # Parse "to unit": VERY IMPORTANT, BECAUSE DOING SELF.UNIT = TO_UNIT WILL OTHERWISE REPARSE AND WILL BE OFTEN INCORRECT!! (NO DENSITY OR BRIGHTNESS INFO)
        to_unit = PhotometricUnit(to_unit, density=density, brightness=brightness, brightness_strict=brightness_strict, density_strict=density_strict)

        # Already in the correct unit
        if to_unit == self.unit:
            log.debug("Image is already in the desired unit")
            return 1.

        # Inform the user
        log.debug("Converting the unit of the image from " + tostr(self.unit, add_physical_type=True) + " to " + tostr(to_unit, add_physical_type=True) + " ...")

        # Set the distance
        if distance is None: distance = self.distance

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

    def convolved(self, kernel, allow_huge=True, fft=True):

        """
        This function ...
        :param kernel:
        :param allow_huge:
        :param fft:
        :return:
        """

        new = self.copy()
        new.convolve(kernel, allow_huge=allow_huge, fft=fft)
        return new

    # -----------------------------------------------------------------

    def rebin(self, reference_wcs, exact=False, parallel=True, convert=None, add_footprint=True, silent=False):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :param add_footprint:
        :param silent:
        """

        # Check whether the image has a WCS
        if not self.has_wcs: raise RuntimeError("Cannot rebin a frame without coordinate system")

        # Create a copy of the current wcs
        original_wcs = self.wcs.deepcopy()

        footprint = None

        # Loop over all currently selected frames
        for frame_name in self.frames:

            # Inform the user
            if not silent: log.debug("Rebinning the '" + frame_name + "' frame ...")

            # Rebin this frame (the reference wcs is automatically set in the new frame)
            footprint = self.frames[frame_name].rebin(reference_wcs, exact=exact, parallel=parallel, convert=convert)

        # Loop over the masks
        for mask_name in self.masks:

            # Inform the user
            if not silent: log.debug("Rebinning the '" + mask_name + "' mask ...")

            # Rebin
            self.masks[mask_name].rebin(reference_wcs)

        # Loop over the segmentation maps
        for segments_name in self.segments:

            # Infomr the user
            if not silent: log.debug("Rebinning the '" + segments_name + "' segmentation map ...")

            # Rebin
            self.segments[segments_name].rebin(reference_wcs)

        # If there was any frame or mask, we have a footprint
        if footprint is not None and add_footprint:

            # Add mask for padded pixels after rebinning
            # this mask now covers pixels added to the frame after rebinning plus more (radius 10 pixels)
            padded = Mask(footprint < 0.9, wcs=self.wcs, pixelscale=self.pixelscale)
            padded.disk_dilate(radius=10)
            self.add_mask(padded, "padded")

        # Return the original WCS
        return original_wcs

    # -----------------------------------------------------------------

    def rebinned(self, reference_wcs, exact=False, parallel=True, convert=None, add_footprint=False, silent=False):

        """
        This function ...
        :param reference_wcs:
        :param exact:
        :param parallel:
        :param convert:
        :param add_footprint:
        :param silent:
        :return:
        """

        new = self.copy()
        new.rebin(reference_wcs, exact=exact, parallel=parallel, convert=convert, add_footprint=add_footprint, silent=silent)
        return new

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

        # Loop over all segmentation maps
        for segments_name in self.segments:

            # Inform the user
            log.debug("Cropping the " + segments_name + " segmentation map ...")

            # Crop
            self.segments[segments_name] = self.segments[segments_name][y_min:y_max, x_min:x_max]

        # Loop over all regions
        for region_name in self.regions:

            # Inform the user
            log.debug("Cropping the " + region_name + " region ...")

            # Crop the region
            self.regions[region_name] = self.regions[region_name].cropped(x_min, x_max, y_min, y_max)

    # -----------------------------------------------------------------

    def downsample(self, factor, order=3, dilate_masks=False, dilate_nans=True, dilate_infs=True, convert=None):

        """
        This function ...
        :param factor:
        :param order:
        :param dilate_masks:
        :param dilate_nans:
        :param dilate_infs:
        :param convert:
        :return:
        """

        new_wcs = None
        new_shape = None

        # Are there frames?
        if self.has_frames:

            # Debugging
            log.debug("Downsampling the '" + self.primary_name + "' frame [primary] ...")

            # Downsample the primary frame
            new_wcs = self.primary.downsample(factor, order=order, dilate_nans=dilate_nans, dilate_infs=dilate_infs, convert=convert)
            new_shape = self.primary.shape

            # Loop over all other frames
            for frame_name in self.frames:

                # Skip the primary frame
                if self.is_primary(frame_name): continue

                # Inform the user
                log.debug("Downsampling the '" + frame_name + "' frame ...")

                # Rebin this frame to the new coordinate system
                if new_wcs is not None: self.frames[frame_name].rebin(new_wcs, convert=convert)
                else: self.frames[frame_name].downsample(factor, order=order, dilate_nans=dilate_nans, dilate_infs=dilate_infs, convert=convert)

                # Check shape
                if new_shape != self.frames[frame_name].shape: raise RuntimeError("Something went wrong")

        # Are there masks?
        for mask_name in self.masks:

            # Inform the user
            log.debug("Downsampling the '" + mask_name + "' mask ...")

            # Rebin this mask
            if new_wcs is not None: self.masks[mask_name].rebin(new_wcs, dilate=dilate_masks)
            else: new_wcs = self.masks[mask_name].downsample(factor, dilate=dilate_masks)
            if new_shape is None: new_shape = self.masks[mask_name].shape

            # Check shape
            if new_shape != self.masks[mask_name].shape: raise RuntimeError("Something went wrong")

        # Loop over all segmentation maps
        for segments_name in self.segments:

            # Inform the user
            #log.debug("Downsampling the '" + segments_name + "' segmentation map ...")

            raise NotImplementedError("Not implemented")

        # Loop over all regions
        for region_name in self.regions:

            # Inform the user
            #log.debug("Downsampling the " + region_name + " region ...")

            raise NotImplementedError("Not implemented")

    # -----------------------------------------------------------------

    def intersect_masks(self):

        """
        This function ...
        :return:
        """

        # Intersection
        return intersection(*self.masks.values())

    # -----------------------------------------------------------------

    def unite_masks(self):

        """
        This function ...
        :return:
        """

        # Calculate the total mask
        return union(*self.masks.values())

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

    def is_identical(self, other):

        """
        This function ...
        :param other:
        :return:
        """

        # Not an image
        if not isinstance(other, Image): return False

        # First checks
        if self.nframes != other.nframes: return False
        if self.nregions != other.nregions: return False
        if self.nmasks != other.nmasks: return False
        if self.nsegments != other.nsegments: return False

        # Check names (with their order)
        if self.frame_names != other.frame_names: return False
        if self.mask_names != other.mask_names: return False
        if self.regions_names != other.regions_names: return False
        if self.segments_names != other.segments_names: return False

        # Check the frames
        for frame_name in self.frame_names:
            if not self.frames[frame_name].is_identical(other.frames[frame_name]): return False

        # Check the regions
        for regions_name in self.regions_names:
            if self.regions[regions_name] != other.regions[regions_name]: return False

        # Check the masks
        for mask_name in self.mask_names:
            if not self.masks[mask_name].is_identical(other.masks[mask_name]): return False

        # Check the segmentation maps
        for segments_name in self.segments_names:
            if not self.segments[segments_name].is_identical(other.segments[segments_name]): return False

        # All checks passed
        return True

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

    def rename_frame(self, old_name, new_name, keep_position=False, silent=False):

        """
        This function renames a frame
        :param old_name:
        :param new_name:
        :param keep_position:
        :param silent:
        :return:
        """

        # Get the current position
        if keep_position: current_position = self.frames.keys().index(old_name)
        else: current_position = None

        # Remove the frame from the dictionary of frames and re-add it under a different key
        self.frames[new_name] = self.frames.pop(old_name)

        # Move the frame to the original position
        if keep_position: self.move_frame(new_name, current_position, silent=silent)

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

    def load_frames(self, path, index=None, name=None, description=None, always_call_first_primary=True,
                    hdulist_index=0, no_filter=False, silent=False, density=False, brightness=False,
                    density_strict=False, brightness_strict=False, indices=None, absolute_index_names=True):

        """
        This function ...
        :param path:
        :param index:
        :param name:
        :param description:
        :param always_call_first_primary:
        :param hdulist_index:
        :param no_filter:
        :param silent:
        :param density:
        :param brightness:
        :param density_strict:
        :param brightness_strict:
        :param indices:
        :param absolute_index_names:
        :return:
        """

        # Check if the file exists
        if not fs.is_file(path): raise IOError("File '" + path + "' does not exist")

        # Show which image we are importing
        if not silent: log.debug("Reading in file '" + path + "' ...")

        # Load frames
        from . import fits as pts_fits
        try:
            frames, masks, segments, meta = pts_fits.load_frames(path, index, name, description, always_call_first_primary,
                                                           hdulist_index, no_filter, density=density, brightness=brightness,
                                                           density_strict=density_strict, brightness_strict=brightness_strict,
                                                           indices=indices, absolute_index_names=absolute_index_names)
        except pts_fits.DamagedFITSFileError: raise IOError("File is possibly damaged")

        # Set frames, masks and meta information
        for frame_name in frames: self.add_frame(frames[frame_name], frame_name, silent=True)
        for mask_name in masks: self.add_mask(masks[mask_name], mask_name, silent=True)
        for segments_name in segments: self.add_segments(segments[segments_name], segments_name, silent=True)
        for keyword in meta: self.metadata[keyword] = meta[keyword]

    # -----------------------------------------------------------------

    def add_frame(self, frame, name, overwrite=False, silent=False, copy=False):

        """
        This function ...
        :param frame:
        :param name:
        :param overwrite:
        :param silent:
        :param copy:
        :return:
        """

        # Inform the user
        if not silent: log.debug("Adding '" + name + "' to the set of frames ...")

        # Check whether a frame with this name already exists
        if name in self.frames and not overwrite: raise RuntimeError("A frame with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if frame.shape != self.shape: raise ValueError("Frame does not have the correct shape for this image")

        # Create copy?
        if copy: frame = frame.copy()

        # Set the WCS
        if self.wcs is not None: frame.wcs = self.wcs

        # Add the layer to the layers dictionary
        self.frames[name] = frame

    # -----------------------------------------------------------------

    def remove_frames_except(self, *names, **kwargs):

        """
        This function ...
        :param names:
        :param kwargs:
        :return:
        """

        # Get options
        silent = kwargs.pop("silent", False)

        # Loop over all frames
        for frame_name in list(self.frames.keys()):

            # Don't remove the frame with the specified name
            if frame_name in names: continue

            # Remove all other frames
            self.remove_frame(frame_name, silent=silent)

    # -----------------------------------------------------------------

    def remove_frame(self, name, silent=False):

        """
        This function ...
        :param name:
        :param silent:
        :return:
        """

        # Inform the user
        if not silent: log.debug("Removing the '" + name + "' frame ...")

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

    def move_frame(self, name, position, silent=False):

        """
        This function ...
        :param name:
        :param position:
        :param silent:
        :return:
        """

        # Inform the user
        if not silent: log.debug("Moving frame '" + name + "' to position " + str(position) + " in the frame layers ...")

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

    def add_mask(self, mask, name, overwrite=False, silent=False, copy=False):

        """
        This function ...
        :param mask:
        :param name:
        :param overwrite:
        :param silent:
        :param copy:
        :return:
        """

        # Inform the user
        if not silent: log.debug("Adding '" + name + "' to the set of masks ...")

        # Check whether a mask with this name already exists
        if name in self.masks and not overwrite: raise RuntimeError("A mask with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if mask.shape != self.shape: raise ValueError("Mask does not have the correct shape for this image")

        # Copy
        if copy: mask = mask.copy()

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

    def add_segments(self, segments, name, overwrite=False, silent=False, copy=False):

        """
        This function ...
        :param segments:
        :param name:
        :param overwrite:
        :param silent:
        :param copy:
        :return:
        """

        # Inform the user
        if not silent: log.debug("Adding '" + name + "' to the set of segmentation maps ...")

        # Check whether a segmentation map with this name already exists
        if name in self.segments and not overwrite: raise RuntimeError("A segmentation map with this name already exists")

        # Check if the shape matches the shape of this image
        if self.shape is not None:
            if segments.shape != self.shape: raise ValueError("Segmentation map does not have the correct shape for this image")

        # Copy
        if copy: segments = segments.copy()

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
